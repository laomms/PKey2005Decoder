#include "pch.h"
#include "Helper.h"
#include "BigInteger.h"
#include "pKeyCalc.h"
#include "Montgomery.h"
#include "TatePairing.h"
#include "H1Search.h"
#include "PCoreAffinity.h"
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cstring>
#include <cstdlib>
#include <thread>
#include <atomic>
#include <mutex>
#include <cstdint>
#include <Windows.h>
#include <objbase.h>

#ifdef _WIN32
#define DLL_EXPORT __declspec(dllexport)
#else
#define DLL_EXPORT
#endif

// =====================================================================
// CPU tuning notes (i7-13700KF: 8 P-cores + HT = 16 P-threads, 8 E-cores)
//
// The OUTER loop below runs INDEPENDENT groups, so -- unlike the inner
// MITM chunk split -- it is a throughput problem, not a balance problem.
// We therefore:
//   * spread it across all 16 homogeneous P-THREADS (P-cores + their HT
//     siblings), not just the 8 physical P-cores. HT adds ~20-30% here
//     because the 8 MB table probes and the modmul dependency chains
//     leave issue slots the sibling thread can use.
//   * skip the 8 E-cores on purpose: if the *winning* group landed on a
//     Gracemont core it would run ~2x slower and stall the whole result
//     while the P-threads sit idle with nothing left to steal.
//   * dispatch groups from a shared atomic counter (work-stealing), so a
//     P-thread that finishes a doomed group immediately grabs the next
//     one instead of waiting on a static slice.
//   * force the INNER MITM serial in this mode (FreeSlots = 0): with 16
//     groups already in flight there are no spare cores to nest into, and
//     nesting would just oversubscribe. One group == one serial MITM on
//     one pinned P-thread.
//
// The single-group path (targetGroupId set, or a config with one entry)
// does the opposite: 1 outer worker, inner MITM claims all 8 P-cores and
// runs fully parallel -- the warm path that beats a serial validator.
// =====================================================================

// Internal Base64 decoder for Helper::ParseKeyEntriesSimple
static std::vector<unsigned char> InternalBase64Decode(std::string_view b64) {
    static constexpr int kLookup[] = {
        -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
        -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
        -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,62,-1,-1,-1,63,
        52,53,54,55,56,57,58,59,60,61,-1,-1,-1,-1,-1,-1,
        -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,
        15,16,17,18,19,20,21,22,23,24,25,-1,-1,-1,-1,-1,
        -1,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,
        41,42,43,44,45,46,47,48,49,50,51,-1,-1,-1,-1,-1
    };
    int in_len = (int)b64.size();
    int i = 0, in_ = 0;
    unsigned char char_array_4[4], char_array_3[3];
    std::vector<unsigned char> ret;

    while (in_len-- && (b64[in_] != '=') && (std::isalnum(b64[in_]) || b64[in_] == '+' || b64[in_] == '/')) {
        char_array_4[i++] = b64[in_++];
        if (i == 4) {
            for (i = 0; i < 4; i++)
                char_array_4[i] = (unsigned char)kLookup[char_array_4[i]];
            char_array_3[0] = (char_array_4[0] << 2) + ((char_array_4[1] & 0x30) >> 4);
            char_array_3[1] = ((char_array_4[1] & 0xf) << 4) + ((char_array_4[2] & 0x3c) >> 2);
            char_array_3[2] = ((char_array_4[2] & 0x3) << 6) + char_array_4[3];
            for (i = 0; i < 3; i++) ret.push_back(char_array_3[i]);
            i = 0;
        }
    }
    if (i) {
        for (int j = i; j < 4; j++) char_array_4[j] = 0;

        for (int j = 0; j < i; j++) {
            char_array_4[j] = (unsigned char)kLookup[char_array_4[j]];
        }

        char_array_3[0] = (char_array_4[0] << 2) + ((char_array_4[1] & 0x30) >> 4);
        char_array_3[1] = ((char_array_4[1] & 0xf) << 4) + ((char_array_4[2] & 0x3c) >> 2);
        char_array_3[2] = ((char_array_4[2] & 0x3) << 6) + char_array_4[3];

        for (int j = 0; j < i - 1; j++) ret.push_back(char_array_3[j]);
    }
    return ret;
}

#ifdef _WIN32
// One affinity mask per logical P-THREAD (splits each physical P-core mask
// into its individual HT-lane bits). On the 13700KF this yields 16 masks.
// Empty on non-hybrid / detection failure -> caller falls back.
static std::vector<DWORD_PTR> GetPerformanceThreadMasks() {
    std::vector<DWORD_PTR> out;
    for (DWORD_PTR m : GetPerformanceCoreMasks()) {
        while (m) {
            DWORD_PTR lsb = m & (DWORD_PTR)(~m + 1); // lowest set bit
            out.push_back(lsb);
            m &= (m - 1);
        }
    }
    return out;
}
#endif

// =====================================================================
// Shared scan core (used by BOTH the XML export and the bink-list export)
// =====================================================================

// Status codes returned by PkeyVerifyBinks. Only PKEY_OK (1) means "valid":
// test for == 1, not for "non-zero" (the error codes are negative).
#define PKEY_OK              1   // key is valid; outputs are filled
#define PKEY_NO_MATCH        0   // every usable bink rejected the key
#define PKEY_E_ARGS         -1   // null pointer or rawKeySize != 16
#define PKEY_E_GROUP_LIST   -2   // group list: wrong type, bad length, or empty
#define PKEY_E_BINK_LIST    -3   // bink list: wrong type, or length does not fit the group count
#define PKEY_E_NO_USABLE    -4   // no entry passed the bink checks
#define PKEY_E_NO_GROUP     -5   // targetGroupId is not in the group list
#define PKEY_E_INTERNAL     -9   // unexpected exception

// Both lists handed to PkeyVerifyBinks start with the same 8-byte header
// (little-endian), followed directly by the items:
//
//   offset 0   uint32  Type         PKEY_LIST_GROUPS (1) or PKEY_LIST_BINKS (2)
//   offset 4   uint32  TotalLength  size of the WHOLE list in bytes,
//                                   these 8 header bytes included
//   offset 8   items
//       Type 1: int32 group ids, one per bink
//       Type 2: the binks back to back, all the same size
//
// The entry count comes from the group list ((TotalLength - 8) / 4), and the
// size of one bink from the bink list ((TotalLength - 8) / count).
#define PKEY_LIST_GROUPS     1
#define PKEY_LIST_BINKS      2
#define PKEY_LIST_HEADER     8

#if defined(_WIN32) && defined(_MSC_VER)
#define PKEY_CALL __cdecl        // explicit, so /Gz or /Gr cannot change it
#else
#define PKEY_CALL
#endif

struct ParsedEntry {
    int groupId;
    int index;      // position in the caller's list
    PubKey pubKey;
};

// Abort / FreeSlots in H1Search are process-wide, so two scans running at
// the same time would cancel each other and fight over the thread budget.
// One scan at a time; a second caller simply waits.
static std::mutex g_scanGate;

// Runs the key against every entry and reports the first one that accepts it.
// 'entries' must be non-empty and already filtered to the wanted group(s).
static bool ScanEntries(
    const std::vector<ParsedEntry>& parsedEntries,
    const std::vector<unsigned char>& bEncryptArray,
    unsigned char* outUid8Bytes,
    int& outGroupId,
    int& outIndex
) {
    std::lock_guard<std::mutex> scanLock(g_scanGate);

    // Fresh run: clear any leftover abort flag BEFORE any worker starts,
    // so it can never race with the winner raising it below.
    H1Search::Abort.store(0, std::memory_order_relaxed);

    const bool singleGroup = (parsedEntries.size() == 1);

    // Inner-MITM budget:
    //   single group  -> let the ONE search grab every P-core (parallel)
    //   many groups    -> force each search serial; throughput comes from
    //                     running many groups at once (no nested spawn).
    H1Search::FreeSlots.store(singleGroup ? H1Search::TotalSlots() : 0,
        std::memory_order_relaxed);

    std::atomic<bool>   foundValid(false);
    std::atomic<size_t> nextIdx(0);
    int resolvedGroupId = 0;
    int resolvedIndex = -1;
    std::vector<unsigned char> resolvedUid;

    // Worker set.
    unsigned int numWorkers;
#ifdef _WIN32
    std::vector<DWORD_PTR> pinMasks;
    if (!singleGroup) {
        pinMasks = GetPerformanceThreadMasks();          // 16 on 13700KF
        if (pinMasks.empty()) {                          // non-hybrid fallback
            unsigned int hc = std::thread::hardware_concurrency();
            if (hc == 0) hc = 4;
            for (unsigned int i = 0; i < hc && i < 8 * sizeof(DWORD_PTR); ++i)
                pinMasks.push_back((DWORD_PTR)1 << i);
        }
        numWorkers = (unsigned int)pinMasks.size();
    }
    else {
        numWorkers = 1;                                  // inner MITM is the parallel one
    }
#else
    numWorkers = singleGroup ? 1u
        : (std::max)(1u, std::thread::hardware_concurrency());
#endif
    if (numWorkers > (unsigned int)parsedEntries.size())
        numWorkers = (unsigned int)parsedEntries.size();

    std::vector<std::thread> workers;
    workers.reserve(numWorkers);

    for (unsigned int t = 0; t < numWorkers; ++t) {
        workers.emplace_back([&, t]() {
#ifdef _WIN32
            if (!singleGroup && t < pinMasks.size()) {
                SetThreadAffinityMask(GetCurrentThread(), pinMasks[t]);
            }
            SetThreadPriority(GetCurrentThread(), THREAD_PRIORITY_HIGHEST);
#endif
            for (;;) {
                if (foundValid.load(std::memory_order_relaxed)) break;
                size_t i = nextIdx.fetch_add(1, std::memory_order_relaxed);
                if (i >= parsedEntries.size()) break;

                const auto& entry = parsedEntries[i];
                std::string actPkeyConfig;
                std::vector<unsigned char> h1Coeffs;
                std::vector<unsigned char> uid;

                // An exception escaping a std::thread terminates the host
                // process, so a bink that makes the math throw just counts
                // as "this group did not accept the key".
                bool success = false;
                try {
                    success = PKeyCalc::TryParsedPubKey(entry.pubKey, bEncryptArray, actPkeyConfig, h1Coeffs, uid);
                }
                catch (...) {
                    success = false;
                }
                if (success) {
                    bool expected = false;
                    if (foundValid.compare_exchange_strong(expected, true)) {
                        // Only the CAS winner writes these; safe without a lock.
                        resolvedUid = uid;
                        resolvedGroupId = entry.groupId;
                        resolvedIndex = entry.index;
                        // Stop every other in-flight MITM at its next poll,
                        // instead of grinding its full 2^18 scan to a miss.
                        H1Search::Abort.store(1, std::memory_order_relaxed);
                    }
                    break;
                }
            }
            });
    }

    for (auto& worker : workers) {
        if (worker.joinable()) worker.join();
    }

    if (!foundValid.load()) return false;

    std::memset(outUid8Bytes, 0, 8);
    std::memcpy(outUid8Bytes, resolvedUid.data(), (std::min)(resolvedUid.size(), (size_t)8));
    outGroupId = resolvedGroupId;
    outIndex = resolvedIndex;
    return true;
}

// Reads the header of a msft:rm/algorithm/pkey/2005 bink and returns the
// number of bytes the bink says it occupies, or 0 when it is not a bink this
// engine can evaluate. PubKeyParser::Parse trusts the sizes stored in the
// blob, so nothing reaches it without passing here first.
//
// Layout (sM = modulus bytes, sO = order bytes, e1/e2 = extension degrees):
//   0   magic 0x44556677        8  fieldDataSize        12  magic 0x00112233
//   20  0, sM, sO               23 e1 (u32)             27  e2 (u32)
//   60  H1 bases[sM], modulus[sM], order[sO], K3 minpoly[e1+1],
//       K6 minpoly[e2+1], reserved[2*sM], curve a[sM], curve b[sM]
//   fieldDataSize + 12:  sM points (x, y in Fp^e1), then the pairing value
//
// The arithmetic in BigInteger.h / Montgomery.h / H1Search.h is written for
// one shape only: a 14-byte prime, Fp3 = Fp[u]/(u^3+u+4), Fp6 = Fp3[y]/(y^2+2)
// and 14 points. A bink of any other shape is refused here instead of being
// fed to code that would index out of bounds or return a wrong answer.
static size_t BinkRequiredSize(const unsigned char* d, size_t avail) {
    if (!d || avail < 60) return 0;

    auto U16 = [&](size_t off) -> uint32_t {
        return (uint32_t)d[off] | ((uint32_t)d[off + 1] << 8);
        };
    auto U32 = [&](size_t off) -> uint32_t {
        return (uint32_t)d[off] | ((uint32_t)d[off + 1] << 8) |
            ((uint32_t)d[off + 2] << 16) | ((uint32_t)d[off + 3] << 24);
        };

    // Both headers: magic + version word 0x0100, as pidgenx checks them.
    if (U32(0) != 0x44556677u || U16(4) != 0x0100) return 0;
    if (U32(12) != 0x00112233u || U16(16) != 0x0100) return 0;

    // Flags byte and the three unused counts must be zero.
    if (d[20] != 0 || U32(31) != 0 || U32(35) != 0 || U32(39) != 0) return 0;

    // pidgenx takes these from the blob; this engine implements one shape.
    const size_t sM = d[21];        // modulus bytes
    const size_t sO = d[22];        // order bytes
    const size_t e1 = 3, e2 = 2;    // extension degrees
    const size_t n = 14;            // points / H1 digits
    if (sM != 14 || sO == 0 || U32(23) != e1 || U32(27) != e2 || U32(43) != n) return 0;

    // The field data must be exactly as long as its own header describes...
    const size_t fieldSize = n + 50 + sO + e1 + e2 + 10 * sM + 2 * e1 * sM;
    if (U32(8) != fieldSize) return 0;

    // ...and the whole bink exactly header + field data + points + pairing value.
    return 12 + fieldSize + 2 * n * e1 * sM + e1 * e2 * sM;
}

extern "C" {
    // Raw-buffer variant (Zero Disk I/O, no strings, no key decoding).
    //
    // rawKey / rawKeySize:
    //     The already-encoded binary key: exactly the 16 bytes that
    //     Helper::EncodeBinaryKey produces from the key text.
    //     rawKeySize must be 16.
    // xmlData / xmlSize:
    //     The raw bytes of the pkeyconfig XML, exactly as stored on disk.
    // outUid8Bytes:
    //     Caller buffer, 8 bytes.
    // outGroupId:
    //     Optional, may be null.
    // targetGroupId:
    //     0  -> scan every group (dynamic, all P-threads)
    //     >0 -> validate ONLY that group, inner MITM fully parallel
    DLL_EXPORT bool VerifyBinaryKey(
        const unsigned char* rawKey,
        int rawKeySize,
        const unsigned char* xmlData,
        int xmlSize,
        unsigned char* outUid8Bytes,
        int* outGroupId,
        int targetGroupId
    ) {
        try {
            if (!rawKey || rawKeySize != 16 ||
                !xmlData || xmlSize <= 0 ||
                !outUid8Bytes) {
                return false;
            }

            // The key arrives already encoded; just copy the 16 bytes.
            std::vector<unsigned char> bEncryptArray(rawKey, rawKey + 16);

            std::string outerXml(reinterpret_cast<const char*>(xmlData), (size_t)xmlSize);

            std::vector<PublicKeyEntry> pkEntries;
            if (!Helper::ParseKeyEntriesSimple(outerXml, pkEntries, InternalBase64Decode)) {
                return false;
            }

            std::vector<ParsedEntry> parsedEntries;
            parsedEntries.reserve(pkEntries.size());

            for (size_t i = 0; i < pkEntries.size(); ++i) {
                const auto& entry = pkEntries[i];
                try {
                    PubKey k = PubKeyParser::Parse(entry.pubKeyBytes);
                    int gId = 0;
                    try { gId = std::stoi(entry.groupId); }
                    catch (...) { gId = 0; }
                    parsedEntries.push_back({ gId, (int)i, std::move(k) });
                }
                catch (...) {
                    // Skip malformed entries safely
                }
            }

            // Optional single-group fast path.
            if (targetGroupId > 0) {
                std::vector<ParsedEntry> only;
                for (auto& pe : parsedEntries) {
                    if (pe.groupId == targetGroupId) { only.push_back(std::move(pe)); break; }
                }
                parsedEntries.swap(only);
            }

            if (parsedEntries.empty()) {
                return false;
            }

            int resolvedGroupId = 0, resolvedIndex = -1;
            if (!ScanEntries(parsedEntries, bEncryptArray, outUid8Bytes, resolvedGroupId, resolvedIndex)) {
                return false;
            }
            if (outGroupId) *outGroupId = resolvedGroupId;
            return true;
        }
        catch (...) {
            return false;
        }
    }

    // Bink-list variant: no XML, no Base64, no disk. The caller extracts
    // the public keys from <pkc:PublicKeys> itself and hands them over raw,
    // as two lists that each carry their own type and length (see the
    // PKEY_LIST_* layout above).
    //
    // rawKey / rawKeySize:
    //     The 16 encoded key bytes, as for VerifyBinaryKey. rawKeySize must be 16.
    // groupList:
    //     Type 1 list: the <pkc:GroupId> of every entry, as int32.
    // binkList:
    //     Type 2 list: the Base64-decoded <pkc:PublicKeyValue> of every
    //     msft:rm/algorithm/pkey/2005 entry, in the same order as groupList.
    //     All binks must have the same size (1579 for the known 2005 binks);
    //     that size is (TotalLength - 8) / count, and every bink's own header
    //     must describe exactly that many bytes or the entry is skipped.
    // outUid8Bytes:
    //     Caller buffer, 8 bytes. Written only when the result is PKEY_OK.
    // outGroupId / outIndex:
    //     Optional, may be null. Group id and list position of the bink
    //     that accepted the key (0 / -1 when there is no match).
    // targetGroupId:
    //     0  -> scan every entry
    //     >0 -> validate ONLY the first entry with that group id
    //
    // Returns PKEY_OK (1) for a valid key, PKEY_NO_MATCH (0) when no bink
    // accepts it, or a negative PKEY_E_* code when the input is unusable.
    DLL_EXPORT int PKEY_CALL VerifyBinks(
        const unsigned char* rawKey,
        int rawKeySize,
        const unsigned char* groupList,
        const unsigned char* binkList,
        unsigned char* outUid8Bytes,
        int* outGroupId,
        int* outIndex,
        int targetGroupId
    ) {
        try {
            if (outGroupId) *outGroupId = 0;
            if (outIndex) *outIndex = -1;

            if (!rawKey || rawKeySize != 16 ||
                !groupList || !binkList ||
                !outUid8Bytes) {
                return PKEY_E_ARGS;
            }

            // Byte-wise reads: the caller's buffers need not be aligned.
            auto U32 = [](const unsigned char* p) -> uint32_t {
                return (uint32_t)p[0] | ((uint32_t)p[1] << 8) |
                    ((uint32_t)p[2] << 16) | ((uint32_t)p[3] << 24);
                };

            // ---- group list: gives the entry count ----
            const uint32_t groupBytes = U32(groupList + 4);
            if (U32(groupList) != PKEY_LIST_GROUPS ||
                groupBytes <= PKEY_LIST_HEADER ||
                groupBytes > 0x7FFFFFFFu ||
                (groupBytes - PKEY_LIST_HEADER) % 4 != 0) {
                return PKEY_E_GROUP_LIST;
            }
            const size_t count = (groupBytes - PKEY_LIST_HEADER) / 4;

            // ---- bink list: must hold exactly 'count' equal-sized binks ----
            const uint32_t binkBytes = U32(binkList + 4);
            if (U32(binkList) != PKEY_LIST_BINKS ||
                binkBytes <= PKEY_LIST_HEADER ||
                binkBytes > 0x7FFFFFFFu ||
                (binkBytes - PKEY_LIST_HEADER) % count != 0) {
                return PKEY_E_BINK_LIST;
            }
            const size_t binkSize = (binkBytes - PKEY_LIST_HEADER) / count;

            const unsigned char* groupItems = groupList + PKEY_LIST_HEADER;
            const unsigned char* binkItems = binkList + PKEY_LIST_HEADER;

            std::vector<unsigned char> bEncryptArray(rawKey, rawKey + 16);

            std::vector<ParsedEntry> parsedEntries;
            bool sawTarget = false;

            for (size_t i = 0; i < count; ++i) {
                const int groupId = (int)U32(groupItems + i * 4);
                if (targetGroupId > 0) {
                    if (groupId != targetGroupId) continue;
                    if (sawTarget) break;        // first entry of that group only
                    sawTarget = true;
                }

                const unsigned char* bink = binkItems + i * binkSize;
                if (BinkRequiredSize(bink, binkSize) != binkSize) continue;

                try {
                    std::vector<unsigned char> blob(bink, bink + binkSize);
                    PubKey k = PubKeyParser::Parse(blob);
                    parsedEntries.push_back({ groupId, (int)i, std::move(k) });
                }
                catch (...) {
                    // Different base field than the binks seen so far, or
                    // a value the parser refuses: skip this entry.
                }
            }

            if (targetGroupId > 0 && !sawTarget) {
                return PKEY_E_NO_GROUP;
            }
            if (parsedEntries.empty()) {
                return PKEY_E_NO_USABLE;
            }

            int resolvedGroupId = 0, resolvedIndex = -1;
            if (!ScanEntries(parsedEntries, bEncryptArray, outUid8Bytes, resolvedGroupId, resolvedIndex)) {
                return PKEY_NO_MATCH;
            }
            if (outGroupId) *outGroupId = resolvedGroupId;
            if (outIndex) *outIndex = resolvedIndex;
            return PKEY_OK;
        }
        catch (...) {
            return PKEY_E_INTERNAL;
        }
    }

    // Convenience variant: key as text, config as a file path.
    // Encodes the key to its 16 raw bytes, reads the file, scans every group.
    DLL_EXPORT bool VerifyKey(
        const char* cdKeyStr,
        const char* configFilePath,
        unsigned char* outUid8Bytes,
        int* outGroupId
    )
    {
        try {
            if (!cdKeyStr ||
                !configFilePath ||
                !outUid8Bytes ||
                !outGroupId) {
                return false;
            }

            *outGroupId = 0;

            // Key text -> 16 raw bytes (VerifyKeyBytes no longer does this).
            std::vector<unsigned char> rawKey;
            try {
                rawKey = Helper::EncodeBinaryKey(std::string(cdKeyStr));
            }
            catch (const std::exception&) {
                return false;
            }
            if (rawKey.size() != 16) {
                return false;
            }

            std::ifstream file(configFilePath, std::ios::binary);
            if (!file.is_open()) {
                return false;
            }

            std::stringstream buffer;
            buffer << file.rdbuf();

            std::string outerXml = buffer.str();

            unsigned char uid[8] = {};
            int groupId = 0;

            const bool success = VerifyBinaryKey(
                rawKey.data(),
                static_cast<int>(rawKey.size()),
                reinterpret_cast<const unsigned char*>(outerXml.data()),
                static_cast<int>(outerXml.size()),
                uid,
                &groupId,
                0
            );

            if (!success) {
                return false;
            }

            std::memcpy(outUid8Bytes, uid, 8);
            *outGroupId = groupId;

            return true;
        }
        catch (...) {
            return false;
        }
    }
}