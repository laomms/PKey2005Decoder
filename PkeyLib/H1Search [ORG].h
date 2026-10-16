#pragma once
#include <vector>
#include <cstdint>
#include <algorithm>
#include <string>
#include <thread>
#include <atomic>
#include <chrono>
#include <stdexcept>
#include "BigInteger.h"
#include "Montgomery.h"
#include "TatePairing.h"
#include "PCoreAffinity.h"

class H1Search {
public:
    inline static std::atomic<int> Abort{ 0 };
    inline static std::atomic<long> ProbeCalls{ 0 };
    inline static std::atomic<long> ProbeHashHits{ 0 };
    inline static double ProbeProdMms = 0.0;

    // v2: cache the P-core masks once. If the CPU is hybrid we run the MITM
    // on the P-cores only (identical cores -> equal chunks balance perfectly,
    // no HT contention, no E-core stragglers). If not hybrid, masks is empty
    // and everything falls back to the previous hardware_concurrency path.
    inline static std::vector<DWORD_PTR> PCoreMasks = GetPerformanceCoreMasks();

    // ---- Shared global thread budget -------------------------------------
    // Drawn by Solve when threads<=0. The DLL entry points set this to
    // TotalSlots() for a single-key search (inner MITM parallel) or to 0 for
    // the many-groups scan (inner MITM serial; throughput comes from running
    // many groups at once). See Export.cpp.
    inline static std::atomic<int> FreeSlots{ -1 };

    static int TotalSlots() {
        if (!PCoreMasks.empty()) return (int)PCoreMasks.size();
        int hc = (std::max)(1, (int)std::thread::hardware_concurrency());
        return hc;
    }

    static int ClaimSlots(int want) {
        int expected = FreeSlots.load(std::memory_order_relaxed);
        if (expected < 0) { // lazy init on first use
            int init = TotalSlots();
            FreeSlots.compare_exchange_strong(expected, init, std::memory_order_relaxed);
            expected = FreeSlots.load(std::memory_order_relaxed);
        }
        for (;;) {
            int avail = (std::max)(0, expected);
            int take = (std::min)(want, avail);
            if (take <= 1) {
                return 1;
            }
            if (FreeSlots.compare_exchange_weak(expected, expected - take,
                std::memory_order_acq_rel, std::memory_order_relaxed)) {
                return take;
            }
        }
    }

    static void ReleaseSlots(int claimed) {
        if (claimed > 1) FreeSlots.fetch_add(claimed, std::memory_order_acq_rel);
    }

    FORCE_INLINE static uint32_t HashM(const Fp6m& v) {
        uint64_t mix = v.R.C0.A0 ^ v.R.C1.A0 ^ v.R.C2.A0 ^ v.I.C0.A0 ^ v.I.C1.A0 ^ v.I.C2.A0;
        // Apply a fast finalizer (MurmurHash3 style avalanche) to prevent clustering
        mix ^= mix >> 33;
        mix *= 0xff51afd7ed558ccdULL;
        mix ^= mix >> 33;
        return (uint32_t)(mix ^ (mix >> 32));
    }

private:
    FORCE_INLINE static Fp6m ProdM(const std::vector<Fp6m>& g, const std::vector<int>& lo, const std::vector<int>& hi, int from, int to, int64_t s, int* digsOut) {
        Fp6m v = Fp6mOp::One();
        int64_t t = s;
        int count = to - from + 1;
        for (int i = 0; i < count; i++) {
            int radix = hi[from + i] - lo[from + i] + 1;
            int d = (int)(t % radix); t /= radix;
            digsOut[i] = d;
            BigInteger expVal(lo[from + i] + d);
            v = Fp6mOp::Mul6(v, Fp6mOp::PowBig(g[from + i], expVal));
        }
        return v;
    }

    template <typename Body>
    static void RunChunks(int nth, int64_t total, Body body) {
        int nchunk = (std::max)(1, nth);
        int64_t chunk = (total + nchunk - 1) / nchunk;
        if (nchunk == 1) {
            body((int64_t)0, total);
            return;
        }
        std::vector<std::thread> workers;
        workers.reserve(nchunk);
        for (int c = 0; c < nchunk; c++) {
            int64_t begin = (int64_t)c * chunk;
            int64_t end = (std::min)(begin + chunk, total);
            workers.emplace_back([&, c, begin, end]() {
                PinToPerformanceCore(PCoreMasks, (size_t)c);
                body(begin, end);
            });
        }
        for (auto& t : workers) t.join();
    }

public:
    static std::vector<int> Solve(const std::vector<Fp6>& g, const Fp6& target, const std::vector<int>& lo, const std::vector<int>& hi, std::string& note, int threads = 0) {
        // NOTE: Abort is NOT reset here. Resetting per-Solve would let a
        // late-starting search wipe out the abort signal a winning thread
        // raised for its siblings. The DLL entry points reset it once, before
        // any worker starts (see Export.cpp).
        int n = (int)lo.size();
        // 'order' below and the R/L split are hard-wired for 14 bases. Anything
        // else would index out of bounds, so refuse instead of corrupting memory.
        if (n != 14 || (int)hi.size() != n || (int)g.size() != n) {
            note = "unsupported key shape";
            return {};
        }
        int claimed = 1;
        int nth;
        if (threads > 0) {
            nth = threads;
        } else {
            claimed = ClaimSlots(TotalSlots());
            nth = claimed;
        }
        struct SlotGuard { int c; ~SlotGuard() { ReleaseSlots(c); } } slotGuard{ claimed };
        Fpm::Init();

        std::vector<Fp6> ginv(n);
        {
            std::vector<Fp6> pref(n + 1);
            pref[0] = Fp6::One();
            for (int k = 0; k < n; k++) pref[k + 1] = pref[k] * g[k];
            Fp6 acc = pref[n].Inverse();
            for (int k = n - 1; k >= 0; k--) {
                ginv[k] = acc * pref[k];
                acc = acc * g[k];
            }
        }

        std::vector<Fp6m> gm(n), gim(n);
        for (int k = 0; k < n; k++) {
            gm[k] = Fp6m::FromFp6(g[k]);
            gim[k] = Fp6m::FromFp6(ginv[k]);
        }
        Fp6m tgt = Fp6m::FromFp6(target);

        std::vector<int> order = { 0, 1, 2, 3, 4, 5, 13, 6, 7, 8, 9, 10, 11, 12 };
        std::vector<int> loP(n), hiP(n);
        std::vector<Fp6m> gmP(n), gimP(n);
        for (int i = 0; i < n; i++) {
            loP[i] = lo[order[i]];
            hiP[i] = hi[order[i]];
            gmP[i] = gm[order[i]];
            gimP[i] = gim[order[i]];
        }

        const int R0 = 0, R1 = 6;
        int rn = R1 - R0 + 1;
        int64_t rsize = 1;
        for (int k = R0; k <= R1; k++) rsize *= (hiP[k] - loP[k] + 1);

        std::vector<int> actR;
        for (int i = 0; i < rn; i++) {
            if (hiP[R0 + i] > loP[R0 + i]) actR.push_back(i);
        }

        std::vector<Fp6m> rwrap(rn);
        for (int i = 0; i < rn; i++) {
            BigInteger expVal(hiP[R0 + i] - loP[R0 + i]);
            rwrap[i] = Fp6mOp::PowBig(gimP[R0 + i], expVal);
        }

        // ---- Reused build-side key buffer (per calling thread) --------------
        // Every slot [0,rsize) is overwritten by the build loop below, so no
        // zero-fill is needed -- only a size guarantee. rsize is constant for
        // all 2005 configs, so after the first Solve this never reallocates.
        static thread_local std::vector<uint64_t> keysBuf;
        if ((int64_t)keysBuf.size() < rsize) keysBuf.resize(rsize);
        std::vector<uint64_t>& keys = keysBuf;

        RunChunks(nth, rsize, [&](int64_t begin, int64_t end) {
            if (begin >= end) return;
            int d0[16];
            Fp6m prod = ProdM(gmP, loP, hiP, R0, R1, begin, d0);
            int dg[16];
            for (int i = 0; i < rn; i++) dg[i] = d0[i];

            for (int64_t s = begin; s < end; s++) {
                if ((s & 0x3FF) == 0 && Abort.load(std::memory_order_relaxed) != 0) return;
                keys[s] = ((uint64_t)HashM(prod) << 32) | (uint64_t)s;
                for (size_t _j = 0; _j < actR.size(); _j++) {
                    int i = actR[_j];
                    if (dg[i] < hiP[R0 + i] - loP[R0 + i]) {
                        dg[i]++;
                        prod = Fp6mOp::Mul6(prod, gmP[R0 + i]);
                        break;
                    }
                    dg[i] = 0;
                    prod = Fp6mOp::Mul6(prod, rwrap[i]);
                }
            }
            });

        // If a sibling search already won while we were building, bail before
        // paying for the table build + full scan -- both are guaranteed waste.
        if (Abort.load(std::memory_order_relaxed) != 0) { note = "aborted"; return {}; }

        // =========================================================================
        // Flat open-addressing hash table with linear probing.
        //
        // Reused across calls on this thread via a per-call generation stamp:
        // a slot is "occupied for THIS call" iff stampTbl[idx] == gen. That
        // removes the 8 MB zero-fill that the previous (fresh) table paid on
        // every single Solve -- at ~180 groups that was ~1.4 GB of memset
        // traffic and 180 alloc/free cycles per validation.
        //
        // valTbl holds R-side indices, which are < rsize (<= 2^18), so int32
        // is plenty and halves this buffer's footprint vs int64.
        // =========================================================================
        const size_t tableSize = 1 << 20;
        const size_t tableMask = tableSize - 1;

        static thread_local std::vector<int32_t>  valTbl;
        static thread_local std::vector<uint32_t> stampTbl;
        static thread_local uint32_t gen = 0;
        if (valTbl.size() != tableSize) {
            valTbl.assign(tableSize, -1);
            stampTbl.assign(tableSize, 0);
            gen = 0;
        }
        if (++gen == 0) { // wrapped after 4 billion calls on this thread
            std::fill(stampTbl.begin(), stampTbl.end(), 0u);
            gen = 1;
        }

        for (int64_t i = 0; i < rsize; i++) {
            uint32_t full_hash = (uint32_t)(keys[i] >> 32);
            size_t idx = full_hash & tableMask;
            while (stampTbl[idx] == gen) {
                idx = (idx + 1) & tableMask;
            }
            stampTbl[idx] = gen;
            valTbl[idx] = (int32_t)i;
        }

        const int L0 = 7, L1 = 13;
        int ln = L1 - L0 + 1;
        int64_t lsize = 1;
        for (int k = L0; k <= L1; k++) lsize *= (hiP[k] - loP[k] + 1);

        std::vector<int> actL;
        for (int i = 0; i < ln; i++) {
            if (hiP[L0 + i] > loP[L0 + i]) actL.push_back(i);
        }

        std::vector<Fp6m> lwrapi(ln);
        for (int i = 0; i < ln; i++) {
            BigInteger expVal(hiP[L0 + i] - loP[L0 + i]);
            lwrapi[i] = Fp6mOp::PowBig(gmP[L0 + i], expVal);
        }

        Fp6 tgtI_fp6 = Fp6::One() / target;
        Fp6m tgtI = Fp6m::FromFp6(tgtI_fp6);

        // ---- FIX: snapshot the per-thread buffers for the workers ----------
        // valTbl / stampTbl / gen are thread_local (function-static), and such
        // variables are NOT captured by a lambda: naming them inside 'probe'
        // resolves to the instance of whichever thread RUNS the lambda. On a
        // RunChunks worker thread those instances are empty vectors (and
        // gen == 0), so 'stampTbl[idx]' read through a null pointer -> access
        // violation -> host process dies. It only showed up when nth > 1,
        // i.e. single-group configs / targetGroupId / the raw entry point.
        // Plain locals ARE captured (which is why 'keys', a local reference,
        // already worked), so take the calling thread's pointers here and use
        // only these from inside the lambdas.
        const uint32_t* const stampP = stampTbl.data();
        const int32_t*  const valP   = valTbl.data();
        const uint64_t* const keysP  = keys.data();
        const uint32_t        genL   = gen;

        std::string hitNote = "";
        std::vector<int> hitSol;
        bool found = false;
        std::atomic<int> stopFlag{ 0 }; // was a plain int written by several threads

        auto probe = [&](const Fp6m& needed, int64_t s, std::vector<int>& solutionOut) -> bool {
            ProbeCalls++;
            uint32_t fh = HashM(needed);
            size_t idx = fh & tableMask;

            while (stampP[idx] == genL) {
                int64_t pp = (int64_t)valP[idx];
                if ((uint32_t)(keysP[pp] >> 32) == fh) {
                    ProbeHashHits++;
                    int64_t rs = (int64_t)(uint32_t)keysP[pp];
                    int rd[16];
                    Fp6m rv = ProdM(gmP, loP, hiP, R0, R1, rs, rd);
                    if (rv.Eq(needed)) {
                        int ld[16];
                        ProdM(gm, loP, hiP, L0, L1, s, ld);
                        solutionOut.assign(n, 0);
                        for (int i = 0; i < n; i++) solutionOut[order[i]] = loP[i];
                        for (int i = 0; i < ln; i++) solutionOut[order[L0 + i]] = loP[L0 + i] + ld[i];
                        for (int i = 0; i < rn; i++) solutionOut[order[R0 + i]] = loP[R0 + i] + rd[i];
                        return true;
                    }
                }
                idx = (idx + 1) & tableMask;
            }
            return false;
            };

        int64_t chunkL = (lsize + (std::max)(1, nth) - 1) / (std::max)(1, nth);
        RunChunks(nth, lsize, [&](int64_t begin, int64_t end) {
            if (begin >= end) return;
            int ci = (chunkL > 0) ? (int)(begin / chunkL) : 0; // for the diagnostic note only
            int d0[16];
            Fp6m li0 = ProdM(gimP, loP, hiP, L0, L1, begin, d0);
            int dg[16];
            for (int i = 0; i < ln; i++) dg[i] = d0[i];

            Fp6m la = Fp6mOp::Mul6(tgt, li0);
            Fp6m lb = Fp6mOp::Mul6(tgtI, li0);

            std::vector<int> tempSol;

            for (int64_t s = begin; s < end; s++) {
                if ((s & 0x3FF) == 0 && Abort.load(std::memory_order_relaxed) != 0) return;
                if (stopFlag.load(std::memory_order_relaxed) != 0) return;

                bool matched = false;
                if (probe(la, s, tempSol)) matched = true;
                else if (probe(lb, s, tempSol)) matched = true;

                if (matched) {
                    // CAS: exactly one worker may write hitSol/hitNote/found.
                    int expectedStop = 0;
                    if (stopFlag.compare_exchange_strong(expectedStop, 1)) {
                        hitSol = tempSol;
                        hitNote = "左半 #" + std::to_string(s) + " 命中右半（并行块 " + std::to_string(ci) + "）";
                        found = true;
                    }
                    return;
                }

                for (size_t _j = 0; _j < actL.size(); _j++) {
                    int i = actL[_j];
                    if (dg[i] < hiP[L0 + i] - loP[L0 + i]) {
                        dg[i]++;
                        la = Fp6mOp::Mul6(la, gimP[L0 + i]);
                        lb = Fp6mOp::Mul6(lb, gimP[L0 + i]);
                        break;
                    }
                    dg[i] = 0;
                    la = Fp6mOp::Mul6(la, lwrapi[i]);
                    lb = Fp6mOp::Mul6(lb, lwrapi[i]);
                }
            }
            });

        if (found) {
            note = hitNote;
            return hitSol;
        }
        note = "未找到";
        return {};
    }
};
