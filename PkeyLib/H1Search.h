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
            v = Fp6mOp::Mul6(v, Fp6mOp::PowSmall(g[from + i], (uint32_t)(lo[from + i] + d)));
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

    // =====================================================================
    // One coordinate of an Fp6 product as a 6-term dot product.
    //
    // The MITM only ever needs a HASH of each product, never the product
    // itself (the full value is recomputed on the rare hash match). In
    // Fp6 = Fp3[y]/(y^2+2), Fp3 = Fp[u]/(u^3+u+4), the first coordinate of
    // z = x*y is
    //
    //   z.R.C0 =  xR0*yR0 - 4*(xR1*yR2 + xR2*yR1)
    //           - 2*xI0*yI0 + 8*(xI1*yI2 + xI2*yI1)
    //
    // so with the constants folded into y once (DotForm), one table entry
    // costs 6 lazy multiply-accumulates and ONE Montgomery reduction instead
    // of a whole Mul6 (18 multiplies, 6 reductions, and the Karatsuba
    // bookkeeping around them).
    // =====================================================================
    struct Dot6 { Fpx w[6]; };

    static Dot6 DotForm(const Fp6m& y) {
        const Fpx zero = Fpx::Make(0, 0);
        auto dbl = [](const Fpx& a) { return Fpm::Add(a, a); };
        auto neg = [&](const Fpx& a) { return Fpm::Sub(zero, a); };
        Dot6 d;
        d.w[0] = y.R.C0;
        d.w[1] = neg(dbl(dbl(y.R.C2)));     // -4 * yR2
        d.w[2] = neg(dbl(dbl(y.R.C1)));     // -4 * yR1
        d.w[3] = neg(dbl(y.I.C0));          // -2 * yI0
        d.w[4] = dbl(dbl(dbl(y.I.C2)));     //  8 * yI2
        d.w[5] = dbl(dbl(dbl(y.I.C1)));     //  8 * yI1
        return d;
    }

    // (x * y).R.C0, canonical Montgomery form, for d = DotForm(y).
    // Six products < p^2 summed: < 6 * 2^224, fits the 256-bit accumulator
    // and is below p * 2^128, which is what RedcWide requires.
    FORCE_INLINE static Fpx DotCoord(const Fp6m& x, const Dot6& d) {
        Fp6mOp::W4 acc; Fp6mOp::z4(acc);
        Fp6mOp::addmul2(acc, x.R.C0.A0, x.R.C0.A1, d.w[0].A0, d.w[0].A1);
        Fp6mOp::addmul2(acc, x.R.C1.A0, x.R.C1.A1, d.w[1].A0, d.w[1].A1);
        Fp6mOp::addmul2(acc, x.R.C2.A0, x.R.C2.A1, d.w[2].A0, d.w[2].A1);
        Fp6mOp::addmul2(acc, x.I.C0.A0, x.I.C0.A1, d.w[3].A0, d.w[3].A1);
        Fp6mOp::addmul2(acc, x.I.C1.A0, x.I.C1.A1, d.w[4].A0, d.w[4].A1);
        Fp6mOp::addmul2(acc, x.I.C2.A0, x.I.C2.A1, d.w[5].A0, d.w[5].A1);
        return Fp6mOp::RedcWide(acc);
    }

    FORCE_INLINE static uint32_t HashX(const Fpx& v) {
        uint64_t mix = v.A0 ^ (v.A1 * 0x9E3779B97F4A7C15ULL);
        mix ^= mix >> 33;
        mix *= 0xff51afd7ed558ccdULL;
        mix ^= mix >> 33;
        return (uint32_t)(mix ^ (mix >> 32));
    }

    // Number of index values over positions from..to (mixed radix).
    static int64_t SpanSize(const std::vector<int>& lo, const std::vector<int>& hi, int from, int to) {
        int64_t size = 1;
        for (int k = from; k <= to; k++) size *= (hi[k] - lo[k] + 1);
        return size;
    }

    // How many leading positions of from..to go into the "low" table so that
    // lowSize + highSize (the two tables that get built with real Mul6s) is
    // as small as possible.
    static int BalancedSplit(const std::vector<int>& lo, const std::vector<int>& hi, int from, int to) {
        int best = 0;
        int64_t bestCost = -1;
        for (int k = 0; k <= to - from + 1; k++) {
            int64_t cost = SpanSize(lo, hi, from, from + k - 1) + SpanSize(lo, hi, from + k, to);
            if (bestCost < 0 || cost < bestCost) { bestCost = cost; best = k; }
        }
        return best;
    }

    // out[s] = base * prod_i gen[from+i]^(lo[from+i] + digit_i(s)), where the
    // digits of s are mixed radix with position 'from' the fastest -- the
    // same numbering ProdM uses. Walks the index like an odometer, so each
    // entry costs about one Mul6. 'genInv' are the inverses of 'gen'.
    static void FillProducts(std::vector<Fp6m>& out, const Fp6m& base,
        const std::vector<Fp6m>& gen, const std::vector<Fp6m>& genInv,
        const std::vector<int>& lo, const std::vector<int>& hi, int from, int to) {
        const int cnt = to - from + 1;      // may be 0 -> a single entry
        const int64_t size = SpanSize(lo, hi, from, to);
        out.resize((size_t)size);

        Fp6m prod = base;
        Fp6m wrap[16];
        int span[16], dg[16];
        for (int i = 0; i < cnt; i++) {
            span[i] = hi[from + i] - lo[from + i];
            dg[i] = 0;
            prod = Fp6mOp::Mul6(prod, Fp6mOp::PowSmall(gen[from + i], (uint32_t)lo[from + i]));
            wrap[i] = Fp6mOp::PowSmall(genInv[from + i], (uint32_t)span[i]);
        }
        for (int64_t s = 0;;) {
            out[(size_t)s] = prod;
            if (++s == size) break;
            for (int i = 0; i < cnt; i++) {
                if (span[i] == 0) continue;
                if (dg[i] < span[i]) {
                    dg[i]++;
                    prod = Fp6mOp::Mul6(prod, gen[from + i]);
                    break;
                }
                dg[i] = 0;
                prod = Fp6mOp::Mul6(prod, wrap[i]);
            }
        }
    }

    // Writes lo + digit for each position from..to of index s into
    // solution[order[position]].
    static void DecodeDigits(std::vector<int>& solution, const std::vector<int>& order,
        const std::vector<int>& lo, const std::vector<int>& hi, int from, int to, int64_t s) {
        for (int k = from; k <= to; k++) {
            int radix = hi[k] - lo[k] + 1;
            solution[order[k]] = lo[k] + (int)(s % radix);
            s /= radix;
        }
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

        // ---- generators, their inverses, and both targets -------------------
        // Pairing outputs are unitary (the final exponentiation contains the
        // factor p^3 - 1), and for a unitary element the inverse is just the
        // conjugate (c0, -c1). That replaces the BigInteger Fp6 inversions
        // with sign flips. Each conjugate is checked with one Mul6; if any
        // input is NOT unitary the old inversion path is used, so the result
        // is the same for every input.
        std::vector<Fp6m> gm(n), gim(n);
        const Fp6m oneM = Fp6mOp::One();
        bool unitary = true;
        for (int k = 0; k < n; k++) {
            gm[k] = Fp6m::FromFp6(g[k]);
            gim[k] = Fp6mOp::Conj(gm[k]);
            if (!Fp6mOp::Mul6(gm[k], gim[k]).Eq(oneM)) unitary = false;
        }
        Fp6m tgt = Fp6m::FromFp6(target);
        Fp6m tgtI = Fp6mOp::Conj(tgt);
        if (!Fp6mOp::Mul6(tgt, tgtI).Eq(oneM)) unitary = false;

        if (!unitary) {
            std::vector<Fp6> pref(n + 1);
            pref[0] = Fp6::One();
            for (int k = 0; k < n; k++) pref[k + 1] = pref[k] * g[k];
            Fp6 acc = pref[n].Inverse();
            for (int k = n - 1; k >= 0; k--) {
                gim[k] = Fp6m::FromFp6(acc * pref[k]);
                acc = acc * g[k];
            }
            tgtI = Fp6m::FromFp6(Fp6::One() / target);
        }

        std::vector<int> order = { 0, 1, 2, 3, 4, 5, 13, 6, 7, 8, 9, 10, 11, 12 };
        std::vector<int> loP(n), hiP(n);
        std::vector<Fp6m> gmP(n), gimP(n);
        for (int i = 0; i < n; i++) {
            loP[i] = lo[order[i]];
            hiP[i] = hi[order[i]];
            gmP[i] = gm[order[i]];
            gimP[i] = gim[order[i]];
        }

        // ---- R side: rv(rs) = prod gmP^(lo+digit) over positions 0..6 -------
        // Split into low | high positions:  rv(i + rLow*j) = RLo[i] * RHi[j].
        // Only the two small tables are built with real Mul6s; the 2^18
        // products themselves are never formed, only their first coordinate.
        const int R0 = 0, R1 = 6;
        const int64_t rsize = SpanSize(loP, hiP, R0, R1);
        const int rSplit = BalancedSplit(loP, hiP, R0, R1);
        std::vector<Fp6m> RLo, RHi;
        FillProducts(RLo, oneM, gmP, gimP, loP, hiP, R0, R0 + rSplit - 1);
        FillProducts(RHi, oneM, gmP, gimP, loP, hiP, R0 + rSplit, R1);
        const int64_t rLow = (int64_t)RLo.size(), rHigh = (int64_t)RHi.size();
        std::vector<Dot6> RHiDot((size_t)rHigh);
        for (int64_t j = 0; j < rHigh; j++) RHiDot[(size_t)j] = DotForm(RHi[(size_t)j]);

        // ---- Reused build-side key buffer (per calling thread) --------------
        // Every slot [0,rsize) is overwritten by the build loop below, so no
        // zero-fill is needed -- only a size guarantee. rsize is constant for
        // all 2005 configs, so after the first Solve this never reallocates.
        static thread_local std::vector<uint64_t> keysBuf;
        if ((int64_t)keysBuf.size() < rsize) keysBuf.resize(rsize);
        std::vector<uint64_t>& keys = keysBuf;

        RunChunks(nth, rHigh, [&](int64_t begin, int64_t end) {
            for (int64_t j = begin; j < end; j++) {
                if (Abort.load(std::memory_order_relaxed) != 0) return;
                const Dot6& d = RHiDot[(size_t)j];
                const Fp6m* x = RLo.data();
                uint64_t* out = keys.data() + j * rLow;
                const uint64_t base = (uint64_t)(j * rLow);
                for (int64_t i = 0; i < rLow; i++) {
                    out[i] = ((uint64_t)HashX(DotCoord(x[i], d)) << 32) | (base + (uint64_t)i);
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

        // ---- Pre-filter: one bit per hash, 2^22 bits = 512 KB ----------------
        // The table above is 8 MB, so every probe into it is a cache miss, and
        // almost every probe is a miss in the other sense too (only one of the
        // 2^19 L-side values is in the table). This bitmap stays in L2 and
        // answers "definitely not there" for ~94% of probes (2^18 bits set out
        // of 2^22); only the rest go on to the table.
        const size_t filterBits = (size_t)1 << 22;
        static thread_local std::vector<uint64_t> filterBuf;
        filterBuf.assign(filterBits / 64, 0);

        for (int64_t i = 0; i < rsize; i++) {
            uint32_t full_hash = (uint32_t)(keys[i] >> 32);
            const uint32_t bit = full_hash >> 10;                // top 22 bits
            filterBuf[bit >> 6] |= (uint64_t)1 << (bit & 63);
            size_t idx = full_hash & tableMask;
            while (stampTbl[idx] == gen) {
                idx = (idx + 1) & tableMask;
            }
            stampTbl[idx] = gen;
            valTbl[idx] = (int32_t)i;
        }

        // ---- L side: needed(s) = target * prod gimP^(lo+digit), positions 7..13
        // Same split:  needed(i + lLow*j) = LLoA[i] * LHi[j]   (target folded
        // into the low table), and LLoB[i] * LHi[j] for the inverse target.
        const int L0 = 7, L1 = 13;
        const int lSplit = BalancedSplit(loP, hiP, L0, L1);
        std::vector<Fp6m> LLoA, LLoB, LHi;
        FillProducts(LLoA, tgt, gimP, gmP, loP, hiP, L0, L0 + lSplit - 1);
        FillProducts(LLoB, tgtI, gimP, gmP, loP, hiP, L0, L0 + lSplit - 1);
        FillProducts(LHi, oneM, gimP, gmP, loP, hiP, L0 + lSplit, L1);
        const int64_t lLow = (int64_t)LLoA.size(), lHigh = (int64_t)LHi.size();
        std::vector<Dot6> LHiDot((size_t)lHigh);
        for (int64_t j = 0; j < lHigh; j++) LHiDot[(size_t)j] = DotForm(LHi[(size_t)j]);

        // ---- snapshot the per-thread buffers for the workers ----------------
        // valTbl / stampTbl / gen are thread_local (function-static), and such
        // variables are NOT captured by a lambda: naming them inside 'probe'
        // resolves to the instance of whichever thread RUNS the lambda. Plain
        // locals ARE captured, so take the calling thread's pointers here and
        // use only these from inside the lambdas.
        const uint32_t* const stampP = stampTbl.data();
        const int32_t*  const valP   = valTbl.data();
        const uint64_t* const keysP  = keys.data();
        const uint64_t* const filterP = filterBuf.data();
        const uint32_t        genL   = gen;

        std::string hitNote = "";
        std::vector<int> hitSol;
        bool found = false;
        std::atomic<int> stopFlag{ 0 }; // was a plain int written by several threads

        // fh is the hash of the first coordinate of xLow * LHi[j]. Only on a
        // 32-bit hash match are the two full products formed and compared.
        auto probe = [&](uint32_t fh, const Fp6m& xLow, int64_t j, int64_t s, std::vector<int>& solutionOut) -> bool {
            const uint32_t bit = fh >> 10;
            if (((filterP[bit >> 6] >> (bit & 63)) & 1) == 0) return false;

            size_t idx = fh & tableMask;

            while (stampP[idx] == genL) {
                int64_t pp = (int64_t)valP[idx];
                if ((uint32_t)(keysP[pp] >> 32) == fh) {
                    int64_t rs = (int64_t)(uint32_t)keysP[pp];
                    Fp6m rv = Fp6mOp::Mul6(RLo[(size_t)(rs % rLow)], RHi[(size_t)(rs / rLow)]);
                    Fp6m needed = Fp6mOp::Mul6(xLow, LHi[(size_t)j]);
                    if (rv.Eq(needed)) {
                        solutionOut.assign(n, 0);
                        DecodeDigits(solutionOut, order, loP, hiP, R0, R1, rs);
                        DecodeDigits(solutionOut, order, loP, hiP, L0, L1, s);
                        return true;
                    }
                }
                idx = (idx + 1) & tableMask;
            }
            return false;
            };

        int64_t chunkL = (lHigh + (std::max)(1, nth) - 1) / (std::max)(1, nth);
        RunChunks(nth, lHigh, [&](int64_t begin, int64_t end) {
            if (begin >= end) return;
            int ci = (chunkL > 0) ? (int)(begin / chunkL) : 0; // for the diagnostic note only
            std::vector<int> tempSol;

            for (int64_t j = begin; j < end; j++) {
                if (Abort.load(std::memory_order_relaxed) != 0) return;
                if (stopFlag.load(std::memory_order_relaxed) != 0) return;
                const Dot6& d = LHiDot[(size_t)j];

                for (int64_t i = 0; i < lLow; i++) {
                    const int64_t s = i + lLow * j;
                    bool matched = false;
                    if (probe(HashX(DotCoord(LLoA[(size_t)i], d)), LLoA[(size_t)i], j, s, tempSol)) matched = true;
                    else if (probe(HashX(DotCoord(LLoB[(size_t)i], d)), LLoB[(size_t)i], j, s, tempSol)) matched = true;

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
