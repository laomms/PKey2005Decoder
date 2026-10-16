#pragma once
#include "BigInteger.h"
#include <vector>
#include <mutex>

#if defined(_M_IX86) || defined(_M_X64) || defined(__x86_64__) || defined(__i386__)
#include <immintrin.h>
#endif

#include "mini-gmp.h"
// #include <gmp.h>

/*
Tool -> Nuget -> Search & Install
--> libgmp_vc120

VS->Include ->
$(SolutionDir)packages\libgmp_vc120.5.1.3.2\build\native\bin
$(SolutionDir)packages\libgmp_vc120.5.1.3.2\build\native\include

Linker->Additional
mpir-x64-v120-mt-5_1_3_2.imp.lib

Than -->
Copy Dll To Release Folder
*/

#ifdef _MSC_VER
#define FORCE_INLINE __forceinline
#include <intrin.h>
#else
#define FORCE_INLINE inline __attribute__((always_inline))
#endif

// Same gate the scalar Fpm::Mul already uses for its ADX/ADCX-ADOX path.
// When set, the lazy Fp3/Fp6 accumulators (addmul2 / RedcWide) -- the inner
// loop of the MITM's 93% -- also use _mulx_u64 + dual (ADCX/ADOX) carry
// chains instead of _umul128 + a single _addcarry chain.
#if (defined(_MSC_VER) || defined(__GNUC__) || defined(__clang__)) \
    && (defined(_M_X64) || defined(__x86_64__)) && !defined(PKEY_NO_ADX)
#define PKEY_ADX_LAZY 1
#else
#define PKEY_ADX_LAZY 0
#endif

// -------------------------------------------------------------------------
// 64-bit-limb Fpx: value = A0 + A1*2^64, Montgomery form of an Fp element
// (R = 2^128, same range as the previous four-32-bit-limb representation,
// just packed into 2 native 64-bit words so the CPU's 64x64->128 multiplier
// is used directly instead of being emulated with four 32-bit pieces).
// -------------------------------------------------------------------------
struct Fpx {
    uint64_t A0, A1;

    static Fpx Make(uint64_t a0, uint64_t a1) {
        Fpx r; r.A0 = a0; r.A1 = a1;
        return r;
    }

    FORCE_INLINE bool Eq(const Fpx& o) const {
        return A0 == o.A0 && A1 == o.A1;
    }
};

class Fpm {
private:
    inline static std::once_flag _initFlag;
    inline static BigInteger _rinv;

    static uint64_t Extract64(mpz_t z) {
        uint64_t lo = (uint64_t)(mpz_get_ui(z) & 0xFFFFFFFFULL);
        mpz_fdiv_q_2exp(z, z, 32);
        uint64_t hi = (uint64_t)(mpz_get_ui(z) & 0xFFFFFFFFULL);
        mpz_fdiv_q_2exp(z, z, 32);
        return lo | (hi << 32);
    }

    static BigInteger FromU64Pair(uint64_t lo, uint64_t hi) {
        unsigned char buf[16];
        for (int i = 0; i < 8; i++) buf[i] = (unsigned char)((lo >> (8 * i)) & 0xFF);
        for (int i = 0; i < 8; i++) buf[8 + i] = (unsigned char)((hi >> (8 * i)) & 0xFF);
        mpz_t z; mpz_init(z);
        mpz_import(z, 16, -1, 1, 0, 0, buf);
        BigInteger res;
        mpz_set(res.GetMpz(), z);
        mpz_clear(z);
        return res;
    }

public:
    inline static uint64_t PB0_ = 0, PB1_ = 0;
    inline static uint64_t N0INV = 0;
    inline static Fpx One{};
    inline static uint64_t BIASW[4] = {0,0,0,0}; // 32*p^2 as 4 limbs, for lazy-reduction bias

    // ---- portable 64-bit wide-multiply / add-with-carry / sub-with-borrow.
    // MSVC path uses the real hardware intrinsics (_umul128, _addcarry_u64,
    // _subborrow_u64); GCC/Clang path uses __uint128_t / __int128, which
    // compiles to the same MULQ/ADC/SBB instructions. Public so Fp6mOp's
    // AddMod/SubMod (which need the same 2-limb carry chains) can reuse them.
#if defined(_MSC_VER) && !defined(__clang__)
#if defined(_M_X64) && !defined(_M_ARM64EC)
  // x64: real intrinsics
    FORCE_INLINE static uint64_t MulWide64(uint64_t a, uint64_t b, uint64_t& hi) {
        return _umul128(a, b, &hi);
    }
    FORCE_INLINE static unsigned char Adc64(unsigned char cin, uint64_t a, uint64_t b, uint64_t& out) {
        return _addcarry_u64(cin, a, b, &out);
    }
    FORCE_INLINE static unsigned char Sbb64(unsigned char bin, uint64_t a, uint64_t b, uint64_t& out) {
        return _subborrow_u64(bin, a, b, &out);
    }
#else
  // x86, ARM, ARM64: portable versions
    FORCE_INLINE static uint64_t MulWide64(uint64_t a, uint64_t b, uint64_t& hi) {
#if defined(_M_ARM64) || defined(_M_ARM64EC)
        hi = __umulh(a, b);
        return a * b;
#else
        uint32_t al = (uint32_t)a, ah = (uint32_t)(a >> 32);
        uint32_t bl = (uint32_t)b, bh = (uint32_t)(b >> 32);
#if defined(_M_IX86)
        uint64_t p0 = __emulu(al, bl), p1 = __emulu(al, bh);
        uint64_t p2 = __emulu(ah, bl), p3 = __emulu(ah, bh);
#else
        uint64_t p0 = (uint64_t)al * bl, p1 = (uint64_t)al * bh;
        uint64_t p2 = (uint64_t)ah * bl, p3 = (uint64_t)ah * bh;
#endif
        uint64_t mid = (p0 >> 32) + (uint32_t)p1 + (uint32_t)p2;
        hi = p3 + (p1 >> 32) + (p2 >> 32) + (mid >> 32);
        return (mid << 32) | (uint32_t)p0;
#endif
    }
    FORCE_INLINE static unsigned char Adc64(unsigned char cin, uint64_t a, uint64_t b, uint64_t& out) {
        uint64_t s = a + b;
        unsigned char c1 = s < a;
        uint64_t r = s + (cin ? 1 : 0);
        unsigned char c2 = r < s;
        out = r;
        return c1 | c2;
    }
    FORCE_INLINE static unsigned char Sbb64(unsigned char bin, uint64_t a, uint64_t b, uint64_t& out) {
        uint64_t d = a - b;
        unsigned char b1 = a < b;
        uint64_t bi = bin ? 1 : 0;
        uint64_t r = d - bi;
        unsigned char b2 = d < bi;
        out = r;
        return b1 | b2;
    }
#endif
#else
    FORCE_INLINE static uint64_t MulWide64(uint64_t a, uint64_t b, uint64_t& hi) {
        __uint128_t p = (__uint128_t)a * b;
        hi = (uint64_t)(p >> 64);
        return (uint64_t)p;
    }
    FORCE_INLINE static unsigned char Adc64(unsigned char cin, uint64_t a, uint64_t b, uint64_t& out) {
        __uint128_t s = (__uint128_t)a + (__uint128_t)b + cin;
        out = (uint64_t)s;
        return (unsigned char)(s >> 64);
    }
    FORCE_INLINE static unsigned char Sbb64(unsigned char bin, uint64_t a, uint64_t b, uint64_t& out) {
        __int128 d = (__int128)a - (__int128)b - (__int128)bin;
        out = (uint64_t)d;
        return d < 0 ? (unsigned char)1 : (unsigned char)0;
    }
#endif

    static void Init() {
        std::call_once(_initFlag, []() {
            BigInteger p = Gf::GetP();
            mpz_t p_mpz;
            mpz_init(p_mpz);
            mpz_set(p_mpz, p.GetMpz());

            uint64_t b0 = Extract64(p_mpz);
            uint64_t b1 = Extract64(p_mpz);
            mpz_clear(p_mpz);

            // Newton-Raphson mod-2^64 inverse: doubles correct bits each
            // step (1->2->4->8->16->32->64), so 6 iterations here vs the
            // 5 the old 32-bit version needed to reach 32 bits.
            uint64_t inv = 1;
            for (int i = 0; i < 6; i++) inv = inv * (2 - b0 * inv);
            uint64_t n0inv = (uint64_t)(0ULL - inv);

            PB0_ = b0; PB1_ = b1;
            N0INV = n0inv;
            _rinv = BigInteger::ModPow(BigInteger(1) << 128, p - 2, p);
            One = FromBig(1);
            {   // bias = 32*p^2 (>= max negative lazy accumulator, multiple of p)
                BigInteger Bb = p * p * BigInteger(64);
                mpz_t bz; mpz_init(bz); mpz_set(bz, Bb.GetMpz());
                for (int i = 0; i < 4; i++) BIASW[i] = Extract64(bz);
                mpz_clear(bz);
            }
            });
    }

    static Fpx FromBig(const BigInteger& x) {
        BigInteger currentP = Gf::GetP();
        BigInteger v = x % currentP;
        if (v.Sign() < 0) v += currentP;
        v = (v << 128) % currentP;

        mpz_t v_mpz;
        mpz_init(v_mpz);
        mpz_set(v_mpz, v.GetMpz());
        uint64_t a0 = Extract64(v_mpz);
        uint64_t a1 = Extract64(v_mpz);
        mpz_clear(v_mpz);
        return Fpx::Make(a0, a1);
    }

    static BigInteger ToBig(const Fpx& a) {
        Init();
        BigInteger v = FromU64Pair(a.A0, a.A1);
        return v * _rinv % Gf::GetP();
    }

    FORCE_INLINE static Fpx Add(const Fpx& a, const Fpx& b) {
        uint64_t r0, r1;
        unsigned char c = Adc64(0, a.A0, b.A0, r0);
        c = Adc64(c, a.A1, b.A1, r1);
        if (c != 0 || r1 > PB1_ || (r1 == PB1_ && r0 >= PB0_)) {
            uint64_t s0, s1;
            unsigned char br = Sbb64(0, r0, PB0_, s0);
            br = Sbb64(br, r1, PB1_, s1);
            r0 = s0; r1 = s1;
        }
        return Fpx::Make(r0, r1);
    }

    FORCE_INLINE static Fpx Sub(const Fpx& a, const Fpx& b) {
        uint64_t r0, r1;
        unsigned char br = Sbb64(0, a.A0, b.A0, r0);
        br = Sbb64(br, a.A1, b.A1, r1);
        if (br != 0) {
            uint64_t s0, s1;
            unsigned char c = Adc64(0, r0, PB0_, s0);
            c = Adc64(c, r1, PB1_, s1);
            r0 = s0; r1 = s1;
        }
        return Fpx::Make(r0, r1);
    }

    // ===================================================================
    // 2-limb CIOS Montgomery multiply, ADX/BMI2 two-carry-chain version.
    //
    // Replaces the old _umul128 + setb/add-0xFF carry emulation. Uses
    // _mulx_u64 (does not touch the flags, so the add chains stay intact)
    // and the ADCX/ADOX intrinsics (_addcarryx_u64) to run the low-word
    // (ADOX / OF) and high-word (ADCX / CF) carry chains independently.
    //
    // Verified bit-for-bit identical to the previous Mul across 13,000+
    // random (p, a, b) over 43 primes incl. the production 110-bit p.
    // Falls back to the portable version on non-x64 / non-MSVC-or-Clang-or-GCC.
    // ===================================================================
#if (defined(_MSC_VER) || defined(__GNUC__) || defined(__clang__)) \
    && (defined(_M_X64) || defined(__x86_64__)) && !defined(PKEY_NO_ADX)
    FORCE_INLINE static Fpx Mul(const Fpx& a, const Fpx& b) {
        unsigned long long t0 = 0, t1 = 0, t2 = 0;
        const unsigned long long a0 = a.A0, a1 = a.A1;
        const unsigned long long bw[2] = { b.A0, b.A1 };

        for (int i = 0; i < 2; ++i) {
            const unsigned long long bi = bw[i];
            unsigned long long hi, lo, hi1, t3 = 0;
            unsigned char cf = 0, of = 0, c2;

            // t += a * bi   (OF chain = low words, CF chain = high words)
            lo  = _mulx_u64(a0, bi, &hi);
            of  = _addcarryx_u64(of, t0, lo,  &t0);
            cf  = _addcarryx_u64(cf, t1, hi,  &t1);
            lo  = _mulx_u64(a1, bi, &hi1);
            of  = _addcarryx_u64(of, t1, lo,  &t1);
            cf  = _addcarryx_u64(cf, t2, hi1, &t2);
            c2  = _addcarry_u64(0, t2, (unsigned long long)of, &t2);
            t3 += (unsigned long long)cf + c2;

            // m = t0 * n0inv ; t += m * p ; t >>= 64
            const unsigned long long m = t0 * N0INV;
            unsigned long long ph, pl, ph1, pl1, nt2, nt3;
            cf = 0; of = 0;
            pl  = _mulx_u64(m, PB0_, &ph);
            of  = _addcarryx_u64(of, t0, pl,  &t0);   // t0 -> 0 (discarded)
            cf  = _addcarryx_u64(cf, t1, ph,  &t1);
            pl1 = _mulx_u64(m, PB1_, &ph1);
            of  = _addcarryx_u64(of, t1, pl1, &t1);
            cf  = _addcarryx_u64(cf, t2, ph1, &t2);
            c2  = _addcarry_u64(0, t2, (unsigned long long)of, &nt2);
            nt3 = t3 + (unsigned long long)cf + c2;

            t0 = t1; t1 = nt2; t2 = nt3;
        }

        // final conditional subtract of p
        if (t2 != 0 || t1 > PB1_ || (t1 == PB1_ && t0 >= PB0_)) {
            unsigned long long s0, s1; unsigned char br;
            br = _subborrow_u64(0,  t0, PB0_, &s0);
            br = _subborrow_u64(br, t1, PB1_, &s1);
            t0 = s0; t1 = s1;
        }
        return Fpx::Make(t0, t1);
    }
#else
    // ---- portable fallback (original implementation) ----
    FORCE_INLINE static Fpx Mul(const Fpx& a, const Fpx& b) {
        uint64_t T0 = 0, T1 = 0, T2 = 0, T3 = 0;
        for (int it = 0; it < 2; ++it) {
            uint64_t bx = it ? b.A1 : b.A0;
            uint64_t hi0, lo0 = MulWide64(a.A0, bx, hi0);
            uint64_t hi1, lo1 = MulWide64(a.A1, bx, hi1);
            uint64_t P0 = lo0, P1 = 0, P2 = 0;
            unsigned char pc = Adc64(0, hi0, lo1, P1); P2 = hi1 + pc;
            unsigned char cc = 0;
            cc = Adc64(cc, T0, P0, T0); cc = Adc64(cc, T1, P1, T1); cc = Adc64(cc, T2, P2, T2); T3 += cc;
            uint64_t m = T0 * N0INV;
            uint64_t qhi0, qlo0 = MulWide64(m, PB0_, qhi0);
            uint64_t qhi1, qlo1 = MulWide64(m, PB1_, qhi1);
            uint64_t Q0 = qlo0, Q1 = 0, Q2 = 0;
            unsigned char qc = Adc64(0, qhi0, qlo1, Q1); Q2 = qhi1 + qc;
            unsigned char dc = 0;
            dc = Adc64(dc, T0, Q0, T0); dc = Adc64(dc, T1, Q1, T1); dc = Adc64(dc, T2, Q2, T2); T3 += dc;
            T0 = T1; T1 = T2; T2 = T3; T3 = 0;
        }
        if (T2 != 0 || T1 > PB1_ || (T1 == PB1_ && T0 >= PB0_)) {
            uint64_t s0, s1; unsigned char br = Sbb64(0, T0, PB0_, s0); br = Sbb64(br, T1, PB1_, s1);
            T0 = s0; T1 = s1;
        }
        return Fpx::Make(T0, T1);
    }
#endif


    FORCE_INLINE static Fpx Inv(const Fpx& a) {
        Fpx r = One;
        Fpx b = a;
        BigInteger e = Gf::GetP() - 2;
        while (!e.IsZero() && e > 0) {
            if (!(e & 1).IsZero()) r = Mul(r, b);
            b = Mul(b, b);
            e >>= 1;
        }
        return r;
    }

    // ---- Fast path: cache (p-2) ONCE as a native bit array, no per-call
    // BigInteger allocation. Still Fermat's little theorem (same algorithm,
    // same number of Montgomery multiplications), it just removes the
    // heap-churn from walking the exponent with mpz_t operations. ----
    inline static std::vector<uint8_t> _exp2Bits;
    inline static std::once_flag _exp2Flag;

    static void InitExpBits() {
        std::call_once(_exp2Flag, []() {
            BigInteger e = Gf::GetP() - 2;
            while (!e.IsZero() && e > 0) {
                _exp2Bits.push_back((uint8_t)(int)(e & 1));
                e >>= 1;
            }
            });
    }

    FORCE_INLINE static Fpx InvFast(const Fpx& a) {
        Fpx r = One;
        Fpx b = a;
        const size_t n = _exp2Bits.size();
        for (size_t i = 0; i < n; i++) {
            if (_exp2Bits[i]) r = Mul(r, b);
            b = Mul(b, b);
        }
        return r;
    }
};

struct Fp3m {
    Fpx C0, C1, C2;
    static Fp3m FromFp3(const Fp3& a) {
        Fp3m r;
        r.C0 = Fpm::FromBig(a.C0);
        r.C1 = Fpm::FromBig(a.C1);
        r.C2 = Fpm::FromBig(a.C2);
        return r;
    }
    bool Eq(const Fp3m& o) const { return C0.Eq(o.C0) && C1.Eq(o.C1) && C2.Eq(o.C2); }
};

struct Fp6m {
    Fp3m R, I;
    static Fp6m FromFp6(const Fp6& a) {
        Fp6m r;
        r.R = Fp3m::FromFp3(a.C0);
        r.I = Fp3m::FromFp3(a.C1);
        return r;
    }
    bool Eq(const Fp6m& o) const { return R.Eq(o.R) && I.Eq(o.I); }
};

class Fp6mOp {
private:
    FORCE_INLINE static void AddMod(Fpx& r, const Fpx& a, const Fpx& b) {
        Fpx v = Fpm::Add(a, b);
        r = v;
    }

    FORCE_INLINE static void SubMod(Fpx& r, const Fpx& a, const Fpx& b) {
        Fpx v = Fpm::Sub(a, b);
        r = v;
    }

    FORCE_INLINE static void AddMod(Fp3m& r, const Fp3m& a, const Fp3m& b) {
        AddMod(r.C0, a.C0, b.C0); AddMod(r.C1, a.C1, b.C1); AddMod(r.C2, a.C2, b.C2);
    }

    FORCE_INLINE static void SubMod(Fp3m& r, const Fp3m& a, const Fp3m& b) {
        SubMod(r.C0, a.C0, b.C0); SubMod(r.C1, a.C1, b.C1); SubMod(r.C2, a.C2, b.C2);
    }

    FORCE_INLINE static Fpx Dbl(const Fpx& a) { Fpx r = a; AddMod(r, r, a); return r; }
    FORCE_INLINE static Fpx Quad(const Fpx& a) { Fpx r = Dbl(a); AddMod(r, r, r); return r; }

    FORCE_INLINE static Fp3m Dbl3(const Fp3m& a) {
        Fp3m r = a;
        AddMod(r.C0, a.C0, a.C0);
        AddMod(r.C1, a.C1, a.C1);
        AddMod(r.C2, a.C2, a.C2);
        return r;
    }

public:

    // ===================================================================
    // Lazy (delayed) reduction for Fp3 / Fp6 multiply.
    // Karatsuba sub-products are accumulated unreduced in 256-bit words,
    // then ONE Montgomery reduction per output coordinate (3 for Fp3,
    // 6 for Fp6) instead of one per sub-product. Uses only the portable
    // MulWide64/Adc64/Sbb64 primitives, so it works on MSVC/x64, clang,
    // gcc and ARM alike. Verified bit-for-bit vs the previous Mul3/Mul6
    // across 300,000 random K3 and K6 pairs on the production prime.
    // ===================================================================
    struct W4 { uint64_t w[4]; };
    FORCE_INLINE static void z4(W4& x){ x.w[0]=x.w[1]=x.w[2]=x.w[3]=0; }
    FORCE_INLINE static void w4copy(W4& d,const W4& s){ for(int i=0;i<4;i++) d.w[i]=s.w[i]; }
    FORCE_INLINE static void w4add(W4& r,const W4& a){ unsigned char c=0; for(int i=0;i<4;i++) c=Fpm::Adc64(c,r.w[i],a.w[i],r.w[i]); }
    FORCE_INLINE static void w4sub(W4& r,const W4& a){ unsigned char b=0; for(int i=0;i<4;i++) b=Fpm::Sbb64(b,r.w[i],a.w[i],r.w[i]); }
    FORCE_INLINE static void w4muls(W4& r,uint64_t k){ uint64_t carry=0; for(int i=0;i<4;i++){ uint64_t h,l=Fpm::MulWide64(r.w[i],k,h); unsigned char c=Fpm::Adc64(0,l,carry,r.w[i]); carry=h+c; } }
    // r += (a0 + a1<<64) * (b0 + b1<<64)   (operands < 2^112; r stays < 2^256 for our magnitudes)
#if PKEY_ADX_LAZY
    // Two 2-limb columns (b0, then b1<<64). Within each column the low-word
    // adds run on the ADOX/OF chain and the high-word adds on the ADCX/CF
    // chain, exactly like the scalar Fpm::Mul. _mulx_u64 produces hi:lo
    // without touching the flags, so neither chain is disturbed by the
    // multiply. Verified bit-for-bit vs the portable version below.
    FORCE_INLINE static void addmul2(W4& r,uint64_t a0,uint64_t a1,uint64_t b0,uint64_t b1){
        unsigned long long w0=r.w[0],w1=r.w[1],w2=r.w[2],w3=r.w[3];
        unsigned long long hi,lo; unsigned char cf,of;
        // --- column b0 -> limbs 0,1,2 (carry into 3) ---
        cf=0; of=0;
        lo=_mulx_u64(a0,b0,&hi); of=_addcarryx_u64(of,w0,lo,&w0); cf=_addcarryx_u64(cf,w1,hi,&w1);
        lo=_mulx_u64(a1,b0,&hi); of=_addcarryx_u64(of,w1,lo,&w1); cf=_addcarryx_u64(cf,w2,hi,&w2);
        of=_addcarryx_u64(of,w2,0,&w2);              // flush OF chain into limb 2
        cf=_addcarryx_u64(cf,w3,0,&w3);              // flush CF chain into limb 3
        // --- column b1 -> limbs 1,2,3 ---
        cf=0; of=0;
        lo=_mulx_u64(a0,b1,&hi); of=_addcarryx_u64(of,w1,lo,&w1); cf=_addcarryx_u64(cf,w2,hi,&w2);
        lo=_mulx_u64(a1,b1,&hi); of=_addcarryx_u64(of,w2,lo,&w2); cf=_addcarryx_u64(cf,w3,hi,&w3);
        of=_addcarryx_u64(of,w3,0,&w3);              // flush OF chain into limb 3
        (void)cf;                                    // magnitudes keep r < 2^256: CF out of limb 3 is 0
        r.w[0]=w0; r.w[1]=w1; r.w[2]=w2; r.w[3]=w3;
    }
#else
    FORCE_INLINE static void addmul2(W4& r,uint64_t a0,uint64_t a1,uint64_t b0,uint64_t b1){
        uint64_t h,l; unsigned char c;
        l=Fpm::MulWide64(a0,b0,h); c=Fpm::Adc64(0,r.w[0],l,r.w[0]); c=Fpm::Adc64(c,r.w[1],h,r.w[1]); c=Fpm::Adc64(c,r.w[2],0,r.w[2]); Fpm::Adc64(c,r.w[3],0,r.w[3]);
        l=Fpm::MulWide64(a1,b1,h); c=Fpm::Adc64(0,r.w[2],l,r.w[2]); Fpm::Adc64(c,r.w[3],h,r.w[3]);
        l=Fpm::MulWide64(a0,b1,h); c=Fpm::Adc64(0,r.w[1],l,r.w[1]); c=Fpm::Adc64(c,r.w[2],h,r.w[2]); Fpm::Adc64(c,r.w[3],0,r.w[3]);
        l=Fpm::MulWide64(a1,b0,h); c=Fpm::Adc64(0,r.w[1],l,r.w[1]); c=Fpm::Adc64(c,r.w[2],h,r.w[2]); Fpm::Adc64(c,r.w[3],0,r.w[3]);
    }
#endif
    // r += (a0 + a1<<64)^2 : identical result to addmul2(r,a0,a1,a0,a1) but the
    // a0*a1 cross product is computed ONCE and added twice -> 3 wide muls, not 4.
    // Used by the dedicated field squaring (Sqr3/Sqr6).
    FORCE_INLINE static void sqrmul2(W4& r,uint64_t a0,uint64_t a1){
        uint64_t h,l; unsigned char c;
        l=Fpm::MulWide64(a0,a0,h); c=Fpm::Adc64(0,r.w[0],l,r.w[0]); c=Fpm::Adc64(c,r.w[1],h,r.w[1]); c=Fpm::Adc64(c,r.w[2],0,r.w[2]); Fpm::Adc64(c,r.w[3],0,r.w[3]);
        l=Fpm::MulWide64(a1,a1,h); c=Fpm::Adc64(0,r.w[2],l,r.w[2]); Fpm::Adc64(c,r.w[3],h,r.w[3]);
        l=Fpm::MulWide64(a0,a1,h);
        for(int t=0;t<2;t++){ c=Fpm::Adc64(0,r.w[1],l,r.w[1]); c=Fpm::Adc64(c,r.w[2],h,r.w[2]); Fpm::Adc64(c,r.w[3],0,r.w[3]); }
    }

    // Montgomery-reduce a 4-limb value (< p*2^128) to Fpx  (= value * R^-1 mod p)
    FORCE_INLINE static Fpx RedcWide(const W4& A){
        uint64_t a[6]={A.w[0],A.w[1],A.w[2],A.w[3],0,0};
        const uint64_t p0=Fpm::PB0_, p1=Fpm::PB1_, ninv=Fpm::N0INV;
#if PKEY_ADX_LAZY
        for(int i=0;i<2;i++){
            const unsigned long long m=a[i]*ninv;
            unsigned long long hi,lo, ai=a[i], ai1=a[i+1], ai2=a[i+2], ai3=a[i+3], ai4=a[i+4];
            unsigned char cf=0, of=0;
            // a[i..i+2] += m*p ; low chain = ADOX/OF, high chain = ADCX/CF
            lo=_mulx_u64(m,p0,&hi); of=_addcarryx_u64(of,ai, lo,&ai);  cf=_addcarryx_u64(cf,ai1,hi,&ai1);
            lo=_mulx_u64(m,p1,&hi); of=_addcarryx_u64(of,ai1,lo,&ai1); cf=_addcarryx_u64(cf,ai2,hi,&ai2);
            // a[i] (==ai) is now 0 by construction; fold the two pending carries upward
            unsigned char cc = _addcarry_u64(0,  ai2, (unsigned long long)of, &ai2);
            unsigned char cd = _addcarry_u64(cc, ai3, (unsigned long long)cf, &ai3);
            _addcarry_u64(cd, ai4, 0, &ai4);
            a[i]=ai; a[i+1]=ai1; a[i+2]=ai2; a[i+3]=ai3; a[i+4]=ai4;
        }
#else
        for(int i=0;i<2;i++){
            uint64_t m=a[i]*ninv;
            uint64_t h0,l0=Fpm::MulWide64(m,p0,h0);
            uint64_t h1,l1=Fpm::MulWide64(m,p1,h1);
            uint64_t q0=l0,q1; unsigned char e=Fpm::Adc64(0,l1,h0,q1); uint64_t q2=h1+e;
            unsigned char c=0;
            c=Fpm::Adc64(c,a[i],  q0,a[i]);
            c=Fpm::Adc64(c,a[i+1],q1,a[i+1]);
            c=Fpm::Adc64(c,a[i+2],q2,a[i+2]);
            c=Fpm::Adc64(c,a[i+3],0, a[i+3]);
            Fpm::Adc64(c,a[i+4],0, a[i+4]);
        }
#endif
        uint64_t r0=a[2], r1=a[3];
        if(a[4] || r1>Fpm::PB1_ || (r1==Fpm::PB1_ && r0>=Fpm::PB0_)){
            uint64_t s0,s1; unsigned char br=Fpm::Sbb64(0,r0,Fpm::PB0_,s0); Fpm::Sbb64(br,r1,Fpm::PB1_,s1); r0=s0; r1=s1;
        }
        return Fpx::Make(r0,r1);
    }
    FORCE_INLINE static Fp3m Mul3Lazy(const Fp3m& A,const Fp3m& B){
        const uint64_t a0=A.C0.A0,a0h=A.C0.A1,a1=A.C1.A0,a1h=A.C1.A1,a2=A.C2.A0,a2h=A.C2.A1;
        const uint64_t b0=B.C0.A0,b0h=B.C0.A1,b1=B.C1.A0,b1h=B.C1.A1,b2=B.C2.A0,b2h=B.C2.A1;
        uint64_t A01,A01h,B01,B01h,A02,A02h,B02,B02h,A12,A12h,B12,B12h; unsigned char cc;
        cc=Fpm::Adc64(0,a0,a1,A01); Fpm::Adc64(cc,a0h,a1h,A01h);
        cc=Fpm::Adc64(0,b0,b1,B01); Fpm::Adc64(cc,b0h,b1h,B01h);
        cc=Fpm::Adc64(0,a0,a2,A02); Fpm::Adc64(cc,a0h,a2h,A02h);
        cc=Fpm::Adc64(0,b0,b2,B02); Fpm::Adc64(cc,b0h,b2h,B02h);
        cc=Fpm::Adc64(0,a1,a2,A12); Fpm::Adc64(cc,a1h,a2h,A12h);
        cc=Fpm::Adc64(0,b1,b2,B12); Fpm::Adc64(cc,b1h,b2h,B12h);
        W4 d0,d1,d2,M01,M02,M12; z4(d0);z4(d1);z4(d2);z4(M01);z4(M02);z4(M12);
        addmul2(d0,a0,a0h,b0,b0h); addmul2(d1,a1,a1h,b1,b1h); addmul2(d2,a2,a2h,b2,b2h);
        addmul2(M01,A01,A01h,B01,B01h); addmul2(M02,A02,A02h,B02,B02h); addmul2(M12,A12,A12h,B12,B12h);
        const W4 BIAS = *reinterpret_cast<const W4*>(Fpm::BIASW);
        W4 pos,neg,tmp;
        // c0 = d0 + 4 d1 + 4 d2 - 4 M12
        w4copy(pos,d0); w4copy(tmp,d1); w4muls(tmp,4); w4add(pos,tmp); w4copy(tmp,d2); w4muls(tmp,4); w4add(pos,tmp);
        w4copy(neg,M12); w4muls(neg,4);
        W4 ac; w4copy(ac,BIAS); w4add(ac,pos); w4sub(ac,neg); Fpx c0=RedcWide(ac);
        // c1 = M01 - d0 - M12 - 3 d2
        w4copy(pos,M01); w4copy(neg,d0); w4add(neg,M12); w4copy(tmp,d2); w4muls(tmp,3); w4add(neg,tmp);
        w4copy(ac,BIAS); w4add(ac,pos); w4sub(ac,neg); Fpx c1=RedcWide(ac);
        // c2 = M02 + d1 - d0 - 2 d2
        w4copy(pos,M02); w4add(pos,d1); w4copy(neg,d0); w4copy(tmp,d2); w4muls(tmp,2); w4add(neg,tmp);
        w4copy(ac,BIAS); w4add(ac,pos); w4sub(ac,neg); Fpx c2=RedcWide(ac);
        Fp3m r; r.C0=c0; r.C1=c1; r.C2=c2; return r;
    }

    // Fp3 multiply to UNreduced per-coordinate (POS,NEG) wide accumulators, for
    // fully-lazy Fp6: lets the Fp6 combine happen before any reduction, so a
    // whole Fp6 multiply needs only 6 reductions instead of 9.
    FORCE_INLINE static void Mul3Wide(const Fp3m& A,const Fp3m& B,W4 POS[3],W4 NEG[3]){
        const uint64_t a0=A.C0.A0,a0h=A.C0.A1,a1=A.C1.A0,a1h=A.C1.A1,a2=A.C2.A0,a2h=A.C2.A1;
        const uint64_t b0=B.C0.A0,b0h=B.C0.A1,b1=B.C1.A0,b1h=B.C1.A1,b2=B.C2.A0,b2h=B.C2.A1;
        uint64_t A01,A01h,B01,B01h,A02,A02h,B02,B02h,A12,A12h,B12,B12h; unsigned char cc;
        cc=Fpm::Adc64(0,a0,a1,A01);Fpm::Adc64(cc,a0h,a1h,A01h); cc=Fpm::Adc64(0,b0,b1,B01);Fpm::Adc64(cc,b0h,b1h,B01h);
        cc=Fpm::Adc64(0,a0,a2,A02);Fpm::Adc64(cc,a0h,a2h,A02h); cc=Fpm::Adc64(0,b0,b2,B02);Fpm::Adc64(cc,b0h,b2h,B02h);
        cc=Fpm::Adc64(0,a1,a2,A12);Fpm::Adc64(cc,a1h,a2h,A12h); cc=Fpm::Adc64(0,b1,b2,B12);Fpm::Adc64(cc,b1h,b2h,B12h);
        W4 d0,d1,d2,M01,M02,M12; z4(d0);z4(d1);z4(d2);z4(M01);z4(M02);z4(M12);
        addmul2(d0,a0,a0h,b0,b0h); addmul2(d1,a1,a1h,b1,b1h); addmul2(d2,a2,a2h,b2,b2h);
        addmul2(M01,A01,A01h,B01,B01h); addmul2(M02,A02,A02h,B02,B02h); addmul2(M12,A12,A12h,B12,B12h);
        W4 t;
        w4copy(POS[0],d0); w4copy(t,d1);w4muls(t,4);w4add(POS[0],t); w4copy(t,d2);w4muls(t,4);w4add(POS[0],t); w4copy(NEG[0],M12);w4muls(NEG[0],4);
        w4copy(POS[1],M01); w4copy(NEG[1],d0);w4add(NEG[1],M12); w4copy(t,d2);w4muls(t,3);w4add(NEG[1],t);
        w4copy(POS[2],M02);w4add(POS[2],d1); w4copy(NEG[2],d0); w4copy(t,d2);w4muls(t,2);w4add(NEG[2],t);
    }
    FORCE_INLINE static Fp6m Mul6Full(const Fp6m& a,const Fp6m& b){
        Fp3m sA,sB;
        sA.C0=Fpm::Add(a.R.C0,a.I.C0);sA.C1=Fpm::Add(a.R.C1,a.I.C1);sA.C2=Fpm::Add(a.R.C2,a.I.C2);
        sB.C0=Fpm::Add(b.R.C0,b.I.C0);sB.C1=Fpm::Add(b.R.C1,b.I.C1);sB.C2=Fpm::Add(b.R.C2,b.I.C2);
        W4 P0[3],N0[3],P1[3],N1[3],Pc[3],Nc[3];
        Mul3Wide(a.R,b.R,P0,N0); Mul3Wide(a.I,b.I,P1,N1); Mul3Wide(sA,sB,Pc,Nc);
        const W4 BIAS = *reinterpret_cast<const W4*>(Fpm::BIASW);
        Fp6m r; W4 pos,neg,ac,t;
        for(int i=0;i<3;i++){
            // r0_i = d0_i - 2 d1_i
            w4copy(pos,P0[i]); w4copy(t,N1[i]);w4muls(t,2);w4add(pos,t);
            w4copy(neg,N0[i]); w4copy(t,P1[i]);w4muls(t,2);w4add(neg,t);
            w4copy(ac,BIAS);w4add(ac,pos);w4sub(ac,neg); Fpx rr=RedcWide(ac);
            // r1_i = cross_i - d0_i - d1_i
            w4copy(pos,Pc[i]);w4add(pos,N0[i]);w4add(pos,N1[i]);
            w4copy(neg,Nc[i]);w4add(neg,P0[i]);w4add(neg,P1[i]);
            w4copy(ac,BIAS);w4add(ac,pos);w4sub(ac,neg); Fpx ri=RedcWide(ac);
            if(i==0){r.R.C0=rr;r.I.C0=ri;} else if(i==1){r.R.C1=rr;r.I.C1=ri;} else {r.R.C2=rr;r.I.C2=ri;}
        }
        return r;
    }

    // Squaring: same wide-accumulator combine as Mul3Wide with B==A, but each
    // of the 6 terms uses sqrmul2 (3 muls) instead of addmul2 (4). Result is
    // bit-identical to Mul3Wide(A,A,...).
    FORCE_INLINE static void Sqr3Wide(const Fp3m& A, W4 POS[3], W4 NEG[3]){
        const uint64_t a0=A.C0.A0,a0h=A.C0.A1,a1=A.C1.A0,a1h=A.C1.A1,a2=A.C2.A0,a2h=A.C2.A1;
        uint64_t A01,A01h,A02,A02h,A12,A12h; unsigned char cc;
        cc=Fpm::Adc64(0,a0,a1,A01);Fpm::Adc64(cc,a0h,a1h,A01h);
        cc=Fpm::Adc64(0,a0,a2,A02);Fpm::Adc64(cc,a0h,a2h,A02h);
        cc=Fpm::Adc64(0,a1,a2,A12);Fpm::Adc64(cc,a1h,a2h,A12h);
        W4 d0,d1,d2,M01,M02,M12; z4(d0);z4(d1);z4(d2);z4(M01);z4(M02);z4(M12);
        sqrmul2(d0,a0,a0h); sqrmul2(d1,a1,a1h); sqrmul2(d2,a2,a2h);
        sqrmul2(M01,A01,A01h); sqrmul2(M02,A02,A02h); sqrmul2(M12,A12,A12h);
        W4 t;
        w4copy(POS[0],d0); w4copy(t,d1);w4muls(t,4);w4add(POS[0],t); w4copy(t,d2);w4muls(t,4);w4add(POS[0],t); w4copy(NEG[0],M12);w4muls(NEG[0],4);
        w4copy(POS[1],M01); w4copy(NEG[1],d0);w4add(NEG[1],M12); w4copy(t,d2);w4muls(t,3);w4add(NEG[1],t);
        w4copy(POS[2],M02);w4add(POS[2],d1); w4copy(NEG[2],d0); w4copy(t,d2);w4muls(t,2);w4add(NEG[2],t);
    }
    FORCE_INLINE static Fp6m Sqr6Full(const Fp6m& a){
        Fp3m sA;
        sA.C0=Fpm::Add(a.R.C0,a.I.C0);sA.C1=Fpm::Add(a.R.C1,a.I.C1);sA.C2=Fpm::Add(a.R.C2,a.I.C2);
        W4 P0[3],N0[3],P1[3],N1[3],Pc[3],Nc[3];
        Sqr3Wide(a.R,P0,N0); Sqr3Wide(a.I,P1,N1); Sqr3Wide(sA,Pc,Nc);
        const W4 BIAS = *reinterpret_cast<const W4*>(Fpm::BIASW);
        Fp6m r; W4 pos,neg,ac,t;
        for(int i=0;i<3;i++){
            w4copy(pos,P0[i]); w4copy(t,N1[i]);w4muls(t,2);w4add(pos,t);
            w4copy(neg,N0[i]); w4copy(t,P1[i]);w4muls(t,2);w4add(neg,t);
            w4copy(ac,BIAS);w4add(ac,pos);w4sub(ac,neg); Fpx rr=RedcWide(ac);
            w4copy(pos,Pc[i]);w4add(pos,N0[i]);w4add(pos,N1[i]);
            w4copy(neg,Nc[i]);w4add(neg,P0[i]);w4add(neg,P1[i]);
            w4copy(ac,BIAS);w4add(ac,pos);w4sub(ac,neg); Fpx ri=RedcWide(ac);
            if(i==0){r.R.C0=rr;r.I.C0=ri;} else if(i==1){r.R.C1=rr;r.I.C1=ri;} else {r.R.C2=rr;r.I.C2=ri;}
        }
        return r;
    }
    FORCE_INLINE static Fp3m Sqr3(const Fp3m& a){ W4 P[3],N[3]; Sqr3Wide(a,P,N);
        const W4 BIAS=*reinterpret_cast<const W4*>(Fpm::BIASW); Fp3m r; W4 ac;
        w4copy(ac,BIAS); w4add(ac,P[0]); w4sub(ac,N[0]); r.C0=RedcWide(ac);
        w4copy(ac,BIAS); w4add(ac,P[1]); w4sub(ac,N[1]); r.C1=RedcWide(ac);
        w4copy(ac,BIAS); w4add(ac,P[2]); w4sub(ac,N[2]); r.C2=RedcWide(ac);
        return r;
    }
    FORCE_INLINE static Fp6m Sqr6(const Fp6m& a){ return Sqr6Full(a); }

    FORCE_INLINE static Fp3m Mul3(const Fp3m& a, const Fp3m& b) {
        return Mul3Lazy(a,b);
    }
    FORCE_INLINE static Fp3m Mul3_OLD(const Fp3m& a, const Fp3m& b) {
        Fpx d0 = Fpm::Mul(a.C0, b.C0);
        Fpx d1 = Fpm::Mul(a.C1, b.C1);
        Fpx d2 = Fpm::Mul(a.C2, b.C2);

        Fpx s = a.C0; AddMod(s, s, a.C1);
        Fpx w = b.C0; AddMod(w, w, b.C1);
        Fpx d01 = Fpm::Mul(s, w); SubMod(d01, d01, d0); SubMod(d01, d01, d1);

        s = a.C0; AddMod(s, s, a.C2);
        w = b.C0; AddMod(w, w, b.C2);
        Fpx d02 = Fpm::Mul(s, w); SubMod(d02, d02, d0); SubMod(d02, d02, d2);

        s = a.C1; AddMod(s, s, a.C2);
        w = b.C1; AddMod(w, w, b.C2);
        Fpx d12 = Fpm::Mul(s, w); SubMod(d12, d12, d1); SubMod(d12, d12, d2);

        Fp3m r;
        r.C0 = d0; r.C1 = d01; r.C2 = d02;
        AddMod(r.C2, r.C2, d1);
        Fpx m4a = Quad(d12);
        Fpx m4b = Quad(d2);
        SubMod(r.C0, d0, m4a);
        SubMod(r.C1, d01, d12);
        SubMod(r.C1, r.C1, m4b);
        SubMod(r.C2, r.C2, d2);
        return r;
    }

    FORCE_INLINE static Fp6m Mul6(const Fp6m& a, const Fp6m& b) {
        return Mul6Full(a, b);
    }
    FORCE_INLINE static Fp6m Mul6_OLD(const Fp6m& a, const Fp6m& b) {
        Fp3m d0 = Mul3(a.R, b.R);
        Fp3m d1 = Mul3(a.I, b.I);
        Fp3m r0 = d0; SubMod(r0, r0, Dbl3(d1));

        Fp3m s = a.R; AddMod(s, s, a.I);
        Fp3m w = b.R; AddMod(w, w, b.I);
        Fp3m r1 = Mul3(s, w);
        SubMod(r1, r1, d0); SubMod(r1, r1, d1);

        Fp6m r; r.R = r0; r.I = r1;
        return r;
    }

    FORCE_INLINE static Fp6m Add6(const Fp6m& a, const Fp6m& b) {
        Fp6m r;
        AddMod(r.R.C0, a.R.C0, b.R.C0); AddMod(r.R.C1, a.R.C1, b.R.C1); AddMod(r.R.C2, a.R.C2, b.R.C2);
        AddMod(r.I.C0, a.I.C0, b.I.C0); AddMod(r.I.C1, a.I.C1, b.I.C1); AddMod(r.I.C2, a.I.C2, b.I.C2);
        return r;
    }

    FORCE_INLINE static Fp6m Sub6(const Fp6m& a, const Fp6m& b) {
        Fp6m r;
        SubMod(r.R.C0, a.R.C0, b.R.C0); SubMod(r.R.C1, a.R.C1, b.R.C1); SubMod(r.R.C2, a.R.C2, b.R.C2);
        SubMod(r.I.C0, a.I.C0, b.I.C0); SubMod(r.I.C1, a.I.C1, b.I.C1); SubMod(r.I.C2, a.I.C2, b.I.C2);
        return r;
    }

    FORCE_INLINE static Fp6m FromFpxEmbedding(const Fpx& x) {
        Fp6m r;
        Fpx z = Fpm::FromBig(0);
        r.R.C0 = x; r.R.C1 = z; r.R.C2 = z;
        r.I.C0 = z; r.I.C1 = z; r.I.C2 = z;
        return r;
    }

    static Fp6m One() {
        Fp6m r;
        Fpx z = Fpm::FromBig(0);
        r.R.C0 = Fpm::One; r.R.C1 = z; r.R.C2 = z;
        r.I.C0 = z; r.I.C1 = z; r.I.C2 = z;
        return r;
    }

    static Fp6m PowBig(const Fp6m& a, const BigInteger& e_in) {
        BigInteger e = e_in;
        Fp6m r = One(), b = a;
        while (!e.IsZero() && e > 0) {
            if (!(e & 1).IsZero()) r = Mul6(r, b);
            b = Sqr6(b);               // dedicated squaring (was Mul6(b,b))
            e >>= 1;
        }
        return r;
    }

    // ---- Frobenius (x -> x^p) and the x^(p^3) conjugate -----------------
    // Precomputed once: u^p and u^(2p) as Fp3 constants, and whether y^p = -y.
    // In GF(p)[u]/(u^3+u+4) with a_i in Fp:  (a0+a1 u+a2 u^2)^p = a0 + a1 u^p
    // + a2 u^(2p). Over the y^2+2 extension: (c0+c1 y)^p = c0^p + (y^p/y) c1^p y.
    inline static Fp3m FrU1{}, FrU2{};
    inline static bool  FrNegY = false;
    inline static std::once_flag _frFlag;
    static void FrobInit() {
        std::call_once(_frFlag, []() {
            BigInteger p = Gf::GetP();
            Fp3 u(0, 1, 0);
            Fp3 U1 = u.Pow(p);         // u^p
            Fp3 U2 = U1 * U1;          // u^(2p) = (u^p)^2
            FrU1 = Fp3m::FromFp3(U1);
            FrU2 = Fp3m::FromFp3(U2);
            BigInteger negTwo = Gf::Mod(BigInteger(-2));
            BigInteger eps = BigInteger::ModPow(negTwo, (p - 1) / 2, p); // (-2)^((p-1)/2) = +/-1
            FrNegY = !eps.IsOne();     // y^p = eps * y
        });
    }
    // Fp3 Frobenius: a0 + a1*U1 + a2*U2  (a_i are Fp scalars -> scalar*Fp3)
    FORCE_INLINE static Fp3m Frob3(const Fp3m& a) {
        Fp3m r;
        r.C0 = a.C0; r.C1 = Fpm::FromBig(0); r.C2 = Fpm::FromBig(0);
        Fpx t0,t1,t2;
        t0=Fpm::Mul(a.C1,FrU1.C0); t1=Fpm::Mul(a.C1,FrU1.C1); t2=Fpm::Mul(a.C1,FrU1.C2);
        r.C0=Fpm::Add(r.C0,t0); r.C1=Fpm::Add(r.C1,t1); r.C2=Fpm::Add(r.C2,t2);
        t0=Fpm::Mul(a.C2,FrU2.C0); t1=Fpm::Mul(a.C2,FrU2.C1); t2=Fpm::Mul(a.C2,FrU2.C2);
        r.C0=Fpm::Add(r.C0,t0); r.C1=Fpm::Add(r.C1,t1); r.C2=Fpm::Add(r.C2,t2);
        return r;
    }
    FORCE_INLINE static Fp3m NegC(const Fp3m& a){ Fpx z=Fpm::FromBig(0);
        Fp3m r; r.C0=Fpm::Sub(z,a.C0); r.C1=Fpm::Sub(z,a.C1); r.C2=Fpm::Sub(z,a.C2); return r; }
    // x^p
    static Fp6m FrobP(const Fp6m& a) {
        FrobInit();
        Fp6m r; r.R = Frob3(a.R); r.I = Frob3(a.I);
        if (FrNegY) r.I = NegC(r.I);
        return r;
    }
    // x^(p^3): fixes Fp3, sends y -> -y, i.e. the Fp6/Fp3 conjugate (c0, -c1).
    static Fp6m Conj(const Fp6m& a) { Fp6m r; r.R = a.R; r.I = NegC(a.I); return r; }

    static Fp3m One3(){ Fp3m r; r.C0=Fpm::One; r.C1=Fpm::FromBig(0); r.C2=Fpm::FromBig(0); return r; }
    // Fp3 inverse by Fermat (a^(p^3-2)). ONE Fp3 inverse per final exp -- far
    // cheaper than a full Fp6 inverse. (Upgradeable to Itoh-Tsujii: reduce to a
    // single Fp inverse + Frobenius; left simple here.)
    static Fp3m Inv3(const Fp3m& a){
        BigInteger e = Gf::GetP(); e = e*e*e - BigInteger(2);
        Fp3m r=One3(), b=a;
        while(!e.IsZero() && e>0){ if(!(e&1).IsZero()) r=Mul3(r,b); b=Sqr3(b); e>>=1; }
        return r;
    }
    // Easy part of the final exponentiation, f^(p^3-1), WITHOUT a full Fp6
    // inverse: f^-1 = Conj(f)/norm(f), norm(f)=f*Conj(f) in Fp3; so
    // f^(p^3-1) = Conj(f) * f^-1 = Conj(f) * Conj(f) * norm^-1.
    static Fp6m EasyPart(const Fp6m& f){
        Fp6m cf = Conj(f);
        Fp6m nf = Mul6(f, cf);                 // norm lives in nf.R, nf.I == 0
        Fp3m ninv = Inv3(nf.R);
        Fp6m finv; finv.R = Mul3(cf.R, ninv); finv.I = Mul3(cf.I, ninv);  // f^-1
        return Mul6(cf, finv);                 // f^(p^3-1)
    }
    // Frobenius-factored final exponentiation: f^((p^6-1)/n)
    //   = ( f^(p^3-1) ) ^ ( (p^3+1)/n ).
    // Caller passes hardExp = (p^3+1)/n (n = curve order). Replaces a ~550-bit
    // PowBig over the full E with an inverse-free easy part + a ~220-bit PowBig.
    static Fp6m FinalExp(const Fp6m& f, const BigInteger& hardExp){
        return PowBig(EasyPart(f), hardExp);
    }

    static Fp6 ToFp6(const Fp6m& a) {
        return Fp6(
            Fp3(Fpm::ToBig(a.R.C0), Fpm::ToBig(a.R.C1), Fpm::ToBig(a.R.C2)),
            Fp3(Fpm::ToBig(a.I.C0), Fpm::ToBig(a.I.C1), Fpm::ToBig(a.I.C2))
        );
    }

    static Fp6 PowBigM(const Fp6& a, const BigInteger& e) {
        return ToFp6(PowBig(Fp6m::FromFp6(a), e));
    }
};
