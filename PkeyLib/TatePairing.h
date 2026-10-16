#pragma once
#include <cstdint>
#include "BigInteger.h"
#include "Montgomery.h"
#include <vector>

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

class TatePairing {
public:
    inline static BigInteger Order;

    static int BitLength(const BigInteger& v_in) {
        BigInteger v = v_in;
        int n = 0;
        while (!v.IsZero() && v > 0) {
            v >>= 1;
            n++;
        }
        return n;
    }

private:
    static BigInteger IntegerPow(const BigInteger& base, unsigned long exp) {
        BigInteger res;
        mpz_pow_ui(res.GetMpz(), base.GetMpz(), exp);
        return res;
    }

public:
    static Fp6 Pair(const BigInteger& n, const BigInteger& Tx, const BigInteger& Ty, const Pt6& Q, const BigInteger& a) {
        Fpm::Init();
        Fp6m f = Fp6mOp::One();

        Fpx vxm = Fpm::FromBig(Tx);
        Fpx vym = Fpm::FromBig(Ty);
        Fpx txm = vxm, tym = vym; // Keep track of initial T coordinates for addition
        Fpx am = Fpm::FromBig(a);
        Fpx zeroM = Fpm::FromBig(0);
        Fpx threeM = Fpm::FromBig(3);

        bool vInf = false;
        int L = BitLength(n);

        Fp6m Qx = Fp6m::FromFp6(Q.X);
        Fp6m Qy = Fp6m::FromFp6(Q.Y);

        for (int i = L - 2; i >= 0; i--) {
            // ---- Point Doubling V -> 2V ----
            if (!vInf) {
                if (vym.Eq(zeroM)) {
                    vInf = true;
                }
                else {
                    Fpx twoVy = Fpm::Add(vym, vym);
                    Fpx invTwoVy = Fpm::Inv(twoVy);
                    Fpx vxSq = Fpm::Mul(vxm, vxm);
                    Fpx numPart = Fpm::Add(Fpm::Mul(threeM, vxSq), am);
                    Fpx lam = Fpm::Mul(numPart, invTwoVy);

                    Fp6m embVy = Fp6mOp::FromFpxEmbedding(vym);
                    Fp6m embVx = Fp6mOp::FromFpxEmbedding(vxm);
                    Fp6m embLam = Fp6mOp::FromFpxEmbedding(lam);

                    Fp6m term1 = Fp6mOp::Sub6(Qy, embVy);
                    Fp6m term2 = Fp6mOp::Mul6(embLam, Fp6mOp::Sub6(Qx, embVx));
                    Fp6m l = Fp6mOp::Sub6(term1, term2);

                    Fpx lamSq = Fpm::Mul(lam, lam);
                    Fpx x2 = Fpm::Sub(lamSq, Fpm::Add(vxm, vxm));
                    Fpx y2 = Fpm::Sub(Fpm::Mul(lam, Fpm::Sub(vxm, x2)), vym);

                    f = Fp6mOp::Mul6(Fp6mOp::Mul6(f, f), l);
                    vxm = x2; vym = y2;
                }
            }
            else {
                f = Fp6mOp::Mul6(f, f);
            }

            // ---- Addition of T based on current bit of n ----
            // Allocation-free bit test (was: (n >> i) & One(), two mpz_t/bit).
            if (n.TestBit(i)) {
                if (vInf) {
                    vxm = txm; vym = tym; vInf = false;
                }
                else if (vxm.Eq(txm) && vym.Eq(tym)) {
                    Fpx twoVy = Fpm::Add(vym, vym);
                    Fpx invTwoVy = Fpm::Inv(twoVy);
                    Fpx vxSq = Fpm::Mul(vxm, vxm);
                    Fpx numPart = Fpm::Add(Fpm::Mul(threeM, vxSq), am);
                    Fpx lam = Fpm::Mul(numPart, invTwoVy);

                    Fp6m embVy = Fp6mOp::FromFpxEmbedding(vym);
                    Fp6m embVx = Fp6mOp::FromFpxEmbedding(vxm);
                    Fp6m embLam = Fp6mOp::FromFpxEmbedding(lam);

                    Fp6m term1 = Fp6mOp::Sub6(Qy, embVy);
                    Fp6m term2 = Fp6mOp::Mul6(embLam, Fp6mOp::Sub6(Qx, embVx));
                    Fp6m l = Fp6mOp::Sub6(term1, term2);

                    Fpx lamSq = Fpm::Mul(lam, lam);
                    Fpx x2 = Fpm::Sub(lamSq, Fpm::Add(vxm, vxm));
                    Fpx y2 = Fpm::Sub(Fpm::Mul(lam, Fpm::Sub(vxm, x2)), vym);

                    f = Fp6mOp::Mul6(f, l);
                    vxm = x2; vym = y2;
                }
                else if (vxm.Eq(txm)) {
                    f = Fp6mOp::Mul6(f, Fp6mOp::Sub6(Qx, Fp6mOp::FromFpxEmbedding(vxm)));
                    vInf = true;
                }
                else {
                    Fpx diffVxTx = Fpm::Sub(vxm, txm);
                    Fpx invDiff = Fpm::Inv(diffVxTx);
                    Fpx numDiff = Fpm::Sub(vym, tym);
                    Fpx lam = Fpm::Mul(numDiff, invDiff);

                    Fp6m embVy = Fp6mOp::FromFpxEmbedding(vym);
                    Fp6m embVx = Fp6mOp::FromFpxEmbedding(vxm);
                    Fp6m embLam = Fp6mOp::FromFpxEmbedding(lam);

                    Fp6m term1 = Fp6mOp::Sub6(Qy, embVy);
                    Fp6m term2 = Fp6mOp::Mul6(embLam, Fp6mOp::Sub6(Qx, embVx));
                    Fp6m l = Fp6mOp::Sub6(term1, term2);

                    Fpx lamSq = Fpm::Mul(lam, lam);
                    Fpx x3 = Fpm::Sub(Fpm::Sub(lamSq, vxm), txm);
                    Fpx y3 = Fpm::Sub(Fpm::Mul(lam, Fpm::Sub(vxm, x3)), vym);

                    f = Fp6mOp::Mul6(f, l);
                    vxm = x3; vym = y3;
                }
            }
        }

        // ---- Final Exponentiation ----
        BigInteger p = Gf::GetP();
        BigInteger p6 = IntegerPow(p, 6);
        BigInteger e = (p6 - 1) / n;
        return Fp6mOp::PowBigM(Fp6mOp::ToFp6(f), e);
    }

    // s * a for a base-field scalar s: six Fp multiplies. Same value as
    // Mul6(FromFpxEmbedding(s), a), without the full Fp6 multiply.
    FORCE_INLINE static Fp6m ScaleFp(const Fpx& s, const Fp6m& a) {
        Fp6m r;
        r.R.C0 = Fpm::Mul(s, a.R.C0); r.R.C1 = Fpm::Mul(s, a.R.C1); r.R.C2 = Fpm::Mul(s, a.R.C2);
        r.I.C0 = Fpm::Mul(s, a.I.C0); r.I.C1 = Fpm::Mul(s, a.I.C1); r.I.C2 = Fpm::Mul(s, a.I.C2);
        return r;
    }

    // Shared-ladder multi-pairing: runs the T-side Miller ladder ONCE and
    // evaluates the line function against every Q in the batch at each step,
    // instead of re-running the whole ladder (and re-paying every inversion)
    // once per Q. Mathematically this must produce the same f_i per point as
    // calling Pair() separately for each Q_i, since the doubling/addition of
    // V never depends on Q at all.
    static std::vector<Fp6> MultiPair(const BigInteger& n, const BigInteger& Tx, const BigInteger& Ty,
        const std::vector<Pt6>& Qs, const BigInteger& a) {
        Fpm::Init();
        Fpm::InitExpBits();
        int m = (int)Qs.size();
        std::vector<Fp6m> f(m, Fp6mOp::One());
        std::vector<Fp6m> Qx(m), Qy(m);
        for (int j = 0; j < m; j++) { Qx[j] = Fp6m::FromFp6(Qs[j].X); Qy[j] = Fp6m::FromFp6(Qs[j].Y); }

        Fpx vxm = Fpm::FromBig(Tx), vym = Fpm::FromBig(Ty);
        Fpx txm = vxm, tym = vym;
        Fpx am = Fpm::FromBig(a);
        Fpx zeroM = Fpm::FromBig(0);
        Fpx threeM = Fpm::FromBig(3);
        bool vInf = false;
        int L = BitLength(n);

        for (int i = L - 2; i >= 0; i--) {
            // ---- doubling: computed ONCE regardless of m ----
            // Three cases, mirroring the original exactly:
            //   (a) was already infinity  -> square f only
            //   (b) just became infinity (vym == 0) -> touch f not at all
            //   (c) normal                -> square-and-multiply by line
            enum { WAS_INF, BECAME_INF, NORMAL } dcase;
            Fpx vxOld = vxm, vyOld = vym, lam{};

            if (vInf) {
                dcase = WAS_INF;
            }
            else if (vym.Eq(zeroM)) {
                vInf = true;
                dcase = BECAME_INF;
            }
            else {
                dcase = NORMAL;
                Fpx invTwoVy = Fpm::InvFast(Fpm::Add(vym, vym));
                Fpx vxSq = Fpm::Mul(vxm, vxm);
                Fpx numPart = Fpm::Add(Fpm::Mul(threeM, vxSq), am);
                lam = Fpm::Mul(numPart, invTwoVy);

                Fpx lamSq = Fpm::Mul(lam, lam);
                Fpx x2 = Fpm::Sub(lamSq, Fpm::Add(vxm, vxm));
                Fpx y2 = Fpm::Sub(Fpm::Mul(lam, Fpm::Sub(vxm, x2)), vym);
                vxm = x2; vym = y2;
            }
            // V's coordinates are the same for every Q: embed them once.
            Fp6m embVy{}, embVx{};
            if (dcase == NORMAL) {
                embVy = Fp6mOp::FromFpxEmbedding(vyOld);
                embVx = Fp6mOp::FromFpxEmbedding(vxOld);
            }
            for (int j = 0; j < m; j++) {
                if (dcase == NORMAL) {
                    Fp6m term1 = Fp6mOp::Sub6(Qy[j], embVy);
                    Fp6m term2 = ScaleFp(lam, Fp6mOp::Sub6(Qx[j], embVx));
                    Fp6m l = Fp6mOp::Sub6(term1, term2);
                    f[j] = Fp6mOp::Mul6(Fp6mOp::Sqr6(f[j]), l);
                }
                else if (dcase == WAS_INF) {
                    f[j] = Fp6mOp::Sqr6(f[j]);
                }
                // BECAME_INF: untouched, matches original
            }

            // ---- addition of T based on current bit of n: also shared ----
            // Allocation-free bit test (was: (n >> i) & One(), two mpz_t/bit).
            if (n.TestBit(i)) {
                Fpx addVxOld{}, addVyOld{}, addLam{};
                bool haveAddLine = false;
                bool becameInf = false, tookT = false;

                if (vInf) {
                    vxm = txm; vym = tym; vInf = false; tookT = true;
                }
                else if (vxm.Eq(txm) && vym.Eq(tym)) {
                    addVxOld = vxm; addVyOld = vym;
                    Fpx invTwoVy = Fpm::InvFast(Fpm::Add(vym, vym));
                    Fpx vxSq = Fpm::Mul(vxm, vxm);
                    Fpx numPart = Fpm::Add(Fpm::Mul(threeM, vxSq), am);
                    addLam = Fpm::Mul(numPart, invTwoVy);
                    haveAddLine = true;

                    Fpx lamSq = Fpm::Mul(addLam, addLam);
                    Fpx x2 = Fpm::Sub(lamSq, Fpm::Add(vxm, vxm));
                    Fpx y2 = Fpm::Sub(Fpm::Mul(addLam, Fpm::Sub(vxm, x2)), vym);
                    vxm = x2; vym = y2;
                }
                else if (vxm.Eq(txm)) {
                    addVxOld = vxm;
                    becameInf = true;
                    vInf = true;
                }
                else {
                    addVxOld = vxm; addVyOld = vym;
                    Fpx diffVxTx = Fpm::Sub(vxm, txm);
                    Fpx invDiff = Fpm::InvFast(diffVxTx);
                    Fpx numDiff = Fpm::Sub(vym, tym);
                    addLam = Fpm::Mul(numDiff, invDiff);
                    haveAddLine = true;

                    Fpx lamSq = Fpm::Mul(addLam, addLam);
                    Fpx x3 = Fpm::Sub(Fpm::Sub(lamSq, vxm), txm);
                    Fpx y3 = Fpm::Sub(Fpm::Mul(addLam, Fpm::Sub(vxm, x3)), vym);
                    vxm = x3; vym = y3;
                }

                Fp6m embAddVy{}, embAddVx{};
                if (!tookT && (becameInf || haveAddLine)) {
                    embAddVx = Fp6mOp::FromFpxEmbedding(addVxOld);
                    if (haveAddLine) embAddVy = Fp6mOp::FromFpxEmbedding(addVyOld);
                }
                for (int j = 0; j < m; j++) {
                    if (tookT) continue; // pure reset, no line factor multiplied in original code
                    if (becameInf) {
                        f[j] = Fp6mOp::Mul6(f[j], Fp6mOp::Sub6(Qx[j], embAddVx));
                    }
                    else if (haveAddLine) {
                        Fp6m term1 = Fp6mOp::Sub6(Qy[j], embAddVy);
                        Fp6m term2 = ScaleFp(addLam, Fp6mOp::Sub6(Qx[j], embAddVx));
                        Fp6m l = Fp6mOp::Sub6(term1, term2);
                        f[j] = Fp6mOp::Mul6(f[j], l);
                    }
                }
            }
        }

        // Final exponentiation f^((p^6-1)/n).  (p^6-1)/n = (p^3-1) * ((p^3+1)/n),
        // and f^(p^3-1) costs one conjugate and one Fp3 inverse, which leaves a
        // power less than half as long. Only valid when n divides p^3+1 (true
        // for every pairing-friendly 2005 curve); otherwise use the plain power.
        BigInteger p = Gf::GetP();
        BigInteger p3 = IntegerPow(p, 3);
        std::vector<Fp6> out(m);
        if (((p3 + 1) % n).IsZero()) {
            BigInteger hard = (p3 + 1) / n;
            // Expand the hard exponent ONCE, reuse for every point. hard depends
            // on the curve order n (varies per bink), so it is NOT globally
            // cached -- but within this call all m points share it.
            std::vector<uint8_t> hardBits;
            Fp6mOp::ExpandBits(hard, hardBits);
            for (int j = 0; j < m; j++) out[j] = Fp6mOp::ToFp6(Fp6mOp::FinalExpBits(f[j], hardBits));
        }
        else {
            BigInteger p6 = IntegerPow(p, 6);
            BigInteger e = (p6 - 1) / n;
            for (int j = 0; j < m; j++) out[j] = Fp6mOp::PowBigM(Fp6mOp::ToFp6(f[j]), e);
        }
        return out;
    }
};
