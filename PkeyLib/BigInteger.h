#pragma once
#include <vector>
#include <mutex>
#include <stdexcept>
#include <string>
#include <sstream>
#include <iomanip>
#include <array>
#include <cstdlib>
#include <optional>

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

// -------------------------------------------------------------------------
// BigInteger Wrapper using GMP's raw C API (mpz_t) with RAII memory management
// -------------------------------------------------------------------------
class BigInteger {
private:
    mpz_t val;

public:
    BigInteger() {
        mpz_init(val);
    }

    BigInteger(int v) {
        mpz_init_set_si(val, v);
    }

    BigInteger(long long v) {
        mpz_init(val);
        mpz_set_str(val, std::to_string(v).c_str(), 10);
    }

    BigInteger(unsigned long long v) {
        mpz_init(val);
        mpz_set_str(val, std::to_string(v).c_str(), 10);
    }

    BigInteger(const char* s) {
        mpz_init_set_str(val, s ? s : "0", 10);
    }

    BigInteger(const std::string& s) {
        mpz_init_set_str(val, s.c_str(), 10);
    }

    // Copy Constructor
    BigInteger(const BigInteger& other) {
        mpz_init_set(val, other.val);
    }

    // Move Constructor
    BigInteger(BigInteger&& other) noexcept {
        mpz_init(val);
        mpz_swap(val, other.val);
    }

    // Destructor
    ~BigInteger() {
        mpz_clear(val);
    }

    // Copy Assignment
    BigInteger& operator=(const BigInteger& other) {
        if (this != &other) {
            mpz_set(val, other.val);
        }
        return *this;
    }

    // Move Assignment
    BigInteger& operator=(BigInteger&& other) noexcept {
        if (this != &other) {
            mpz_swap(val, other.val);
        }
        return *this;
    }

    // ★ Added Bitwise OR Operator using GMP's mpz_ior
    BigInteger operator|(const BigInteger& o) const {
        BigInteger res;
        mpz_ior(res.val, val, o.val);
        return res;
    }

    // Conversion operator to int
    explicit operator int() const {
        return (int)mpz_get_si(val);
    }

    friend BigInteger operator*(int lhs, const BigInteger& rhs) {
        return BigInteger(lhs) * rhs;
    }

    friend BigInteger operator*(const BigInteger& lhs, int rhs) {
        return lhs * BigInteger(rhs);
    }

    // Access underlying mpz_t
    const mpz_t& GetMpz() const { return val; }
    mpz_t& GetMpz() { return val; }

    // C# Style Properties & Methods
    bool IsZero() const { return mpz_sgn(val) == 0; }
    bool IsOne() const { return mpz_cmp_si(val, 1) == 0; }
    bool IsEven() const { return mpz_even_p(val) != 0; }
    int Sign() const { return mpz_sgn(val); }

    // Allocation-free single-bit read. Replaces the (x >> i) & 1 idiom in the
    // hot exponent/ladder walks, which constructs two temporary mpz_t per bit.
    bool TestBit(int i) const { return mpz_tstbit(val, (mp_bitcnt_t)i) != 0; }

    static BigInteger Zero() { return BigInteger(0); }
    static BigInteger One() { return BigInteger(1); }

    // One truncating division yielding BOTH quotient and remainder (like MS's
    // mp_divrem_knuth), instead of a separate B/p and B%p. For non-negative
    // B and p the remainder is already in [0,p).
    static BigInteger DivRem(const BigInteger& n, const BigInteger& d, BigInteger& rem) {
        BigInteger q;
        mpz_tdiv_qr(q.val, rem.val, n.val, d.val);
        return q;
    }

    // GMP 模幂加速
    static BigInteger ModPow(const BigInteger& base, const BigInteger& exp, const BigInteger& mod) {
        BigInteger res;
        mpz_powm(res.val, base.val, exp.val, mod.val);
        return res;
    }

    // Arithmetic Operators
    BigInteger operator+(const BigInteger& o) const {
        BigInteger res;
        mpz_add(res.val, val, o.val);
        return res;
    }

    BigInteger operator-(const BigInteger& o) const {
        BigInteger res;
        mpz_sub(res.val, val, o.val);
        return res;
    }

    BigInteger operator-() const {
        BigInteger res;
        mpz_neg(res.val, val);
        return res;
    }

    BigInteger operator*(const BigInteger& o) const {
        BigInteger res;
        mpz_mul(res.val, val, o.val);
        return res;
    }

    BigInteger operator/(const BigInteger& o) const {
        BigInteger res;
        mpz_tdiv_q(res.val, val, o.val);
        return res;
    }

    BigInteger operator%(const BigInteger& o) const {
        BigInteger res;
        mpz_tdiv_r(res.val, val, o.val);
        return res;
    }

    BigInteger& operator+=(const BigInteger& o) { mpz_add(val, val, o.val); return *this; }
    BigInteger& operator-=(const BigInteger& o) { mpz_sub(val, val, o.val); return *this; }
    BigInteger& operator*=(const BigInteger& o) { mpz_mul(val, val, o.val); return *this; }
    BigInteger& operator/=(const BigInteger& o) { mpz_tdiv_q(val, val, o.val); return *this; }
    BigInteger& operator%=(const BigInteger& o) { mpz_tdiv_r(val, val, o.val); return *this; }

    // Shift & Bitwise Operators
    BigInteger operator>>(int s) const {
        BigInteger res;
        mpz_fdiv_q_2exp(res.val, val, s);
        return res;
    }

    BigInteger operator<<(int s) const {
        BigInteger res;
        mpz_mul_2exp(res.val, val, s);
        return res;
    }

    BigInteger& operator>>=(int s) { mpz_fdiv_q_2exp(val, val, s); return *this; }
    BigInteger& operator<<=(int s) { mpz_mul_2exp(val, val, s); return *this; }

    BigInteger operator&(const BigInteger& o) const {
        BigInteger res;
        mpz_and(res.val, val, o.val);
        return res;
    }

    // Comparison Operators
    bool operator==(const BigInteger& o) const { return mpz_cmp(val, o.val) == 0; }
    bool operator!=(const BigInteger& o) const { return mpz_cmp(val, o.val) != 0; }
    bool operator<(const BigInteger& o) const { return mpz_cmp(val, o.val) < 0; }
    bool operator>(const BigInteger& o) const { return mpz_cmp(val, o.val) > 0; }
    bool operator<=(const BigInteger& o) const { return mpz_cmp(val, o.val) <= 0; }
    bool operator>=(const BigInteger& o) const { return mpz_cmp(val, o.val) >= 0; }

    friend std::ostream& operator<<(std::ostream& os, const BigInteger& b) {
        char* str = mpz_get_str(nullptr, 10, b.val);
        if (str) {
            os << str;
            free(str);
        }
        return os;
    }
};

// -------------------------------------------------------------------------
// 异常与基域 GF(p) 管理
// -------------------------------------------------------------------------
class BaseFieldMismatchException : public std::runtime_error {
public:
    explicit BaseFieldMismatchException(const std::string& message)
        : std::runtime_error(message) {}
};

class Gf {
private:
    inline static BigInteger _p = 0;
    inline static std::mutex _pGate;   // guards WRITES only (parse phase)

public:
    // Lock-free read. Safe because all writes (SetP) happen during the serial
    // parse phase, before any worker thread starts -- thread creation provides
    // the happens-before. This mirrors MS, which keeps p in static storage and
    // never locks on read; it removes a global-mutex acquire that every Fp3/
    // Curve/lift_x op was paying, 16x contended under the parallel scan.
    static BigInteger GetP() {
        return _p;
    }

    static void SetP(const BigInteger& value) {
        std::lock_guard<std::mutex> lock(_pGate);
        if (_p.IsZero()) {
            _p = value;
            return;
        }
        if (_p != value) {
            std::stringstream ss1, ss2;
            ss1 << _p; ss2 << value;
            throw BaseFieldMismatchException(
                "公钥基域不一致：此前已用 p=" + ss1.str() +
                "，本组公钥是 p=" + ss2.str() + "。"
            );
        }
    }

    static BigInteger Mod(const BigInteger& x) {
        BigInteger currentP = GetP();
        if (x.Sign() >= 0 && x < currentP) return x;
        BigInteger r = x % currentP;
        return (r.Sign() < 0) ? (r + currentP) : r;
    }

    // 1/(-2) mod p in closed form. p is odd, so -2 * (p-1)/2 = -(p-1) = 1 (mod p),
    // i.e. the inverse of -2 is exactly (p-1)/2. No extended-GCD loop at all.
    static BigInteger InvNeg2() {
        BigInteger half = GetP() - BigInteger(1);
        half >>= 1;
        return half;
    }

    static BigInteger Inv(const BigInteger& a) {
        BigInteger currentP = GetP();
        BigInteger modA = Mod(a);
        if (modA.IsZero()) throw std::runtime_error("Gf.Inv(0)");
        if (modA.IsOne()) return BigInteger::One();

        BigInteger u = modA, v = currentP, x1 = BigInteger::One(), x2 = BigInteger::Zero();
        while (!u.IsOne() && !v.IsOne()) {
            while (u.IsEven()) {
                u >>= 1;
                x1 = x1.IsEven() ? (x1 >> 1) : ((x1 + currentP) >> 1);
            }
            while (v.IsEven()) {
                v >>= 1;
                x2 = x2.IsEven() ? (x2 >> 1) : ((x2 + currentP) >> 1);
            }
            if (u >= v) { u -= v; x1 -= x2; }
            else { v -= u; x2 -= x1; }
        }
        BigInteger r = u.IsOne() ? x1 : x2;
        r %= currentP;
        return (r.Sign() < 0) ? (r + currentP) : r;
    }
};

// -------------------------------------------------------------------------
// GF(p^3) 域扩展
// -------------------------------------------------------------------------
struct Fp3 {
    BigInteger C0, C1, C2;

    Fp3() : C0(0), C1(0), C2(0) {}
    Fp3(const BigInteger& c0, const BigInteger& c1, const BigInteger& c2) {
        C0 = Gf::Mod(c0); C1 = Gf::Mod(c1); C2 = Gf::Mod(c2);
    }

    static Fp3 Zero() { return Fp3(0, 0, 0); }
    static Fp3 One() { return Fp3(1, 0, 0); }

    bool IsZero() const { return C0.IsZero() && C1.IsZero() && C2.IsZero(); }

    friend Fp3 operator+(const Fp3& a, const Fp3& b) {
        return Fp3(a.C0 + b.C0, a.C1 + b.C1, a.C2 + b.C2);
    }
    friend Fp3 operator-(const Fp3& a, const Fp3& b) {
        return Fp3(a.C0 - b.C0, a.C1 - b.C1, a.C2 - b.C2);
    }
    Fp3 operator-() const { return Fp3(-C0, -C1, -C2); }

    inline static thread_local std::array<BigInteger, 5> _scratch;

    friend Fp3 operator*(const Fp3& a, const Fp3& b) {
        auto& t = _scratch;
        t[0] = 0; t[1] = 0; t[2] = 0; t[3] = 0; t[4] = 0;
        for (int i = 0; i < 3; i++) {
            BigInteger ai = a.Get(i);
            if (ai.IsZero()) continue;
            for (int j = 0; j < 3; j++) t[i + j] += ai * b.Get(j);
        }
        for (int k = 4; k >= 3; k--) {
            BigInteger c = t[k];
            if (c.IsZero()) continue;
            t[k] = 0;
            t[k - 3] -= 4 * c;
            t[k - 2] -= c;
        }
        return Fp3(t[0], t[1], t[2]);
    }

    friend Fp3 operator*(const Fp3& a, const BigInteger& s) {
        return Fp3(a.C0 * s, a.C1 * s, a.C2 * s);
    }

    BigInteger Get(int i) const { return (i == 0) ? C0 : ((i == 1) ? C1 : C2); }

    Fp3 Pow(const BigInteger& e_in) const {
        BigInteger e = e_in;
        Fp3 r = One(), b = *this;
        while (!e.IsZero() && e > 0) {
            if (!(e & 1).IsZero()) r = r * b;
            b = b * b;
            e >>= 1;
        }
        return r;
    }

    Fp3 Inverse() const {
        BigInteger p = Gf::GetP();
        BigInteger p3 = p * p * p;
        return Pow(p3 - 2);
    }

    bool operator==(const Fp3& o) const { return C0 == o.C0 && C1 == o.C1 && C2 == o.C2; }
    bool operator!=(const Fp3& o) const { return !(*this == o); }

    std::string ToString() const {
        std::stringstream ss;
        ss << "(" << C0 << ", " << C1 << ", " << C2 << ")";
        return ss.str();
    }
};

// -------------------------------------------------------------------------
// GF(p^6) 域扩展
// -------------------------------------------------------------------------
struct Fp6 {
    Fp3 C0, C1;

    Fp6() : C0(Fp3::Zero()), C1(Fp3::Zero()) {}
    Fp6(const Fp3& c0, const Fp3& c1) : C0(c0), C1(c1) {}

    static Fp6 Zero() { return Fp6(Fp3::Zero(), Fp3::Zero()); }
    static Fp6 One() { return Fp6(Fp3::One(), Fp3::Zero()); }

    static Fp6 FromFp3(const Fp3& a) { return Fp6(a, Fp3::Zero()); }

    friend Fp6 operator+(const Fp6& a, const Fp6& b) {
        return Fp6(a.C0 + b.C0, a.C1 + b.C1);
    }

    friend Fp6 operator-(const Fp6& a, const Fp6& b) {
        return Fp6(a.C0 - b.C0, a.C1 - b.C1);
    }

    Fp6 operator-() const {
        return Fp6(-C0, -C1);
    }

    friend Fp6 operator*(const Fp6& a, const Fp6& b) {
        Fp3 p1 = a.C1 * b.C1;
        Fp3 two_p1 = Fp3(p1.C0 + p1.C0, p1.C1 + p1.C1, p1.C2 + p1.C2);
        Fp3 r0 = a.C0 * b.C0 - two_p1;
        Fp3 r1 = a.C0 * b.C1 + a.C1 * b.C0;
        return Fp6(r0, r1);
    }

    Fp6 Pow(const BigInteger& e_in) const {
        BigInteger e = e_in;
        Fp6 r = One(), b = *this;
        while (!e.IsZero() && e > 0) {
            if (!(e & 1).IsZero()) r = r * b;
            b = b * b;
            e >>= 1;
        }
        return r;
    }

    bool operator==(const Fp6& o) const { return C0 == o.C0 && C1 == o.C1; }
    bool operator!=(const Fp6& o) const { return !(*this == o); }

    std::string ToString() const {
        return C0.ToString() + " + " + C1.ToString() + "*y";
    }

    Fp6 Inverse() const {
        BigInteger p = Gf::GetP();
        BigInteger p2 = p * p;
        BigInteger p3 = p2 * p;
        BigInteger p6 = p3 * p3;
        return Pow(p6 - 2);
    }

    friend Fp6 operator/(const Fp6& a, const Fp6& b) {
        return a * b.Inverse();
    }
};

// -------------------------------------------------------------------------
// GF(p) 上的开平方（Tonelli-Shanks 算法）
// -------------------------------------------------------------------------
class Fp {
public:
    static bool IsQR(const BigInteger& a_in) {
        BigInteger a = Gf::Mod(a_in);
        if (a.IsZero()) return true;
        BigInteger p = Gf::GetP();
        return BigInteger::ModPow(a, (p - 1) / 2, p).IsOne();
    }

private:
    static BigInteger NonResidue() {
        BigInteger p = Gf::GetP();
        for (BigInteger z = 2; ; z = z + 1) {
            if (!BigInteger::ModPow(z, (p - 1) / 2, p).IsOne())
                return z;
        }
    }

public:
    static std::optional<BigInteger> Sqrt(const BigInteger& a_in) {
        BigInteger a = Gf::Mod(a_in);
        BigInteger p = Gf::GetP();
        if (a.IsZero()) return BigInteger::Zero();
        if (!IsQR(a)) return std::nullopt;

        BigInteger q = p - 1;
        int s = 0;
        while (q.IsEven()) {
            q >>= 1;
            s++;
        }

        BigInteger z = NonResidue();
        BigInteger c = BigInteger::ModPow(z, q, p);
        BigInteger x = BigInteger::ModPow(a, (q + 1) / 2, p);
        BigInteger t = BigInteger::ModPow(a, q, p);
        int m = s;

        while (!t.IsOne()) {
            int i = 0;
            BigInteger tt = t;
            while (!tt.IsOne()) {
                tt = (tt * tt) % p;
                i++;
                if (i >= m) return std::nullopt;
            }

            BigInteger b = c;
            for (int k = 0; k < m - i - 1; k++) {
                b = (b * b) % p;
            }
            x = (x * b) % p;
            BigInteger b2 = (b * b) % p;
            t = (t * b2) % p;
            c = b2;
            m = i;
        }

        BigInteger hi = p - x;
        return (x <= hi) ? x : hi;
    }
};

// -------------------------------------------------------------------------
// 椭圆曲线工具类 Curve
// -------------------------------------------------------------------------
class Curve {
public:
    static BigInteger Rhs(const BigInteger& x_in, const BigInteger& a, const BigInteger& b) {
        BigInteger x = Gf::Mod(x_in);
        BigInteger p = Gf::GetP();
        // x * x % p * x + a * x + b
        BigInteger term1 = (x * x) % p;
        term1 = (term1 * x) % p;
        BigInteger term2 = (a * x) % p;
        BigInteger res = term1 + term2 + b;
        return Gf::Mod(res);
    }

    /// <summary>返回曲线上的两个点 (x, ±y)，无解时返回 std::nullopt。</summary>
    static std::optional<std::vector<BigInteger>> LiftX(const BigInteger& x_in, const BigInteger& a, const BigInteger& b) {
        BigInteger x = Gf::Mod(x_in);
        std::optional<BigInteger> y = Fp::Sqrt(Rhs(x, a, b));
        if (!y.has_value()) return std::nullopt;

        BigInteger yv = y.value();
        BigInteger negYv = Gf::Mod(-yv);
        return std::vector<BigInteger>{ yv, negYv };
    }
};

// -------------------------------------------------------------------------
// E6 上的点（坐标在 Fp6 里）
// -------------------------------------------------------------------------
struct Pt6 {
    Fp6 X, Y;
    bool Inf;

    Pt6() : X(Fp6::Zero()), Y(Fp6::Zero()), Inf(true) {}

    static Pt6 Infinity() {
        Pt6 p;
        p.X = Fp6::Zero();
        p.Y = Fp6::Zero();
        p.Inf = true;
        return p;
    }

    static Pt6 From(const Fp6& x, const Fp6& y) {
        Pt6 p;
        p.X = x;
        p.Y = y;
        p.Inf = false;
        return p;
    }

    /// <summary>倍点：λ = (3x² + a) / 2y</summary>
    Pt6 Dbl(const BigInteger& a) const {
        if (Inf || Y == Fp6::Zero()) return Infinity();
        Fp3 two(2, 0, 0);
        Fp6 num = (X * X) * Fp6::FromFp3(Fp3(3, 0, 0)) + Fp6::FromFp3(Fp3(a, 0, 0));
        Fp6 den = Y * Fp6::FromFp3(two);
        Fp6 lam = num / den;
        Fp6 x3 = lam * lam - X - X;
        Fp6 y3 = lam * (X - x3) - Y;
        return From(x3, y3);
    }

    static Pt6 Add(const Pt6& P, const Pt6& Q) {
        if (P.Inf) return Q;
        if (Q.Inf) return P;
        if (P.X == Q.X) {
            if (P.Y == Q.Y) {
                throw std::runtime_error("Add: 需要倍点，请调用 Dbl");
            }
            return Infinity(); // 互为逆元
        }
        Fp6 lam = (Q.Y - P.Y) / (Q.X - P.X);
        Fp6 x3 = lam * lam - P.X - Q.X;
        Fp6 y3 = lam * (P.X - x3) - P.Y;
        return From(x3, y3);
    }

    Pt6 Mul(const BigInteger& k_in, const BigInteger& a) const {
        BigInteger k = k_in;
        Pt6 r = Infinity(), b = *this;
        while (!k.IsZero() && k > 0) {
            if (!(k & 1).IsZero()) r = r.Inf ? b : Add(r, b);
            b = b.Dbl(a);
            k >>= 1;
        }
        return r;
    }
};
