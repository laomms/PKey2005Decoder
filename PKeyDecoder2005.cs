using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Numerics;
using System.Runtime.CompilerServices;


namespace PidKeyPlugIn
{
    // -------------------------------------------------------------------------
    // 基域 GF(p)：直接用 BigInteger，只做模约简
    // -------------------------------------------------------------------------
    /// <summary>
    /// 公钥基域不一致。只在「多组公钥的 p 不同」时抛出 —— 见 <see cref="Gf.P"/>。
    /// 单独一个类型是为了让 TryPubKey 能精确地把它放行，不跟公钥解析失败混在一起。
    /// </summary>
    public class BaseFieldMismatchException : InvalidOperationException
    {
        public BaseFieldMismatchException(string message) : base(message) { }
    }

    public static class Gf
    {
        static BigInteger _p;
        static readonly object _pGate = new object();

        /// <summary>
        /// 基域素数 p，由 PubKeyParser.Parse 写入。
        ///
        /// ★ 带一致性校验：首次写入生效；之后再写**必须相同**，否则抛
        ///   <see cref="BaseFieldMismatchException"/>。
        ///
        ///   为什么需要：Parse 会写这个静态字段，而 Parallel.ForEach 里每个线程
        ///   都会 Parse。如果配置里混了不同基域的公钥，p 会被互相覆盖、算到一半
        ///   被换掉，结果是**错的而且不报错**。宁可在这里炸掉。
        /// </summary>
        public static BigInteger P
        {
            get { return _p; }
            set
            {
                lock (_pGate)
                {
                    if (_p.IsZero) { _p = value; return; }
                    if (_p != value)
                        throw new BaseFieldMismatchException(
                            "公钥基域不一致：此前已用 p=" + _p.ToString("X") +
                            "，本组公钥是 p=" + value.ToString("X") +
                            "。多组公钥必须共用同一个基域；若配置里确实混了不同基域的公钥，" +
                            "请先在循环外筛掉，否则并发解码会静默算错。");
                }
            }
        }

        public static BigInteger Mod(BigInteger x)
        {
            // 热路径里绝大多数输入已经是规范形式，先做一次廉价的范围判断，
            // 命中就不必走 `%`（BigInteger 的除法比比较贵得多）。
            if (x.Sign >= 0 && x < P) return x;
            BigInteger r = x % P;
            return r.Sign < 0 ? r + P : r;
        }

        /// <summary>
        /// 模逆。用**二进制扩展欧几里得**而不是费马小定理的 ModPow：
        /// 110 位下前者约 220 轮「位移/加减」，后者要做约 165 次模乘 ——
        /// 而 Miller 循环每一比特都要调用一次，是配对里最热的一段。
        /// </summary>
        public static BigInteger Inv(BigInteger a)
        {
            a = Mod(a);
            if (a.IsZero) throw new DivideByZeroException("Gf.Inv(0)");
            if (a.IsOne) return BigInteger.One;

            BigInteger u = a, v = P, x1 = BigInteger.One, x2 = BigInteger.Zero;
            while (!u.IsOne && !v.IsOne)
            {
                while (u.IsEven) { u >>= 1; x1 = x1.IsEven ? (x1 >> 1) : ((x1 + P) >> 1); }
                while (v.IsEven) { v >>= 1; x2 = x2.IsEven ? (x2 >> 1) : ((x2 + P) >> 1); }
                if (u >= v) { u -= v; x1 -= x2; } else { v -= u; x2 -= x1; }
            }
            BigInteger r = u.IsOne ? x1 : x2;
            r %= P;
            return r.Sign < 0 ? r + P : r;
        }
    }

    // -------------------------------------------------------------------------
    // GF(p^3) = GF(p)[u] / (u^3 + u + 4)
    // 元素 c0 + c1*u + c2*u^2
    // 归约关系：u^3 = -(u + 4) = -4 - u
    // -------------------------------------------------------------------------
    public struct Fp3 : IEquatable<Fp3>
    {
        public BigInteger C0, C1, C2;

        public Fp3(BigInteger c0, BigInteger c1, BigInteger c2)
        {
            C0 = Gf.Mod(c0); C1 = Gf.Mod(c1); C2 = Gf.Mod(c2);
        }

        // 注意：必须是属性而不是 static readonly 字段 —— 字段初始化会在 Gf.P
        // 赋值之前跑（静态构造顺序），Fp3 的构造函数要取模，会除零。
        public static Fp3 Zero { get { return new Fp3(0, 0, 0); } }
        public static Fp3 One { get { return new Fp3(1, 0, 0); } }

        public bool IsZero { get { return C0.IsZero && C1.IsZero && C2.IsZero; } }

        public static Fp3 operator +(Fp3 a, Fp3 b)
        {
            return new Fp3(a.C0 + b.C0, a.C1 + b.C1, a.C2 + b.C2);
        }

        public static Fp3 operator -(Fp3 a, Fp3 b)
        {
            return new Fp3(a.C0 - b.C0, a.C1 - b.C1, a.C2 - b.C2);
        }

        public static Fp3 operator -(Fp3 a)
        {
            return new Fp3(-a.C0, -a.C1, -a.C2);
        }

        // 线程静态暂存：Fp3 乘法在热路径里被调用几百万次，每次 new 一个数组
        // 会让 GC 变成瓶颈。operator* 不会重入，所以复用是安全的。
        [ThreadStatic] static BigInteger[] _scratch;

        public static Fp3 operator *(Fp3 a, Fp3 b)
        {
            BigInteger[] t = _scratch;
            if (t == null) { t = new BigInteger[5]; _scratch = t; }
            t[0] = BigInteger.Zero; t[1] = BigInteger.Zero; t[2] = BigInteger.Zero;
            t[3] = BigInteger.Zero; t[4] = BigInteger.Zero;
            for (int i = 0; i < 3; i++)
            {
                if (a.Get(i).IsZero) continue;
                for (int j = 0; j < 3; j++)
                    t[i + j] += a.Get(i) * b.Get(j);
            }
            // 从高次往低次归约：u^k = u^(k-3) * u^3 = u^(k-3) * (-4 - u)
            for (int k = 4; k >= 3; k--)
            {
                BigInteger c = t[k];
                if (c.IsZero) continue;
                t[k] = 0;
                t[k - 3] -= 4 * c;
                t[k - 2] -= c;
            }
            return new Fp3(t[0], t[1], t[2]);
        }

        public static Fp3 operator *(Fp3 a, BigInteger s)
        {
            return new Fp3(a.C0 * s, a.C1 * s, a.C2 * s);
        }

        public BigInteger Get(int i)
        {
            return i == 0 ? C0 : (i == 1 ? C1 : C2);
        }

        public Fp3 Pow(BigInteger e)
        {
            Fp3 r = One, b = this;
            while (e > 0)
            {
                if (!(e & 1).IsZero) r = r * b;
                b = b * b;
                e >>= 1;
            }
            return r;
        }

        /// <summary>逆元：a^(p^3 - 2)</summary>
        public Fp3 Inverse()
        {
            return Pow(BigInteger.Pow(Gf.P, 3) - 2);
        }

        public bool Equals(Fp3 o)
        {
            return C0 == o.C0 && C1 == o.C1 && C2 == o.C2;
        }

        public override bool Equals(object o) { return o is Fp3 && Equals((Fp3)o); }

        public override int GetHashCode()
        {
            return C0.GetHashCode() ^ (C1.GetHashCode() << 1) ^ (C2.GetHashCode() << 2);
        }

        public override string ToString()
        {
            return "(" + C0 + ", " + C1 + ", " + C2 + ")";
        }
    }

    // -------------------------------------------------------------------------
    // GF(p^6) = GF(p^3)[y] / (y^2 + 2)
    // 元素 c0 + c1*y，其中 c0/c1 是 Fp3；y^2 = -2
    // -------------------------------------------------------------------------
    public struct Fp6 : IEquatable<Fp6>
    {
        public Fp3 C0, C1;

        public Fp6(Fp3 c0, Fp3 c1) { C0 = c0; C1 = c1; }

        // 同上：用属性，避免静态构造顺序问题
        public static Fp6 Zero { get { return new Fp6(Fp3.Zero, Fp3.Zero); } }
        public static Fp6 One { get { return new Fp6(Fp3.One, Fp3.Zero); } }

        /// <summary>把 K3 元素嵌进 K6（虚部为 0）</summary>
        public static Fp6 FromFp3(Fp3 a) { return new Fp6(a, Fp3.Zero); }

        public static Fp6 operator +(Fp6 a, Fp6 b)
        {
            return new Fp6(a.C0 + b.C0, a.C1 + b.C1);
        }

        public static Fp6 operator -(Fp6 a, Fp6 b)
        {
            return new Fp6(a.C0 - b.C0, a.C1 - b.C1);
        }

        public static Fp6 operator -(Fp6 a)
        {
            return new Fp6(-a.C0, -a.C1);
        }

        public static Fp6 operator *(Fp6 a, Fp6 b)
        {
            //  y^2 = -2
            //  (a0 + a1 y)(b0 + b1 y) = (a0 b0 - 2 a1 b1) + (a0 b1 + a1 b0) y
            //  ★ 原来写成 `new Fp3(2,0,0) * (a1*b1)` —— 为了乘个常数 2
            //    白做一整次 Fp3 乘法（9 次 BigInteger 乘法）+ 一次构造（3 次取模）。
            //    直接对系数加倍即可。
            Fp3 p1 = a.C1 * b.C1;
            Fp3 two_p1 = new Fp3(p1.C0 + p1.C0, p1.C1 + p1.C1, p1.C2 + p1.C2);
            Fp3 r0 = a.C0 * b.C0 - two_p1;
            Fp3 r1 = a.C0 * b.C1 + a.C1 * b.C0;
            return new Fp6(r0, r1);
        }

        public Fp6 Pow(BigInteger e)
        {
            Fp6 r = One, b = this;
            while (e > 0)
            {
                if (!(e & 1).IsZero) r = r * b;
                b = b * b;
                e >>= 1;
            }
            return r;
        }

        public bool Equals(Fp6 o) { return C0.Equals(o.C0) && C1.Equals(o.C1); }

        public override bool Equals(object o) { return o is Fp6 && Equals((Fp6)o); }

        public override int GetHashCode() { return C0.GetHashCode() ^ C1.GetHashCode(); }

        public override string ToString() { return C0 + " + " + C1 + "*y"; }

        /// <summary>逆元：a^(p^6 - 2)</summary>
        public Fp6 Inverse()
        {
            return Pow(BigInteger.Pow(Gf.P, 6) - 2);
        }

        public static Fp6 operator /(Fp6 a, Fp6 b) { return a * b.Inverse(); }
    }

    // -------------------------------------------------------------------------
    // GF(p) 上的开平方（Tonelli-Shanks）
    //
    // 注意 p ≡ 1 (mod 4)，所以**不能**用 a^((p+1)/4) 那个简化公式，必须走完整
    // Tonelli-Shanks：
    //     p - 1 = q * 2^s   (q 为奇数)
    //     找非二次剩余 z，令 c = z^q, x = a^((q+1)/2), t = a^q, m = s
    //     while t != 1: 找最小的 i 使 t^(2^i) == 1; b = c^(2^(m-i-1));
    //                   x *= b; t *= b^2; c = b^2; m = i
    // -------------------------------------------------------------------------
    public static class Fp
    {
        /// <summary>勒让德符号是否为 1（即是否为非零二次剩余）</summary>
        public static bool IsQR(BigInteger a)
        {
            a = Gf.Mod(a);
            if (a.IsZero) return true;
            return BigInteger.ModPow(a, (Gf.P - 1) / 2, Gf.P).IsOne;
        }

        /// <summary>返回一个非二次剩余（用于 Tonelli-Shanks 的 z）</summary>
        static BigInteger NonResidue()
        {
            for (BigInteger z = 2; ; z++)
                if (!BigInteger.ModPow(z, (Gf.P - 1) / 2, Gf.P).IsOne)
                    return z;
        }

        /// <summary>模 p 开平方。不是二次剩余时返回 null；返回的是两根中较小的那个。</summary>
        public static BigInteger? Sqrt(BigInteger a)
        {
            a = Gf.Mod(a);
            if (a.IsZero) return BigInteger.Zero;
            if (!IsQR(a)) return null;

            // p - 1 = q * 2^s
            BigInteger q = Gf.P - 1;
            int s = 0;
            while (q.IsEven) { q >>= 1; s++; }

            BigInteger z = NonResidue();
            BigInteger c = BigInteger.ModPow(z, q, Gf.P);
            BigInteger x = BigInteger.ModPow(a, (q + 1) / 2, Gf.P);
            BigInteger t = BigInteger.ModPow(a, q, Gf.P);
            int m = s;

            while (!t.IsOne)
            {
                // 找最小的 i (0 < i < m) 使 t^(2^i) == 1
                int i = 0;
                BigInteger tt = t;
                while (!tt.IsOne)
                {
                    tt = tt * tt % Gf.P;
                    i++;
                    if (i >= m) return null;    // 理论上不会发生
                }
                // b = c^(2^(m-i-1))
                BigInteger b = c;
                for (int k = 0; k < m - i - 1; k++) b = b * b % Gf.P;
                x = x * b % Gf.P;
                BigInteger b2 = b * b % Gf.P;
                t = t * b2 % Gf.P;
                c = b2;
                m = i;
            }

            BigInteger hi = Gf.P - x;
            return x <= hi ? x : hi;
        }
    }

    // -------------------------------------------------------------------------
    // 曲线 y^2 = x^3 + a*x + b over GF(p)
    // 公钥解析出来的 B 是 bEncryptArray（产品密钥的整数形式），
    // 这里 x = B mod p，T = lift_x(x)。
    // -------------------------------------------------------------------------
    public static class Curve
    {
        public static BigInteger Rhs(BigInteger x, BigInteger a, BigInteger b)
        {
            x = Gf.Mod(x);
            return Gf.Mod(x * x % Gf.P * x + a * x + b);
        }

        /// <summary>返回曲线上的两个点 (x, ±y)，无解时返回 null。</summary>
        public static BigInteger[] LiftX(BigInteger x, BigInteger a, BigInteger b)
        {
            x = Gf.Mod(x);
            BigInteger? y = Fp.Sqrt(Rhs(x, a, b));
            if (y == null) return null;
            BigInteger yv = y.Value;
            return new BigInteger[] { yv, Gf.Mod(-yv) };
        }

    }

    // -------------------------------------------------------------------------
    // E6 上的点（坐标在 K6 里）
    // -------------------------------------------------------------------------
    public struct Pt6
    {
        public Fp6 X, Y;
        public bool Inf;

        public static Pt6 Infinity { get { return new Pt6 { X = Fp6.Zero, Y = Fp6.Zero, Inf = true }; } }
        public static Pt6 From(Fp6 x, Fp6 y) { return new Pt6 { X = x, Y = y, Inf = false }; }

        /// <summary>倍点：λ = (3x² + a) / 2y</summary>
        public Pt6 Dbl(BigInteger a)
        {
            if (Inf || Y.Equals(Fp6.Zero)) return Infinity;
            Fp3 two = new Fp3(2, 0, 0);
            Fp6 num = (X * X) * Fp6.FromFp3(new Fp3(3, 0, 0)) + Fp6.FromFp3(new Fp3(a, 0, 0));
            Fp6 den = Y * Fp6.FromFp3(two);
            Fp6 lam = num / den;
            Fp6 x3 = lam * lam - X - X;
            Fp6 y3 = lam * (X - x3) - Y;
            return From(x3, y3);
        }

        public static Pt6 Add(Pt6 P, Pt6 Q)
        {
            if (P.Inf) return Q;
            if (Q.Inf) return P;
            if (P.X.Equals(Q.X))
                return P.Y.Equals(Q.Y) ? throw new InvalidOperationException("Add: 需要倍点，请调用 Dbl")
                                       : Infinity;   // 互为逆元
            Fp6 lam = (Q.Y - P.Y) / (Q.X - P.X);
            Fp6 x3 = lam * lam - P.X - Q.X;
            Fp6 y3 = lam * (P.X - x3) - P.Y;
            return From(x3, y3);
        }

        public Pt6 Mul(BigInteger k, BigInteger a)
        {
            Pt6 r = Infinity, b = this;
            while (k > 0)
            {
                if (!(k & 1).IsZero) r = r.Inf ? b : Add(r, b);
                b = b.Dbl(a);
                k >>= 1;
            }
            return r;
        }
    }

    // -------------------------------------------------------------------------
    // Tate 配对  e(T, Q) = f_{n,T}(Q) ^ ((p^6 - 1) / n)
    //
    // Miller 循环里把累加点 V 保持在 **E(K)** 上（坐标是普通 BigInteger，运算快），
    // 只把直线函数在 **Q ∈ E6** 上求值。直线与垂直线都取完整形式（不做分母消除）——
    // 分母能否消除取决于群阶与子域的关系，先不赌，两种都算、用结果对答案。
    // -------------------------------------------------------------------------
    //  ★ 分母消除（denominator elimination）—— 本实现里最要紧的一处优化。
    //
    //  Miller 递推是  f_{m+m'} = f_m · f_{m'} · l / v，其中 v 是垂直线。
    //  在 Q 上求值时 v(Q) = Q.X − x([m]T)。本例中
    //      Q.X ∈ Fp3  （扭曲线映射 x' = x·t⁻² 把 K3 坐标嵌进 K6 的实部）
    //      x([m]T) ∈ GF(p)  （累加点 V 始终在 E(K) 上）
    //  所以所有分母都在真子域 Fp3* 里。而最终幂指数
    //      E = (p⁶−1)/n = (p³−1)·(p³+1)/n      （n | p²−p+1，p³+1 = (p+1)(p²−p+1)）
    //  对任意 d ∈ Fp3* 给出  d^E = (d^(p³−1))^((p³+1)/n) = 1。
    //  ⇒ 分母被最终幂整个吃掉，可以**直接不除**。
    //
    //  这一步的收益是数量级的：原先每处理一位都要做一次 Fp6 除法，而 Fp6.Inverse()
    //  是 Pow(p⁶−2)，约 1000 次 Fp6 乘法 —— 110 位的循环就是 11 万次乘法，
    //  而循环本身只需约 550 次。也就是说**求逆占了配对耗时的 99% 以上**。
    public static class TatePairing
    {
        public static int BitLength(BigInteger v)
        {
            int n = 0;
            while (v > 0) { v >>= 1; n++; }
            return n;
        }

        static Fp6 Emb(BigInteger v) { return Fp6.FromFp3(new Fp3(v, 0, 0)); }

        /// <summary>直线 l_{P1,P2}(Q)，P1/P2 的坐标在 K 里（BigInteger），Q 在 E6 里</summary>
        static Fp6 Line(BigInteger x1, BigInteger y1, Fp6 slopeNumFp6, Fp6 slopeDenFp6, Pt6 Q)
        {
            Fp6 lam = slopeNumFp6 / slopeDenFp6;
            return (Q.Y - Emb(y1)) - lam * (Q.X - Emb(x1));
        }

        public static Fp6 Pair(BigInteger n, BigInteger Tx, BigInteger Ty, Pt6 Q, BigInteger a)
        {
            Fp6 f = Fp6.One;
            BigInteger vx = Tx, vy = Ty;
            bool vInf = false;
            int L = BitLength(n);

            for (int i = L - 2; i >= 0; i--)
            {
                // ---- 倍点 V -> 2V ----
                if (!vInf)
                {
                    if (vy == 0) { vInf = true; }                       // 2-挠点
                    else
                    {
                        BigInteger lam = Gf.Mod(Gf.Mod(3 * vx * vx + a) * Fpm.InvBig(Gf.Mod(2 * vy)));
                        Fp6 l = (Q.Y - Emb(vy)) - Emb(lam) * (Q.X - Emb(vx));
                        BigInteger x2 = Gf.Mod(lam * lam - 2 * vx);
                        BigInteger y2 = Gf.Mod(lam * (vx - x2) - vy);
                        f = (f * f) * l;              // 分母消除，见上
                        vx = x2; vy = y2;
                    }
                }
                else f = f * f;                                         // V = ∞：直线平凡

                // ---- 按 n 的当前位决定是否加 T ----
                if (!((n >> i) & BigInteger.One).IsZero)
                {
                    if (vInf)
                    {
                        // ∞ + T = T，直线平凡
                        vx = Tx; vy = Ty; vInf = false;
                    }
                    else if (vx == Tx && vy == Ty)
                    {
                        // V == T：退化为倍点
                        BigInteger lam = Gf.Mod(Gf.Mod(3 * vx * vx + a) * Fpm.InvBig(Gf.Mod(2 * vy)));
                        Fp6 l = (Q.Y - Emb(vy)) - Emb(lam) * (Q.X - Emb(vx));
                        BigInteger x2 = Gf.Mod(lam * lam - 2 * vx);
                        BigInteger y2 = Gf.Mod(lam * (vx - x2) - vy);
                        f = f * l;                    // 分母消除
                        vx = x2; vy = y2;
                    }
                    else if (vx == Tx)
                    {
                        // V == -T：直线是垂直线，V + T = ∞，除数取 1
                        // （这是**正常**情形 —— 循环最后一轮必然出现：
                        //   V = [(n-1)/2]T 倍点后 = [(n-1)]T = -T）
                        f = f * (Q.X - Emb(vx));
                        vInf = true;
                    }
                    else
                    {
                        BigInteger lam = Gf.Mod(Gf.Mod(vy - Ty) * Fpm.InvBig(Gf.Mod(vx - Tx)));
                        Fp6 l = (Q.Y - Emb(vy)) - Emb(lam) * (Q.X - Emb(vx));
                        BigInteger x3 = Gf.Mod(lam * lam - vx - Tx);
                        BigInteger y3 = Gf.Mod(lam * (vx - x3) - vy);
                        f = f * l;                    // 分母消除
                        vx = x3; vy = y3;
                    }
                }
            }

            // ---- 最终幂 ----
            BigInteger e = (BigInteger.Pow(Gf.P, 6) - 1) / n;
            // 最终幂占了单次配对约 3/4 的时间，而这里还在用 BigInteger 版 Fp6。
            // 接到 Montgomery 层后 14 次配对从 ~0.50s 降到 ~0.1s。
            return Fp6mOp.PowBigM(f, e);
        }

        public static BigInteger Order;   // 由调用方设置
    }

    // -------------------------------------------------------------------------
    // 公钥
    // -------------------------------------------------------------------------
    public sealed class PubKey
    {
        public int SizeModulus, SizeOrder, ExtDeg1, ExtDeg2;
        public byte[] H1Bases;
        public BigInteger Modulus, Order;
        public byte[] K3Minpoly, K6Minpoly;
        public BigInteger CurveA, CurveB;
        public List<Fp3[]> Points = new List<Fp3[]>();   // 每项 = [x, y]
        public Fp3[] PairingVal;                          // K6 = [c0, c1]
    }

    public static class PubKeyParser
    {
        public static PubKey Parse(byte[] data)
        {
            int pos = 0;
            Func<int, uint> U32 = off => (uint)(data[off] | (data[off + 1] << 8) |
                                                (data[off + 2] << 16) | (data[off + 3] << 24));
            Func<int, int, BigInteger> Big = (off, n) =>
            {
                byte[] b = new byte[n + 1];
                Array.Copy(data, off, b, 0, n);
                b[n] = 0;                       // 补 0 保证正数（小端）
                return new BigInteger(b);
            };

            if (U32(0) != 0x44556677u) throw new InvalidDataException("bad magic1");
            uint fieldDataSize = U32(8);
            if (U32(12) != 0x00112233u) throw new InvalidDataException("bad magic2");

            pos = 20;
            byte mustBe0 = data[pos];
            int sizeModulus = data[pos + 1];
            int sizeOrder = data[pos + 2];
            if (mustBe0 != 0) throw new InvalidDataException("bad field data");
            uint[] dataList = new uint[9];
            for (int i = 0; i < 9; i++)
                dataList[i] = U32(pos + 3 + 4 * i);
            int extDeg1 = (int)dataList[0];
            int extDeg2 = (int)dataList[1];
            pos = 20 + 3 + 36;

            PubKey k = new PubKey();
            k.SizeModulus = sizeModulus;
            k.SizeOrder = sizeOrder;
            k.ExtDeg1 = extDeg1;
            k.ExtDeg2 = extDeg2;
            pos += 1;                                                // 跳过 1 字节

            k.H1Bases = new byte[sizeModulus];
            Array.Copy(data, pos, k.H1Bases, 0, sizeModulus); pos += sizeModulus;
            k.Modulus = Big(pos, sizeModulus); pos += sizeModulus;
            // Fp3 的构造函数要用 Gf.P 取模，所以必须在这里就设好
            Gf.P = k.Modulus;
            k.Order = Big(pos, sizeOrder); pos += sizeOrder;
            k.K3Minpoly = new byte[extDeg1 + 1];
            Array.Copy(data, pos, k.K3Minpoly, 0, extDeg1 + 1); pos += extDeg1 + 1;
            k.K6Minpoly = new byte[extDeg2 + 1];
            Array.Copy(data, pos, k.K6Minpoly, 0, extDeg2 + 1); pos += extDeg2 + 1;
            pos += sizeModulus * 2;                                  // 跳过一段
            k.CurveA = Big(pos, sizeModulus); pos += sizeModulus;
            k.CurveB = Big(pos, sizeModulus); pos += sizeModulus;

            // 点与 pairing_val 从 fieldDataSize + 12 开始
            pos = (int)fieldDataSize + 12;
            for (int i = 0; i < sizeModulus; i++)
            {
                Fp3 x = new Fp3(Big(pos, sizeModulus), Big(pos + sizeModulus, sizeModulus),
                                Big(pos + 2 * sizeModulus, sizeModulus));
                pos += sizeModulus * extDeg1;
                Fp3 y = new Fp3(Big(pos, sizeModulus), Big(pos + sizeModulus, sizeModulus),
                                Big(pos + 2 * sizeModulus, sizeModulus));
                pos += sizeModulus * extDeg1;
                k.Points.Add(new Fp3[] { x, y });
            }
            Fp3 pv0 = new Fp3(Big(pos, sizeModulus), Big(pos + sizeModulus, sizeModulus),
                              Big(pos + 2 * sizeModulus, sizeModulus));
            pos += sizeModulus * extDeg1;
            Fp3 pv1 = new Fp3(Big(pos, sizeModulus), Big(pos + sizeModulus, sizeModulus),
                              Big(pos + 2 * sizeModulus, sizeModulus));
            pos += sizeModulus * extDeg1;
            k.PairingVal = new Fp3[] { pv0, pv1 };
            return k;
        }
    }

    // =========================================================================
    //  ExtractM 的 C# 复刻
    //
    //  C++ 侧真正参与的就两段：
    //    CalculateH1:  h1Coeffs[0] = keyByte;  h1Coeffs[1] = 1;
    //                  h1Coeffs[2] = keyByte mod (H1Bases[1]+1)   ← mp_extract_byte_digits 第 1 位
    //                  h1Coeffs[i4+1] = bufB[i4]  (i4 = dw(a1,40)+1 .. n-1)
    //    ExtractM:     extract_scalar_digits(h1Coeffs + 2, dataOff + 2, sizeModulus - 1, ...)
    //
    //  extract_scalar_digits 是个 Horner 循环（k 从 n-1 递减到 0）：
    //        acc = acc * (radix[k] + 1) + digit[k]
    //  其中
    //        digit[k] = h1Coeffs[2 + k]
    //        radix[k] = *(dataOff + 2 + k) = H1Bases[1 + k]
    //
    //  为什么是 H1Bases[1+k]：dataOff 指向公钥里那个被跳过的字节（偏移 59），
    //  所以 dataOff+2 = 偏移 61 = H1Bases[1]。这一条已用三组黄金 M 验证。
    // =========================================================================
    public static class PKeyCalc
    {
        /// <summary>
        /// ★ 单个候选公钥的完整判定：公钥 + 密钥 → h1Coeffs / uid。
        ///
        /// 调用方自己用 Parallel.ForEach 遍历公钥字典，groupId 由字典键取得，例如：
        ///
        ///     int groupId = 0; byte[] h1 = null, uid = null; object gate = new object();
        ///     Parallel.ForEach(ConfigData2005.PublicKeyPart2, (item, state) =&gt;
        ///     {
        ///         byte[] pk = ConfigData2005.PublicKeyPart1.Concat(item.Value).ToArray();
        ///         int seq; string cfg; byte[] h, u;
        ///         if (!PKeyCalc.TryPubKey(pk, bEncryptArray, out seq, out cfg, out h, out u)) return;
        ///         lock (gate) { if (groupId != 0) return; groupId = item.Key; h1 = h; uid = u; }
        ///         state.Stop();
        ///     });
        ///
        /// 说明：
        ///  · 内部**不持有任何跨调用状态**，可安全并发调用；不需要任何 Prepare/初始化，
        ///    Fpm.Init() 是线程安全的懒初始化，第一次用到时自动完成。
        ///  · 内层 MITM 固定单线程（外层已经按候选铺满了核），两层都并行会超订。
        ///  · 公钥无效、或这把密钥不属于这组公钥，返回 false。
        /// </summary>
        /// <param name="publicKey">完整公钥字节 = PublicKeyPart1 + 该组的 Part2</param>
        /// <param name="bEncryptArray">产品密钥的 16 字节整数形式</param>
        /// <param name="channelSeq">输出的渠道序号</param>
        /// <param name="actPkeyConfig">输出的 Base64 配置串</param>
        /// <param name="h1Coeffs">输出的 15 字节 h1 系数</param>
        /// <param name="uid">输出的 8 字节 uid</param>
        public static bool TryPubKey(byte[] publicKey, byte[] bEncryptArray,
                                     out int channelSeq, out string actPkeyConfig,
                                     out byte[] h1Coeffs, out byte[] uid)
        {
            channelSeq = 0; actPkeyConfig = null; h1Coeffs = null; uid = null;
            if (publicKey == null || publicKey.Length == 0) return false;
            if (bEncryptArray == null || bEncryptArray.Length != 16) return false;

            PubKey k;
            try { k = PubKeyParser.Parse(publicKey); }   // 会写静态 Gf.P；各组公钥的 p 相同，并发写同值是安全的
            catch (BaseFieldMismatchException) { throw; }   // ★ 基域不一致必须让调用方看到，不能吞
            catch { return false; }
            if (k == null) return false;

            byte[] h1, u;
            if (!TryCandidate(k, bEncryptArray, out h1, out u)) return false;

            h1Coeffs = h1; uid = u;
            byte[] cfg = new byte[12];
            Buffer.BlockCopy(u, 0, cfg, 0, 8);           // 尾部 4 字节保持 0
            actPkeyConfig = Convert.ToBase64String(cfg);
            int ev = u[0] | ((u[1] | ((u[2] | ((u[3] & 0x7F) << 8)) << 8)) << 8);
            channelSeq = ((ev >> 1) % 1000000) + (1000000 * ((ev >> 1) / 1000000));
            return true;
        }

        // ---- 分段计时（单线程测量用；并发调用时互相覆盖，只在 bench 里读）----
        public static double TLastPairMs, TLastMitmMs;

        /// <summary>候选判定主体：由密钥推出 T，做 14 次配对，再用 MITM 搜 h1Coeffs[4..14]。</summary>
        static bool TryCandidate(PubKey k, byte[] enc, out byte[] h1, out byte[] uid)
        {
            h1 = null; uid = null;

            byte[] eb = new byte[enc.Length + 1];          // 补 0 保证 BigInteger 为正
            Buffer.BlockCopy(enc, 0, eb, 0, enc.Length);
            BigInteger B = new BigInteger(eb);

            BigInteger xT = Gf.Mod(B);                     // 余数 = x(T)
            int keyByte = (int)((B / Gf.P) & 0xFF);        // 商 = h1Coeffs[0]

            BigInteger[] T = Curve.LiftX(xT, k.CurveA, k.CurveB);
            if (T == null) return false;

            // h1Coeffs[0..3]：keyByte、常量 1、keyByte 的两个混合进制位
            h1 = new byte[15];
            int b1 = k.H1Bases[1] + 1, b2 = k.H1Bases[2] + 1;
            h1[0] = (byte)keyByte;
            h1[1] = 1;
            h1[2] = (byte)(keyByte % b1);
            h1[3] = (byte)((keyByte / b1) % b2);

            // 14 个配对 g_i = e(T, Q_i)
            var swPair = System.Diagnostics.Stopwatch.StartNew();
            Fp3 inv2 = new Fp3(-2, 0, 0).Inverse();
            Fp3 inv2sq = inv2 * inv2;
            int np = k.Points.Count;
            Fp6[] g = new Fp6[np];
            for (int i = 0; i < np; i++)
            {
                Pt6 Qi = Pt6.From(Fp6.FromFp3(k.Points[i][0] * inv2),
                                  new Fp6(Fp3.Zero, k.Points[i][1] * inv2sq));
                g[i] = TatePairing.Pair(k.Order, xT, T[0], Qi, k.CurveA);
            }
            swPair.Stop();
            TLastPairMs = swPair.Elapsed.TotalMilliseconds;

            // MITM 搜 h1Coeffs[4..14]
            // h_k = h1Coeffs[1+k]，所以搜索范围 lo[0..2] 取的是 h1Coeffs[1..3]
            int[] lo = new int[np], hi = new int[np];
            lo[0] = h1[1]; lo[1] = h1[2]; lo[2] = h1[3];
            for (int i = 0; i < 3; i++) hi[i] = lo[i];
            for (int i = 3; i < np; i++) hi[i] = k.H1Bases[i];

            Fp6 pv = new Fp6(k.PairingVal[0], k.PairingVal[1]);
            string note;
            var swMitm = System.Diagnostics.Stopwatch.StartNew();
            int[] sol = H1Search.Solve(g, Fp6.One / pv, lo, hi, out note, 1);
            swMitm.Stop();
            TLastMitmMs = swMitm.Elapsed.TotalMilliseconds;
            if (sol == null) return false;                 // Solve 内部已同时查过 pv 与其逆元

            for (int i = 3; i < np; i++) h1[4 + i - 3] = (byte)sol[i];

            uid = ExtractMBytes(k, h1, 8);
            return true;
        }
        /// <summary>公钥里被 `pos += 1` 跳过的那 1 字节 = C++ 的 dataOff。</summary>
        public const int DataOff = 59;

        /// <summary>
        /// M = extract_scalar_digits(h1Coeffs+2, dataOff+2, sizeModulus-1)。
        /// 每一位必须严格小于自己的进制，否则 C++ 侧会置 error_code 6 并返回 0。
        /// </summary>
        public static BigInteger ExtractM(PubKey k, byte[] h1Coeffs)
        {
            int n = k.SizeModulus - 1;                  // v6 - 1
            BigInteger acc = 0;
            for (int i = n - 1; i >= 0; i--)
            {
                int radix = k.H1Bases[1 + i] + 1;
                int digit = h1Coeffs[2 + i];
                if (digit >= radix)
                    throw new InvalidDataException(
                        "h1Coeffs[" + (2 + i) + "] = " + digit + "  >= 进制 " + radix);
                acc = acc * radix + digit;
            }
            return acc;
        }

        /// <summary>同 ExtractM，但按 DLL 的输出约定写成 nBytes 字节（小端）。</summary>
        public static byte[] ExtractMBytes(PubKey k, byte[] h1Coeffs, int nBytes)
        {
            byte[] src = ExtractM(k, h1Coeffs).ToByteArray();   // 小端，可能多带一个 0x00
            byte[] r = new byte[nBytes];
            Array.Copy(src, r, Math.Min(src.Length, nBytes));
            return r;
        }
    }

    // =========================================================================
    //  sa_match_entry 的摘要内核
    //
    //  sa_checksum4(a1, a2, a3, a4)  —— PKeyKernelBackend.cpp:2298
    //  这是 sa_match_entry 整个闭包里**唯一的摘要函数**：
    //
    //      v4 = 0
    //      for i in 0..3:
    //          if (i < *a2) v4 ^= a1[i];      // 最多 4 个 dword（16 字节）
    //          a3[i] = (byte)v4;              // 每轮先写低字节
    //          v4 = 123456789 * (v4 >> 8);    // mod 2^32
    //
    //  实测输入是 4 个 limb（n=12 但只有前 4 个 dword 参与），最高 limb 很小
    //  （d3 最大 0x26D4），量级正好是 110 位的 p —— 即被摘要的是一个 GF(p) 元素。
    //
    //  排序键（两个独立位置算出，互为印证）：
    //      sa_count_sort2 :  256 * rec[0] + rec[1]        （记录前两字节）
    //      sa_match_block :  ((ck >> 8) & 0xFF) + ((ck & 0xFF) << 8)
    //  两者都等于 256*ck[0] + ck[1]，即摘要前两字节按大端拼成的 16 位值 → 65536 个桶。
    // =========================================================================
    public static class SaChecksum
    {
    }

    // =========================================================================
    //  加速层：Montgomery 形式的 GF(p) / GF(p^3) / GF(p^6)
    //
    //  动机：MITM 要算约 65 万次 Fp6 乘法 = 约 1700 万次模乘。用 BigInteger 时
    //  **每次取模都要分配堆对象**，GC 成了瓶颈 —— 实测把建表并行到 8 核只拿到
    //  约 2× 加速，正是这个原因。换成定长 limb + Montgomery 归约后零堆分配。
    //
    //  p < 2^110 ⇒ 4 个 32 位 limb；R = 2^128。
    //  只在 H1Search 内部使用，边界处与 BigInteger 版互转，不影响其它代码。
    // =========================================================================
    public struct Fpx
    {
        public ulong A0, A1, A2, A3;
        public static Fpx Make(ulong a0, ulong a1, ulong a2, ulong a3)
        { Fpx r; r.A0 = a0; r.A1 = a1; r.A2 = a2; r.A3 = a3; return r; }
        public bool Eq(in Fpx o)
        { return A0 == o.A0 && A1 == o.A1 && A2 == o.A2 && A3 == o.A3; }
    }

    public static class Fpm
    {
        const ulong M = 0xFFFFFFFFUL;
        public static ulong PB0_, PB1_, PB2_, PB3_;
        static ulong N0INV;
        public static Fpx One;
        static volatile bool _init;
        static readonly object _initGate = new object();

        /// <summary>
        /// 建立 p 的 Montgomery 常数。★ 线程安全 —— 任何线程第一次用到时自动完成，
        /// 你在外面不用做任何初始化，Parallel.ForEach 直接开跑即可。
        /// 原来是把 _init 先置 true 再计算，多线程下别的线程会立刻返回、
        /// 读到 PB0_..PB3_ 还是 0，算出来的全是垃圾。这里改成全部算完再发布。
        /// </summary>
        public static void Init()
        {
            if (_init) return;
            lock (_initGate)
            {
                if (_init) return;
                BigInteger p = Gf.P;
                ulong b0 = (ulong)(p & M), b1 = (ulong)((p >> 32) & M);
                ulong b2 = (ulong)((p >> 64) & M), b3 = (ulong)((p >> 96) & M);
                // N0INV = -p0^{-1} mod 2^32；牛顿迭代 5 次，每步精度翻倍
                ulong inv = 1;
                for (int i = 0; i < 5; i++) inv = (inv * (2 - b0 * inv)) & M;
                Fpx one = FromBig(1);
                BigInteger rinv = BigInteger.ModPow(BigInteger.One << 128, p - 2, p);
                PB0_ = b0; PB1_ = b1; PB2_ = b2; PB3_ = b3;
                N0INV = (0x100000000UL - inv) & M;
                One = one;
                _rinv = rinv;
                _init = true;                       // ★ 最后才发布
            }
        }

        public static Fpx FromBig(BigInteger x)
        {
            BigInteger v = x % Gf.P;
            if (v.Sign < 0) v += Gf.P;
            v = (v << 128) % Gf.P;                       // 转 Montgomery 形式
            return Fpx.Make((ulong)(v & M), (ulong)((v >> 32) & M),
                            (ulong)((v >> 64) & M), (ulong)((v >> 96) & M));
        }

        static BigInteger _rinv;    // R^{-1}，在 Init() 里算好

        public static BigInteger ToBig(in Fpx a)
        {
            BigInteger v = a.A0 | ((BigInteger)a.A1 << 32) | ((BigInteger)a.A2 << 64)
                         | ((BigInteger)a.A3 << 96);
            // R^{-1} 只跟 p 有关，已在 Init() 里算好 —— 原来这里是懒加载，
            // 多线程并发时会重复算，而且 _rinv 的写入不原子。
            Init();
            return v * _rinv % Gf.P;
        }

        // 一次 CIOS 迭代：t[0..4] += a * bi，再做一次 Montgomery 归约并右移 32 位
        /// <summary>
        /// 求逆。a 已是 Montgomery 形式 aR，而
        ///     (aR)^(p−2) = a^(p−2) · R^(p−2) = a^(−1) · R · R^(p−1) = a^(−1) R
        /// （R ∈ GF(p)* ⇒ R^(p−1) = 1），所以直接幂就得到「逆的 Montgomery 形式」。
        /// 110 位指数 ≈ 165 次 Fpm.Mul，约 18µs；而 BigInteger 版要 70µs。
        /// </summary>
        public static Fpx Inv(in Fpx a)
        {
            Fpx r = One, b = a;
            BigInteger e = Gf.P - 2;
            while (e > 0)
            {
                if (!(e & 1).IsZero) r = Mul(r, b);
                b = Mul(b, b);
                e >>= 1;
            }
            return r;
        }

        /// <summary>BigInteger 进出的求逆，供 TatePairing 这种还在用 BigInteger 的地方调用。</summary>
        public static BigInteger InvBig(BigInteger a)
        {
            Init();
            return ToBig(Inv(FromBig(a)));
        }

        public static Fpx Mul(in Fpx a, in Fpx b)
        {
            // ★ 这里刻意**不手工展开**。展开版同时有 5 个 ulong 累加器 + 4 个 limb
            //   ≈ 10 个活跃的 64 位值 = 需要 20 个寄存器：x64 有 16 个还凑合，
            //   x86 只有 8 个（esp/ebp 还占掉 2 个）必然疯狂 spill。而 .NET Framework
            //   的 x86 用的是老 JIT32，寄存器分配比 x64 的 RyuJIT 差得多 ——
            //   为 x64 展开的代码在 x86 上是负优化。循环版 IL 体积小得多，
            //   JIT32 才能好好分配寄存器。实测见 PKeyDecoder.cs 顶部说明。
            ulong a0 = a.A0, a1 = a.A1, a2 = a.A2, a3 = a.A3;
            ulong b0 = b.A0, b1 = b.A1, b2 = b.A2, b3 = b.A3;
            ulong t0 = 0, t1 = 0, t2 = 0, t3 = 0, t4 = 0;
            for (int i = 0; i < 4; i++)
            {
                ulong bi = i == 0 ? b0 : i == 1 ? b1 : i == 2 ? b2 : b3;
                ulong carry, q;
                // t += a * bi
                q = a0 * bi + t0; t0 = q & M; carry = q >> 32;
                q = a1 * bi + t1 + carry; t1 = q & M; carry = q >> 32;
                q = a2 * bi + t2 + carry; t2 = q & M; carry = q >> 32;
                q = a3 * bi + t3 + carry; t3 = q & M; carry = q >> 32;
                ulong s = t4 + carry; t4 = s & M; ulong t5 = s >> 32;
                // Montgomery 归约：低 32 位必归零，直接丢弃
                ulong m = (t0 * N0INV) & M;
                carry = 0;
                q = m * PB0_ + t0; carry = q >> 32;
                q = m * PB1_ + t1 + carry; t0 = q & M; carry = q >> 32;
                q = m * PB2_ + t2 + carry; t1 = q & M; carry = q >> 32;
                q = m * PB3_ + t3 + carry; t2 = q & M; carry = q >> 32;
                q = t4 + carry; t3 = q & M; carry = q >> 32;
                t4 = carry + t5;
            }
            // 结果 < 2p，最多减一次 p
            if (t4 != 0 || t3 > PB3_
                || (t3 == PB3_ && (t2 > PB2_
                    || (t2 == PB2_ && (t1 > PB1_ || (t1 == PB1_ && t0 >= PB0_))))))
            {
                long c = 0, qq;
                qq = (long)t0 - (long)PB0_; t0 = (ulong)qq & M; c = qq >> 32;
                qq = (long)t1 - (long)PB1_ + c; t1 = (ulong)qq & M; c = qq >> 32;
                qq = (long)t2 - (long)PB2_ + c; t2 = (ulong)qq & M; c = qq >> 32;
                qq = (long)t3 - (long)PB3_ + c; t3 = (ulong)qq & M;
            }
            return Fpx.Make(t0, t1, t2, t3);
        }
    }

    public struct Fp3m
    {
        public Fpx C0, C1, C2;
        public static Fp3m FromFp3(Fp3 a)
        { Fp3m r; r.C0 = Fpm.FromBig(a.C0); r.C1 = Fpm.FromBig(a.C1); r.C2 = Fpm.FromBig(a.C2); return r; }
        public bool Eq(in Fp3m o) { return C0.Eq(o.C0) && C1.Eq(o.C1) && C2.Eq(o.C2); }
    }

    public struct Fp6m
    {
        public Fp3m R, I;                       // 实部 + 虚部*y，y^2 = -2
        public static Fp6m FromFp6(Fp6 a)
        { Fp6m r; r.R = Fp3m.FromFp3(a.C0); r.I = Fp3m.FromFp3(a.C1); return r; }
        public bool Eq(in Fp6m o) { return R.Eq(o.R) && I.Eq(o.I); }
    }

    public static class Fp6mOp
    {
        // ---- 模加减：limb 表示要求每位 < 2^32，且域里加减必须环绕，不能裸做 ----
        const ulong M32 = 0xFFFFFFFFUL;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        static bool GeP(in Fpx a)                       // a >= p ?
        {
            if (a.A3 != Fpm.PB3_) return a.A3 > Fpm.PB3_;
            if (a.A2 != Fpm.PB2_) return a.A2 > Fpm.PB2_;
            if (a.A1 != Fpm.PB1_) return a.A1 > Fpm.PB1_;
            return a.A0 >= Fpm.PB0_;
        }
        static void SubP(ref Fpx a)
        {
            long c = 0, t;
            t = (long)a.A0 - (long)Fpm.PB0_; a.A0 = (ulong)(t & (long)M32); c = t >> 32;
            t = (long)a.A1 - (long)Fpm.PB1_ + c; a.A1 = (ulong)(t & (long)M32); c = t >> 32;
            t = (long)a.A2 - (long)Fpm.PB2_ + c; a.A2 = (ulong)(t & (long)M32); c = t >> 32;
            t = (long)a.A3 - (long)Fpm.PB3_ + c; a.A3 = (ulong)(t & (long)M32);
        }
        static void AddP(ref Fpx a)
        {
            ulong c = 0, t;
            t = a.A0 + Fpm.PB0_; a.A0 = t & M32; c = t >> 32;
            t = a.A1 + Fpm.PB1_ + c; a.A1 = t & M32; c = t >> 32;
            t = a.A2 + Fpm.PB2_ + c; a.A2 = t & M32; c = t >> 32;
            t = a.A3 + Fpm.PB3_ + c; a.A3 = t & M32;
        }

        // ★ 热路径上最频繁的两个操作（每个 Mul3 要跑 16 次）。原来它们内部调用
        //   GeP/SubP/AddP 三个独立方法，IL 体积把 JIT 的内联预算撑爆，每次都是
        //   真调用 + Fpx 结构体进出内存。把条件减/加手工内联成叶子方法后：
        //   实测 Mul3 1355 -> 1138 ns（-16%），Mul6 4364 -> 3898 ns（-11%）。
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        static void AddMod(ref Fpx r, in Fpx a, in Fpx b)
        {
            ulong c = 0, t;
            t = a.A0 + b.A0; r.A0 = t & M32; c = t >> 32;
            t = a.A1 + b.A1 + c; r.A1 = t & M32; c = t >> 32;
            t = a.A2 + b.A2 + c; r.A2 = t & M32; c = t >> 32;
            t = a.A3 + b.A3 + c; r.A3 = t & M32; c = t >> 32;   // ★ 别漏最高位进位
            if (c != 0
                || r.A3 > Fpm.PB3_
                || (r.A3 == Fpm.PB3_ && (r.A2 > Fpm.PB2_
                || (r.A2 == Fpm.PB2_ && (r.A1 > Fpm.PB1_
                || (r.A1 == Fpm.PB1_ && r.A0 >= Fpm.PB0_))))))
            {
                long k = 0, u;
                u = (long)r.A0 - (long)Fpm.PB0_; r.A0 = (ulong)(u & (long)M32); k = u >> 32;
                u = (long)r.A1 - (long)Fpm.PB1_ + k; r.A1 = (ulong)(u & (long)M32); k = u >> 32;
                u = (long)r.A2 - (long)Fpm.PB2_ + k; r.A2 = (ulong)(u & (long)M32); k = u >> 32;
                u = (long)r.A3 - (long)Fpm.PB3_ + k; r.A3 = (ulong)(u & (long)M32);
            }
        }

        /// <summary>(a - b) mod p</summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        static void SubMod(ref Fpx r, in Fpx a, in Fpx b)
        {
            long c = 0, t;
            t = (long)a.A0 - (long)b.A0; r.A0 = (ulong)(t & (long)M32); c = t >> 32;
            t = (long)a.A1 - (long)b.A1 + c; r.A1 = (ulong)(t & (long)M32); c = t >> 32;
            t = (long)a.A2 - (long)b.A2 + c; r.A2 = (ulong)(t & (long)M32); c = t >> 32;
            t = (long)a.A3 - (long)b.A3 + c; r.A3 = (ulong)(t & (long)M32); c = t >> 32;
            if (c != 0)                                   // 借位 ⇒ a < b，加 p 环绕
            {
                ulong k = 0, u;
                u = r.A0 + Fpm.PB0_; r.A0 = u & M32; k = u >> 32;
                u = r.A1 + Fpm.PB1_ + k; r.A1 = u & M32; k = u >> 32;
                u = r.A2 + Fpm.PB2_ + k; r.A2 = u & M32; k = u >> 32;
                u = r.A3 + Fpm.PB3_ + k; r.A3 = u & M32;
            }
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        static void AddMod(ref Fp3m r, in Fp3m a, in Fp3m b)
        { AddMod(ref r.C0, a.C0, b.C0); AddMod(ref r.C1, a.C1, b.C1); AddMod(ref r.C2, a.C2, b.C2); }
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        static void SubMod(ref Fp3m r, in Fp3m a, in Fp3m b)
        { SubMod(ref r.C0, a.C0, b.C0); SubMod(ref r.C1, a.C1, b.C1); SubMod(ref r.C2, a.C2, b.C2); }


        // ---- 常数乘法用「加倍」做，比走一次完整 Montgomery 乘法便宜约 5 倍 ----
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        static Fpx Dbl(in Fpx a) { Fpx r = a; AddMod(ref r, r, a); return r; }
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        static Fpx Quad(in Fpx a) { Fpx r = Dbl(a); AddMod(ref r, r, r); return r; }
        static Fp3m Dbl3(in Fp3m a)
        {
            Fp3m r = a;                       // 先定值，ref 传参要求 definite assignment
            AddMod(ref r.C0, a.C0, a.C0);
            AddMod(ref r.C1, a.C1, a.C1);
            AddMod(ref r.C2, a.C2, a.C2);
            return r;
        }

        /// <summary>K3 乘法，u^3 = -4 - u。用 Karatsuba：6 次 Fpm.Mul（朴素要 9 次）。</summary>
        public static Fp3m Mul3(in Fp3m a, in Fp3m b)
        {
            Fpx d0 = Fpm.Mul(a.C0, b.C0);
            Fpx d1 = Fpm.Mul(a.C1, b.C1);
            Fpx d2 = Fpm.Mul(a.C2, b.C2);

            Fpx s = a.C0; AddMod(ref s, s, a.C1);
            Fpx w = b.C0; AddMod(ref w, w, b.C1);
            Fpx d01 = Fpm.Mul(s, w); SubMod(ref d01, d01, d0); SubMod(ref d01, d01, d1);

            s = a.C0; AddMod(ref s, s, a.C2);
            w = b.C0; AddMod(ref w, w, b.C2);
            Fpx d02 = Fpm.Mul(s, w); SubMod(ref d02, d02, d0); SubMod(ref d02, d02, d2);

            s = a.C1; AddMod(ref s, s, a.C2);
            w = b.C1; AddMod(ref w, w, b.C2);
            Fpx d12 = Fpm.Mul(s, w); SubMod(ref d12, d12, d1); SubMod(ref d12, d12, d2);

            // 卷积：t0=d0  t1=d01  t2=d02+d1  t3=d12  t4=d2
            // u^3 = -4-u, u^4 = -4u-u²  ⇒  r0 = t0-4t3, r1 = t1-t3-4t4, r2 = t2-t4
            Fp3m r;
            r.C0 = d0; r.C1 = d01; r.C2 = d02;
            AddMod(ref r.C2, r.C2, d1);                         // t2 = d02 + d1
            Fpx m4a = Quad(d12);
            Fpx m4b = Quad(d2);
            SubMod(ref r.C0, d0, m4a);
            SubMod(ref r.C1, d01, d12);
            SubMod(ref r.C1, r.C1, m4b);
            SubMod(ref r.C2, r.C2, d2);
            return r;
        }

        /// <summary>K6 乘法，(a0+a1y)(b0+b1y) = (a0b0 - 2a1b1) + (a0b1+a1b0)y。
        /// 用 Karatsuba：3 次 Mul3（朴素要 4 次），且 -2a1b1 用加倍而非乘法。</summary>
        public static Fp6m Mul6(in Fp6m a, in Fp6m b)
        {
            Fp3m d0 = Mul3(a.R, b.R);
            Fp3m d1 = Mul3(a.I, b.I);
            Fp3m r0 = d0; SubMod(ref r0, r0, Dbl3(d1));         // d0 - 2 d1

            Fp3m s = a.R; AddMod(ref s, s, a.I);
            Fp3m w = b.R; AddMod(ref w, w, b.I);
            Fp3m r1 = Mul3(s, w);                                // (aR+aI)(bR+bI)
            SubMod(ref r1, r1, d0); SubMod(ref r1, r1, d1);      // = aR*bI + aI*bR

            Fp6m r; r.R = r0; r.I = r1;
            return r;
        }

        public static Fp6m One()
        {
            Fp6m r;
            Fpx z = Fpm.FromBig(0);
            r.R.C0 = Fpm.One; r.R.C1 = z; r.R.C2 = z;
            r.I.C0 = z; r.I.C1 = z; r.I.C2 = z;
            return r;
        }

        public static Fp6m Pow(in Fp6m a, int e)
        {
            Fp6m r = One(), b = a;
            while (e > 0) { if ((e & 1) != 0) r = Mul6(r, b); b = Mul6(b, b); e >>= 1; }
            return r;
        }

        /// <summary>大指数幂（配对的最终幂要用 658 位的指数）。</summary>
        public static Fp6m PowBig(in Fp6m a, BigInteger e)
        {
            Fp6m r = One(), b = a;
            while (e > 0)
            {
                if (!(e & 1).IsZero) r = Mul6(r, b);
                b = Mul6(b, b);
                e >>= 1;
            }
            return r;
        }

        public static Fp6 ToFp6(in Fp6m a)
        {
            return new Fp6(new Fp3(Fpm.ToBig(a.R.C0), Fpm.ToBig(a.R.C1), Fpm.ToBig(a.R.C2)),
                           new Fp3(Fpm.ToBig(a.I.C0), Fpm.ToBig(a.I.C1), Fpm.ToBig(a.I.C2)));
        }

        /// <summary>把 BigInteger 版 Fp6 提到 Montgomery 层做幂，再转回来。</summary>
        public static Fp6 PowBigM(Fp6 a, BigInteger e)
        {
            return ToFp6(PowBig(Fp6m.FromFp6(a), e));
        }

    }

    // =========================================================================
    //  h1Coeffs[4..14] 的 MITM 搜索
    //
    //  约束（检验 [12] 实证）：
    //      Π_{k=0..13} g_k^{h_k}  ==  pairing_val（或其逆元）
    //      g_k = e(T, Q_k)，Q_k = 公钥第 k 个点
    //
    //  索引 k 与 h1Coeffs 的对应：h_k = h1Coeffs[1 + k]（k = 0..13）
    //  搜索范围（由 compute_h1_core 的 buf / bufA 定）：
    //      k = 0,1,2  → 固定为 buf[k]（v42[k]=0，无自由度）
    //      k = 3..13  → [0, H1Bases[k+1]]
    //  空间 = 16·16·16·32 · 4·4·4·4·4·8·64 = 2^17 · 2^19 = 2^36
    //
    //  劈半与 C++ 的 ia_1 = ia_2 = 7 一致：
    //      左半 k = 0..6   → 2^17 个候选，枚举
    //      右半 k = 7..13  → 2^19 个候选，建表（C++ 里 v14 = 524288 就是这个数）
    //
    //  因为群阶 ~2^110 远大于 2^36，解唯一 ⇒ 用哪种哈希都一样，不必复刻
    //  DLL 内部（Montgomery 形式）的字节表示。
    // =========================================================================
    public static class H1Search
    {
        /// <summary>摘要：对 Fp6m 第一个 GF(p) 分量的 4 个 limb 跑 sa_checksum4。
        /// 用的就是 Montgomery 形式的 limb —— 一致性只要求两侧相同，不必等于规范形式。</summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static uint HashM(in Fp6m v)
        {
            // 展开 sa_checksum4 且**不做任何分配** —— 建表要跑 52 万次，
            // 每次 new 一个 uint[] 加一个 byte[] 就是 100 万次堆分配，GC 会吃掉大半时间。
            uint w0 = (uint)v.R.C0.A0, w1 = (uint)v.R.C0.A1,
                 w2 = (uint)v.R.C0.A2, w3 = (uint)v.R.C0.A3;
            uint v4 = 0, b0, b1, b2, b3;
            v4 ^= w0; b0 = (byte)v4; v4 = unchecked(123456789u * (v4 >> 8));
            v4 ^= w1; b1 = (byte)v4; v4 = unchecked(123456789u * (v4 >> 8));
            v4 ^= w2; b2 = (byte)v4; v4 = unchecked(123456789u * (v4 >> 8));
            v4 ^= w3; b3 = (byte)v4;
            return b0 | (b1 << 8) | (b2 << 16) | (b3 << 24);
        }

        static Fp6m ProdM(Fp6m[] g, int[] lo, int[] hi, int from, int to, long s, out int[] digs)
        {
            digs = new int[to - from + 1];
            Fp6m v = Fp6mOp.One();
            long t = s;
            for (int i = 0; i <= to - from; i++)
            {
                int radix = hi[from + i] - lo[from + i] + 1;
                int d = (int)(t % radix); t /= radix;
                digs[i] = d;
                v = Fp6mOp.Mul6(v, Fp6mOp.Pow(g[from + i], lo[from + i] + d));
            }
            return v;
        }

        /// <summary>返回 h[0..n-1]，失败返回 null。</summary>
        /// <param name="threads">0 = 自动（核数）；1 = 单线程。
        /// 外层已经在按候选并行时传 1 —— 两层并行会在 4 物理核上超订，实测慢 3 倍。</param>
        /// <summary>置 1 后所有进行中的 Solve 会尽快退出（供多候选并行时用：
        /// 某个候选命中后，其它候选不必把整个 MITM 跑完）。</summary>
        public static int Abort;

        // ---- 诊断计数器：只在 bench 里读。单候选跑的时候是单线程（nth=1），
        //      所以普通静态字段就够，不加锁。
        public static long ProbeCalls, ProbeHashHits;
        public static double ProbeProdMms;

        public static int[] Solve(Fp6[] g, Fp6 target, int[] lo, int[] hi, out string note,
                                  int threads = 0)
        {
            Abort = 0;
            int n = lo.Length;
            int nth = threads > 0 ? threads : Environment.ProcessorCount;
            Fpm.Init();

            // 批量求逆（Montgomery 技巧）：1 次求逆 + 3n 次乘法，代替 n 次求逆。
            // 这一小段仍在 BigInteger 侧做 —— 只 14 个，代价可忽略。
            var ginv = new Fp6[n];
            {
                var pref = new Fp6[n + 1];
                pref[0] = Fp6.One;
                for (int k = 0; k < n; k++) pref[k + 1] = pref[k] * g[k];
                Fp6 acc = pref[n].Inverse();
                for (int k = n - 1; k >= 0; k--)
                {
                    ginv[k] = acc * pref[k];
                    acc = acc * g[k];
                }
            }

            // 转 Montgomery 形式，之后整个 MITM 都在定长 limb 上跑（零堆分配）
            var gm = new Fp6m[n];
            var gim = new Fp6m[n];
            for (int k = 0; k < n; k++) { gm[k] = Fp6m.FromFp6(g[k]); gim[k] = Fp6m.FromFp6(ginv[k]); }
            Fp6m tgt = Fp6m.FromFp6(target);

            // ---------- 分裂点：按**数字值**而不是下标来分 ----------
            // 原来按下标分：建表侧 4*4*4*4*4*8*64 = 2^19、扫描侧 2^17，
            // 而 k=0,1,2 的进制是 1（固定值，没有自由度，白占位置）。
            // 重排成 16*16*16*64 = 2^18 与 32*4*4*4*4*4*8 = 2^18 —— 两边平衡。
            // 失败的候选要跑满整个 MITM，所以这一步直接省掉 20% 的工作量，
            // 而且省下来的大头是**不能并行**的建表那一半。
            int[] order = { 0, 1, 2, 3, 4, 5, 13, 6, 7, 8, 9, 10, 11, 12 };
            var loP = new int[n];
            var hiP = new int[n];
            // ★ 生成元也一起重排 —— 位置 i 对应的是 g[order[i]]，不是 g[i]。
            //   漏了这一步会一直算不出解（踩过）。
            var gmP = new Fp6m[n];
            var gimP = new Fp6m[n];
            for (int i = 0; i < n; i++)
            {
                loP[i] = lo[order[i]]; hiP[i] = hi[order[i]];
                gmP[i] = gm[order[i]]; gimP[i] = gim[order[i]];
            }

            // ---------- 建表侧（重排后的 0..6）----------
            const int R0 = 0, R1 = 6;
            int rn = R1 - R0 + 1;
            long rsize = 1;
            for (int k = R0; k <= R1; k++) rsize *= (hiP[k] - loP[k] + 1);
            // ★ 固定位（hi == lo，范围 1）没有自由度，但里程计每步都会对它走一次
            //   "回绕"分支，乘一次 rwrap = g^(hi-lo) = g^0 = **单位元** —— 纯浪费。
            //   建表侧 7 位里有 3 位是固定位 ⇒ 每步 4 次 Mul6 里 3 次白做。
            //   只遍历真正有自由度的位，每步恰好 1 次 Mul6。
            var _actR = new System.Collections.Generic.List<int>();
            for (int i = 0; i < rn; i++) if (hiP[R0 + i] > loP[R0 + i]) _actR.Add(i);
            int[] actR = _actR.ToArray();

            var rwrap = new Fp6m[rn];
            for (int i = 0; i < rn; i++) rwrap[i] = Fp6mOp.Pow(gimP[R0 + i], hiP[R0 + i] - loP[R0 + i]);

            var keys = new ulong[rsize];
            var swBuild = System.Diagnostics.Stopwatch.StartNew();
            // 并行建表：里程计本身是串行的（每步依赖上一步），但可以按**块**切开，
            // 每块用 ProdM 直接算出自己的起始乘积，再在块内跑里程计。块数 = 核数。
            {
                int nchunk = Math.Max(1, nth);
                long chunk = (rsize + nchunk - 1) / nchunk;
                System.Threading.Tasks.Parallel.For(0, nchunk, c =>
                {
                    long begin = (long)c * chunk;
                    long end = Math.Min(begin + chunk, rsize);
                    if (begin >= end) return;
                    int[] d0;
                    Fp6m prod = ProdM(gmP, loP, hiP, R0, R1, begin, out d0);
                    var dg = (int[])d0.Clone();
                    for (long s = begin; s < end; s++)
                    {
                        if ((s & 0x3FF) == 0 && System.Threading.Volatile.Read(ref Abort) != 0) return;
                        keys[s] = ((ulong)HashM(prod) << 32) | (ulong)s;
                        for (int _j = 0; _j < actR.Length; _j++)
                        {
                            int i = actR[_j];
                            if (dg[i] < hiP[R0 + i] - loP[R0 + i])
                            { dg[i]++; prod = Fp6mOp.Mul6(prod, gmP[R0 + i]); break; }
                            dg[i] = 0; prod = Fp6mOp.Mul6(prod, rwrap[i]);
                        }
                    }
                });
            }
            Array.Sort(keys);

            // ---------- 扫描侧（重排后的 7..13）----------
            const int L0 = 7, L1 = 13;
            int ln = L1 - L0 + 1;
            long lsize = 1;
            for (int k = L0; k <= L1; k++) lsize *= (hiP[k] - loP[k] + 1);
            var _actL = new System.Collections.Generic.List<int>();
            for (int i = 0; i < ln; i++) if (hiP[L0 + i] > loP[L0 + i]) _actL.Add(i);
            int[] actL = _actL.ToArray();

            var lwrap = new Fp6m[ln];
            var lwrapi = new Fp6m[ln];
            for (int i = 0; i < ln; i++)
            {
                // 回绕：该位从 (hi-lo) 归零，乘积要**除掉** g^(hi-lo) ⇒ 乘它的逆
                lwrap[i] = Fp6mOp.Pow(gimP[L0 + i], hiP[L0 + i] - loP[L0 + i]);
                lwrapi[i] = Fp6mOp.Pow(gmP[L0 + i], hiP[L0 + i] - loP[L0 + i]);
            }

            var ldig = new int[ln];
            Fp6m lprod = Fp6mOp.One(), linv = Fp6mOp.One();
            for (int i = 0; i < ln; i++)
            {
                lprod = Fp6mOp.Mul6(lprod, Fp6mOp.Pow(gmP[L0 + i], loP[L0 + i]));
                linv = Fp6mOp.Mul6(linv, Fp6mOp.Pow(gimP[L0 + i], loP[L0 + i]));
            }

            long tBuild = swBuild.ElapsedMilliseconds;
            var swScan = System.Diagnostics.Stopwatch.StartNew();
            // 真实约束是 Π g_k^{h_k} == pairing_val，但**正负号取决于 lift_x 取了哪个 y 根**
            // （三组样本里 J6999 落在逆元、另两组直接相等）。两个目标放在同一趟里查 ——
            // 失败后重跑整个 Solve 会白扫一遍 131072 次，那才是真正的时间大头。
            Fp6m tgtI = Fp6m.FromFp6(Fp6.One / target);
            int[] h = new int[n];
            int hits = 0;
            string hitNote = null;
            int[] hitSol = null;

            // 只读查找：命中就返回解，否则 null。无共享写入 ⇒ 可并行。
            int[] Probe(in Fp6m needed, long s)
            {
                ProbeCalls++;
                uint hh = HashM(needed);
                long pp = LowerBound(keys, hh);
                while (pp < keys.Length && (uint)(keys[pp] >> 32) == hh)
                {
                    ProbeHashHits++;
                    long rs = (long)(uint)keys[pp];
                    int[] rd;
                    var swp = System.Diagnostics.Stopwatch.StartNew();
                    Fp6m rv = ProdM(gmP, loP, hiP, R0, R1, rs, out rd);
                    swp.Stop();
                    ProbeProdMms += swp.Elapsed.TotalMilliseconds;
                    if (rv.Eq(needed))
                    {
                        int[] ld;
                        ProdM(gm, loP, hiP, L0, L1, s, out ld);
                        var res = new int[n];
                        for (int i = 0; i < n; i++) res[order[i]] = loP[i];
                        for (int i = 0; i < ln; i++) res[order[L0 + i]] = loP[L0 + i] + ld[i];
                        for (int i = 0; i < rn; i++) res[order[R0 + i]] = loP[R0 + i] + rd[i];
                        return res;
                    }
                    pp++;
                }
                return null;
            }

            // 扫描同样按块并行（和建表一个套路）：每块用 ProdM(gim, ...) 算出自己的
            // 起始逆乘积，块内跑里程计。解唯一，所以谁先命中都行。
            // ★ 这里只维护 linv —— needed = tgt * linv，lprod 根本用不到，
            //   原先每轮白做一次 Mul6。
            int[] stop = new int[1];
            int nchunkL = Math.Max(1, nth);
            long chunkL = (lsize + nchunkL - 1) / nchunkL;
            System.Threading.Tasks.Parallel.For(0, nchunkL, ci =>
            {
                long begin = (long)ci * chunkL;
                long end = Math.Min(begin + chunkL, lsize);
                if (begin >= end) return;
                int[] d0;
                Fp6m li0 = ProdM(gimP, loP, hiP, L0, L1, begin, out d0);
                var dg = (int[])d0.Clone();
                // ★ 原来每步是 3 次 Mul6：算 tgt*li、算 tgtI*li、再推进 li。
                //   但这两个探测值每步乘的都是同一个因子，各自增量维护即可 ——
                //   探测变成零成本，每步只剩 2 次 Mul6（两个累加器各推进一次）。
                Fp6m la = Fp6mOp.Mul6(tgt, li0);
                Fp6m lb = Fp6mOp.Mul6(tgtI, li0);
                for (long s = begin; s < end; s++)
                {
                    if ((s & 0x3FF) == 0 && System.Threading.Volatile.Read(ref Abort) != 0) return;
                    if (System.Threading.Volatile.Read(ref stop[0]) != 0) return;
                    int[] r = Probe(la, s);
                    if (r == null) r = Probe(lb, s);
                    if (r != null)
                    {
                        if (System.Threading.Interlocked.CompareExchange(ref stop[0], 1, 0) == 0)
                        {
                            hitSol = r;
                            hits++;
                            hitNote = "左半 #" + s + " 命中右半（并行块 " + ci + "）";
                        }
                        return;
                    }
                    for (int _j = 0; _j < actL.Length; _j++)
                    {
                        int i = actL[_j];
                        if (dg[i] < hiP[L0 + i] - loP[L0 + i])
                        { dg[i]++; la = Fp6mOp.Mul6(la, gimP[L0 + i]); lb = Fp6mOp.Mul6(lb, gimP[L0 + i]); break; }
                        dg[i] = 0;
                        la = Fp6mOp.Mul6(la, lwrapi[i]);
                        lb = Fp6mOp.Mul6(lb, lwrapi[i]);
                    }
                }
            });
            if (hitSol != null)
            {
                note = hitNote;
                return hitSol;
            }
            note = "未找到（命中 " + hits + " 次）";
            return null;
        }

        static long LowerBound(ulong[] a, uint hh)
        {
            ulong key = (ulong)hh << 32;
            int lo = 0, hi = a.Length;
            while (lo < hi) { int m = (lo + hi) >> 1; if (a[m] < key) lo = m + 1; else hi = m; }
            return lo;
        }
    }

    // =========================================================================
    //  加速层自检 —— 动 Montgomery 核心之前**必须**先跑它。
    //  上次改 2x64 位 limb 时因为没自检，只能靠「整个解码失败」判断，
    //  根本定位不到是哪一层错了。
    // =========================================================================
    public static class PKeySelfTest
    {
        /// <param name="bench">true 时额外跑一段吞吐量测量（只在调优时用，生产路径不要开）</param>
        public static bool Run(Action<string> log, BigInteger modulus, bool bench = false)
        {
            if (log == null) log = s => Console.WriteLine(s);
            if (modulus != BigInteger.Zero) Gf.P = modulus;
            Fpm.Init();
            var rnd = new Random(12345);
            BigInteger R()
            {
                byte[] b = new byte[15]; rnd.NextBytes(b); b[14] = 0;
                BigInteger v = new BigInteger(b) % Gf.P;
                return v.IsZero ? BigInteger.One : v;
            }
            const int N = 200;
            int ok1 = 0, ok2 = 0, ok3 = 0, ok4 = 0;
            BigInteger x0 = R(), y0 = R();

            for (int i = 0; i < N; i++)
            {
                BigInteger x = R(), y = R(), z = R();

                // 1) GF(p) 乘法
                if (Fpm.ToBig(Fpm.Mul(Fpm.FromBig(x), Fpm.FromBig(y))) == x * y % Gf.P) ok1++;
                // 2) GF(p) 求逆
                if (Gf.Mod(x * Gf.Inv(x)) == BigInteger.One) ok2++;

                // 3) K3 乘法
                Fp3 a3 = new Fp3(x, y, z), b3 = new Fp3(R(), R(), R());
                Fp3 w3 = a3 * b3;
                Fp3m g3 = Fp6mOp.Mul3(Fp3m.FromFp3(a3), Fp3m.FromFp3(b3));
                if (Fpm.ToBig(g3.C0) == w3.C0 && Fpm.ToBig(g3.C1) == w3.C1
                 && Fpm.ToBig(g3.C2) == w3.C2) ok3++;

                // 4) K6 乘法
                Fp6 a6 = new Fp6(a3, b3), b6 = new Fp6(new Fp3(R(), R(), R()), new Fp3(R(), R(), R()));
                Fp6 w6 = a6 * b6;
                Fp6m g6 = Fp6mOp.Mul6(Fp6m.FromFp6(a6), Fp6m.FromFp6(b6));
                if (Fpm.ToBig(g6.R.C0) == w6.C0.C0 && Fpm.ToBig(g6.R.C1) == w6.C0.C1
                 && Fpm.ToBig(g6.R.C2) == w6.C0.C2 && Fpm.ToBig(g6.I.C0) == w6.C1.C0
                 && Fpm.ToBig(g6.I.C1) == w6.C1.C1 && Fpm.ToBig(g6.I.C2) == w6.C1.C2) ok4++;
            }
            if (bench)
            {
                // 吞吐量测量：只在调优时开。生产路径不要跑，白费时间。
                var swb = System.Diagnostics.Stopwatch.StartNew();
                Fpx accb = Fpm.FromBig(x0), bb = Fpm.FromBig(y0);
                const int MN = 400000;
                for (int i = 0; i < MN; i++) accb = Fpm.Mul(accb, bb);
                swb.Stop();
                log(string.Format("  Fpm.Mul {0,7:F1} ns/op", swb.Elapsed.TotalMilliseconds * 1e6 / MN));

                Fp6m a6b = Fp6m.FromFp6(new Fp6(new Fp3(x0, y0, x0), new Fp3(y0, x0, y0)));
                Fp6m b6b = Fp6m.FromFp6(new Fp6(new Fp3(y0, x0, y0), new Fp3(x0, y0, x0)));
                Fp6m acc6b = a6b;
                swb.Restart();
                for (int i = 0; i < MN / 12; i++) acc6b = Fp6mOp.Mul6(acc6b, b6b);
                swb.Stop();
                log(string.Format("  Mul6    {0,7:F1} ns/op", swb.Elapsed.TotalMilliseconds * 1e6 / (MN / 12)));
            }

            // 5) 基域守卫：重复写入同一个 p 必须放行；写入不同的 p 必须抛异常
            int ok5 = 0;
            try { Gf.P = Gf.P; ok5++; } catch { }
            try { Gf.P = Gf.P + 1; } catch (BaseFieldMismatchException) { ok5++; }

            bool all = ok1 == N && ok2 == N && ok3 == N && ok4 == N && ok5 == 2;
            log("[自检] Fpm.Mul " + ok1 + "/" + N + "   Gf.Inv " + ok2 + "/" + N
                + "   K3 " + ok3 + "/" + N + "   K6 " + ok4 + "/" + N
                + "   基域守卫 " + ok5 + "/2"
                + "   -> " + (all ? "PASS" : "*** FAIL ***"));
            return all;
        }
    }
}
