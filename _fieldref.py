# -*- coding: utf-8 -*-
"""
GF(p) -> GF(p^3) -> GF(p^6) 参考实现 + 自洽性验证
用规范里给出的断言 (pairing_val^order == 1) 检验实现是否正确。

pubkey: Windows 7 Retail/GVLK, group 170
"""
import io

P = 886368969471450739924935101400677
ORDER = 886368969471450710152985728350703


def is_prime(n):
    if n < 2:
        return False
    for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % q == 0:
            return n == q
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


# ---------------- K3 = GF(p)[u]/(u^3 + u + 4) ----------------
# 元素 = [c0, c1, c2] 表示 c0 + c1*u + c2*u^2
K3_DEG = 3
# minpoly [4,1,0,1] -> u^3 = -(4 + u)  =>  系数表: u^3 = -4 - u
K3_RED = [(-4) % P, (-1) % P]          # u^3 用到的低次系数
K3_DEGREE_TERM = 3


def k3_zero():
    return [0, 0, 0]


def k3_one():
    return [1, 0, 0]


def k3_add(a, b):
    return [(a[i] + b[i]) % P for i in range(3)]


def k3_sub(a, b):
    return [(a[i] - b[i]) % P for i in range(3)]


def k3_mul(a, b):
    """两级乘法：先卷积到 5 次，再按 u^3 = -4 - u 归约"""
    t = [0] * 5
    for i in range(3):
        if a[i]:
            for j in range(3):
                t[i + j] = (t[i + j] + a[i] * b[j]) % P
    # 从高次往低次归约
    for k in range(4, 2, -1):
        c = t[k] % P
        if c:
            t[k] = 0
            # u^k = u^(k-3) * u^3 = u^(k-3) * (-4 - u)
            t[k - 3] = (t[k - 3] - 4 * c) % P
            t[k - 2] = (t[k - 2] - c) % P
    return t[:3]


def k3_pow(a, e):
    r = k3_one()
    while e:
        if e & 1:
            r = k3_mul(r, a)
        a = k3_mul(a, a)
        e >>= 1
    return r


def k3_eq(a, b):
    return a == b


def k3_inv(a):
    return k3_pow(a, P ** 3 - 2)


# ---------------- K6 = K3[y]/(y^2 + 2) ----------------
# 元素 = [d0, d1]，d0/d1 是 K3 元素，表示 d0 + d1*y
def k6_zero():
    return [k3_zero(), k3_zero()]


def k6_one():
    return [k3_one(), k3_zero()]


def k6_add(a, b):
    return [k3_add(a[0], b[0]), k3_add(a[1], b[1])]


def k6_mul(a, b):
    """(a0 + a1 y)(b0 + b1 y) = (a0b0 - 2 a1b1) + (a0b1 + a1b0) y   [y^2 = -2]"""
    a0, a1 = a
    b0, b1 = b
    c0 = k3_sub(k3_mul(a0, b0), k3_mul([2, 0, 0], k3_mul(a1, b1)))
    c1 = k3_add(k3_mul(a0, b1), k3_mul(a1, b0))
    return [c0, c1]


def k6_pow(a, e):
    r = k6_one()
    while e:
        if e & 1:
            r = k6_mul(r, a)
        a = k6_mul(a, a)
        e >>= 1
    return r


def k6_eq(a, b):
    return k3_eq(a[0], b[0]) and k3_eq(a[1], b[1])


# ---------------- 从 pubkey 取数据 ----------------
d = io.open('E:/Project/pubkey.bin', 'rb').read()
SM, SO, E1, E2 = 14, 14, 3, 2


def rd(buf, off, n):
    return int.from_bytes(buf[off:off + n] + bytes([0]), "little")


modulus = rd(d, 74, SM)
order = rd(d, 88, SO)
assert modulus == P, "modulus 与规范不一致"
assert order == ORDER, "order 与规范不一致"

pos = 307 + 12
pts = []
for _ in range(SM):
    x = [rd(d, pos + SM * i, SM) for i in range(E1)]
    pos += SM * E1
    y = [rd(d, pos + SM * i, SM) for i in range(E1)]
    pos += SM * E1
    pts.append((x, y))
pval = []
for _ in range(E2):
    pval.append([rd(d, pos + SM * i, SM) for i in range(E1)])
    pos += SM * E1

print("P 是素数:", is_prime(P))
print("modulus/order 与规范一致: True")
print()

# ---- 检验 1: pairing_val 落在 K6 里，且 pairing_val^order == 1 ----
pv = [pval[0], pval[1]]
print("pairing_val = K3(%s) + K3(%s)*y" % (pval[0], pval[1]))
r = k6_pow(pv, ORDER)
print("pairing_val^order == 1 ?", k6_eq(r, k6_one()))
if not k6_eq(r, k6_one()):
    print("  实得:", r)
print()

# ---- 检验 2: 曲线方程 y^2 = x^3 + a x + b（扭曲线映射后）----
a_curve = rd(d, 137, SM)
b_curve = rd(d, 151, SM)
print("curve.a = %d" % a_curve)
print("curve.b = %d" % b_curve)

inv2 = k3_inv([(-2) % P, 0, 0])            # t^2 = -2  =>  t^-2 = 1/(-2)
inv2sq = k3_mul(inv2, inv2)                 # t^-3 的 K3 部分 = inv(-2)^2


def k6_from_k3(x):
    return [x, k3_zero()]


ok_all = True
for idx, (px, py) in enumerate(pts):
    xp = k3_mul(px, inv2)                       # x' = x * t^-2   (K3 部分)
    yp = [k3_zero(), k3_mul(py, inv2sq)]        # y' = y * t^-3   (t 部分)
    lhs = k6_mul(yp, yp)
    x3 = k6_mul(k6_mul(k6_from_k3(xp), k6_from_k3(xp)), k6_from_k3(xp))
    rhs = k6_add(k6_add(x3, k6_from_k3(k3_mul(xp, [a_curve % P, 0, 0]))),
                 k6_from_k3([b_curve % P, 0, 0]))
    if not k6_eq(lhs, rhs):
        ok_all = False
        print("  points[%d] 不满足曲线方程" % idx)
        break

print("14 个点全部满足 y^2 = x^3 + a x + b ?", ok_all)
