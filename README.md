# PKey2005 Pairing-Based Cryptography

## 这个 DLL 是干什么的

`PidKeyData.dll` 是微软 PKEY2005 产品密钥验证算法的核心库。它负责验证一个 Windows 产品密钥的“数字签名”是否合法——具体手段是用 **Squared Tate 配对**（一种基于超奇异椭圆曲线的配对密码学）对产品密钥的数据做一次配对运算，得到校验值 `H1`，再和密钥里内嵌的签名字段比对。
# PidKeyData.dll — 算法分析

> 对 `PidKeyData.dll`（Windows / Office 产品密钥校验模块，2005 版配置体系）的
> 逆向分析与等价重建。
>
> 本文只描述**校验**方向的计算过程：给定一个产品密钥和一组公钥，
> 如何算出校验所需的各个中间量。

---

## 目录

- [1. 概览](#1-概览)
- [2. 数学基础](#2-数学基础)
- [3. 算法流程](#3-算法流程)
- [4. MITM 搜索](#4-mitm-搜索)
- [5. 常量与偏移](#5-常量与偏移)
- [6. 性能](#6-性能)

---

## 1. 概览

产品密钥的整数形式被映射到一条椭圆曲线上的点，通过与一组公钥做配对运算，
反解出一串 15 字节的系数 `h1Coeffs`，再由它按混合进制还原出目标值 `M`。

```
                    产品密钥字符串
              2XVT3-BM23Y-XY67K-YTBF7-8MB87
                          │
                          │ base24 解码
                          ▼
              B = bEncryptArray（16 字节小端）
                          │
              ┌───────────┴───────────┐
              │                       │
          B ÷ p                   B mod p        ← ★ 同一次除法
              │                       │             （mp_divrem_knuth
              ▼                       ▼              同时给出商和余数）
          keyByte                   x(T)
              │                       │
              │                       │ lift_x（Tonelli–Shanks）
              │                       ▼
              │              T = (x, ±y)  ∈ K3
              │                       │
              │                       │ 14 次 Tate 配对  e(T, Q_k)
              │                       ▼
              │              g₀ … g₁₃ ∈ K6
              │                       │
              │                       │ MITM 搜索 2³⁶
              │                       ▼
              │              h₃ … h₁₃
              │                       │
              └───────────┬───────────┘
                          ▼
                 h1Coeffs[0 … 14]  （15 字节）
                          │
                          │ 混合进制提取（基数来自 H1Bases）
                          ▼
                    M （8 字节）
```

### 三个导出函数

| 函数 | 职责 |
|---|---|
| `fPubkeyParser` | 解析公钥 → 参数对象（`p`, `n`, `Q_k`, 目标配对值, 曲线参数） |
| `fCalculateH1` | `B` + 公钥参数 → `h1Coeffs[0..14]` |
| `fExtractM` | `h1Coeffs[2..14]` → `M` |

**公钥对象的内容：**

| 字段 | 说明 |
|---|---|
| `Modulus` | 基域素数 `p`（110 位） |
| `Order` | 曲线阶 `n` |
| `H1Bases[14]` | 各位的进制上界（格式常量，见 [§5](#5-常量与偏移)） |
| `Points[14][2]` | 14 个 K3 上的点 `Q_k` |
| `PairingVal` | 目标配对值，供 MITM 比对 |
| `CurveA` / `CurveB` | 曲线参数 `a` / `b` |

> **注意**：`p` 存放在静态存储中。同一批公钥必须共用同一个基域 ——
> 并发解析时若混入不同 `p` 的公钥，静态变量会被互相覆盖，
> 算出的结果是错的且不会报错。等价实现中对此加了显式校验。

---

## 2. 数学基础

### 2.1 域塔

```
        GF(p)
          │  u³ + u + 4
          ▼
      K3 = GF(p³)
          │  y² + 2
          ▼
      K6 = GF(p³)[y] = GF(p⁶)
```

- `p = 886368969471450739924935101400677`（110 位素数）
- **`p ≡ 1 (mod 4)`** —— 开平方必须用完整的 Tonelli–Shanks，
  不能用 `p ≡ 3 (mod 4)` 的快速公式

### 2.2 曲线与扭曲线映射

曲线方程：

```
y² = x³ + a·x + b
```

公钥点 `Q_k` 定义在 **K3** 上。要参与 K6 上的配对，需经扭曲线映射：

```
x' = x · t⁻²
y' = y · t⁻³          其中  t² = −2
```

### 2.3 Tate 配对

`squared_tate_pairing` 计算 `e(T, Q) ∈ K6`，内部走 Miller 循环。

**分母消除**：Miller 循环中的分母都落在真子域 `Fp3*` 里，而

```
E = (p⁶ − 1) / n = (p³ − 1)(p³ + 1) / n
```

对任意 `d ∈ Fp3*` 都有 `d^E = 1`，所以**分母可以整体丢掉**，
不影响最终幂的结果。

> 最终幂 `f^E` 占单次配对约 3/4 的时间，是整个配对运算里最值得优化的单点。

---

## 3. 算法流程

### 3.1 `h1Coeffs` 的结构

理解整个算法的关键，是先看清这 15 个字节分别从哪来：

```
 下标:   0     1     2     3   │   4     5     6    …    13    14
      ┌─────┬─────┬─────┬─────┬─────┬─────┬─────┬────┬─────┬─────┐
      │ 1C  │ 01  │ 00  │ 07  │ 01  │ 04  │ 0F  │ …  │ 00  │ 13  │
      └─────┴─────┴─────┴─────┴─────┴─────┴─────┴────┴─────┴─────┘
         │     │     └──┬──┘     └────────────┬────────────┘
         │     │        │                     │
         │     │   keyByte 的拆位          MITM 搜索得到
         │     │   基数 4 → 9              基数 16,16,16,32,
         │     │                            4,4,4,4,4,8,64
         │     │
         │     └── 常量 1
         │
         └── B ÷ p

      ┌──────────────────────────────────────────────────────┐
      │  ▼ 进入 M 的部分：[2] … [14]，共 13 位                │
      │    [0] keyByte 和 [1] 常量 1 都【不】参与              │
      └──────────────────────────────────────────────────────┘
```

| 下标 | 来源 | 含义 |
|---|---|---|
| `[0]` | 一次除法 | `keyByte = B / p` |
| `[1]` | 常量 | 固定为 `1` |
| `[2]` | `keyByte` 拆位 | 基数 4，`M` 的最低有效位 |
| `[3]` | `keyByte` 拆位 | 基数 9 |
| `[4] … [14]` | MITM | `h₃ … h₁₃`，基数见表 |

**两条关键的对应关系：**

**① `h1Coeffs[1 + k] = h_k`**

配对基元是 `g₀ … g₁₃`，对应指数 `h₀ … h₁₃`，存放在 `h1Coeffs[1] … h1Coeffs[14]`
—— **整体偏移 1**。这就是为什么搜索下标从 `k = 3` 开始，却写进 `h1Coeffs[4]`：

```
h1Coeffs[4 + i − 3] = sol[i]          // i = 3 … 13
```

**② `M` 只用 `h1Coeffs[2..14]`**

前两个字节不参与编码。详见 [§3.6](#36-m-的提取)。

---

### 3.2 第 0 位 —— 一次除法同时产生两个东西

```
B        = 产品密钥整数（16 字节小端）
keyByte  = B / p          ← 商，赋值给 h1Coeffs[0]
x(T)     = B mod p        ← 余数，作为曲线点的 x 坐标
```

> **这一步只做一次除法。** `mp_divrem_knuth` 同时输出商和余数，
> 而不是"先算 `B / p` 再算 `B % p`"。

`B` 是 128 位、`p` 是 110 位，所以 `keyByte` 只有 **1 个字节**。

由余数得到曲线点：

```
T = lift_x(x(T))          // 在 K3 上开平方求 y
```

`lift_x` 有两个根 `±y`，返回哪一个会影响后面 —— 见 [§4.1](#41-约束与解的唯一性)。

---

### 3.3 第 1 位 —— 常量

```
h1Coeffs[1] = 1
```

固定值，无自由度。作用是让 `h₀ = 1`，即 `g₀` 的指数恒为 1。

---

### 3.4 第 2、3 位 —— `keyByte` 的混合进制拆位

```
h1Coeffs[2] = keyByte mod (H1Bases[1] + 1)                        = keyByte mod 4
h1Coeffs[3] = (keyByte / (H1Bases[1] + 1)) mod (H1Bases[2] + 1)   = (keyByte / 4) mod 9
```

即把 `keyByte` 按基数 `4 → 9` 逐级拆开。

`keyByte` 是 8 位，而 `4 × 9 = 36 < 256`，所以**只够拆出 2 位**，
剩下的高位不参与 —— 这也是搜索从 `k = 3` 才开始的原因。

> 这一步和 [§3.6](#36-m-的提取) 的 `M` 提取用的是**同一套基数**。
> 也就是说：**`M` 的最低两位就是 `keyByte` 本身**（在混合进制下）。

---

### 3.5 第 4–14 位 —— MITM

剩余 11 个字节（`h₁ … h₃` 之外的 `h₃ … h₁₃`）没有闭式解，只能搜索：

```
求 h₃ … h₁₃  使得   Π_{k=0..13} g_k^{h_k}  ==  pairing_val^{±1}
```

其中 `g_k = e(T, Q_k)`。搜索范围：

| k | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 上界 | 15 | 15 | 15 | 31 | 3 | 3 | 3 | 3 | 3 | 7 | 63 |
| 基数 | 16 | 16 | 16 | 32 | 4 | 4 | 4 | 4 | 4 | 8 | 64 |

连乘 = `16·16·16·32·4·4·4·4·4·8·64 = 2³⁶`。详细算法见 [§4](#4-mitm-搜索)。

---

### 3.6 `M` 的提取

```
M = extract_scalar_digits(h1Coeffs + 2,  dataOff + 2,  sizeModulus − 1)
                         └─ 从下标 2 开始读  └─ 基数表起点  └─ 位数
```

**三个参数：**

| 参数 | 值 | 含义 |
|---|---|---|
| `h1Coeffs + 2` | — | 从 `h1Coeffs[2]` 开始读，跳过 `[0]` 和 `[1]` |
| `dataOff + 2` | `61` | 基数表在公钥中的起点（`dataOff = 59`，故 `61 = H1Bases[1]`） |
| `sizeModulus − 1` | `13` | 要提取的位数 |

**算法：**

```
n   = sizeModulus − 1                      // = 13
acc = 0
for i = n−1 downto 0:                      // 从高位到低位
    radix = H1Bases[1 + i] + 1
    digit = h1Coeffs[2 + i]
    acc   = acc * radix + digit
return acc
```

等价的数学写法：

```
M = Σᵢ h1Coeffs[2+i] · Πⱼ<ᵢ (H1Bases[1+j] + 1)        i = 0 … 12
```

**对应表：**

| `i` | 数字 | 基数 |
|---|---|---|
| 0 | `h1Coeffs[2]` | `H1Bases[1] + 1` = 4 |
| 1 | `h1Coeffs[3]` | `H1Bases[2] + 1` = 9 |
| 2 | `h1Coeffs[4]` | `H1Bases[3] + 1` = 16 |
| … | … | … |
| 12 | `h1Coeffs[14]` | `H1Bases[13] + 1` = 64 |

`h1Coeffs[2]` 是**最低位**，`h1Coeffs[14]` 是**最高位**。

> 各位的基数，就是它在 [§3.5](#35-第-4-14-位--mitm) 里的搜索上界加一 ——
> **搜索范围和 `M` 的编码共用同一套基数**，这是整个设计里最容易看漏的一环。

**提取过程示意：**

```
  数字            00   07   01   04   0F   09   01   01   00   02   01   00   13
  基数             4    9   16   16   16   32    4    4    4    4    4    8   64
                  └────┴────┴────┴────┴────┴────┴────┴────┴────┴────┴────┴────┘
                                 Horner 累乘（从右往左）

  acc = 0
  acc = acc × 64 + 19    = 19                    ← 最高位
  acc = acc ×  8 +  0    = 152
  acc = acc ×  4 +  1    = 609
  acc = acc ×  4 +  2    = 2438
  acc = acc ×  4 +  0    = 9752
  acc = acc ×  4 +  1    = 39009
  acc = acc ×  4 +  1    = 156037
  acc = acc × 32 +  9    = 4993193
  acc = acc × 16 + 15    = 79891103
  acc = acc × 16 +  4    = 1278257652
  acc = acc × 16 +  1    = 20452122433
  acc = acc ×  9 +  7    = 184069101904
  acc = acc ×  4 +  0    = 736276407616          ← 最低位
                         = 0xAB6D7E6540
```

---

### 3.7 `buf` / `bufA` / `bufB`

`compute_h1_core` 内部有三个长度为 `3n` 的缓冲：

```
buf[0]     = 1
buf[1..3]  = keyByte 的 2 个混合进制位
buf[3..]   = 0

bufA[k]    = H1Bases[k]              ★ k ≥ 3 时是【同下标】，不是 k + 1
bufB[k]    = buf[k] + ((v42[k] + v35[k]) >> 1)
              其中 v42[k] = bufA[k] − buf[k]
                   v35[k] = 2·h_k − hi_k − lo_k     （带符号的中点偏移）
```

`v35` 是每个搜索位的中点偏移量，`bufB` 因此给出该位的中心值。

---

## 4. MITM 搜索

### 4.1 约束与解的唯一性

```
Π_{k=0..13} g_k^{h_k}  ==  pairing_val        或其逆元
```

取逆元与否，取决于 `lift_x` 选了哪个 `y` 根 —— 因为 `e(−T, Q) = e(T, Q)⁻¹`，
错误符号会让 14 个 `g_k` 全部变成逆元，乘积随之取逆。

搜索空间 `2³⁶` 远小于群阶 `2¹¹⁰`，所以**解是唯一的**。
因此哈希函数可以自选，不必复刻原实现内部的 Montgomery 字节表示。

### 4.2 劈半

```
  k:     0      1       2      │    3    4    5    6    │   7    8    9   10   11   12   13
      ┌──────┬───────┬───────┬───────┬────┬────┬────┬────┬────┬────┬────┬────┬────┬────┐
  h_k │  =1  │keyByte│keyByte│   ?   │ ?  │ ?  │ ?  │ ?  │ ?  │ ?  │ ?  │ ?  │ ?  │ ?  │
      └──────┴───────┴───────┴───────┴────┴────┴────┴────┴────┴────┴────┴────┴────┴────┘
        └──── 固定，不搜 ────┘     └──────────── 2³⁶ 种组合，需要搜 ────────────┘

                                       │
                    ┌──────────────────┴──────────────────┐
                    │           劈成两半做 MITM            │
                    └──────────────────┬──────────────────┘
                                       │
                    ┌──────────────────┴──────────────────┐
                    │                                     │
              ┌─────▼──────┐                       ┌──────▼─────┐
              │  建表侧     │                       │  扫描侧     │
              │  k = 7…13  │                       │  k = 0…6   │
              │   2¹⁹       │                       │   2¹⁷      │
              └─────┬──────┘                       └──────┬─────┘
                    │                                     │
              枚举 ∏_R                              枚举 ∏_L
                    │                                     │
                    ▼                                     ▼
            hash(∏_R) ──► 排序成表              hash(tgt · ∏_L)
                    │                                     │
                    └──────────► 查表 ◄───────────────────┘
                                  │
                        哈希相同 → 重建表项 → 全量比对
                                  │
                                  ▼
                          命中 → 还原 h₃ … h₁₃
```

原实现取 `ia₁ = ia₂ = 7`：左半 `k = 0..6`（`2¹⁷`），右半 `k = 7..13`（`2¹⁹`）。

### 4.3 哈希与分桶

```
sa_checksum4(a1, a2, a3):
    v4 = 0
    for i in 0..3:
        if i < *a2:  v4 ^= a1[i]
        a3[i] = (uint8_t)v4
        v4 = 123456789u * (v4 >> 8)
```

只吃前 4 个 limb（16 字节）。

- **桶键** = `256 * ck[0] + ck[1]`，共 65536 个桶
- **记录格式** = `[4 字节摘要][v6 额外字节]`，步长 `a1[17]`

### 4.4 累加器

```
累加器 = Π_k g_k^{h_k}  ∈ K6（24 limb）
```

- `mp_b_sub_exp_loop` —— **Montgomery ladder 标量幂**
- `mp_signed_block_loop` —— 遍历数字数组，在
  `mp_slot_limb_loop(a5, Src, Src, ...)` 那一行累加
- `mp_slot_limb_loop` 取 **vtable 槽 8**，是逐 limb 的域运算

---

## 5. 常量与偏移

### `H1Bases` —— 格式常量

```
[1, 3, 8, 15, 15, 15, 31, 3, 3, 3, 3, 3, 7, 63]
```

**三组公钥完全相同**，说明它是格式常量而非公钥数据。

### `dataOff = 59`

公钥中偏移 **59** 的那个字节（`pos += 1;  // 跳过 1 字节`）决定了数据区起点。

```
dataOff + 2 = 61 = H1Bases[1]
```

**这一条是解开 `M` 计算的关键。**

---

## 6. 性能

### 6.1 单候选耗时

| 实现 | 语言 / 平台 | 单候选 |
|---|---|---|
| 原版 `PidKeyData.dll` | x86 原生 | ~0.94 s |
| 等价重建 | C++（x64） | 0.82 ~ 0.87 s |
| 等价重建 | C#（.NET Framework 4.7.2） | 约 1.9 × 原生 |

### 6.2 单候选内的开销分布

| 阶段 | 占比 |
|---|---|
| MITM（有限域乘法） | **93 %** |
| 14 次 Tate 配对 | 7 % |

**瓶颈几乎全部集中在 MITM 的乘法上。** 配对虽然看着复杂，
但由于分母消除 + 最终幂接入 Montgomery 层，实际开销占比很小。

### 6.3 与原生实现的差距

```
1.9 倍  =  1.2 倍（乘法次数）  ×  1.6 倍（单次乘法成本）
```

**1.2 倍 —— 乘法次数**

| | 次数 |
|---|---|
| 原实现 | `2¹⁷ + 2¹⁹ = 655360` |
| 优化后的等价实现 | `2¹⁸ + 2¹⁹ = 786432` |

> 建表与扫描的乘法次数都是 `B`、`S`，在 `B·S = 2³⁶` 约束下 `B = S = 2¹⁸`
> 才是最优。原实现的劈半是不平衡的。
>
> 而等价实现因为需要同时比对 `pv` 与 `pv⁻¹` 两个目标，扫描侧每步至少要 2 次乘法 ——
> 把常量折进建表侧代价相同，`786432` 已是该形式的**下界**。

**1.6 倍 —— 单次乘法成本**

托管代码与手写原生汇编在大数乘法上的固有差距。
完全展开、无分支化、更宽的 limb 等手段均已验证无法消除。

### 6.4 有效的优化点

等价实现相对朴素写法累计提速约 **5 倍**，主要来自：

| 手段 | 效果 |
|---|---|
| **分母消除**（Miller 循环分母落在真子域） | 约 29× |
| **里程计跳过固定位** —— 范围 1 的位每步只会乘一次单位元 | 建表 6.5 → 1.1 s |
| **扫描侧双累加器增量维护** —— 两个目标共用推进因子 | 每步 3 次乘法 → 2 次 |
| 热路径辅助运算手工内联 + `AggressiveInlining` | −10 % |
| 大结构体参数按值传递改为 `in` | −7 % |
| Montgomery 域表示（避免 BigInteger 堆分配） | 消除 GC 压力 |
| Tate 配对最终幂接入 Montgomery 层 | 14 次配对 0.50 → 0.10 s |

---

## 声明

本文档为**互操作性研究**目的，描述产品密钥的**校验**过程。

- 不包含、也不提供任何密钥生成方法
- 数学描述基于对二进制行为的观察与等价重建
- 与 Microsoft 无关联



## How to use
#### c++版调用
```c
typedef int(__fastcall* PubkeyParserDelegate)(int* pDstMem, unsigned char* PublicKeyBytes, unsigned int dwSize, int* retValue);
typedef int(__fastcall* CalculateH1Delegate)(unsigned char* pMem1, unsigned char* pMem2, unsigned char* PID3Array, unsigned char* isValid, unsigned char* h1Coeffs, int* retValue);
typedef int(__fastcall* ExtractMDelegate)(unsigned char* pMem1, unsigned char* pMem2, unsigned char* M, int* retValue);

HMODULE hModule = LoadLibrary(L"PidKeyData.dll");
PubkeyParserDelegate pPubkeyParser = (PubkeyParserDelegate)GetProcAddress(hModule, "PubkeyParser");
CalculateH1Delegate pCalculateH1 = (CalculateH1Delegate)GetProcAddress(hModule, "CalculateH1");
ExtractMDelegate pExtractM = (ExtractMDelegate)GetProcAddress(hModule, "ExtractM");

// Parser pkeyconfig data to pMem
int pMem[8] = { 0 };
int retValue[5] = { 0 };
int result = pPubkeyParser(pMem, bPublicKey, 0x62b, retValue);

// Calculate h1Coeffs from pid3 key array and pkeyconfig data
typedef struct _PubkeyData {
    unsigned char header[44];
    unsigned char bytes1[44];
    unsigned char bytes2[36];
    int end_marker;
} PubkeyData;

PubkeyData* pData = (PubkeyData*)pMem[6];
unsigned char* bytes1 = pData->bytes1;
unsigned char* bytes2 = pData->bytes2;
unsigned char ifTrue[4] = { 0 };
unsigned char h1Coeffs[15] = { 0 };
result = pCalculateH1(bytes1, bytes2, KeyArray, ifTrue, h1Coeffs, retValue);

// Extract M value from h1Coeffs
unsigned char M[8] = { 0 };
if (ifTrue[0] == 1) {
    retValue[4] = 1;
    result = pExtractM(bytes1, h1Coeffs, M, retValue);
}
```

#### c#版调用
```c#
var sw = Stopwatch.StartNew();

int groupId = 0, keyId = 0;
string actPkeyConfig = null;
byte[] h1Out = null, uidOut = null;
object gate = new object();
int done = 0;

var popts = new ParallelOptions { MaxDegreeOfParallelism = dop };
Parallel.ForEach(ConfigData2005.PublicKeyPart2, popts, (item, state) =>
{
    byte[] pk = ConfigData2005.PublicKeyPart1.Concat(item.Value).ToArray();

    int seq; string cfg; byte[] h, u;
    if (!PKeyCalc.TryPubKey(pk, bEncryptArray, out seq, out cfg, out h, out u))
    {
        if (!quiet)
            Console.WriteLine("  ✗ groupId=" + item.Key + "   累计 "
                + sw.Elapsed.TotalSeconds.ToString("F1") + "s");
        return;
    }

    lock (gate)
    {
        if (groupId != 0) return;      
        groupId = item.Key; keyId = seq; actPkeyConfig = cfg;
        h1Out = h; uidOut = u;
    }
    Interlocked.Exchange(ref done, 1);
    if (!quiet) Console.WriteLine("命中 groupId=" + item.Key);
    if (Array.IndexOf(args, "all") < 0) state.Stop();
});

sw.Stop();

if (groupId == 0)
{
    Console.WriteLine("未找到匹配的公钥（耗时 " + sw.Elapsed.TotalSeconds.ToString("F2") + "s）");
    return 1;
}

long Secret = BitConverter.ToInt64(h1Out, 0);
int value1 = (uidOut[3] >> 7) | (uidOut[4] * 2);
int value2 = (uidOut[5] * 2 | uidOut[4] >> 7) & 3;
byte b1 = (byte)(value1 & 0xFF);
byte b2 = (byte)(value2 & 0xFF);
byte[] hashData = new byte[] { b1, b2, 0, 0, 0, 0, 0, 0 };

Console.WriteLine();
Console.WriteLine("groupId       = " + groupId);
Console.WriteLine("keyId         = " + keyId);
Console.WriteLine("actPkeyConfig = " + actPkeyConfig);
Console.WriteLine("h1Coeffs      = " + BitConverter.ToString(h1Out).Replace("-", ""));
Console.WriteLine("uniqueId      = " + BitConverter.ToString(uidOut).Replace("-", ""));
Console.WriteLine("hashData      = " + BitConverter.ToString(hashData).Replace("-", ""));
Console.WriteLine("耗时           = " + sw.Elapsed.TotalSeconds.ToString("F2") + "s");
```

c++调用单个公钥对象大约0.7秒, C#版本release编译的计算时间是c++原版的2倍以上.极力推荐用c++版.
