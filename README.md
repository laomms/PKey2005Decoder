# PKey2005 Pairing-Based Cryptography

## 这个 DLL 是干什么的

`PidKeyData.dll` 是微软 PKEY2005 产品密钥验证算法的核心库。它负责验证一个 Windows 产品密钥的“数字签名”是否合法——具体手段是用 **Squared Tate 配对**（一种基于超奇异椭圆曲线的配对密码学）对产品密钥的数据做一次配对运算，得到校验值 `H1`，再和密钥里内嵌的签名字段比对。
# PidKeyData.dll — 算法分析

对 `PidKeyData.dll`（Windows / Office 产品密钥校验模块，2005 版配置体系）的逆向分析与等价重建。

本文只描述**校验**方向：给定一个产品密钥和一组公钥，如何算出校验所需的中间量。

---

## 目录

- [1. 概览](#1-概览)
- [2. 三个导出函数](#2-三个导出函数)
- [3. 数学基础](#3-数学基础)
- [4. 完整算法流程](#4-完整算法流程)
- [5. MITM 搜索](#5-mitm-搜索)
- [6. 常量与偏移](#6-常量与偏移)
- [7. 验证数据](#7-验证数据)
- [8. 走过的弯路（否定结论）](#8-走过的弯路否定结论)
- [9. 性能](#9-性能)

---

## 1. 概览

这套体系的核心是：**把产品密钥的整数形式映射到一条椭圆曲线上的点，与一组公钥做配对运算，反解出一串 15 字节的系数 `h1Coeffs`，再由它导出机器相关的标识 `M`。**

```
产品密钥字符串
    │  base24 解码
    ▼
bEncryptArray（16 字节小端整数 B）
    │
    ▼
┌──────────────────────────────────────────────┐
│ fCalculateH1：B ──► h1Coeffs[15]             │
│   · 一次除法同时得到商和余数                  │
│   · 余数 → 曲线点 T                          │
│   · 14 次 Tate 配对                          │
│   · MITM 搜索剩余 11 个系数                   │
└──────────────────────────────────────────────┘
    │
    ▼
┌──────────────────────────────────────────────┐
│ fExtractM：h1Coeffs ──► M（8 字节）          │
│   · 混合进制位提取                            │
└──────────────────────────────────────────────┘
```

公钥侧由 `fPubkeyParser` 解析，产出后续两步所需的全部参数。

---

## 2. 三个导出函数

### 2.1 `fPubkeyParser`

```c
int fPubkeyParser(void** pMem, const uint8_t* pubkey, uint32_t flags, int32_t* retVal);
```

解析一组公钥（约 1579 字节）到内存对象。校验头部魔数，再按内嵌的长度字段切出各个域元素。

**输出对象包含：**

| 字段 | 说明 |
|---|---|
| `Modulus` | 基域素数 `p`（110 位） |
| `Order` | 曲线阶 `n` |
| `H1Bases[14]` | 各位的进制上界（**格式常量，见 §6**） |
| `Points[14][2]` | 14 个 K3 上的点 `Q_k` |
| `PairingVal` | 目标配对值，供 MITM 比对 |
| `CurveA`, `CurveB` | 曲线参数 `a`, `b` |

> **注意：`p` 存放在静态存储中。** 同一批公钥必须共用同一个基域 —— 并发解析时若混入不同 `p` 的公钥，静态变量会被互相覆盖，算出的结果是错的且不会报错。等价实现中对此加了显式校验。

---

### 2.2 `fCalculateH1`

```c
int fCalculateH1(const uint8_t* bytes1, const uint8_t* bytes2,
                 const uint8_t* bEncryptArray,
                 uint8_t* ifTrue, uint8_t* h1Coeffs, int32_t* retVal);
```

整个算法的核心。输入产品密钥整数 `B` 与公钥参数，输出 15 字节 `h1Coeffs` 和一个"是否命中"标志 `ifTrue`。

算法细节见 [§4](#4-完整算法流程)。

---

### 2.3 `fExtractM`

```c
int fExtractM(const uint8_t* bytes1, const uint8_t* h1Coeffs,
              uint8_t* uid, int32_t* retVal);
```

把 `h1Coeffs` 的若干个系数按混合进制还原成一个整数 `M`：

```
acc = 0
for i = n-1 downto 0:
    radix = H1Bases[1+i] + 1
    acc   = acc * radix + h1Coeffs[2+i]
```

其中 `n = sizeModulus - 1`。即

```
M = Σ  h1Coeffs[2+i] · Π (H1Bases[1+j] + 1)
```

---

## 3. 数学基础

### 3.1 域塔

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
- **`p ≡ 1 (mod 4)`** —— 开平方必须用完整的 Tonelli–Shanks，不能用 `p ≡ 3 (mod 4)` 的快速公式

### 3.2 曲线与扭曲线映射

曲线方程 `y² = x³ + ax + b`。

公钥点 `Q_k` 定义在 **K3** 上。要参与 K6 上的配对，需经扭曲线映射：

```
x' = x · t⁻²
y' = y · t⁻³        其中 t² = −2
```

### 3.3 Tate 配对

`squared_tate_pairing` 计算 `e(T, Q) ∈ K6`，中间用到 Miller 循环。

**分母消除**：Miller 循环中的分母都落在真子域 `Fp3*` 里，而

```
E = (p⁶ − 1) / n = (p³ − 1)(p³ + 1) / n
```

对 `d ∈ Fp3*` 有 `d^E = 1`，所以**分母可以整体丢掉**，不影响最终幂的结果。实测这一步带来约 29× 加速。

**最终幂** `f^E` 占单次配对约 3/4 的时间 —— 它是整个算法里最值得优化的单点。

---

## 4. 完整算法流程

```
1.  PubKeyParser.Parse(pubkey)                        → 公钥对象
2.  B = 产品密钥整数（bEncryptArray，16 字节小端）
3.  keyByte = B / p                                   → h1Coeffs[0]
4.  x(T) = B mod p ;  T = lift_x(x(T))
       ★ 3 和 4 是【同一次除法】：mp_divrem_knuth 同时给出商和余数
5.  h1Coeffs[1] = 1
    h1Coeffs[2] = keyByte mod (H1Bases[1] + 1)
    h1Coeffs[3] = (keyByte / (H1Bases[1] + 1)) mod (H1Bases[2] + 1)
       ★ 即 keyByte 的 2 个混合进制位
6.  g_k = e(T, Q_k)          k = 0…13     —— 14 次 Tate 配对
7.  MITM 求 h[3..13]                      → h1Coeffs[4..14]     见 §5
8.  M = extract_scalar_digits(h1Coeffs + 2, dataOff + 2, sizeModulus − 1)
```

### 4.1 `buf` / `bufA` / `bufB`

`compute_h1_core` 里有三个长度为 `3n` 的缓冲：

```
buf[0]     = 1
buf[1..3]  = keyByte 的 2 个混合进制位
buf[3..]   = 0

bufA[k]    = H1Bases[k]              ★ k ≥ 3 时是【同下标】，不是 k+1
bufB[k]    = buf[k] + ((v42[k] + v35[k]) >> 1)
              其中 v42[k] = bufA[k] − buf[k]
                   v35[k] = 2h_k − hi_k − lo_k     （带符号的中点偏移）
```

`v35` 是每个搜索位的中点偏移量，`bufB` 因此给出该位的中心值。

---

## 5. MITM 搜索

### 5.1 约束

```
Π_{k=0..13} g_k^{h_k}  ==  pairing_val        或其逆元
```

取逆元与否，取决于 `lift_x` 选了哪个 `y` 根（`e(−T, Q) = e(T, Q)⁻¹`）。

搜索空间 `2^36`，远小于群阶 `2^110` ⇒ **解唯一**。因此哈希函数可以自选，不必复刻 DLL 内部的 Montgomery 字节表示。

### 5.2 劈半

DLL 用 `ia_1 = ia_2 = 7`：

- 左半 `k = 0..6`  —— 枚举，`2^17`
- 右半 `k = 7..13` —— 建表，`2^19`

> **这个劈半是不平衡的。** 建表与枚举各需 `B`、`S` 次乘法，在 `B·S = 2^36` 约束下 `B = S = 2^18` 才是最优（`2^19 = 524288` 次乘法）。DLL 的 `2^17 + 2^19 = 655360` 比理论最优多 25%。

### 5.3 哈希与摘要

```
sa_checksum4(a1, a2, a3):
    v4 = 0
    for i in 0..3:
        if i < *a2: v4 ^= a1[i]
        a3[i] = (uint8_t)v4
        v4 = 123456789u * (v4 >> 8)
```

只用前 4 个 limb（16 字节）。

- **桶键** = `256 * ck[0] + ck[1]`（65536 个桶）
- **记录格式** = `[4 字节摘要][v6 额外字节]`，步长 `a1[17]`

### 5.4 累加器

```
累加器 = Π_k g_k^{h_k}  ∈ K6（24 limb）
```

`mp_b_sub_exp_loop` 是 **Montgomery ladder 标量幂**；
`mp_signed_block_loop` 遍历数字数组，在 `mp_slot_limb_loop(a5, Src, Src, ...)` 那一行累加。

`mp_slot_limb_loop` 取 **vtable 槽 8**，是逐 limb 的域运算。

---

## 6. 常量与偏移

### `H1Bases` —— 格式常量

```
[1, 3, 8, 15, 15, 15, 31, 3, 3, 3, 3, 3, 7, 63]
```

**三组公钥完全相同**，因此它是格式常量而非公钥数据。

### `dataOff = 59`

公钥中偏移 **59** 的那个字节（`pos += 1; // 跳过 1 字节`）决定了数据区起点。

```
dataOff + 2 = 61 = H1Bases[1]
```

**这一条是解开 `M` 计算的关键。**

---

## 7. 验证数据

### 7.1 完整向量

| 密钥 | `keyByte` | `h1Coeffs[4..14]` | `M` |
|---|---|---|---|
| `J6999-VHTDQ-KYQWG-4YV9C-YDWHX` | `0x0A` | 2, 8, 11, 8, 1, 1, 0, 0, 3, 6, 41 | `0x178997B9E52` |

### 7.2 黄金向量

重建实现的输出必须逐字节匹配：

```
h1Coeffs = 1C 01 00 07 01 04 0F 09 01 01 00 02 01 00 13
M        = 40 65 7E 6D AB 00 00 00
```

对应 `33PXH` 那组公钥。

### 7.3 其他实测

| 密钥（整数形式，小端） | groupId | `h1Coeffs` |
|---|---|---|
| `259294F8BA2325B1D9E7E7B088880400` | 180 | `1A0102060B07091D03010002020034` |
| `070CA7C00CD065130FBEA27F16C40100` | 98 | `0A01020202080B0801010000030629` |

---

## 8. 走过的弯路（否定结论）

以下都**实测过、确认不成立**，记录在此以免重复尝试。

### 8.1 摘要输入不是点坐标

被 `sa_checksum4` 摘要的候选值曾经长得像点坐标，实测**全部不是**：

| 候选 | 命中 |
|---|---|
| GF(p) 的 x 坐标 | 6/12 |
| K3 的 x 坐标 | 3/6 |
| 扭曲线 x 坐标 | 3/6 |
| `X^order == 1` | 0/6 |
| 范数 1 | 0/6 |
| 相邻差 / 比 / 和 / 二阶差 | 全部互异 |

两次出现的"正好一半"是随机域元素的期望命中率。原因：指数是小数字，乘积是一般元素，当然不满阶。

### 8.2 性能上的三个错误假设

在一个纯 C# 等价实现上做优化时，以下假设**全部被测量否定**：

| 假设 | 实测 |
|---|---|
| MITM 中表项重建是瓶颈 | 占单候选 **0.1%** |
| 14 次 Tate 配对（BigInteger 版）是大头 | 占单候选 **6.8%** |
| MITM 中存在隐藏开销 | 推算的 Mul6 时间与实测吻合（差 12%） |

**真正的瓶颈只有一处：MITM 里的乘法次数与单次乘法成本。**

### 8.3 优化手法的效果

| 手法 | 结果 |
|---|---|
| 里程计跳过"固定位"（`hi == lo`，每步白乘一次单位元） | ✅ 建表 6.5 → 1.1 s |
| 扫描侧双累加器改为增量维护 | ✅ 每步 3 次乘法 → 2 次 |
| 热路径辅助运算手工内联 + `AggressiveInlining` | ✅ −10% |
| 大结构体参数按值传递 → `in` | ✅ −7% |
| **完全手工展开 CIOS** | ❌ **x86 上寄存器溢出，慢 2.5~3 倍** |
| 2×64 位 limb | ❌ .NET 无 `UInt128`，乘法次数相同还多了合并开销 |
| 条件减改无分支 select | ❌ 指令数翻倍，更慢且结果错误 |

> **教训**：同一份大数代码在 x64（RyuJIT）上适合展开，在 x86（JIT32）上必须用循环版 —— 两套 JIT 的寄存器分配能力差一个档次。

---

## 9. 性能

### 9.1 单候选耗时

| 实现 | 单候选 |
|---|---|
| 原版 DLL（x86 原生） | ~0.94 s |
| 本项目 C++ 等价实现（x64） | 0.82 ~ 0.87 s |
| 纯 C# 等价实现（x86 Release） | 约 1.9 × 原生 |

### 9.2 差距的构成

```
1.9 倍 = 1.2 倍（乘法次数 786432 vs 655360） × 1.6 倍（单次乘法：托管代码 vs 原生汇编）
```

- **1.2 倍**：`786432 = 2^18.5 × 2`，已证明是本 MITM 形式的**下界** —— 因需同时比对 `pv` 与 `pv⁻¹` 两个目标，每步至少要 2 次乘法；把常量折进建表侧代价相同
- **1.6 倍**：语言层面的代码生成差距

### 9.3 测量方法（重要）

这类算法的基准测量**极易被环境噪声误导**。本项目踩过的坑：

- 同一个二进制在同一台机器上跑出过 **39 s / 44 s / 52 s / 62 s / 82 s**
- Debug 配置比 Release 慢 **2~4 倍**，且所有内联优化失效
- `AnyCPU` 不带 `<Prefer32Bit>false</Prefer32Bit>` 时，.NET Framework 的 exe 默认以 **32 位**运行（`dumpbin -headers` 显示 `14C machine (x86)`）

**可靠做法**：交替 A/B（各 2~3 轮）、同一次运行内分段计时、只用 Release 配置。

---

## 许可与声明

本文档为**互操作性研究**目的，描述产品密钥的**校验**过程。

- 不包含、也不提供任何密钥生成方法
- 测试向量均来自公开可得的产品密钥
- 数学描述基于对二进制行为的观察与等价重建


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
