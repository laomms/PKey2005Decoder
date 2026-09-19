# PKey2005 Pairing-Based Cryptography (Curve-based Cryptosystem)

## 这个 DLL 是干什么的

`PidKeyData.dll` 是微软 PKEY2005 产品密钥验证算法的核心库。它负责验证一个 Windows 产品密钥的“数字签名”是否合法——具体手段是用 **Squared Tate 配对**（一种基于超奇异椭圆曲线的配对密码学）对产品密钥的数据做一次配对运算，得到校验值 `H1`，再和密钥里内嵌的签名字段比对。

它只导出 3 个函数，形成一条“公钥 → 校验上下文 → 签名 → 标量”的流水线：

```text
1579 字节公钥 + 产品密钥数字
        │  PubkeyParser（解析）
        ▼
124 字节 PubkeyData（曲线参数 + 配对/Montgomery 上下文）
        │  CalculateH1（配对运算）
        ▼
H1 系数（15 字节签名）
        │  ExtractM（提取标量）
        ▼
M 值（8 字节）
```

---

## 三个导出函数的算法过程

### ① `PubkeyParser(pMem, keyBytes, dwSize, retValue)` —— 解析公钥

**输入：** 1579 字节的原始公钥（BINK 格式，内含域参数 + 生成元点）。

**过程：**

1. `ParseFieldData` 解析公钥的三段字段（域大小、模数、以及两组“生成元”点 `p`/`q` 的坐标数组）。
2. `build_montgomery_ctx` 构建 Montgomery 归约上下文：算出 `n' = −modinv(modulus)`（模逆）和 `R = 2^128 mod N`，供后续所有模乘/模约减使用。
3. `build_pairing_ctx` 构建配对上下文：把曲线参数、域大小（limb 数）、点坐标打包成 44 字节的上下文块，并分配工作缓冲区。
4. 组装成 124 字节 `PubkeyData`。

**输出：**

- `pMem[6]` → 指向 124 字节 `PubkeyData`。
- `pMem[0..5]` = 域/点尺寸参数（本例 `73 29 0A 0F 06 14`）。

### ② `CalculateH1(b1, b2, enc, ok, h1Coeffs, retValue)` —— 计算 H1 签名

**输入：** 配对上下文 `b1`（44 字节）+ 公钥数据 `b2`（36 字节）+ 加密数据 `enc`（16 字节）+ 产品密钥的序列号数字。

**过程：** 这就是 Squared Tate 配对的完整流程，也是整个 DLL 最重、最核心的部分。

1. **差分/进位解析**：把序列号数字 `buf` 和 `a6` 合成系数向量 `a7`（`a7[k]=buf[k]+a6[k]`，同时存差值 `v42[k]=a6[k]−buf[k]`）。
2. **第一轮 Miller 环**（`mp_signed_block_loop`）：对 14 个元素逐个做窗口幂运算（`mid^exp`），底层是 Montgomery 乘法 `mont_mul_redc`。
3. **`sa_match_entry` 后缀数组匹配**：这是 PKEY2005 独有的“系数匹配”——把系数向量和公钥生成元点做匹配，经 `sa_descend` 递归下降得到带符号系数（决定配对里每个点的系数）。
4. **第二轮 Miller 环 + 合并回写**：再跑一次 `mp_signed_block_loop`，把两半结果合并回系数数组。
5. **输出 H1 系数**（15 字节）。

**输出：**

- `h1Coeffs = 1C 01 00 07 01 04 0F 09 01 01 00 02 01 00 13`（15 字节）。
- `ifTrue[0] = 1`（校验通过标志）。

**本质：** 把产品密钥的数字映射成椭圆曲线上的点，用 Squared Tate 配对 `e(P,Q)` 算出一个“哈希”，这个哈希就是 H1 签名。

### ③ `ExtractM(bytes1, h1Coeffs, resultBytes, retValue)` —— 提取 M 标量

**输入：** 配对上下文 `bytes1` + H1 系数 `h1Coeffs`（15 字节，`CalculateH1` 的输出）。

**过程：**

1. `mp_extract_byte_digits` 从 H1 系数里提取出“字节数字”（相当于把 15 字节签名拆成数字序列）。
2. `extract_scalar_digits` 用“逐字节 ×(digit+1) + 进位”的循环，把数字序列归约成一个标量——这就是 M 值。

**输出：**

- `resultBytes = 40 65 7E 6D AB 00 00 00`（8 字节的 M 标量）。

---

## 一张图看懂整体流水线

```text
产品密钥（25 位 base-24 数字）
   │
   │  序列号数字（buf / a6）
   ▼
┌─────────────────────────────────────────────┐
│ CalculateH1（Squared Tate 配对）             │
│  数字 → 系数向量 → Miller 环(窗口幂)          │
│  → sa_descend 匹配 → 带符号系数              │
│  → 再配对 → H1 系数(15B)                     │
└─────────────────────────────────────────────┘
   │
   ▼
H1 系数 = 1C 01 00 07 01 04 0F 09 01 01 00 02 01 00 13
   │
   │ ExtractM（字节数字 → 标量）
   ▼
M 值 = 40 65 7E 6D AB 00 00 00
```

---

## 核心数学/实现要点

| 层 | 函数族 | 作用 |
|---|---|---|
| 大数算术 | `mp_*`（已换 mini-gmp） | 多精度加/减/乘/除/比较/移位 |
| 域归约 | `mont_*` | Montgomery 模乘（`n'=−modinv`，CIOS） |
| 椭圆曲线 | `gf2m2_*`（实为派发胶水）+ `bn32_*` | 点运算/域元素批处理编排 |
| 窗口幂 | `mp_b_sub_exp_loop` / `mp_b_pow_chain` | 标量乘法的窗口法 |
| 系数匹配 | `sa_*` | PKEY2005 独有的后缀数组式匹配 |
| 配对内核 | `squared_tate_pairing` | Squared Tate 配对主循环 |

> **一句话总结：** 这个 DLL 用“椭圆曲线配对密码学”来验证产品密钥的签名——公钥先被解析成曲线上下文，产品密钥的数字被映射成点、经 Squared Tate 配对算出 H1 签名，最后从 H1 提取出 M 标量，作为整个验证流程的最终结果。

---

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

Console.WriteLine();
Console.WriteLine("groupId       = " + groupId);
Console.WriteLine("keyId         = " + keyId);
Console.WriteLine("actPkeyConfig = " + actPkeyConfig);
Console.WriteLine("h1Coeffs      = " + BitConverter.ToString(h1Out).Replace("-", ""));
Console.WriteLine("uid           = " + BitConverter.ToString(uidOut).Replace("-", ""));
Console.WriteLine("耗时          = " + sw.Elapsed.TotalSeconds.ToString("F2") + "s");
```

c++调用单个公钥大约0.7秒, C#版本release编译的计算时间是c++原版的两倍以上.极力推荐用c++版.
