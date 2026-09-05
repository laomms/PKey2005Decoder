## PKey2005 Pairing-Based Cryptography
### How to use

```c
typedef int(__fastcall* PubkeyParserDelegate)(int* pDstMem, unsigned char* PublicKeyBytes, unsigned int dwSize, int* retValue);
typedef int(__fastcall* CalculateH1Delegate)(unsigned char* pMem1, unsigned char* pMem2, unsigned char* PID3Array, unsigned char* isValid, unsigned char* h1Coeffs, int* retValue);
typedef int(__fastcall* ExtractMDelegate)(unsigned char* pMem1, unsigned char* pMem2, unsigned char* M, int* retValue);

HMODULE hModule = LoadLibrary(L"PidKeyData.dll");
PubkeyParserDelegate pPubkeyParser = (PubkeyParserDelegate)GetProcAddress(hModule, "PubkeyParser");
CalculateH1Delegate	pCalculateH1 = (CalculateH1Delegate)GetProcAddress(hModule, "CalculateH1");
ExtractMDelegate pExtractM = (ExtractMDelegate)GetProcAddress(hModule, "ExtractM");

//Parser pkeyconfig data to pMem
int pMem[8] = { 0 };
int retValue[5] = { 0 };
int result = pPubkeyParser(pMem, bPublicKey, 0x62b, retValue);

//Calculate h1Coeffs from pid3 key array and and pkeyconfig data
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

//Extract M value from h1Coeffs
unsigned char M[8] = { 0 };
if (ifTrue[0] == 1) {
	retValue[4] = 1;	
	result = pExtractM(bytes1, h1Coeffs, M, retValue);
	
}
```
