using namespace System
using namespace System.Diagnostics
using namespace System.Management.Automation
using namespace System.Runtime.InteropServices

param(
    [string]$BinkFile,
    [string]$KeyFile,
    [switch]$FromParent
)

if (-not $BinkFile -or -not $KeyFile) {
  throw "Error parse the files"
}
if ((!(Test-Path $BinkFile)) -or (!(Test-Path $KeyFile)) ) {
  throw "Error parse the files"
}
if ([Environment]::Is64BitProcess) {
    $ps32 = "$($env:Windir)\SysWOW64\WindowsPowerShell\v1.0\powershell.exe"
    $fwd = @('-ExecutionPolicy', 'Bypass', '-File', $PSCommandPath, '-FromParent')
    if ($BinkFile) { $fwd += @('-BinkFile', $BinkFile) }
    if ($KeyFile)  { $fwd += @('-KeyFile',  $KeyFile) }

    $raw  = & $ps32 @fwd
    try {return ($raw | ConvertFrom-Json)} catch {}
    return $raw
}

#region Misc
function Restore-Util {
    $asm = $null
    try { $asm = [Reflection.Assembly]::LoadWithPartialName('Microsoft.PowerShell.Commands.Utility') } catch { }
    if (-not $asm) {
        foreach ($ver in '3.0.0.0', '3.1.0.0') {
            try {
                $asm = [Reflection.Assembly]::Load(
                    "Microsoft.PowerShell.Commands.Utility, Version=$ver, Culture=neutral, PublicKeyToken=31bf3856ad364e35")
                if ($asm) { break }
            } catch { }
        }
    }

    if (-not $asm) {
        $dll = [IO.Path]::Combine($PSHOME, 'Microsoft.PowerShell.Commands.Utility.dll')
        if ([IO.File]::Exists($dll)) {
            try { $asm = [Reflection.Assembly]::LoadFile($dll) } catch { }
        }
    }

    if (-not $asm) {
        [Console]::Error.WriteLine('Could not locate Microsoft.PowerShell.Commands.Utility. The GAC copy may be missing too.')
        return
    }

    $token = ($asm.GetName().GetPublicKeyToken() | ForEach-Object { $_.ToString('x2') }) -join ''
    if ($token -ne '31bf3856ad364e35') {
        [Console]::Error.WriteLine("Refusing to import: unexpected publisher token '$token' (expected 31bf3856ad364e35).")
        return
    }

    Import-Module -Assembly $asm -ErrorAction Stop
}
function Connect-PidKeyData {
    param([string]$DllPath)

    if ([IntPtr]::Size -ne 4) { throw 'PidKeyData.dll is 32-bit: use 32-bit PowerShell.' }

    if (-not $DllPath) { $DllPath = [IO.Path]::Combine($PSScriptRoot, 'PidKeyData.dll') }
    if (-not [IO.Path]::IsPathRooted($DllPath)) { $DllPath = [IO.Path]::Combine($PSScriptRoot, $DllPath) }
    if (-not [IO.File]::Exists($DllPath)) { throw "DLL not found: $DllPath" }

    $dll = [Pk.K32]::LoadLibraryW($DllPath)
    if ($dll -eq [IntPtr]::Zero) { throw "LoadLibrary failed: $DllPath" }

    # executable memory for the stdcall -> fastcall stubs (16 bytes each)
    $mem = [Pk.K32]::VirtualAlloc([IntPtr]::Zero, [uint32]64, [uint32]0x3000, [uint32]0x40)
    if ($mem -eq [IntPtr]::Zero) { throw 'VirtualAlloc failed.' }

    $exports = @(
        @{ Name = 'PubkeyParser'; Type = [Pk.Fn4] }
        @{ Name = 'CalculateH1';  Type = [Pk.Fn6] }
        @{ Name = 'ExtractM';     Type = [Pk.Fn4] }
    )
    $fn = @{}
    for ($i = 0; $i -lt $exports.Count; $i++) {
        $target = [Pk.K32]::GetProcAddress($dll, $exports[$i].Name)
        if ($target -eq [IntPtr]::Zero) { throw "$($exports[$i].Name) is not exported." }

        # pop eax / pop ecx / pop edx / push eax / mov eax,target / jmp eax
        [byte[]]$stub = (0x58, 0x59, 0x5A, 0x50, 0xB8) + [BitConverter]::GetBytes($target.ToInt32()) + (0xFF, 0xE0)
        $addr = [IntPtr]::Add($mem, 16 * $i)
        [Marshal]::Copy($stub, 0, $addr, $stub.Length)
        $fn[$exports[$i].Name] = [Marshal]::GetDelegateForFunctionPointer($addr, $exports[$i].Type)
    }
    $script:Fn = $fn
}

Restore-Util
if (-not (Get-Command Add-Type -ErrorAction SilentlyContinue)) {
    throw 'Add-Type is still unavailable - Restore-Util could not import the Utility assembly.'
}
if (-not ([PSTypeName]'Pk.K32').Type) {
    Add-Type -TypeDefinition @'
using System;
using System.Runtime.InteropServices;
namespace Pk {
    public static class K32 {
        [DllImport("kernel32", CharSet = CharSet.Unicode, SetLastError = true)]
        public static extern IntPtr LoadLibraryW(string path);
        [DllImport("kernel32", CharSet = CharSet.Ansi, ExactSpelling = true, BestFitMapping = false, SetLastError = true)]
        public static extern IntPtr GetProcAddress(IntPtr module, string name);
        [DllImport("kernel32", SetLastError = true)]
        public static extern IntPtr VirtualAlloc(IntPtr address, uint size, uint type, uint protect);
    }
    // __fastcall exports are called via a stub, so the managed side is StdCall.
    [UnmanagedFunctionPointer(CallingConvention.StdCall)] public delegate int Fn4(IntPtr a, IntPtr b, IntPtr c, IntPtr d);
    [UnmanagedFunctionPointer(CallingConvention.StdCall)] public delegate int Fn6(IntPtr a, IntPtr b, IntPtr c, IntPtr d, IntPtr e, IntPtr f);
}
'@
}
#endregion

$Bink = [IO.File]::ReadAllBytes($BinkFile)
$Key  = [IO.File]::ReadAllBytes($KeyFile)

if ($Key.Length -ne 16) { throw 'Key must be 16 bytes.' }
if (-not $script:Fn) { Connect-PidKeyData }

# --- run ---------------------------------------------------------------------
# one allocation per C++ variable
$blob   = [Marshal]::AllocHGlobal($Bink.Length)   # BINK bytes (input)
$pMem   = [Marshal]::AllocHGlobal(32)             # intptr_t pMem[8]   (pMem[6] = struct ptr)
$pRet   = [Marshal]::AllocHGlobal(20)             # int      retValue[5]
$pValid = [Marshal]::AllocHGlobal(4)              # byte     isValid[4]
$pH1    = [Marshal]::AllocHGlobal(16)             # byte     h1Coeffs[15]
$pM     = [Marshal]::AllocHGlobal(8)              # byte     M[8]
$pKey   = [Marshal]::AllocHGlobal(16)             # 16-byte binary key (input)

$out = $null
try {
    [Marshal]::Copy($Bink, 0, $blob, $Bink.Length)
    [Marshal]::Copy([byte[]]::new(32), 0, $pMem,   32)
    [Marshal]::Copy([byte[]]::new(20), 0, $pRet,   20)
    [Marshal]::Copy([byte[]]::new(4),  0, $pValid, 4)
    [Marshal]::Copy([byte[]]::new(16), 0, $pH1,    16)
    [Marshal]::Copy([byte[]]::new(8),  0, $pM,     8)
    [Marshal]::Copy($Key, 0, $pKey, 16)

    if ($script:Fn.PubkeyParser.Invoke($pMem, $blob, [IntPtr]$Bink.Length, $pRet) -eq 0) {
        throw 'PubkeyParser rejected the BINK.'
    }
    $pData = [Marshal]::ReadIntPtr($pMem, 24)  # pMem[6]
    if ($pData -eq [IntPtr]::Zero) { throw 'PubkeyParser returned no data.' }
    $bytes1 = [IntPtr]::Add($pData, 44)
    $bytes2 = [IntPtr]::Add($pData, 88)

    [void]$script:Fn.CalculateH1.Invoke($bytes1, $bytes2, $pKey, $pValid, $pH1, $pRet)
    if ([Marshal]::ReadByte($pValid, 0) -ne 1) {
        $out = [pscustomobject]@{ Valid = $false; Upgrade = $null; Serial = $null; Security = $null }
    }
    else {
        [Marshal]::WriteInt32($pRet, 16, 1)  # retValue[4] = 1
        [void]$script:Fn.ExtractM.Invoke($bytes1, $pH1, $pM, $pRet)
        $value = [Marshal]::ReadInt64($pM, 0)

        $out = [pscustomobject]@{
            Valid    = $true
            Upgrade  = [bool]($value -band 1)
            Serial   = [uint32](($value -shr 1)  -band 0x3FFFFFFF)
            Security = [uint32](($value -shr 31) -band 0x3FF)
        }
    }
}
finally {
    foreach ($p in $blob, $pMem, $pRet, $pValid, $pH1, $pM, $pKey) {
        [Marshal]::FreeHGlobal($p)
    }
}

return $(if ($FromParent) {
  $out | ConvertTo-Json -Compress 
} else { 
  $out 
})