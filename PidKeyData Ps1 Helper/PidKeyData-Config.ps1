using namespace System
using namespace System.Diagnostics
using namespace System.Management.Automation
using namespace System.Runtime.InteropServices

param(
    [string]$CDKey,
    [string]$Config,
    [int]$Threads = 0,
    [switch]$FromParent
)

function Resolve-InputAbs([string]$p) {
    if (-not $p) { return $p }
    if ([IO.Path]::IsPathRooted($p)) { return $p }
    [IO.Path]::Combine($PWD.Path, $p)
}
if ([Environment]::Is64BitProcess) {
    $ps32 = "$($env:Windir)\SysWOW64\WindowsPowerShell\v1.0\powershell.exe"
    $fwd  = @('-ExecutionPolicy', 'Bypass', '-File', $PSCommandPath, '-FromParent')
    if ($CDKey)        { $fwd += @('-CDKey',   $CDKey) }
    if ($Config)       { $fwd += @('-Config',  (Resolve-InputAbs $Config)) }
    if ($Threads -gt 0){ $fwd += @('-Threads', $Threads) }

    $raw  = & $ps32 @fwd
    $text = ($raw -join [Environment]::NewLine)
    $a = $text.IndexOf('{'); $b = $text.LastIndexOf('}')
    if ($a -ge 0 -and $b -gt $a) { return ($text.Substring($a, $b - $a + 1) | ConvertFrom-Json) }
    return $text
}

#region Setup
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
        if ([IO.File]::Exists($dll)) { try { $asm = [Reflection.Assembly]::LoadFile($dll) } catch { } }
    }
    if (-not $asm) {
        [Console]::Error.WriteLine('Could not locate Microsoft.PowerShell.Commands.Utility.')
        return
    }
    $token = ($asm.GetName().GetPublicKeyToken() | ForEach-Object { $_.ToString('x2') }) -join ''
    if ($token -ne '31bf3856ad364e35') {
        [Console]::Error.WriteLine("Refusing to import: unexpected publisher token '$token'.")
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
    $fn
}

# 25-char CD-key -> 16-byte binary form.
function ConvertTo-KeyBytes([string]$ProductKey) {
    $alphabet = 'BCDFGHJKMPQRTVWXY2346789'
    $raw = $ProductKey.Replace('-', '').Trim().ToUpperInvariant()
    if ($raw.Length -ne 25) { throw 'Key must be 25 characters.' }

    $digits = [byte[]]::new(25); $isN = $false; $count = 0
    foreach ($ch in $raw.ToCharArray()) {
        if ($ch -eq 'N' -and -not $isN) {
            $isN = $true
            for ($i = $count; $i -gt 0; $i--) { $digits[$i] = $digits[$i - 1] }
            $digits[0] = [byte]$count; $count++; continue
        }
        $val = $alphabet.IndexOf($ch)
        if ($val -lt 0) { throw 'Invalid character in key.' }
        $digits[$count++] = [byte]$val
    }
    $binary = [byte[]]::new(16)
    foreach ($d in $digits) {
        $carry = [uint32]$d
        for ($i = 0; $i -lt 16; $i++) {
            $res = [uint32]$binary[$i] * 24 + $carry
            $binary[$i] = [byte]($res -band 0xFF); $carry = $res -shr 8
        }
    }
    if ($isN) { $binary[14] = $binary[14] -bor 0x08 }
    , $binary
}

# pkeyconfig.xrm-ms -> list of { GroupId; Bytes(1579) }.
function Get-PublicKeys([string]$configPath) {
    $outer = [IO.File]::ReadAllText($configPath)
    $i = $outer.IndexOf('pkeyConfigData'); if ($i -lt 0) { throw 'pkeyConfigData not found in config.' }
    $cs = $outer.IndexOf('>', $i) + 1
    $ce = $outer.IndexOf('</', $cs)
    $b64 = ($outer.Substring($cs, $ce - $cs) -replace '\s', '')
    $innerBytes = [Convert]::FromBase64String($b64)

    $xml = [System.Xml.XmlDocument]::new()
    $xml.Load([System.IO.MemoryStream]::new($innerBytes))   # auto-detects UTF-8/UTF-16

    $list = [System.Collections.Generic.List[object]]::new()
    foreach ($pk in $xml.GetElementsByTagName('PublicKey', '*')) {
        $kv = $pk.GetElementsByTagName('PublicKeyValue', '*')
        if ($kv.Count -eq 0) { continue }
        $gid = $pk.GetElementsByTagName('GroupId', '*')
        $groupId = if ($gid.Count) { $gid[0].InnerText } else { '' }
        $bytes = [Convert]::FromBase64String(($kv[0].InnerText -replace '\s', ''))
        if ($bytes.Length -eq 1579) { $list.Add([pscustomobject]@{ GroupId = $groupId; Bytes = $bytes }) }
    }
    $list
}

# Parallel search: first matching public key wins.
function Search-Keys($Fn, $entries, [byte[]]$key16, [int]$threads) {
    $total = $entries.Count
    if ($threads -lt 1) { $threads = 1 }
    if ($threads -gt $total) { $threads = $total }
    $sync = [hashtable]::Synchronized(@{ Found = $false; Result = $null })

    $worker = {
        param($Fn, $entries, $key16, $start, $end, $sync)
        $M = [Runtime.InteropServices.Marshal]
        $pMem = $M::AllocHGlobal(32); $pRet = $M::AllocHGlobal(20); $pValid = $M::AllocHGlobal(4)
        $pH1 = $M::AllocHGlobal(16); $pM = $M::AllocHGlobal(8); $pKey = $M::AllocHGlobal(16); $blob = $M::AllocHGlobal(1579)
        try {
            $M::Copy($key16, 0, $pKey, 16)
            for ($i = $start; $i -lt $end; $i++) {
                if ($sync.Found) { break }
                $e = $entries[$i]; $b = $e.Bytes
                $M::Copy([byte[]]::new(32), 0, $pMem, 32)
                $M::Copy([byte[]]::new(20), 0, $pRet, 20)
                $M::Copy([byte[]]::new(4), 0, $pValid, 4)
                $M::Copy($b, 0, $blob, $b.Length)

                if ($Fn.PubkeyParser.Invoke($pMem, $blob, [IntPtr]$b.Length, $pRet) -eq 0) { continue }
                $pData = $M::ReadIntPtr($pMem, 24)
                if ($pData -eq [IntPtr]::Zero) { continue }
                $bytes1 = [IntPtr]::Add($pData, 44); $bytes2 = [IntPtr]::Add($pData, 88)

                [void]$Fn.CalculateH1.Invoke($bytes1, $bytes2, $pKey, $pValid, $pH1, $pRet)
                if ($M::ReadByte($pValid, 0) -ne 1) { continue }

                $M::WriteInt32($pRet, 16, 1)
                [void]$Fn.ExtractM.Invoke($bytes1, $pH1, $pM, $pRet)
                $v = $M::ReadInt64($pM, 0)

                [System.Threading.Monitor]::Enter($sync)
                try {
                    if (-not $sync.Found) {
                        $sync.Found = $true
                        $sync.Result = [pscustomobject]@{
                            Valid    = $true
                            Upgrade  = [bool]($v -band 1)
                            Serial   = [uint32](($v -shr 1)  -band 0x3FFFFFFF)
                            Security = [uint32](($v -shr 31) -band 0x3FF)
                            GroupId  = $e.GroupId
                        }
                    }
                } finally { [System.Threading.Monitor]::Exit($sync) }
                break
            }
        }
        finally {
            foreach ($p in $pMem, $pRet, $pValid, $pH1, $pM, $pKey, $blob) { $M::FreeHGlobal($p) }
        }
    }

    $iss  = [System.Management.Automation.Runspaces.InitialSessionState]::CreateDefault2()
    $pool = [RunspaceFactory]::CreateRunspacePool(1, $threads, $iss, $Host)
    $pool.Open()
    try {
        $chunk = [Math]::Ceiling($total / $threads)
        $jobs = @()
        for ($t = 0; $t -lt $threads; $t++) {
            $s = $t * $chunk; $en = [Math]::Min($s + $chunk, $total)
            if ($s -ge $en) { continue }
            $ps = [PowerShell]::Create(); $ps.RunspacePool = $pool
            [void]$ps.AddScript($worker).
                AddArgument($Fn).AddArgument($entries).AddArgument($key16).
                AddArgument($s).AddArgument($en).AddArgument($sync)
            $jobs += [pscustomobject]@{ PS = $ps; Handle = $ps.BeginInvoke() }
        }
        foreach ($j in $jobs) { $j.PS.EndInvoke($j.Handle); $j.PS.Dispose() }
    }
    finally { $pool.Close(); $pool.Dispose() }

    if ($sync.Result) { $sync.Result }
    else { [pscustomobject]@{ Valid = $false; Upgrade = $null; Serial = $null; Security = $null; GroupId = $null } }
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

# --- run ---------------------------------------------------------------------
if (-not $CDKey)  { throw 'Provide -CDKey (25 characters).' }
if (-not $Config) { throw 'Provide -Config (path to pkeyconfig.xrm-ms).' }
if (-not [IO.Path]::IsPathRooted($Config)) { $Config = [IO.Path]::Combine($PWD.Path, $Config) }
if (-not [IO.File]::Exists($Config)) { throw "Config not found: $Config" }

$key16   = ConvertTo-KeyBytes $CDKey
$entries = Get-PublicKeys $Config
if ($entries.Count -eq 0) { throw 'No public keys found in the config.' }

$Fn = Connect-PidKeyData
if ($Threads -le 0) { $Threads = [Environment]::ProcessorCount }

$sw  = [System.Diagnostics.Stopwatch]::StartNew()
$out = Search-Keys $Fn $entries $key16 $Threads
$sw.Stop()

$out | Add-Member -NotePropertyName Groups  -NotePropertyValue $entries.Count -PassThru |
       Add-Member -NotePropertyName Elapsed -NotePropertyValue ([math]::Round($sw.Elapsed.TotalSeconds, 4)) -PassThru |
       Out-Null

return $(if ($FromParent) {
  $out | ConvertTo-Json -Compress 
} else { 
  $out 
})