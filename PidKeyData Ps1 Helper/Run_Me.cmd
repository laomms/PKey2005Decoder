@echo off
@cls

cd /d "%~dp0"
powershell -nop -f .\PidKeyData-Direct.ps1 -BinkFile .\Bink.bin -KeyFile .\KeyData.bin
powershell -nop -f .\PidKeyData-Config.ps1 -CDKey RHTBY-VWY6D-QJRJ9-JGQ3X-Q2289 -Config .\pkeyconfig.xrm-ms