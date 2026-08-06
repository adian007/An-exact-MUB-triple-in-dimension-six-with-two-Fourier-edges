# Run Macaulay2 on exported .m2 scripts with correct host mount.
# Usage (from project root):
#   .\scripts\docker\run_m2.ps1 pool_F6.m2
#   .\scripts\docker\run_m2.ps1 fourth_mub_w_elimination_theta0.m2

param(
    [Parameter(Mandatory = $true, Position = 0)]
    [string]$ScriptName
)

$ErrorActionPreference = "Stop"

$ProjectRoot = Resolve-Path (Join-Path $PSScriptRoot "..\..")
$ExportDir = Join-Path $ProjectRoot "symbolic_export"
$ScriptPath = Join-Path $ExportDir $ScriptName

if (-not (Test-Path $ScriptPath)) {
    Write-Error "Script not found: $ScriptPath"
}

$LogDir = Join-Path $ProjectRoot "results\track_c_elimination"
New-Item -ItemType Directory -Force -Path $LogDir | Out-Null
$LogFile = Join-Path $LogDir ("m2_" + [IO.Path]::GetFileNameWithoutExtension($ScriptName) + ".log")

Write-Host "Project root : $ProjectRoot"
Write-Host "M2 script    : $ScriptPath"
Write-Host "Log file     : $LogFile"

docker run --rm `
    -v "${ProjectRoot}:/work" `
    -w /work/symbolic_export `
    unlhcc/macaulay2:1.24.05 `
    M2 --script $ScriptName 2>&1 | Tee-Object -FilePath $LogFile

Write-Host "Done. Log: $LogFile"
