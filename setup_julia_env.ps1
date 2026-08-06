# One-time / per-session Julia environment for this project (Windows).
# The user env var JULIA_DEPOT_PATH may still point at a removed LearningHub tree;
# this script uses a project-local depot under .julia-depot/ unless you override.
#
# Usage (PowerShell, from project root):
#   . .\setup_julia_env.ps1
#   julia --project=. scripts/julia/validate_perH_pipeline.jl

$ErrorActionPreference = "Stop"

$JuliaBin = "C:\Users\adian\AppData\Local\Programs\Julia-1.12.6\bin"
if (-not (Test-Path "$JuliaBin\julia.exe")) {
    Write-Error "Julia not found at $JuliaBin — edit setup_julia_env.ps1 if installed elsewhere."
}

$ProjectDepot = Join-Path $PSScriptRoot ".julia-depot"
New-Item -ItemType Directory -Force -Path $ProjectDepot | Out-Null

# Prefer project-local depot (portable, not tied to LearningHub).
$env:JULIA_DEPOT_PATH = $ProjectDepot
$env:PATH = "$JuliaBin;" + $env:PATH

Write-Host "JULIA_DEPOT_PATH = $env:JULIA_DEPOT_PATH"
Write-Host "julia on PATH    = $JuliaBin"
& "$JuliaBin\julia.exe" -e "println(\"Julia \", VERSION); println(\"DEPOT_PATH = \", DEPOT_PATH)"
