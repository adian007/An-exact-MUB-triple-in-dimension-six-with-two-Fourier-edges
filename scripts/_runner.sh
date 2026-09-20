#!/bin/bash
# Runner script for MUB project Julia verification scripts
set -e
export JULIA_DEPOT_PATH=/home/adian/mub-depot
export PATH=$PATH:/home/adian/julia-1.12.6/bin
cd /home/adian/mub6

echo "=== Environment ==="
echo "Julia: $(julia --version)"
echo "Project: $(cat Project.toml | head -5)"
echo ""

echo "=== Checking packages ==="
julia --project=. -e '
try
    using Nemo; println("Nemo: OK")
catch e; println("Nemo: MISSING"); end
try
    using LinearAlgebra; println("LinearAlgebra: OK")
catch e; println("LinearAlgebra: MISSING"); end
'
echo ""

echo "=== Running check_anchors.jl ==="
julia --project=. scripts/julia/check_anchors.jl 2>&1
echo ""

echo "=== Running verify_dita_construction.jl ==="
julia --project=. scripts/julia/verify_dita_construction.jl 2>&1
echo ""

echo "=== Running reconstruct_b3_algebraic.jl ==="
julia --project=. scripts/julia/reconstruct_b3_algebraic.jl --use-cache --method lll --basis-degree 8 2>&1
echo ""

echo "=== Results directory ==="
ls -la results/
echo ""

echo "=== Copying results to Windows FS ==="
mkdir -p "/mnt/d/MUBs in 6-dimension/results"
cp -v results/verify_dita_construction.txt "/mnt/d/MUBs in 6-dimension/results/" 2>&1 || echo "verify_dita_construction.txt not found"
cp -v results/reconstruct_b3_algebraic.txt "/mnt/d/MUBs in 6-dimension/results/" 2>&1 || echo "reconstruct_b3_algebraic.txt not found"
cp -v results/reconstruct_b3_algebraic.meta.txt "/mnt/d/MUBs in 6-dimension/results/" 2>&1 || echo "reconstruct_b3_algebraic.meta.txt not found"
cp -rv results/benchmarks "/mnt/d/MUBs in 6-dimension/results/" 2>&1 || echo "benchmarks not found"

echo ""
echo "=== Done ==="
