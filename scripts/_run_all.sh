#!/bin/bash
export JULIA_DEPOT_PATH=/home/adian/mub-depot
export PATH=$PATH:/home/adian/julia-1.12.6/bin
cd /home/adian/mub6

# Script 1: check_anchors.jl
echo "=== check_anchors.jl START ===" > /tmp/check_anchors.log
julia --project=. scripts/julia/check_anchors.jl >> /tmp/check_anchors.log 2>&1
echo "=== check_anchors.jl EXIT=$? ===" >> /tmp/check_anchors.log

# Script 2: verify_dita_construction.jl
echo "=== verify_dita_construction.jl START ===" > /tmp/verify_dita.log
julia --project=. scripts/julia/verify_dita_construction.jl >> /tmp/verify_dita.log 2>&1
echo "=== verify_dita_construction.jl EXIT=$? ===" >> /tmp/verify_dita.log

# Script 3: reconstruct_b3_algebraic.jl
echo "=== reconstruct_b3_algebraic.jl START ===" > /tmp/reconstruct_b3.log
julia --project=. scripts/julia/reconstruct_b3_algebraic.jl --use-cache --method lll --basis-degree 8 >> /tmp/reconstruct_b3.log 2>&1
echo "=== reconstruct_b3_algebraic.jl EXIT=$? ===" >> /tmp/reconstruct_b3.log

# Copy results to Windows FS
mkdir -p /mnt/d/MUBs\ in\ 6-dimension/results
cp -v /home/adian/mub6/results/verify_dita_construction.txt /mnt/d/MUBs\ in\ 6-dimension/results/ 2>&1 || echo "verify_dita: copy failed"
cp -v /home/adian/mub6/results/reconstruct_b3_algebraic.txt /mnt/d/MUBs\ in\ 6-dimension/results/ 2>&1 || echo "reconstruct_b3: copy failed"
cp -v /home/adian/mub6/results/reconstruct_b3_algebraic.meta.txt /mnt/d/MUBs\ in\ 6-dimension/results/ 2>&1 || echo "reconstruct_b3.meta: copy failed"
cp -v /home/adian/mub6/results/check_anchors.txt /mnt/d/MUBs\ in\ 6-dimension/results/ 2>&1 || true

echo "ALL_DONE" > /tmp/all_done.flag
