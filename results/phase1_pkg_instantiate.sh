#!/usr/bin/env bash
# Phase 1 gate: instantiate Julia environment (WSL), log to results/.
set -x
cd "/mnt/d/MUBs in 6-dimension" || exit 1
export PATH=/home/adian/julia-1.12.6/bin:$PATH
julia --project=. -e 'using Pkg; Pkg.instantiate(); println("INSTANTIATE_OK")'
