#!/usr/bin/env bash
set -eu
ROOT='/mnt/d/MUBs in 6-dimension'
BASE="$ROOT/results/precision_f6_bounded_20260917"
mkdir -p "$BASE"
label="$1"
lam="$2"
seed="${3:-20260917}"
test ! -e "$BASE/$label.log"
export JULIA_DEPOT_PATH=/home/adian/mub-depot
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
(
  set +e
  timeout --signal=TERM --kill-after=5s 120s /home/adian/julia-1.12.6/bin/julia --startup-file=no --project=/home/adian/mub6 "$ROOT/scripts/julia/precision_f6_bounded.jl" "$lam" "$BASE/$label" "$seed"
  rc=$?
  echo "PROCESS_EXIT=$rc"
  echo "$rc" > "$BASE/$label.exit"
) > "$BASE/$label.log" 2>&1 &
echo "launched $label pid=$! (120s wall cap + 5s kill grace)"
