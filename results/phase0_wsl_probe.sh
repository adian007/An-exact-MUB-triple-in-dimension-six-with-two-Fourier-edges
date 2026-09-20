#!/usr/bin/env bash
# Phase 0 environment probe for the MUB-in-dimension-6 audit.
export PATH=/home/adian/julia-1.12.6/bin:$PATH
echo "== julia =="
julia --version
echo "== julia depot =="
ls -d ~/.julia 2>/dev/null && ls ~/.julia/packages 2>/dev/null | head -30 || echo "no ~/.julia depot"
echo "== macaulay2 =="
M2 --version 2>&1 | head -3
echo "== singular =="
Singular -v 2>&1 | head -2
echo "== sage/mathematica/magma =="
which sage mathematica magma 2>/dev/null || echo "not installed"
