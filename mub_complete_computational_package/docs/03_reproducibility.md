# Reproducibility Guide

## Python symbolic environment

Recommended Windows setup:

    conda create -n mub-symbolic python=3.11 sympy numpy
    conda activate mub-symbolic

From the main repository root, run the canonical project implementation:

    python scripts/python/dita_i3.py
    python -m unittest discover -s test -p "test_dita_i3.py" -v

The original package path remains a compatibility wrapper:

    python mub_complete_computational_package/01_symbolic_I3/build_i3.py

Generated artifacts are written to `symbolic_export/I3/` by default. Pass
`--out-dir PATH` to preserve a separate comparison run.

If the project's existing environment already contains SymPy, use it rather than creating a duplicate environment.

## Julia numerical algebraic geometry

Install Julia and add HomotopyContinuation.jl in the Julia REPL:

    using Pkg
    Pkg.add("HomotopyContinuation")

The official documentation describes polynomial systems, parameter homotopies, mixed-volume starts, continuation, and certification: https://www.juliahomotopycontinuation.org/HomotopyContinuation.jl/stable/

For this project, record Julia version, package versions, random seeds, parameter values, tolerances, and whether certification was used. This is especially important for the 252-path witness computations.

## Exact CAS

The project may use Macaulay2, Mathematica, SageMath, or SymPy depending on the subproblem. Exact coefficient fields must be written explicitly. Never silently replace an algebraic coefficient by a floating-point approximation when claiming an exact ideal computation.

## Reproduction checklist

1. Verify the exact H_D(z) matrix.
2. Regenerate `symbolic_export/I3/I3_single_vector.txt`.
3. Check the single-vector equations and the branch counts (97 generators,
   62 variables before reduction).
4. Run the automated exact conjugation, parameter-shift, and fixed-z
   generator identity tests.
5. Re-run the fixed-z monomial search at the documented generic sample.
6. Verify every discovered monomial automorphism symbolically.
7. Recompute the A4 closure/order.
8. Convert pool cliques to the exact `H_D(z)` gauge with a verified
   equivalence map, then run `a4_clique_orbits` before interpreting singular
   roots. The orbit helper intentionally rejects missing images.
9. Re-run fold detection with logged precision/tolerances.
10. Only then attach the fourth-vector ideal and run certification.
