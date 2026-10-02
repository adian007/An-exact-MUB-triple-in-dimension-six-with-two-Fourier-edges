# Completed Fold + A4 Computation

This package completes the next computation in the dimension-6 MUB project: refinement of the critical MU-pool singular roots and their A4 orbit decomposition.

## Main result

At lambda* = 0.1114802243779665542913031975274172717685818097497045..., the corrected column-wise MU system shows 96 distinct numerical roots, with 24 refined singular roots and 72 regular roots. The 24 singular roots decompose into two A4 orbits of size 12, with maximum numerical matching error below 1e-9.

The refined roots also satisfy numerical fold diagnostics in a five-equation chart: both the parameter transversality projection and the quadratic null-direction projection are nonzero at all 24 roots.

## Scientific status

These are high-precision-in-principle numerical results based on double-precision finite-difference Jacobians and nonlinear refinement. They are not yet exact or interval-certified singular-root results. HomotopyContinuation.jl certification is designed for nonsingular isolated solutions; singular points require a separate deflation/interval treatment. See the project documentation and official HomotopyContinuation certification documentation.

## Reproduce

From the package root:

    python scripts/refine_singular_augmented.py
    python scripts/a4_orbit_refined.py
    python scripts/fold_diagnostics.py

The scripts expect SciPy and NumPy.

## Research interpretation

The computation supports the geometric picture

    120 roots  ->  96 distinct roots at the critical parameter  ->  72 roots
                       |
                       +-- 24 singular roots = 12 + 12 under A4

The next object is the A4 quotient of complete third-MUB cliques B3, followed by the fourth-vector incidence ideal I4 for inequivalent representatives.
