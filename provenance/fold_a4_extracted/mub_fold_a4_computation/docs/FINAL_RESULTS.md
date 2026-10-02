# Fold + A4 Computation: Final Results

## Critical parameter

lambda* = 0.1114802243779665542913031975274172717685818097497045...

## Pool counts

Using the corrected column-wise MU equations:

- lambda* - 2e-6: 120 distinct numerical MU-vector roots recovered.
- lambda*: 96 distinct numerical roots recovered.
- lambda* + 2e-6: 72 distinct numerical MU-vector roots recovered.

Thus the observed transition is 120 -> 72, a loss of 48 regular solutions.

## Refined singular system

For each candidate, the augmented numerical system was solved:

F(a,lambda*) = 0,
J(a,lambda*) u = 0,
||u||^2 - 1 = 0,

where a are five phase coordinates and J is the 6 x 5 Jacobian of the six real MU equations.

A reduced set of 24 candidates satisfied the augmented system to numerical precision. Across these 24 roots:

- maximum ||F|| < 1e-9,
- maximum ||J u|| < 5e-9,
- ||u||^2 = 1 to the solver tolerance,
- smallest Jacobian singular values are approximately 1e-10 to 1e-9 in double precision.

These are numerical singular-root results, not interval-certified singularity proofs.

## A4 orbit decomposition

The previously verified fixed-z monomial generators G2 and G3 were applied to the 24 refined roots. The induced permutations close on the 24-root set with maximum matching error below 1e-9.

The 24 roots decompose into exactly two numerical A4 orbits:

12 + 12.

Therefore the previous hypothesis 24 = 12 + 12 is supported by the refined computation.

The generators have order/group structure consistent with the previously identified A4 action of order 12. The result should still be described as a numerical orbit decomposition until the singular points and group action are certified algebraically or by interval enclosures.

## Fold diagnostics

For the first five MU equations, the 5 x 5 Jacobian has rank 4 at each refined singular root. Let w span its numerical left nullspace, u the right null vector, and F_lambda the parameter derivative. The diagnostics give:

3.018 <= |w^T F_lambda| <= 5.390

1.078 <= |w^T F_{aa}[u,u]| <= 4.660

across the 24 roots (finite-difference evaluation in double precision).

Thus every refined root passes the numerical transversality/nonzero-quadratic tests expected for a simple fold in the reduced five-equation chart. This is evidence for a fold normal form, not a certified theorem.

## Important status distinction

Established computationally in this package:

1. The corrected column-wise MU system reproduces the 120/96/72 transition.
2. Twenty-four singular candidates refine to solutions of the augmented null-vector system.
3. Those 24 roots split numerically into two A4 orbits of size 12.
4. All 24 pass nonzero first-parameter and quadratic fold diagnostics numerically.

Not established by this computation:

1. A rigorous interval certificate for the 24 singular roots.
2. A proof that the A4 stabilizer/action is complete over the exact algebraic singular locus.
3. A global theorem describing every singularity of the Diita-circle MU-pool.
4. Any fourth-MUB nonexistence result at lambda*.
