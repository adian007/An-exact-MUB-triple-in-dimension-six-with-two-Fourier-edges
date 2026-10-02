# Singular-system formulation

At fixed `lambda*`, use phase coordinates `a=(a1,...,a5)` and `xj=exp(i aj)`. Let `F(a,lambda)` denote the six real MU residuals. A singular solution satisfies

`F(a,lambda*) = 0`

and

`rank(J_a F) < 5`.

A robust augmented formulation introduces a nonzero null vector `u` and imposes

`F = 0`,
`J_a F u = 0`,
`||u||^2 = 1`.

Because this produces a redundant/overdetermined representation when all six MU equations are retained, a practical implementation should choose a verified independent 5-equation chart, solve the square deflated system, and then check the sixth MU equation afterward. Alternatively use a projective null-vector chart (`u_j=1`) and a square subsystem, cycling charts to avoid coordinate singularities.

For exact algebra, use the Laurent-polynomial `I3_single_vector.txt` representation with inverse variables and `z t - 1 = 0`. For numerical fold detection, phase coordinates avoid inverse-variable conditioning problems.
