# MUB Solution Attack Package v2

This package is the current research handoff for attacking the dimension-6 MUB problem through Karlsson's K_6^(3) family.

## What is established computationally in this project

- Exact Laurent representation of the Dita circle H_D(z).
- Corrected column-based MU-vector polynomial system I3.
- A verified 12-element fixed-z monomial symmetry group whose row permutations form A4.
- At the critical λ*, numerical refinement of 24 singular MU-vector candidates and a 12+12 A4 orbit decomposition, with residual/orbit diagnostics recorded separately as numerical evidence.
- At λ=0.4, π/3, 2π/3, numerical recovery of 72 MU vectors and four complete B3 cliques; the four cliques form one A4 orbit within each recovered pool.
- Separate finite-permutation numerical fits for selected third-basis transitions at λ=0.4 and 2π/3 had residuals of approximately 2.5e-15 and 5.4e-15. These direct fits are independent of, and not established by, the legacy one-column diagnostic.

## What remains unproved

- Exhaustive completeness of every numerical MU pool.
- Generic completeness of the 12-element A4 stabilizer from one sample.
- Certified singular-root isolation at λ*.
- Exact family assignment for every B3 branch and every relevant orientation remains unproved.
- Nonexistence of a fourth vector for the whole Dita circle.
- The full statement N(6)=3.

## Strategic conclusion
The next structural target is to establish Theorem 1's three-distinct-columns condition for every relevant third-MUB component, then resolve its transposed-Fourier or 2-circulant alternative. The Fourier-family quartet obstruction applies only if the Fourier branch is separately established.

This package intentionally distinguishes numerical evidence from proof.
