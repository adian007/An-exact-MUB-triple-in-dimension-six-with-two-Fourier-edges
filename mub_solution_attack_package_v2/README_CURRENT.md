# MUB Solution Attack Package v2

This package is the current research handoff for attacking the dimension-6 MUB problem through Karlsson's K_6^(3) family.

## What is established computationally in this project

- Exact Laurent representation of the Dita circle H_D(z).
- Corrected column-based MU-vector polynomial system I3.
- A verified 12-element fixed-z monomial symmetry group whose row permutations form A4.
- At the critical λ*, numerical refinement of 24 singular MU-vector candidates and a 12+12 A4 orbit decomposition, with residual/orbit diagnostics recorded separately as numerical evidence.
- At λ=0.4, π/3, 2π/3, numerical recovery of 72 MU vectors and four complete B3 cliques; the four cliques form one A4 orbit within each recovered pool.
- Numerical Fourier-family fits for selected third-basis transitions at λ=0.4 and 2π/3 with errors at approximately 10^-15 scale in the tested chart.

## What remains unproved

- Exhaustive completeness of every numerical MU pool.
- Generic completeness of the 12-element A4 stabilizer from one sample.
- Certified singular-root isolation at λ*.
- Exact Fourier-family membership for every B3 branch and every relevant orientation.
- Nonexistence of a fourth vector for the whole Dita circle.
- The full statement N(6)=3.

## Strategic conclusion
The next theorem to attack is not "Dita has no fourth MUB" by itself. It is the K3 structural theorem that every relevant third-MUB component carries a Fourier-family transition matrix (or the stronger X/F/F^T triplet structure). Combined with the rigorous whole-Fourier-family obstruction, that would produce a genuine quartet obstruction for the K3 component.

This package intentionally distinguishes numerical evidence from proof.
