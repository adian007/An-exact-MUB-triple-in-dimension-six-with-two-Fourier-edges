# Research status

## Established computational inputs

- Exact Laurent representation of the Diţă circle is available in the parent I3 package.
- The single-vector symbolic ideal `I3^(1)` is available.
- The fixed-z search found 12 compatible monomial automorphisms at a generic numerical sample; each found automorphism was symbolically verified.
- The row permutation parts form an order-12 group consistent with `A4`.
- The fold parameter is numerically identified as `lambda* = 0.1114802243779665542913031975...`.
- Numerical continuation reports a 120-to-72 pool transition around this parameter.

## Not yet a theorem

- Exact completeness of the 12-element generic stabilizer.
- Exact existence and multiplicity of all 24 singular roots.
- Exact decomposition of the 24 singular roots into A4 orbits.
- A statement that the fold is globally the only singular event on the relevant interval.
- Fourth-MUB non-extendability for the entire Diţă circle or all of `K_6^(3)`.

## Certification note

HomotopyContinuation.jl documents interval/Krawczyk certification for nonsingular isolated solutions. Singular roots need deflation, a suitable singular certification method, or an exact algebraic argument. Therefore numerical convergence alone must remain labeled numerical evidence.
