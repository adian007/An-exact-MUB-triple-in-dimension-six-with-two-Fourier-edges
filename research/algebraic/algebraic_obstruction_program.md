# Algebraic obstruction program

## Fixed-point formulation

For a flat vector `v = z/sqrt(6)` with gauge `z_0 = 1`, introduce formal inverses `w_i` and impose `z_i w_i - 1 = 0`. MU to a CHM column `h_k` is encoded by

`(sum_j conj(h_jk) z_j)(sum_j h_jk w_j) - 6 = 0`.

For a fixed third basis add the analogous six equations. The resulting `W1` has 10 variables and 17 equations in the current encoding.

## Exact route

At an algebraic anchor, work over a number field containing every coefficient. The strongest fixed-point certificate is a Groebner basis containing `1`, or an equivalent exact ideal-membership result. A computation over floating `CC` is not enough to support that claim; it can diagnose, but not certify, emptiness.

The local D0-equivalent calculation over `Q(zeta_24,sqrt(5))` is therefore a model for the exact route, but it covers one fixed pair and does not prove the whole Dita circle.

## Family route

A family theorem would need elimination in both vector variables and parameters, followed by a real/unit-circle analysis. Direct elimination is likely to suffer from:

- denominators in Mobius parameterisations;
- components at infinity introduced by inverse variables;
- positive-dimensional components at Dita degeneracies;
- gauge copies and permutation orbits;
- coefficient-field growth.

A better decomposition is stratified: generic parameter stratum, D0-equivalent special stratum, and boundary/pole strata. Eliminate only after saturating by the known denominator and gauge factors.

## Concrete next algebraic tests

1. Reproduce the exact D0 unit ideal with an independently generated polynomial file.
2. Compute saturation by all Mobius denominators before interpreting dimension.
3. Compare elimination ideals from two distinct Dita lambda representatives.
4. Use modular Groebner calculations to estimate degree and detect accidental coefficient errors before number-field computation.
5. Record the field, monomial order, saturation factors, and certificate hash.

The central research conclusion is negative but useful: a failed inexact elimination must be labelled inconclusive, not converted into an absence claim.
