# Research Gap Report

## Smallest meaningful next theorem

Prove exact fourth-MUB nonextendability at one algebraic Diţă/Karlsson point that is not CHM-equivalent to `D0`, preferably `lambda=0` or `lambda=pi/3` on the Karlsson Diţă slice.

A strong target is:

> For an explicitly reconstructed exact third basis `B3` at the chosen non-D0-equivalent point, the single-vector witness ideal `W1(H,B3)` is the unit ideal over a named number field.

A stronger, more literature-comparable target is:

> Enumerate the complete MU pool for `{I,H}` exactly or with a complete certificate, enumerate every third basis, and prove no two third bases are MU.

## Why this target is properly scoped

- It changes the quantifier from a floating sample to one exact point.
- It is outside the already-settled D0 case.
- It does not claim all Diţă parameters, all Karlsson matrices, or `N(6)=3`.
- The repository already has the exact-W1 workflow at `(D_bc,F_D)` and exact H reconstructions at algebraic lambda values.

## Required evidence

1. Exact H entries and a named coefficient field.
2. An independently reconstructed exact B3; do not fit a numerical B3 into an unjustified small field.
3. Explicit witness equations with conjugates encoded correctly.
4. Exact Gröbner output over a genuine number field, preferably `GB={1}`.
5. Independent CHM, ONB, and MUB identities for H and B3.
6. A source/equivalence check showing the chosen point is not D0-equivalent.
7. A clear statement whether the result covers one B3 or every third basis.

## Do not count as the target

- More dense lambda sampling.
- `dim I=-1` over `CC` with float coefficients.
- A certified numerical result at a D0-equivalent point.
- A finite search described as a theorem on a continuous family.
- A candidate B3 with only partial LLL reconstruction.
