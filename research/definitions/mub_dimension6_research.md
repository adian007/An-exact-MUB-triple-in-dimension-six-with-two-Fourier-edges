# Definitions and logical reductions

## Core objects

A basis `B` of `C^6` is mutually unbiased to `C` when `|<b,c>|^2 = 1/6` for every pair of vectors. A complex Hadamard matrix `H` is an unnormalised matrix with unimodular entries and `H H^* = 6 I`; its normalised columns form a basis unbiased to the computational basis `I`.

For a fixed pair `{I,H}`, the MU-vector variety is

`P(H) = { [v] : |v_i| = 1/sqrt(6), |<v,H_:k/sqrt(6)>|^2 = 1/6 for all k }`.

A third MUB is a six-element orthonormal clique in `P(H)`. Once a third basis `B3` is fixed, define the one-vector witness

`W1(H,B3) = { [v] : v is MU to I, H, and B3 }`.

If a fourth MUB exists, all six of its vectors lie in `W1(H,B3)`. Therefore `W1 = empty` is sufficient to exclude a fourth MUB for that fixed triple. The converse is false.

## Gauge invariance

Left and right diagonal unitaries, row and column permutations, and global column phases preserve the MUB question after transporting all bases. Consequently, the computational object is a CHM equivalence class plus the transported third basis, not a raw parameter tuple. A numerical comparison must dephase before declaring inequivalence.

## Critical scope distinction

`No fourth MUB for the tested triples` does not imply:

- no fourth MUB for every third basis not recovered by the pool;
- no fourth MUB for every `H` in `K_6^(3)`;
- `N(6) = 3` globally.

A family theorem needs a completeness argument for both the MU pool and the third-basis clique enumeration, plus a uniform witness exclusion.

## Research consequence

The cleanest theorem route is local-to-global:

1. exact algebraic anchors;
2. certified fixed-parameter witness emptiness;
3. certified continuation/covering over a compact parameter region;
4. explicit treatment of singular and boundary strata.

Reference: Grassl, arXiv:quant-ph/0406175; Brierley--Weigert, arXiv:0901.4051; local `docs/METHODS.md`.
