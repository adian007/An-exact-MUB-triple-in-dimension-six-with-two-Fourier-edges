# Numerical strategy for the dimension-six problem

## Why ordinary sweeps are insufficient

A grid in `(theta,phi,lambda)` can miss a one-dimensional locus, as the local 512-grid result demonstrates. Conversely, a near-clique caused by rounded parameters can look like a genuine third MUB. The correct numerical object is a per-H certified MU pool followed by high-precision graph verification.

## Recommended pipeline

For each fixed H:

1. validate the CHM residual;
2. solve the 10-variable pool system with polyhedral homotopy;
3. certify tracked roots and filter to the conjugate locus;
4. deduplicate projectively at several tolerances;
5. require an independent reseeded count agreement;
6. enumerate all size-six cliques, not just one;
7. verify each basis at high precision;
8. run the W1 witness test for every basis.

## Failure modes

A failed path is not automatically a missing physical vector. At degenerate points, roots can escape to infinity or collide. Therefore pool completeness must be a compound claim: path accounting, certification, conjugacy, independent reseed, and physical-count stability.

Similarly, `found_fourth=false` from a single chosen clique is not enough. The all-clique result is the relevant fixed-H statement.

## Continuation opportunity

The Dita circle is ideal for monodromy/continuation, but the current local notes show that the third-basis clique selection can vary between runs. Continue labelled roots and use projective matching rather than raw pool indices. At constellation changes, branch to a new stratum rather than forcing a single smooth labelling.

## Numerical research recommendation

Prioritize certified parameter boxes near the two observed third-MUB loci, then test whether W1 emptiness has a positive residual margin. A margin is the quantity needed to turn isolated certified instances into interval-covered neighborhoods.
