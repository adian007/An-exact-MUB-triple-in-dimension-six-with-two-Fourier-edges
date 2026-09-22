# Independent reproduction protocol

## Goal

Reproduce the strongest local claims without relying on cached pool vectors, cached clique indices, or the same random seed.

## Protocol

1. Run in the documented WSL2 Julia environment.
2. Record Julia, HomotopyContinuation, git commit, CPU, and thread count.
3. Construct H from exact formulas, not serialized floating matrices.
4. Solve the MU pool twice with unrelated seeds.
5. Compare physical counts after projective deduplication at `1e-6`, `1e-8`, `1e-10`, and `1e-12`.
6. Enumerate every size-six clique from the independently reconstructed pool.
7. Verify each candidate B3 at higher precision.
8. Run W1 certification on every clique, saving the complete path accounting.
9. Hash the source, parameter manifest, and result JSON.
10. Compare semantic claims, not pool indices: counts, residual margins, clique cardinalities, and verdicts.

## Independence tests

A reproduction is weak if it shares cached B3 data or copied random matrices. It is stronger when it regenerates the polynomial system and uses a different solver seed, precision, and implementation path.

## Acceptance

A fixed-point result can be called reproduced only when the physical pool count, all-clique count, and W1 verdict agree, with any differences explained by equivalent basis permutations or projective ordering.
