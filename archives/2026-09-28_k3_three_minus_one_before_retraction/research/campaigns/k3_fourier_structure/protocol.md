# Campaign: K3 third-MUB transition structure

## Objective

Test, and ultimately prove algebraically, whether every relevant third-MUB
component in the Diţă slice of Karlsson's `K_6^(3)` family has a transition
matrix in the Fourier family, its transpose, or the stronger `X/F/F^T`
Hadamard-cube structure.

This campaign consumes complete, gauge-reconciled B3 cliques. It is a
parallel structural route to the existing fourth-vector (`W1`) obstruction;
it does not replace direct W1 certification.

## Inputs

- Canonical I3 equations and gauge checks.
- `src/Cliques.jl` clique enumeration.
- `results/campaigns/i3_singular_locus/third_mub_cliques.json`.
- The sampled pools in
  `mub_solution_attack_package_v2/01_results/`.

## Stages

1. Recompute all B3 cliques from the canonical pool representation.
2. Verify the row/column/dephasing map from the Karlsson representative to
   the displayed `H_D(z)` representative.
3. Form `T = B3^* H_D` for each clique, alongside the other two transitions
   of the triplet.
4. Search pivot-row/pivot-column/target-column dephasing charts for a column
   containing three -1 entries; apply the analogous row test to the transpose.
5. Record flatness/unitarity checks, numerical residuals, and the exact chart.
6. Export exact incidence equations for detected charts; pursue elimination
   only after the numerical chart audit is complete.

## Success criteria

The numerical stage succeeds only when every supplied clique has a recorded
chart, residual, and gauge conversion. A theorem-level result requires exact
component containment or a certified continuation argument over every relevant
component. Finite sample fits are not sufficient.

## Stop conditions

Stop and report the component as unresolved if pool completeness, gauge
conversion, transpose coverage, or Fourier-chart coverage fails.
