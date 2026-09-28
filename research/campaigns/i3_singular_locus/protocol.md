# Campaign: I3 singular locus and A4 reconciliation

## Objective

Reconcile the fold/A4 handoff with the repository's canonical column-based
I3 equations and the existing Karlsson pool/clique workflow.

The target parameter is

```text
lambda* = 0.1114802243779665542913031975274...
```

## Inputs

- `scripts/python/dita_i3.py`
- `scripts/julia/i3_singular_locus.jl`
- `mub_fold_a4_computation/` (provenance and independent reproduction)
- `mub_solution_attack_package_v2/04_provenance/` (archived provenance)

## Required checks

1. Use the same displayed `H_D(z)` representative and column orientation.
2. Reproduce the regular transition counts `120, 72, 72`.
3. Recover and validate the 24 physical augmented singular roots.
4. Verify the `12 + 12` A4 decomposition in the canonical gauge.
5. Reconcile the four-clique results at `0.4`, `pi/3`, and `2*pi/3`.

## Evidence boundary

The campaign can establish numerical or interval-certified properties of the
supplied roots and supplied pools. It does not, by itself, prove generic
completeness of the stabilizer, completeness of every pool, or
fourth-MUB nonexistence over a parameter family.

## Outputs

Outputs belong under `results/campaigns/i3_singular_locus/`. Large root tables
are raw evidence; the interpretation must remain in the claim ledger and
campaign summary.
