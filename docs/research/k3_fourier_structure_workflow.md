# K3 Fourier-structure workflow

This campaign is the structural companion to the direct fourth-vector
obstruction workflow. It consumes the same canonical I3 pools and B3 cliques;
it does not bypass pool completeness or gauge reconciliation.

The verified result is Theorem 1 of Matszangosz–Szöllősi,
<https://link.springer.com/article/10.1007/s10623-024-01503-w>. For a
normalized order-six complex Hadamard matrix, its condition is that three
distinct columns each contain at least one \(-1\). The conclusion is that the
matrix belongs to the transposed Fourier family **or** the 2-circulant family.
Applying the theorem to the transpose gives the row formulation by inference.
This is not a Fourier-family-only characterization. Floating-point detections
remain numerical candidates.

## Workflow

```text
canonical H_D / I3
    -> complete MU pool
    -> all B3 cliques
       -> W1 fourth-vector certification
       -> Fourier/X/F^T transition audit
          -> exact incidence/elimination target
```

The source material is preserved in
[`mub_solution_attack_package_v2`](../../mub_solution_attack_package_v2/README_CURRENT.md).
Its Fourier fits are sample-level numerical evidence. The repository claim
ledger must not treat them as a proof of component containment.

## Reproduction inputs

- Protocol: [`research/campaigns/k3_fourier_structure/protocol.md`](../../research/campaigns/k3_fourier_structure/protocol.md)
- Parameters: [`research/campaigns/k3_fourier_structure/parameters.json`](../../research/campaigns/k3_fourier_structure/parameters.json)
- Source pools: `mub_solution_attack_package_v2/01_results/`
- Existing clique implementation: [`src/Cliques.jl`](../../src/Cliques.jl)
- Existing singular-locus campaign: [`i3_singular_locus`](../../research/campaigns/i3_singular_locus/protocol.md)

## Run the implemented stage

From the repository root:

```powershell
python scripts/python/analyze_k3_fourier_structure.py
python scripts/python/export_k3_fourier_incidence.py
```

The legacy diagnostic `three_minus_one_test` searches for three \(-1\)-like
entries in one dephased column. That is not Theorem 1 and does not establish
family membership. The new theorem predicate searches all 36 pivot-row and
pivot-column dephasings and checks three distinct columns, each with a
\(-1\)-like entry. The dephased entry condition is

```text
T[i,j]*T[r,c] + T[i,c]*T[r,j] = 0
```

The campaign applies the theorem predicate to both a matrix and its transpose;
the latter is labelled as the inferred row-form application. The exporter
writes exact I3 incidence equations and the three exact equations for each
detected theorem chart. These are elimination inputs, not an elimination
result; the sample-specific numerical B3 solutions have not been replaced by
exact branch parametrizations. The theorem conclusion is the disjunction
“transposed Fourier or 2-circulant,” not Fourier membership alone.

The current sample result is summarized in
[`k3_fourier_structure_diagnostic.md`](../results/k3_fourier_structure_diagnostic.md).

## Evidence boundary

At the diagnostic stage, a small residual means only a numerical candidate
for the stated theorem predicate on the supplied CHM and chart. It does not
identify which of the two theorem families applies, nor prove the condition
for every clique, branch, or parameter value. Component-wise algebraic
containment or certified continuation remains the theorem target.
