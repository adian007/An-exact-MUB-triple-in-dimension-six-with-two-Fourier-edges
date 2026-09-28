# K3 Fourier-structure workflow

This campaign is the structural companion to the direct fourth-vector
obstruction workflow. It consumes the same canonical I3 pools and B3 cliques;
it does not bypass pool completeness or gauge reconciliation.

The first-stage invariant follows the published order-six characterization
referenced at
<https://link.springer.com/article/10.1007/s10623-024-01503-w>: after
normalization, Fourier-family membership is characterized by an equivalent
matrix with a column containing three entries equal to -1; the transposed
family has the corresponding row criterion. Floating-point detections remain
numerical candidates.

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

The diagnostic searches every pivot-row/pivot-column/target-column dephasing
chart for `T = B3^* H_D` and the other transitions. In a chart, a dephased
entry equals -1 precisely when

```text
T[i,j]*T[r,c] + T[i,c]*T[r,j] = 0
```

The exporter writes the exact I3 branch incidence equations together with the
three exact equations for each detected witness chart. These are elimination
inputs, not an elimination result; the sample-specific numerical B3 solutions
have not been replaced by exact branch parametrizations.

The current sample result is summarized in
[`k3_fourier_structure_diagnostic.md`](../results/k3_fourier_structure_diagnostic.md).

## Evidence boundary

At the diagnostic stage, a small three-minus-one residual means “numerical
candidate for this supplied clique and chart.” It does not mean that every
clique, branch, or parameter value is Fourier-family equivalent. The theorem
target is exact algebraic containment or certified continuation over each
relevant component.
