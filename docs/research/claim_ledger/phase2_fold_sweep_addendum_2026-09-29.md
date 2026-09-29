# Phase 2 fold-sweep claim-ledger addendum

Evidence bundle:
[`fold_sweep_phase2_2026-09-29.md`](../../results/fold_sweep_phase2_2026-09-29.md),
SHA-256 `7E6F9F1E56EB14A7DD4E2D009528B3F0EC738D477DFBF22B5099DE2F253A1F36`.
Authoritative per-point data:
[`fold_sweep_postnorm_20260929_reconciled.json`](../../results/campaigns/i3_singular_locus/fold_sweep_postnorm_20260929_reconciled.json),
SHA-256 `546991EA6AB8AAEF2810148B0D653C943768F6913E621148AFFBCE896DD2198B`.

| Finding | Tier | Evidence SHA-256 | Scope / limitation |
|---|---|---|---|
| Four tested points below \(\lambda_*\) each returned 120 verified distinct pool vectors | NUMERICAL | `546991EA6AB8AAEF2810148B0D653C943768F6913E621148AFFBCE896DD2198B` | Multi-start recovered pools; completeness not certified |
| At \(\lambda_*\), the current Julia pool builder returned 72 vectors rather than the expected 96 | NUMERICAL | `546991EA6AB8AAEF2810148B0D653C943768F6913E621148AFFBCE896DD2198B` | The 24-count difference is not a matched list of missing singular vectors |
| At \(0.05,0.10,\lambda_*-2\cdot10^{-6}\), each 120-vector pool contains only four verified six-cliques | NUMERICAL | `546991EA6AB8AAEF2810148B0D653C943768F6913E621148AFFBCE896DD2198B` | Contradicts the requested qualitative expectation of more than four below-fold cliques at these tested values |
| The just-above-fold and archived Phase 1 comparison points reproduce their expected pool/clique counts | NUMERICAL | `546991EA6AB8AAEF2810148B0D653C943768F6913E621148AFFBCE896DD2198B` | Counts for recovered pools only; not a completeness certificate |
| Phase 2.2 family-fit status at the tested fold points | OPEN | `7E6F9F1E56EB14A7DD4E2D009528B3F0EC738D477DFBF22B5099DE2F253A1F36` | Fits were intentionally not run because the Phase 2 count gate failed |

No row certifies root-pool completeness, the critical singular roots, or
family-level membership.
