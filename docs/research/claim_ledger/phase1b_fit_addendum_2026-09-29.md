# Phase 1b claim-ledger addendum

Evidence bundle:
[`phase1b_certified_clique_audit_2026-09-29.md`](../../results/phase1b_certified_clique_audit_2026-09-29.md),
SHA-256 `127B9B6CD1F8C27593DDBBFA48E9F8CA044876B2C4DA095D1547EE65A7C333EC`.

| Finding | Tier | Evidence SHA-256 | Scope / limitation |
|---|---|---|---|
| Current post-normalization pool generation reproduces the expected pool and clique counts at \(0,\pi,4\pi/3,5\pi/3\) | NUMERICAL | `127B9B6CD1F8C27593DDBBFA48E9F8CA044876B2C4DA095D1547EE65A7C333EC` | Recovered pools from finite multi-start solves; completeness not certified |
| Both specified transition matrices for all 40 represented clique bases fit direct \(F(a,b)\) at \(10^{-6},10^{-8},10^{-10}\) | NUMERICAL | `127B9B6CD1F8C27593DDBBFA48E9F8CA044876B2C4DA095D1547EE65A7C333EC` | 12 stored plus 28 regenerated bases; not an exact membership proof or historical clique-index identity |
| Fitted family parameters agree by monomial equivalence across \((0,\pi)\), \((\pi/3,4\pi/3)\), and \((2\pi/3,5\pi/3)\), with explicit clique maps | NUMERICAL | `127B9B6CD1F8C27593DDBBFA48E9F8CA044876B2C4DA095D1547EE65A7C333EC` | Floating-point residuals; canonical fundamental-triangle coordinates were not computed |
| The sampled \(\lambda=0.4\) Dita point remains a separate class among the four tested Dita classes | NUMERICAL | `127B9B6CD1F8C27593DDBBFA48E9F8CA044876B2C4DA095D1547EE65A7C333EC` | Finite class sample only |
| The exact Dita-to-\(X_6\) chart applies at \(z=1\), hence \(H_D(0)\) has the stated \(X_6\) representation | EXACT | `69ED3DFD16A44936FC9EF13B141DD5CCE02A1DB484250A39B04543D2281756E6` | Exact symbolic identity and chart are in the earlier Phase 1 report; fitter parameterization itself is non-unique |

No row supports pool completeness, exact membership of the fitted transition
matrices, fold-wide behavior, or a global MUB nonexistence claim.
