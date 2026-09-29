# Phase 1b: regenerated pools, Fourier fits, and class equivalence

Date: 2026-09-29

This is a new addendum to the earlier Phase 1 audit. It supersedes that
report's statement that the four exceptional-parameter matrix sets were
unavailable. It does not overwrite or revise the earlier report or any
package files.

## Files and provenance

**Created:** this report,
`docs/research/claim_ledger/phase1b_fit_addendum_2026-09-29.md`, and
`results/campaigns/i3_singular_locus/family_fits_postnorm_20260929_new4.json`
(SHA-256
`4C876B63478EEFD3F7E07C7265B37ADCEDFA7460B3E8763B66594993EE8CB5DB`).

**Modified:** none.

**Untouched evidence artifacts (SHA-256):**

| File | SHA-256 | Note |
|---|---|---|
| `results/campaigns/i3_singular_locus/third_mub_cliques.json` | `EA8176FCDB853A15A4B1D710B8954E5A59CB59C55C7270165A55A53651C2A2A5` | Existing post-normalization export at \(0.4,\pi/3,2\pi/3\) |
| `results/campaigns/i3_singular_locus/pool_rebuild_postnorm_20260929_lambda_0.json` | `54D314D73FD9FF11C195E8A51D97AA2350FFFD799FD84FE4B268BEA7A878151F` | Regenerated pool, count gate passed |
| `results/campaigns/i3_singular_locus/pool_rebuild_postnorm_20260929_lambda_pi.json` | `9E99E233D5643AC70E84AD1E5034045E8232551F0C57BBB6A1A4F9D0C3B678B` | Superseded metadata: pool is valid, expected-count metadata was wrong |
| `results/campaigns/i3_singular_locus/pool_rebuild_postnorm_20260929_lambda_pi_count_gate_corrected.json` | `53999CAA88AC63D5F2783C9508EC070918885D6A91657ECCC471B649F8054728` | Corrected metadata; pool and cliques copied unchanged from preceding file |
| `results/campaigns/i3_singular_locus/pool_rebuild_postnorm_20260929_lambda_4pi_over_3.json` | `9FD6A7F8B27663C41684B82F971147B01AF76B02C84CA3254C26C72FDEF2C2A0` | Regenerated pool, count gate passed |
| `results/campaigns/i3_singular_locus/pool_rebuild_postnorm_20260929_lambda_5pi_over_3.json` | `B48F4F34A88285E825244F3AEA7064C8F37472ABA7F1B42707E445545A1D181F` | Regenerated pool, count gate passed |
| `docs/results/phase1_certified_clique_audit_2026-09-29.md` | `69ED3DFD16A44936FC9EF13B141DD5CCE02A1DB484250A39B04543D2281756E6` | Earlier report, unchanged |
| `readme.md` | not a result artifact | Pre-existing user modification, left untouched |

The final check found no diff in `Project.toml` or `Manifest.toml`. The
verified Julia environment used Julia 1.12.6, Pkg 1.12.1, and the pre-existing
`/home/adian/mub-depot` on ext4. Its eight pinned artifact trees matched their
expected hashes. The fast regression suite passed 101/101 tests (exit 0).

## Commands and resource use

All science runs were sequential. The four pool solves used
`MubSearch.generate_candidate_pool_fresh(H; tol=1e-8)` with HomotopyContinuation
2.22.1, one Julia thread, and 252 start paths. Pool filtering used \(10^{-8}\),
deduplication \(10^{-6}\), clique orthogonality \(10^{-8}\), and 128-bit
verification at \(10^{-10}\). Each solve and clique reconstruction returned
exit 0. The pool artifact's `peak_rss_bytes` is reported below.

The family-fit runner processed one matrix at a time using
`fit_fourier_family`, `fit_two_circulant`, and the theorem-predicate helper
from `scripts/python/k3_family_membership.py`; it returned exit 0 after
426.83 seconds, with peak working set 85,319,680 bytes. The exceptional-point
and monomial-equivalence checks used
`analyze_k3_fourier_structure.dita_hadamard`,
`chm_equivalence_residual_full`, and `fit_fourier_family`; both sequential
Python runs returned exit 0. The latter run's peak working set was
81,457,152 bytes.

| \(\lambda\) | Seed | Raw solver solutions | Pool | 6-cliques | Peak RSS |
|---:|---:|---:|---:|---:|---:|
| \(0\) | 20260929 | 192 | 120 | 10 | 1,306,402,816 bytes |
| \(\pi\) | 20260930 | 192 | 120 | 10 | 1,308,979,200 bytes |
| \(4\pi/3\) | 20260931 | 240 | 72 | 4 | 1,315,278,848 bytes |
| \(5\pi/3\) | 20260932 | 240 | 72 | 4 | 1,314,230,272 bytes |

The first \(\pi\) output was marked `count_mismatch` because its expected
metadata incorrectly said 72 vectors / 4 cliques. Its observed pool was
120 / 10. The named corrected artifact changes only the expected-count
metadata and correctly passes the gate. No pool data was overwritten.

## Phase 1.4: count gate

The current post-2026-09-17 solver and unit-vector normalization reproduced
the expected counts at all four points. The pool builder's results are
multi-start numerical outputs; pool completeness is not certified.

## Phase 1.5: fits

For each represented clique, both \(B_3\) and \(H_D^\dagger B_3\) were tested
at \(10^{-6},10^{-8},10^{-10}\). The direct \(F(a,b)\) call, a separate
transpose call, \(X_6\) in both orientations, and the three-distinct-columns
predicate on \(M\) and \(M^T\) were evaluated. The result JSON contains the
per-matrix residuals, permutations, fitted parameters, and tolerance statuses.

| \(\lambda\) | Matrices | Direct \(F\), each tolerance | \(F\) transpose call, each tolerance | \(X_6\), either orientation, each tolerance | Predicate \(M/M^T\) | Maximum direct-\(F\) score |
|---:|---:|---:|---:|---:|---:|---:|
| \(0.4\) | 8 | 8 FIT | 0 FIT | 0 FIT | 0 / 8 | \(<5\times10^{-15}\) |
| \(\pi/3\) | 8 | 8 FIT | 0 FIT | 0 FIT | 0 / 8 | \(<5\times10^{-15}\) |
| \(2\pi/3\) | 8 | 8 FIT | 0 FIT | 0 FIT | 0 / 8 | \(<5\times10^{-15}\) |
| \(0\) | 20 | 20 FIT | 0 FIT | 0 FIT | 0 / 20 | \(8.456\times10^{-16}\) |
| \(\pi\) | 20 | 20 FIT | 0 FIT | 0 FIT | 0 / 20 | \(9.695\times10^{-16}\) |
| \(4\pi/3\) | 8 | 8 FIT | 0 FIT | 0 FIT | 0 / 8 | \(2.785\times10^{-15}\) |
| \(5\pi/3\) | 8 | 8 FIT | 0 FIT | 0 FIT | 0 / 8 | \(2.234\times10^{-15}\) |
| **Total** | **80** | **80 FIT at all three tolerances** | **0** | **0** | **0 / 80** | — |

The 80 matrices are the two specified transitions for 40 represented clique
bases: 12 from the stored \(0.4,\pi/3,2\pi/3\) export and 28 from the four
newly regenerated pools. The largest flatness or unitarity discrepancy in
the 56-matrix new fit artifact was below \(1.56\times10^{-15}\). Across all
80 matrices the result is consistent with direct \(F\) fits, with the
three-distinct-columns predicate appearing only on the transpose. There were
no matrices in none of \(F,F^T,X_6\), and no predicate-positive,
family-negative anomaly.

**Special point \(H_D(0)\).** Direct and transpose Fourier fits both returned
NONE (score \(\sqrt2\)). The \(X_6\) fitter returned FIT in both orientations;
the matrix residual was \(2.45\times10^{-16}\) and the constraint residual
was \(2.59\times10^{-16}\). The exact chart identity from the preceding
Phase 1 report applies at \(z=1\), with
\((\beta,\gamma,\epsilon,\phi)=(1,-1,i,i)\). The numerical fitter selected
another equivalent parameterization, so its particular output parameters
are not unique.

## Phase 1.6: paired-parameter equivalence

The source's Section 3 fundamental triangle is cited in
`docs/results/k3_family_membership_2026-09-29.md` (arXiv:0902.0882v2,
printed p. 9); the paper refers to [5] for the parameter equivalences. Rather
than infer coordinates from raw fitter phases, each fitted family matrix was
compared under the full row/column monomial equivalence search. The source
Hadamard pairs all use the same row permutation
\((0,1,2,3,5,4)\), column permutation \((0,1,3,2,4,5)\), and no
conjugation or transpose.

| Pair | \(H_D\) equivalence residual | Target-to-source clique map (0-based) | Maximum projective vector residual | Fitted-family equivalence residual range |
|---|---:|---|---:|---:|
| \(\pi/3 \to 4\pi/3\) | \(4.003\times10^{-16}\) | \(0\to3,\ 1\to2,\ 2\to0,\ 3\to1\) | \(<10^{-12}\) | \(1.271\times10^{-15}\)–\(3.231\times10^{-15}\) |
| \(2\pi/3 \to 5\pi/3\) | \(4.003\times10^{-16}\) | \(0\to0,\ 1\to2,\ 2\to1,\ 3\to3\) | \(<10^{-12}\) | \(2.116\times10^{-15}\)–\(3.388\times10^{-15}\) |
| \(0 \to \pi\) | \(1.155\times10^{-16}\) | \(0\to3,\ 1\to5,\ 2\to0,\ 3\to9,\ 4\to2,\ 5\to8,\ 6\to6,\ 7\to7,\ 8\to1,\ 9\to4\) | \(1.490\times10^{-8}\) | \(1.221\times10^{-15}\)–\(2.907\times10^{-15}\) |

The clique maps are permutations of all cliques at each paired point. The
small \(1.49\times10^{-8}\) discrepancy at \(0\leftrightarrow\pi\) is below
the \(10^{-7}\) projective matching tolerance; the other reported matches
were at floating-point zero. For every mapped pair, the raw fitted complex
\((a,b)\) phases differ, but their \(F(a,b)\) matrices are monomial
equivalent at residuals of a few \(10^{-15}\). Thus they represent the same
numerical equivalence class after reduction to the source's fundamental
region. Explicit canonical triangle coordinates were not computed.
Tier: **NUMERICAL**.

## Phase 1.7: fitted clique coverage and class labels

| \(\lambda\) | Represented cliques | Source of matrices | \(B_3\) and \(H_D^\dagger B_3\) direct-\(F\) fits |
|---:|---:|---|---|
| \(0\) | 10 | Fresh post-normalization pool | 20 / 20 |
| \(0.4\) | 4 | Existing post-normalization export | 8 / 8 |
| \(\pi/3\) | 4 | Existing post-normalization export | 8 / 8 |
| \(2\pi/3\) | 4 | Existing post-normalization export | 8 / 8 |
| \(\pi\) | 10 | Fresh post-normalization pool; corrected metadata artifact | 20 / 20 |
| \(4\pi/3\) | 4 | Fresh post-normalization pool | 8 / 8 |
| \(5\pi/3\) | 4 | Fresh post-normalization pool | 8 / 8 |
| **Total** | **40** | **12 stored + 28 regenerated** | **80 / 80** |

Every represented clique basis has both requested transition matrices fitted;
none is missing among these 40 current/exported or regenerated bases.
However, the regenerated clique index sets need not be identical to the
August native witness-log cliques. Therefore this table does not claim an
identity match to all 40 historical witness objects, nor does it upgrade
their completeness or certification status. The sampled class label at
\(\lambda=0.4\) remains its own class among these four verified Dita classes;
it is not merged with the three paired classes above.

## Tiers, anomalies, and limitations

- **EXACT:** the \(H_D(z)\)-to-\(X_6\) chart and constraint-(3) identity, as
  symbolically verified in the earlier Phase 1 report. The \(z=1\) case is
  included.
- **NUMERICAL:** regenerated counts; all family fits; exceptional-point
  fitter outputs; and all monomial/clique/parameter-class equivalences.
- **Anomalies:** none observed in the tested matrices.
- Pool completeness was not certified. Numerical fits are not exact or
  family-level membership certificates. The finite sampled class labels do
  not establish a full-circle result. Nothing here addresses the fold or
  \(\lambda<\lambda_*\), and no global MUB nonexistence statement follows.
- Theorem 1.4 remains a cited published computer-assisted theorem, not a
  computation independently reproduced in this phase.

## Commit scope

The Phase 1b commit contains only this report, the separate ledger addendum,
and the generated fit artifact
`results/campaigns/i3_singular_locus/family_fits_postnorm_20260929_new4.json`.
Existing pool artifacts and the user's `readme.md` modification are not part
of that commit.
