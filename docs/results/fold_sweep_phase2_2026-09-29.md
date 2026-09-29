# Phase 2: fold-sweep count gate

Date: 2026-09-29

This report stops at the Phase 2 count gate. The recovered pool count at the
critical value did not match the expected total, and below-fold clique counts
also contradicted the requested qualitative expectation. No family fits were
run.

## Files and provenance

**Created:**

- `results/campaigns/i3_singular_locus/fold_sweep_postnorm_20260929_reconciled.json`
  — authoritative full sweep, SHA-256
  `546991EA6AB8AAEF2810148B0D653C943768F6913E621148AFFBCE896DD2198B`.
- This report.
- `docs/research/claim_ledger/phase2_fold_sweep_addendum_2026-09-29.md`.

**Intermediate artifact:** the first collector created
`results/campaigns/i3_singular_locus/fold_sweep_postnorm_20260929.json`,
SHA-256 `6AB8FF00C194D48A55085C97D1D2532C92FA9E8DD2E76C0C85FB0A2A48CA6BD0`.
It contains only the first two points, incorrectly extrapolates 10 expected
cliques from \(\lambda=0\) to \(\lambda=0.05\), and is superseded by the
reconciled artifact. It was removed as a scratch output after its hash and
failure were recorded here; the reconciled artifact is not a rewrite of it.

**Modified:** no pre-existing source, package, or result files. The user's
pre-existing `readme.md` change was left untouched. `Project.toml` and
`Manifest.toml` have no diff.

**Untouched reference artifacts (SHA-256):**

| File | SHA-256 |
|---|---|
| `mub_fold_a4_computation/results/RESULTS.md` | `EB7BF1EA7142550A07730BF906A64079D51A40FD40986D165026A328646152B6` |
| `mub_fold_a4_computation/results/below_summary.json` | `0CF976874B267CBAE6EB3BD84DB48B3B7307E12B7FBBC4B5286175BFA1D98958` |
| `mub_fold_a4_computation/results/critical_summary.json` | `B2CF3220BC784851991FE315F811A0328248D443357CBD1BE781B884466C1C6D` |
| `mub_fold_a4_computation/results/critical_1200_summary.json` | `1C539CE8115FD80744A59DF6F230AD1DDDAA715DC6BF1CC222CB4BB5B3B4B78F` |
| `mub_fold_a4_computation/results/above_summary.json` | `4D9BE36D8EE30DE19A81D53CB4E41F44625FF256E25E24F822275BACF38709A9` |
| `mub_fold_a4_computation/run_fold_a4.py` | `365D200E0A3E37A554624029302D449091A2F8B21D7576117EA60D97DFF83FBC` |
| `results/campaigns/i3_singular_locus/pool_rebuild_postnorm_20260929_lambda_0.json` | `54D314D73FD9FF11C195E8A51D97AA2350FFFD799FD84FE4B268BEA7A878151F` |
| `results/campaigns/i3_singular_locus/pool_rebuild_postnorm_20260929_lambda_pi_count_gate_corrected.json` | `53999CAA88AC63D5F2783C9508EC070918885D6A91657ECCC471B649F8054728` |
| `results/campaigns/i3_singular_locus/pool_rebuild_postnorm_20260929_lambda_4pi_over_3.json` | `9FD6A7F8B27663C41684B82F971147B01AF76B02C84CA3254C26C72FDEF2C2A0` |
| `results/campaigns/i3_singular_locus/pool_rebuild_postnorm_20260929_lambda_5pi_over_3.json` | `B48F4F34A88285E825244F3AEA7064C8F37472ABA7F1B42707E445545A1D181F` |

No files under `mub_solution_attack_package_v2` or legacy JSONs were edited.

## Archived evidence and expected root-count direction

The direction was checked from both the fold records and the regenerated
Phase 1 pools rather than inferred from parameter order:

| Regime / archived point | Archived evidence | Expected distinct physical vectors |
|---|---|---:|
| Below fold, \(\lambda_* - 2\times10^{-6}\) | `below_summary.json`: 114 found in the incomplete 400-start search; status notes a 120-to-72 transition | 120 |
| At \(\lambda_*\) | 400-start `critical_summary.json`: 93; 1200-start `critical_1200_summary.json`: 96, described as 72 regular plus 24 singular candidates | 96 |
| Above fold, \(\lambda_* + 2\times10^{-6}\) | `above_summary.json`: 70 found in the incomplete 400-start search; status notes a 120-to-72 transition | 72 |

Thus the local fold expectation is 120 below, 96 at the critical value, and
72 above. The 96 includes numerical singular-root candidates and is not a
singular-root certificate. The archived circle points also show that this is
a local fold regime, not a monotone rule for all \(\lambda\): the regenerated
\(\lambda=\pi\) pool has 120 vectors, while \(0.4,\pi/3,2\pi/3,4\pi/3\), and
\(5\pi/3\) have 72.

## Phase 2.1: sequential fresh-pool sweep

The chosen values were
\[
0,\ 0.05,\ 0.10,\ \lambda_*-2\times10^{-6},\ \lambda_*,
\lambda_*+2\times10^{-6},\ \lambda_*+0.01,\ 0.4,\ \pi/3,\ 2\pi/3,\
\pi,\ 4\pi/3,\ 5\pi/3,
\]
with \(\lambda_* =
0.1114802243779665542913031975274172717685818097497045\).
This includes four values strictly below the fold, the critical value, two
just above, and all six nonzero archived Phase 1 comparison points.

Every solve ran sequentially in WSL with Julia 1.12.6, the verified depot
`/home/adian/mub-depot`, and the pinned project manifest. The command template
was:

```text
JULIA_DEPOT_PATH=/home/adian/mub-depot \
  /home/adian/julia-1.12.6/bin/julia --project=. --compiled-modules=existing \
  - <lambda> <seed>
```

The inline Julia runner called
`MubSearch.generate_candidate_pool_fresh(H; tol=1e-8)`, deduplicated at
\(10^{-6}\), then called
`MubSearch.enumerate_all_third_mub_bases(pool,H; ortho_tol=1e-8,mu_tol=1e-8,
hp_bits=128)`. Each solve tracked 252 homotopy paths. All 13 solves and
clique reconstructions exited 0. The sum of recorded per-point elapsed times
was 997.283 seconds; peak RSS ranged from 1,305,038,848 to 1,421,021,184
bytes. Every listed clique passed the double-precision verifier and the
128-bit orthogonality/MU verification.

| \(\lambda\) | Expected pool | Observed pool | Expected verified cliques | Observed verified cliques | Peak RSS (bytes) | Seconds |
|---:|---:|---:|---:|---:|---:|---:|
| \(0\) | 120 | 120 | 10 | 10 | 1,319,084,032 | 95.832 |
| \(0.05\) | 120 | 120 | \(>4\) (qualitative below-fold expectation) | **4** | 1,313,808,384 | 62.936 |
| \(0.10\) | 120 | 120 | \(>4\) (qualitative below-fold expectation) | **4** | 1,354,272,768 | 64.115 |
| \(\lambda_*-2\times10^{-6}\) | 120 | 120 | \(>4\) (qualitative below-fold expectation) | **4** | 1,343,647,744 | 68.006 |
| \(\lambda_*\) | **96** | **72** | Not prespecified | 4 | 1,305,038,848 | 62.930 |
| \(\lambda_*+2\times10^{-6}\) | 72 | 72 | 4 | 4 | 1,324,617,728 | 68.495 |
| \(\lambda_*+0.01\) | 72 | 72 | 4 | 4 | 1,421,021,184 | 77.734 |
| \(0.4\) | 72 | 72 | 4 | 4 | 1,352,601,600 | 85.850 |
| \(\pi/3\) | 72 | 72 | 4 | 4 | 1,322,319,872 | 70.484 |
| \(2\pi/3\) | 72 | 72 | 4 | 4 | 1,387,294,720 | 85.945 |
| \(\pi\) | 120 | 120 | 10 | 10 | 1,378,947,072 | 84.255 |
| \(4\pi/3\) | 72 | 72 | 4 | 4 | 1,337,761,792 | 87.455 |
| \(5\pi/3\) | 72 | 72 | 4 | 4 | 1,347,563,520 | 83.246 |

The per-sample JSON also records the seed, raw/conjugate/verified solution
counts, Hadamard defect, individual clique verification errors, exit status,
pool vectors, and clique indices. The current run's critical pool contains 72 verified distinct vectors rather
than the expected 96. The 24-count difference matches the number of
singular candidates reported by the separate 1200-start phase-coordinate
search, but no vector matching was performed, so the missing candidates have
not been individually identified.

## Count-gate result and stop

The count gate **failed** for two reasons:

1. At \(\lambda_*\), the pinned Julia pool builder returned 72 vectors against
   the expected 96. This differs by the 24 numerical singular candidates
   from the earlier Python 1200-start result.
2. The requested below-fold expectation was that clique counts exceed 4.
   At \(0.05,\ 0.10,\ \lambda_*-2\times10^{-6}\), the regenerated pools have
   120 vectors but only 4 verified six-cliques. Only \(\lambda=0\) among the
   below-fold samples has 10.

The pool contains numerical recovered vectors, not a certified-complete
root set. The archived 400-start counts 114/93/70 remain incomplete and are
not treated as the true counts. Since the count gate failed, Phase 2.2 family
fits were **not run**, and no fit-based classification or family conclusion
is made for these fold samples.

### Initial collector failure (verbatim)

The first collector run emitted two rows, then exited 1 while decoding the
third WSL subprocess output:

```text
UnicodeDecodeError: 'charmap' codec can't decode byte 0x8f in position 374: character maps to <undefined>
AttributeError: 'NoneType' object has no attribute 'splitlines'
```

That attempt also incorrectly assigned an expected clique count of 10 at
\(\lambda=0.05\) based only on \(\lambda=0\). The corrected sweep used the
repository's canonical clique enumerator, UTF-8-safe output collection, and
left the exact clique expectation unspecified at previously unsampled
values. The below-fold \(>4\) qualitative expectation was separately tested
and failed at three points. No fit was run in either attempt.

## Tier and what is NOT established

- All sweep counts and fits (none were run) are **NUMERICAL**.
- No completeness certificate was produced for any recovered pool.
- The 24 critical singular candidates remain unverified singular roots.
- The below-fold clique count is not uniformly \(>4\) at the tested points.
- There is no family-fit result at these sweep parameters because the gate
  failed; no statement about all cliques fitting \(F,F^T\), or \(X_6\) follows.
- Nothing here establishes a full-interval or full-circle result, a
  fourth-MUB exclusion, or global MUB nonexistence.

## Gate status

Stop at **Gate 2**. Proposed next action, not performed: resolve whether the
critical 24 are part of the intended physical-root count using an approved
singular-root treatment, and establish the expected clique counts at the
newly sampled below-fold values before resuming any Phase 2.2 fits.
