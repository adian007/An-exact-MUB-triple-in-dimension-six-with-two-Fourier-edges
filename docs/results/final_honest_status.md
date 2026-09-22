# Final Honest Status — Claim ledger (2026-08-05)

Search-set definitions and reproduction: [`../overview/REPRODUCE.md`](../overview/REPRODUCE.md). Scientific interpretation: [`../SCIENTIFIC_STATUS.md`](../SCIENTIFIC_STATUS.md). Finish audit: [`../../results/project_finish_audit.txt`](../../results/project_finish_audit.txt).

## Ledger key

| Ledger | Meaning |
|--------|---------|
| **Proved** | Analytic or logical proof in-repo |
| **HP-verified** | 400-bit (or certified HC) numerical verification |
| **Conjectured** | Strong numerical/structural evidence; not a theorem |
| **Open** | Not established; gap or out of scope |

---

## Claims

| Claim | Ledger | Evidence |
|-------|--------|----------|
| Planted-clique pipeline recovers 6-cliques | **HP-verified** | `audit_clique_pipeline.jl`, T1/T2 tests |
| **T1 Gauge structure (L1–L4)** | **Proved** | `docs/paper/proofs/gauge_structure.tex`, `scripts/julia/formalize_gauge_lemmas.jl` (2026-08-05 ALL PASS) |
| **T2 Third MUB on sampled Dita λ** | **Finding (HP 126/126, 628/628); not a theorem** | `dita_third_mub.tex`; all-λ extension informal |
| **T3 No fourth MUB at seven Dita-λ points** | **Certified (2026-08-14)** — native all-clique W₁ empty **40/40** at Λ_cert; Bertini square + HC certify(); see `fourth_mub_obstruction.tex`, `certify_nwit1_all_cliques_*.txt` |
| **T4 Exact W₁ unit ideal at D₀-equivalent pair** | **Proved (exact Groebner, 2026-08-15)** — independent re-verification of BW 2009 (not an extension): GB={1} for single pair (D_bc, F_D); BW already had stronger complete-pool Np=0 (Nv=120, Nt=10); T4 weaker in coverage, stronger only as exact W₁ certificate form; λ∈{π/2,3π/2} via E0 numerical match; `track_c_elimination/w1_D0_groebner.log`, FINDINGS §17 |
| L5 1D Dita locus | **Finding (HP)** | `locus_geometry_probes.txt`, `locus_phi_sweep_full_circle.txt` |
| L6 F6 bounded λ-arc | **HP-verified** | `f6_boundary_reverify.txt` |
| L7 Per-H pool dim=0 at Dita | **HP-verified** | `pool_Dita_exact.m2`, `witness_decomposition.jl` (240 roots) |
| 512-grid max_clique=2 everywhere (sampled) | **HP-verified** | `karlsson_perH_sweep.csv` — **not** family-wide |
| Third MUB at F6, Dita, θ=0 | **HP-verified** | anchors, HP 400-bit |
| φ is gauge at θ=0 (H φ-independent) | **Proved** | analytic + `gauge_analysis.txt`, `chm_equivalence.txt` |
| λ at Dita: CHM-inequivalent 1D family | **HP-verified** | `chm_equivalence.txt`, Item 1 |
| Dita λ boundary | **HP-verified** | **FULL CIRCLE** [0,2π]; `lambda_periodicity_dita.txt`: 126/126 + 7/7 spot checks |
| F6 θ=0 λ boundary | **HP-verified** | arc width **0.117643**; `f6_boundary_reverify.txt` |
| 1D curve at Dita (9 λ values, φ=π/4) | **HP-verified** | `locus_phi_sweep_full_circle.txt`: 9/9 clique=6; 9/9 drop at \|Δφ\|≥10⁻³ |
| Dita vs F6_theta0 topology contrast | **HP-supported structural** | A-block geometry; `arc_circle_asymmetry.txt` |
| 2D grid clique-6 at Dita θ | **HP-verified** | **0/400** (grid misses φ=π/4) |
| φ drop off Dita curve | **HP-verified** | **\|Δφ\| ≥ 10⁻³** (9/9 λ) |
| Fourth MUB on locus interior (22 pts) | **HP-verified** | not found; `fourth_mub_locus_sweep.txt` |
| Fourth MUB at anchors / dense Dita circle | **HP-verified (628/628); not part of T3 certificate** | `dita_lambda_fourth_dense.csv` |
| Fixed-H NID dim=0 on locus λ values | **HP-verified** | `nid_probe.txt` — per-H only |
| degen_circulant_match (π/2) third MUB | **Refuted (tolerance artifact)** | `circulant_match_anomaly.txt` |
| Liang/Chen cross-family (Path C) | **Closed gap** | `docs/LIANG_CHEN_FAMILY_ASSESSMENT.md` |
| Certified search CSV | **HP-verified** | **928 rows** (899 pool-complete); design target $\mathcal{S}_{588}$; 0/928 fourth |
| **S\* batch 1 (degen200)** | **HP-verified** | **219** primary (subset of degen500), 0 fourth; `special_loci_degen200.csv` |
| **S\* batch 2 (degen500)** | **HP-verified** | **519** primary (subset of degen1845), 0 fourth; `special_loci_degen500.csv` |
| **S\* full primary (degen1845)** | **HP-verified** | **1863/1865** distinct primary, 0 fourth; `special_loci_degen1845.csv` regenerated 2026-08-13 (1863 rows, 21 clique-6, 1852 pool-complete) |
| **S\* mislabeled run (degen500_run2)** | **HP-verified duplicate** | Requested 1845, actual cap 500; identical to degen500; `special_loci_degen500_run2.csv` |
| HP audit of CSV clique-6 rows | **HP-verified** | **36/45** survive HP (9 θ=π/2 circulant_match fail); `third_mub_audit.csv` |
| Fourth MUB per basis (45 clique-6 rows) | **HP-verified** | not found (270/270 basis tests); `fourth_mub_per_basis.csv` |
| Macaulay2 Groebner elimination | **Done at one point (T4); open elsewhere** | exact n_wit=1 unit ideal at (D_bc, F_D) over Q(ζ₂₄,√5), 16 s (`w1_D0_groebner.log`); n_wit=2 (35 eq, dim=-1 CC artifact) and non-D₀ λ not run |
| B3 algebraic reconstruction Q(i,√2,√3) | **Inconclusive → LLL partial** | Nemo.lll: 6/36 trivial row-1 only (bound≤8); non-row-1 no candidate; `reconstruct_b3_algebraic.meta.txt` |
| Fourth-MUB witness HC certify (Dita) | **Skipped — not a certificate** | Reduced n_wit=2 already has mv=119210; `witness_mixed_volume_profile.txt` |
| Homotopy drop-detection (generic K6) | **HP-verified** | tracked=83/156, drop_suspected=true; `drop_detection_stress_test.txt` |
| Track D locus classification (HP components) | **HP-verified** | `locus_classification.json` — F6 arc (27 HP) + Dita circle (9 HP) only |
| Per-H pool M2 export (10 eq) | **HP-verified** | `symbolic_export/pool_F6.m2`, `pool_Dita.m2` |
| Fourth-MUB witness M2 export | **HP-verified export; elimination open** | 35-eq reduced witness; M2 dim=-1 inconclusive over CC |
| Dense Dita λ fourth test (Track A) | **HP-verified** | **628/628** clique≥6, fourth=0 |
| Reproducible CSV third hits | **HP-verified** | **2 primary** (F6_theta0, Dita) |
| No fourth MUB on region R | **Open** | S\* primary **1863/1865** checked (0 fourth); witness ideal + refinement remain |
| No fourth MUB in all of K₆⁽³⁾ | **Open** | — |
| N(6)=3 globally | **Open** | explicitly out of scope |

---

## Summary by ledger

| Ledger | Count | Notes |
|--------|-------|-------|
| Proved | 3 | T1, L1 φ-gauge, T4 (exact Groebner, D₀-equivalent pair, single clique) |
| HP-verified | 28 | incl. T2 (126/126, 628/628), S\* degen1845 primary (1863/1865) |
| Certified (T3) | 1 | seven Dita-λ points, 40/40 all-clique W₁ |
| Conjectured | 1 | analytic κ(H(λ)) formula |
| Open | 3 | Groebner certificate at non-D₀ λ, R-theorem, global claims |

---

## STOP status

**Pipeline NOT stopped.** `found_fourth=false` at all completed tests (928 + 1863 S\* primary + 628 audited).

## Phase branch

**Phase 5 (2026-08-05)** — finish audit, extended HP validation (36/45), S\* batch 1 (degen200), Julia depot repair, publication prep.

## Key outputs (Phases 0–5)

- `results/project_finish_audit.txt` — consolidated finish-line stats
- `results/validation_summary.txt` — 928 rows, 36/45 HP clique-6
- `results/locus_classification.json` — Track D (27 F6 + 9 Dita HP rows)
- `results/special_loci_degen200.csv` — S\* batch 1 (219 primary, subset of degen500)
- `results/special_loci_degen1845.meta.txt` — S\* full primary (1863/1865; log-backed; CSV must be regenerated)
- `results/special_loci_degen1845_CORRUPTED592.csv` — accidental partial overwrite; do not cite
- `results/special_loci_degen500_run2.csv` — mislabeled duplicate (actual cap 500)
- `results/formalize_gauge_lemmas.txt` — T1 re-verified 2026-08-05

---

## 18. Methodology hardening and regression benchmarks (2026-09-14/15)

Audit and fixes: [`docs/METHODOLOGY_AUDIT.md`](../docs/METHODOLOGY_AUDIT.md) (call graph, bug ledger F1–F9). New module: `src/MubSearch.jl` (Karlsson variants + validator, PoolAudit, projective dedup, exhaustive clique enumeration, W1 witness, claim tiers, provenance). Test suite: `test/runtests.jl` (101 fast tests + long benchmarks, gated by `MUB_LONG_TESTS`).

| Claim | Ledger | Evidence |
|-------|--------|----------|
| Karlsson original vs McNulty–Weigert review transcription (C1) | **HP-verified (regression test)** | MW literal transcription fails AA†=2I at 40/40 random points; original passes everywhere; `test/runtests.jl` |
| Fourier-seam resolution (θ=0 exact algebraic slice; z2² = α_A/β_A) | **Proved (algebraic) + regression-tested** | seam identity derivation in `build_karlsson_family`; backward-compatible ~1e-16 with historical anchors |
| Möbius identities in denominator-cleared form; Dita-anchor 0/0 flagged | **Proved (algebraic)** | `karlsson_mobius_audit`, `mobius_identity_residual`; legacy z4_dev ≈ 1.84 at Dita is a conditioning artifact, not a family inconsistency |
| NaN/pole fail-fast in family construction | **HP-verified** | pole point (1.0, 2.0, 0.7) now throws with diagnosis instead of returning NaN-H |
| **Benchmark A: F6 pool = 48 physical MU vectors** | **CERTIFIED_NUMERICAL** | `results/benchmarks/benchmarks.json`; mv 252 → 156 sols → 156 certified → 48 physical → 48 dedup; reseed crosscheck agrees |
| **Benchmark B: D₀-equivalent λ=π/2 pool = 120 vectors, 10 third bases, N_p = 0, W₁ empty 10/10** | **CERTIFIED_NUMERICAL** | independent pipeline reproduction of BW 2009 (N_v, N_t, N_p, per-clique W₁); all audit counts explained; crosscheck agrees |
| Constellation-size variation on the Dita circle (120 near λ∈{0,π/2,π,3π/2} mod 2π; 72 elsewhere; 4 vs 10 third bases) | **NUMERICALLY_SUPPORTED** | new runs + historical `dita_lambda_fourth_dense.csv` (89/628 rows at 120, all in four ≈0.22-rad intervals) |
| Benchmark C: Tao S₆ (isolated) pool = 90 vectors, no third MUB | **CERTIFIED_NUMERICAL** | first pool count for the isolated matrix; transcription self-validating (Wuttig–Tindall Eq. 10) |
| Dedup stability under tolerance sweep 1e-6…1e-12 | **HP-verified** | `dedup_tolerance_diagnostics`: counts 48/120/72/90 invariant at all four benchmark points |
| Claim-tier vocabulary + provenance on all new artifacts | **Enforced** | `Certification.jl` tiers; `write_result` stamps git commit/seed/tolerances |

Environment note (F4): Windows Smart App Control (On since ~2026-09) blocks unsigned Julia core DLLs; all Julia computation now runs in WSL2 Ubuntu (documented in `../overview/REPRODUCE.md`). No security policy was modified.
