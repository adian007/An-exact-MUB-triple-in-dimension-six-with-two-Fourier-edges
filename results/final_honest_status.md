# Final Honest Status — Claim ledger (2026-08-05)

Search-set definitions: [`docs/SEARCH_SETS.md`](../docs/SEARCH_SETS.md). Reproduction: [`REPRODUCE.md`](../REPRODUCE.md). Finish audit: [`project_finish_audit.txt`](project_finish_audit.txt).

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
| **T1 Gauge structure (L1–L4)** | **Proved** | `paper/proofs/gauge_and_locus.tex`, `formalize_gauge_lemmas.jl` (2026-08-05 ALL PASS) |
| **T2 Third MUB on Dita λ-circle** | **HP-supported (126/126, 628/628); all-λ extension informal, not analytic** | `paper/proofs/dita_third_mub.tex`, 126/126 + 628/628; Route B continuity not rigorous |
| **T3 No fourth MUB on Dita circle** | **HP-verified (628/628); Claim not theorem** | `fourth_mub_obstruction.tex`; M2 CC probe + eliminate inconclusive |
| L5 1D Dita locus | **Proved (outline)** | `locus_geometry_probes.txt`, `locus_phi_sweep_full_circle.txt` |
| L6 F6 bounded λ-arc | **HP-verified** | `f6_boundary_reverify.txt` |
| L7 Per-H pool dim=0 at Dita | **HP-verified** | `pool_Dita_exact.m2`, `witness_decomposition.jl` (240 roots) |
| 512-grid max_clique=2 everywhere (sampled) | **HP-verified** | `karlsson_perH_sweep.csv` — **not** family-wide |
| Third MUB at F6, Dita, θ=0 | **HP-verified** | anchors, HP 400-bit |
| φ is gauge at θ=0 (H φ-independent) | **Proved** | analytic + `gauge_analysis.txt`, `chm_equivalence.txt` |
| λ at Dita: CHM-inequivalent 1D family | **HP-verified** | `chm_equivalence.txt`, Item 1 |
| Dita λ boundary | **HP-verified** | **FULL CIRCLE** [0,2π]; `lambda_periodicity_dita.txt`: 126/126 + 7/7 spot checks |
| F6 θ=0 λ boundary | **HP-verified** | arc width **0.117643**; `f6_boundary_reverify.txt` |
| 1D curve at Dita (9 λ values, φ=π/4) | **HP-verified** | `locus_phi_sweep_full_circle.txt`: 9/9 clique=6; 9/9 drop at \|Δφ\|≥10⁻³ |
| Dita vs F6_theta0 topology contrast | **Conjectured** | full λ-circle vs bounded arc |
| 2D grid clique-6 at Dita θ | **HP-verified** | **0/400** (grid misses φ=π/4) |
| φ drop off Dita curve | **HP-verified** | **\|Δφ\| ≥ 10⁻³** (9/9 λ) |
| Fourth MUB on locus interior (22 pts) | **HP-verified** | not found; `fourth_mub_locus_sweep.txt` |
| Fourth MUB at anchors | **HP-verified** | not found |
| Fixed-H NID dim=0 on locus λ values | **HP-verified** | `nid_probe.txt` — per-H only |
| degen_circulant_match (π/2) third MUB | **HP-verified** | **Refuted** under HP |
| Certified search CSV | **HP-verified** | **928 rows** (899 pool-complete); design target $\mathcal{S}_{588}$; 0/928 fourth |
| **S\* batch 1 (degen200)** | **HP-verified** | **219** primary (subset of degen500), 0 fourth; `special_loci_degen200.csv` |
| **S\* batch 2 (degen500)** | **HP-verified** | **519** primary (subset of degen1845), 0 fourth; `special_loci_degen500.csv` |
| **S\* full primary (degen1845)** | **HP-verified (log/meta)** | **1863/1865** distinct primary, 0 fourth; `special_loci_degen1845.meta.txt` + run log. **CSV corrupted** (592-row partial rerun); regenerate before submission. |
| **S\* mislabeled run (degen500_run2)** | **HP-verified duplicate** | Requested 1845, actual cap 500; identical to degen500; `special_loci_degen500_run2.csv` |
| HP audit of CSV clique-6 rows | **HP-verified** | **36/45** survive HP (9 θ=π/2 circulant_match fail); `third_mub_audit.csv` |
| Fourth MUB per basis (45 clique-6 rows) | **HP-verified** | not found (270/270 basis tests); `fourth_mub_per_basis.csv` |
| Macaulay2 Groebner elimination | **Open** | reduced witness M2-probed (35 eq, dim=-1 CC artifact); elimination not run |
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
| Proved | 2 | T1, L1 φ-gauge |
| HP-verified | 29 | incl. T2 (126/126, 628/628), T3 numerical (628/628), S\* degen1845 primary (1863/1865) |
| Conjectured | 2 | Dita vs F6 topology; witness ideal tractability |
| Open | 3 | Groebner certificate, R-theorem, global claims |

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
