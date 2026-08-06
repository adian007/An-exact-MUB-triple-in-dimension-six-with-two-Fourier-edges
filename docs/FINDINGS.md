# Rigorous audit findings (2026-08-03, theorem path 2026-08-05)

**Theorem-path artifacts (T1–T3):**
- Proofs: [`paper/proofs/gauge_and_locus.tex`](../paper/proofs/gauge_and_locus.tex), [`dita_third_mub.tex`](../paper/proofs/dita_third_mub.tex), [`fourth_mub_obstruction.tex`](../paper/proofs/fourth_mub_obstruction.tex)
- Scripts: `formalize_gauge_lemmas.jl`, `verify_dita_construction.jl`, `src/dita_third_mub_construction.jl`
- Results: `formalize_gauge_lemmas.txt`, `verify_dita_construction.txt`, `locus_geometry_probes.txt`

**Search sets and reproduction (Phase 0 baseline):**

- Search-set definitions (**S588**, **S\***, covering region **R**): [`docs/SEARCH_SETS.md`](SEARCH_SETS.md)
- One-command regeneration and environment setup: [`REPRODUCE.md`](../REPRODUCE.md) (Julia 1.12.6, project-local `.julia-depot`)
- Claim ledger (Proved / HP-verified / Conjectured / Open): [`results/final_honest_status.md`](../results/final_honest_status.md)

This document supersedes the conclusion in `paper/mub6_karlsson_ieee.tex` and the
512-point uniform sweep. **Do not cite "max clique 2 everywhere" as a Karlsson-family
impossibility result.**

## 1. Pipeline bug (fixed)

**Julia 1.12 `argmax(f, collection)` returns the maximizing *element*, not its index.**

In `check_four_mub_extension`, the branch for `found_third && !found_fourth` used
`cliques[argmax(...)]`, which crashed when a 6-clique was detected. The prior 512-point
sweep never hit this path because `found_third=false` at every grid point — but anchor
points (F6, Dita, θ=0) *do* produce 6-cliques.

**Fix:** use `best_c = argmax(c -> length(c), cliques)` directly (no indexing).

## 2. Planted-clique tests (PASS)

| Test | Description | Result |
|------|-------------|--------|
| T1 | Synthetic ONB + decoys | **PASS** at ε_ortho ∈ {10⁻⁶…10⁻¹²} |
| T2 | F6 pool 6-clique planted + decoys | **PASS** at ε_ortho, ε_mu ∈ {10⁻⁶,10⁻⁸,10⁻¹⁰} |
| T2b | Single clique → found_fourth=false | **PASS** (negative control) |
| T2c | Random two-ONB (not MU to H) → found_fourth=false | **PASS** |

**Conclusion:** clique-finding is not systematically missing 6-cliques.

## 3. 512-point sweep re-check (read-only CSV parse)

Parsed `results/karlsson_perH_sweep.csv` without re-running the sweep:

| Statistic | Value |
|-----------|-------|
| Points sampled | **512/512** |
| max_clique | **2 at all 512 points** |
| found_third | **false at all 512** |
| found_fourth | **false at all 512** |
| n_pool | 48, 52, 56, or 72 |

**Grid limitations:** θ, φ ≥ 0.02 only — missed θ=0 Fourier subfamily, Dita point
(θ ≈ arccos(1/√3) ≈ 0.955), and exact symmetry planes.

**Verdict:** 512/512 sampled points show max_clique=2; **NOT family-wide.**

### 3.1 Corrected pool completeness re-tag (Task 2, 2026-08-03)

Re-certified all 512 grid points via `scripts/julia/retag_512_completeness.jl`
(fresh `pool_completeness_report` per point; output `results/karlsson_perH_sweep_retagged.csv`):

| Metric | Count |
|--------|-------|
| OLD complete (tracked ≥ 0.9×mixed_volume) | **512/512** |
| NEW complete (`n_certified == n_tracked`) | **512/512** |
| Flip complete→incomplete | **0** |
| Flip incomplete→complete | **0** |
| max_clique ≥ 3 (original sweep clique) | **0** |

Fresh polyhedral solves track 244 or 252 roots (mixed_volume=252); all certify fully.
The corrected `pool_complete_flag` (`n_certified==n_tracked`) agrees with the old
90%-of-mixed-volume gate on this generic grid — no trust downgrade required.

## 3b. Mobius degeneracy candidate (Task 1, 2026-08-03)

From `results/degeneracy_candidates.csv`: **exactly one** point with `z4_dev > 10⁻⁶`:

| Field | Value |
|-------|-------|
| θ | **0.9553166181** |
| φ | **0.7853981634** (= π/4) |
| λ | **0.4** |
| z4_dev | **1.345×10⁻⁵** |
| flags | `mobius_z4` |

This is the **Dita anchor** (θ = arccos(1/√3) within 2.5×10⁻¹¹). `investigate_candidate.jl`
had **not** been run on this point; investigation executed via `scripts/julia/run_task1_investigate.jl`
(`investigate_candidate.jl` blocked by CSV precompile policy on this host).

| Test | CSV rounded coords | Exact Dita coords |
|------|-------------------|-------------------|
| pool_complete | **true** (240/240 certified) | true |
| max_clique | **2** | **6** |
| found_third | **false** | **true** |
| found_fourth | **false** | **false** |
| HP verify (400-bit BigFloat) | N/A (no 6-clique) | third clique certifiable |

**Conclusion:** Rounded degeneracy-scan coordinates sit on a **sharp third-MUB locus** —
10-digit θ truncation destroys the 6-clique. No fourth MUB at either coordinate set.
PSLQ: θ ≈ arccos(1/√3), φ = π/4.

## 3c. Witness decomposition at F6 and Dita (Task 3, 2026-08-03)

`scripts/julia/witness_decomposition.jl` on the per-H MU pool system (10 eqs, 10 vars):

| Anchor | mixed_volume | tracked | certified | Dimension |
|--------|-------------|---------|-----------|-----------|
| F6 (Fourier matrix) | 252 | 156 | 156 | **0** (isolated roots) |
| Dita Karlsson | 252 | 240 | 240 | **0** (isolated roots) |

HomotopyContinuation `witness_set` / full NID hit API issues on this host; finite
solve counts establish **zero-dimensional** components. Degrees ≥156 (F6) and ≥240 (Dita);
BKK bound 252. **Recommend:** local Groebner at F6/Dita subsystems (manageable degree);
deprioritize full-family symbolic elimination.

## 4. Pool completeness (F6 anchor)

`pool_completeness_report(F6)` at the F6 anchor:

| Metric | Value |
|--------|-------|
| mixed_volume | 252 |
| n_tracked / n_certified | 156 / 156 |
| n_distinct_certified | 156 |
| n_conjugate → n_verified → n_dedup_pool | 156 → 156 → **48** |
| pool_complete_flag | **true** |

The 156→48 gap is **expected** (conjugate locus + MU filter + dedup), not incomplete tracking.

## 5. Certified search set S

**Definition:** S = parameter points with `pool_complete_flag == true` and
`n_certified == n_tracked` from `search_special_loci.jl`.

Search covers:
- Anchor loci: F6, θ=0, Dita, Mobius audit-fail, λ=kπ/6 grid
- Degeneracy candidates from `scripts/python/degeneracy_scan.py` (Mobius z4_dev, A/B eigenvalue degeneracy, circulant match)
- Adaptive 5×5×5 refinement when max_clique≥6 or pool count jumps

Output: `results/special_loci_search.csv`

### 5.1 CSV pipeline fix (Step 0, 2026-08-03)

| Item | Status |
|------|--------|
| UTF-16 corruption from PowerShell `Tee-Object` | **Fixed** — Julia writes UTF-8 directly (`open(...,"w")`, no Tee-Object) |
| `CSV.jl` dependency | Optional; manual line verify in `search_special_loci.jl` |
| Current file | **848 rows** in `special_loci_search_backup848.csv` (UTF-8; includes nested `ref_ref_*` from pre-guard run) |
| Current file | **491 rows** (471 baseline + partial resume; verify-only vs 588 target **FAIL**) |
| Interim validation set | `special_loci_search_interim.csv` (347 rows, 16 third-MUB hits) |

Note: 848-row file includes one-level `ref_*` refinements from an earlier run without the `from_primary` guard; exclude `ref_ref_*` for analysis.

### 5.2 Third-MUB candidate audit (Steps 0–5, 2026-08-03)

**Step 0 — full search re-run:** `scripts/julia/search_special_loci.jl` (UTF-8 direct write; no `Tee-Object`). Fixed Julia 1.12 bugs (`open` encoding kw, `_row_to_dict` NamedTuple access). Baseline **471 data rows** (99 primary + 372 `ref_*`; `--verify-only` vs 588 target **FAIL** — 117-row gap documented in `results/csv_reconciliation.txt`). Current `special_loci_search.csv`: **491 rows** after partial resume append. Original 588-point Tee-Object log is **not recoverable** as clean CSV.

**Steps 1–3 — `validate_third_mub_candidates.jl` on clean CSV:**

| Metric | Result |
|--------|--------|
| CSV rows loaded | **491** |
| `max_clique≥6` candidates | **20** |
| HP clique survival (256-bit BigFloat) | **15/20 true** |
| HP failures | **5/5** `degen_circulant_match` / `ref_degen_circulant_match` at θ=π/2 (max_clique=6 in sweep, HP re-verify **false**) |
| Verified third-MUB loci | **F6_theta0** (θ=0 slice) + **Dita** (θ=arccos(1/√3), φ=π/4) + refinements thereof |
| `ortho_max_hp` on survivors | ≤ 4×10⁻¹⁶ |
| Exact params (symbolic) | θ=0, θ=arccos(1/√3), φ=π/4 |
| Fourth MUB per basis | **n_filtered=0**, `found_fourth=false` for **every** basis at **every** candidate (20/20) |
| Aggregate `found_fourth` | **false** (20/20) |
| Parameter-space SVD rank (20 candidates) | **3** (σ = 11.74, 1.76, **0.012** — **aggregate** envelope over F6+Dita loci; not per-locus dimension) |
| Tangent perturbation (ε=10⁻⁴, 5 samples) | **0/5** preserve third-MUB (max_clique→2) |

Outputs: `results/third_mub_audit.csv`, `results/fourth_mub_per_basis.csv`, `results/validation_summary.txt`

**Step 4 — advanced (`investigate_candidate.jl` on HP-confirmed anchors):**

| Anchor | PSLQ | Monodromy | Fourth |
|--------|------|-----------|--------|
| F6_theta0 (0, 0.5, 0.3) | θ=0 [Fourier slice] | **Failed** — `MonodromyOptions` API mismatch (HC v2.22) | false |
| Dita (0.9553…, π/4, 0.4) | θ=arccos(1/√3), φ=π/4 | **Failed** (same) | false |

Logs: `results/investigate_F6_theta0.log`, `results/investigate_Dita.log`. Oscar.jl: **not installed**. Macaulay2 export: stub only (§7).

**Step 5 — verified vs heuristic:**

| Status | Points |
|--------|--------|
| **Verified** third MUB (HP) | F6_theta0 family (12), Dita family (4) = **15** |
| **Heuristic only** (sweep max_clique=6, HP fail) | circulant_match at θ=π/2 (**5**) — treat as **non-reproducing** (consistent with §5.3 post-audit note) |

Per-basis fourth test confirms aggregate `found_fourth=false` is not an aggregate-only artifact.

### 5.3 Quick-run summary (20 anchors) — 18-vs-17 resolved (2026-08-03)

| Result | Count |
|--------|-------|
| Points searched | **20** |
| Hadamard-valid (`status=ok`) | **18** |
| `pool_complete_flag=true` (corrected gate) | **18/18** |
| `pool_complete_flag=true` (legacy gate, quick run log) | **17** |
| `found_third=true` | **2** (F6_theta0, Dita) |
| `found_fourth=true` | **0** |
| `max_clique=6` | **2** |
| `max_clique=2` | **16** |

Note: F6_cyclic (0,0,0) and theta_0 fail Hadamard check; Dita and Mobius audit-fail deduplicate.

**The 17th-vs-18th discrepancy:** exactly one Hadamard-valid point had `pool_complete=false`
under the **legacy** gate used at quick-run time (10:52 AM, `special_loci_run_log.txt`):

| Field | Value |
|-------|-------|
| Point | **F6_theta0** (θ=0 subfamily) |
| (θ, φ, λ) | **(0.0, 0.5, 0.3)** |
| n_tracked / n_certified | **204 / 204** |
| mixed_volume / mv_gap | 252 / 48 |
| max_clique / found_third | **6 / true** |
| Legacy gate fail | 204 < 0.9×252 = 226.8 |
| Corrected gate | **true** (cert_match only) |

Fresh polyhedral solve at exact point: 204/204 certified, max_clique=6 — **not incomplete
tracking**, legacy mixed-volume threshold bug. Fix in `src/mub_zauner_6d_liang_chen.jl` at
10:54 AM. Re-certification: `scripts/julia/verify_quick_run_task1.jl`.

**Post-audit CSV note (2026-08-03 12:32):** `special_loci_search.csv` grew beyond the
20-anchor quick run (degeneracy-candidate append). One appended row reported
`max_clique=6` at (θ,φ,λ)=(π/2, 0.05, 6.2331853072). Fresh post-fix re-solve
(`scripts/julia/verify_extra_third_hit.jl`): **max_clique=2**, found_third=**false**,
192/192 certified — **does not reproduce**; treat quick-run third-MUB positives as
**F6_theta0 + Dita only**.

### 5.4 Pre/post argmax fix audit (Task 2, 2026-08-03)

`scripts/julia/audit_argmax_task2.jl` → `results/task2_argmax_audit.txt`

| Event | Timestamp |
|-------|-----------|
| Quick run (`special_loci_run_log.txt`) | 2026-08-03 **10:52 AM** (pre-fix) |
| argmax fix (`src/mub_zauner_6d_liang_chen.jl`) | 2026-08-03 **10:54 AM** |

`check_four_mub_extension` has **no try/catch**. Julia 1.12 `argmax` bug affects only the
`ortho_defect` branch after `found_fourth` is set — **not a silent false-negative**.

| Point | found_third | found_fourth (pre-fix CSV) | found_fourth (post-fix) | Re-run? |
|-------|-------------|---------------------------|------------------------|---------|
| F6_theta0 | true | false | false | No |
| Dita | true | false | false | No |
| F6 (Fourier) | true | — | false | No |

**Post-audit check:** appended CSV row (π/2, 0.05, 6.2331853072) reported
max_clique=6 pre-fix; post-fix re-solve → max_clique=**2**, found_third=**false**
(does not reproduce; not a trusted third-MUB point).

### 5.5 Third-MUB locus dimension — Part A sampling (2026-08-03)

Scripts: `third_mub_locus_sampling.jl`, `extended_axis_sweep.jl`, `nid_probe.jl`
Outputs: `results/third_mub_locus_sampling.txt`, `results/extended_axis_sweep.txt`, `results/nid_probe.txt`, `results/priority_hits_verify.txt`

#### A1 — Multi-radius isotropic Gaussian (100 samples/radius)

| Anchor | radius | n_clique6 / 100 | n_clique2 / 100 |
|--------|--------|-----------------|-----------------|
| Dita | 1e-6, 1e-5, 1e-4, 1e-3, 1e-2 | **0** each | **100** each |
| F6_theta0 | 1e-6, 1e-5, 1e-4, 1e-3, 1e-2 | **0** each (1 other at 1e-6) | **99–100** each |

Isotropic perturbations **miss** axis-aligned third-MUB persistence (thin/curved components).

#### A2 — Directional sweeps (RESOLVED by Phase 1 gauge analysis)

**High-precision re-verify** (`verify_priority_hits.jl`, 400-bit HP): axis hits at |δ|=0.01 are genuine clique-6 where H differs.

| Anchor | Direction | max_clique=6 persists to | Gauge? |
|--------|-----------|---------------------------|--------|
| Dita | **λ** (θ,φ fixed) | **full circle** — 126/126 on [0,2π], 7/7 spot checks | **No** — CHM-inequivalent λ |
| Dita | θ | drops at **≤1e-6** | No |
| Dita | φ | drops at **≤1e-6** | No |
| F6_theta0 | **φ** | N/A (same H entrywise) | **Yes** — θ=0 makes A φ-independent |
| F6_theta0 | **λ** | **\|Δλ\|=0.055044803035 (+), 0.0625984193874 (−)**; arc **0.117643** | No |
| F6_theta0 | θ | drops at **≤1e-6** | No |

**Phase 1 resolution (2026-08-03, Items 1–2; closeout Items A–B):** NOT "positive-dimensional surface" in (φ,λ) at F6_theta0 — φ is redundant.
At Dita, **genuine 1D λ-curve** at φ=π/4 exact: full-circle λ periodicity (Item 1); 20×20 grid 0/400 clique≥6 (misses φ=π/4 by 2.6×10⁻³); φ-offset sweep **9/9** at λ∈{0.3,…,6.0} — clique=6 only at φ=π/4, drop at |Δφ|≥10⁻³ (400-bit HP). **Topology contrast:** Dita supports third MUB on a **full 2π λ-circle**; F6_theta0 on a **bounded arc** of width 0.117643 rad — symmetry-specific, not universal (see paper Discussion).

Scripts: `lambda_periodicity_dita.jl`, `locus_dimension_confirm.jl`, `locus_phi_sweep_full_circle.jl`, `f6_boundary_reverify.jl`, `gauge_equivalence_analysis.jl`, `chm_equivalence.py`, `extended_axis_sweep.jl`.
Outputs: `results/lambda_periodicity_dita.txt`, `results/locus_dimension_confirm.txt`, `results/locus_phi_sweep_full_circle.txt`, `results/f6_boundary_reverify.txt`, `results/locus_2d_grid.csv`, `results/gauge_analysis.txt`, `results/chm_equivalence.txt`, `results/phase1_resolution.txt`.

#### A3 — NID / witness sets (HomotopyContinuation v2.22.1)

| System | Result |
|--------|--------|
| Fixed-H pool at Dita | **NID succeeded:** ~240 components, each **dim=0, degree=1** (certified isolated roots at fixed H) |
| `witness_set(dim=0)` | **Fails:** `ArgumentError: codim has to be between 0 and n` (API/usage issue, not host crash) |
| Parametric (hc,h) pool | **Fails:** parameters not specialized — `Missing: h₁₋₁, …` |
| (θ,φ,λ) third-MUB NID | **Not implemented** — clique≥6 is combinatorial, not polynomial |

**Certified:** fixed-H MU-pool variety is **0-dimensional** at each H.
**Not certified:** global dimension of third-MUB locus in (θ,φ,λ) by NID; combinatorial probes confirm a **1D λ-curve** at Dita (φ=π/4 exact) and bounded λ interval at F6_theta0 (Items 1–2).

**Do not conflate** these two claims.

## 6. Fourth MUB search (no counterexample)

At all points tested (anchors, degeneracy candidates, refinements):

- `found_fourth = false` everywhere
- Third MUB (`max_clique=6`, `found_third=true`) at F6, Dita, θ=0 subfamily
- Generic Karlsson points: `max_clique=2`

**This is not a proof** that no Karlsson matrix extends to four MUBs.

## 7. Symbolic elimination (Step 4)

| Subsystem | Status |
|-----------|--------|
| θ=0 slice numerical scan (36 points) | **No found_fourth** (6/36 third) |
| Local fourth-MUB witness at F6/Dita/θ0 | Sized (60 vars, 105 eqs); **no Groebner run** |
| Macaulay2 export | `symbolic_export/fourth_mub_w_elimination_*.m2` — **stub only** (header, no equations) |
| Macaulay2 Docker (unlhcc/macaulay2:1.24.05) | **Attempted** — script exits 0; no computation (stub) |
| Oscar.jl | **Not installed** |
| Monodromy / PSLQ / Clifford | PSLQ **run** on F6_theta0 + Dita; monodromy **API fail** (HC v2.22 MonodromyOptions) |

**No algebraic proof** of empty fourth-MUB ideal was obtained. See `results/symbolic_elimination_attempt.txt`.

## 8. What is actually established

| Claim | Status |
|-------|--------|
| Clique code recovers planted 6-cliques (synthetic + F6 pool) | **Verified** |
| Third MUB exists for some Karlsson H (F6, Dita, θ=0) | **Verified numerically** |
| Dita λ-curve: CHM-inequivalent family with third MUB | **Verified** — **full circle** on [0,2π] (126/126) |
| 2D (φ,λ) patch at Dita θ | **1D curve** — 0/400 grid; φ drop at 10⁻³ (**9/9** λ on full circle) |
| Dita vs F6 λ-topology | **Full circle vs bounded arc** (0.117643 rad) | `locus_phi_sweep_full_circle.txt`, `f6_boundary_reverify.txt`, paper Discussion |
| φ at θ=0 is gauge (parametrization artifact) | **Verified analytically + CHM test** |
| Fourth MUB at any tested point | **Not found** |
| degen_circulant_match (π/2, 0.05, 6.233…) | **Non-reproducing** (max_clique=2) |
| Macaulay2 Groebner elimination | **Not performed** (Docker stub export only) |
| No fourth MUB in entire Karlsson family | **NOT proved** |
| 512-point sweep rules out all third MUBs | **FALSE** (missed special loci) |
| Certified root count at F6 | **156 distinct, pool_complete** |
| CSV complete at 588 rows | **FALSE** (471 baseline → **491** after partial resume; gap=117 — see `csv_reconciliation.txt`) |
| HP audit of CSV clique-6 hits | **15/20 verified** (5 θ=π/2 circulant_match fail HP; `third_mub_audit.csv`) |

## 9. Final honest statement

> No fourth MUB was found on certified set S ⊂ K₆⁽³⁾ (**491** CSV rows after partial resume from 471 baseline: anchors +
> partial degeneracy load + one-level refinements; 588 target not met). HP audit: **15/20** clique-6 CSV rows verified;
> **5/5** θ=π/2 circulant_match rows fail HP. Dita third-MUB locus at φ=π/4 is a
> **full-circle λ family** (126/126 on [0,2π], 400-bit HP). Reproducible third-MUB loci:
> F6 at θ=0, Dita, and a **genuine 1D λ-family** of CHM-inequivalent matrices at
> Dita (φ at θ=0 is a gauge coordinate). **λ-topology contrast:** Dita third-MUB locus is a
> **full 2π circle**; F6_theta0 locus is a **bounded arc** (width 0.117643 rad). The π/2 circulant-match CSV row does not
> reproduce. This does **not** prove N(6)=3 for all Karlsson matrices.

## Scripts

| File | Purpose |
|------|---------|
| `scripts/julia/audit_clique_pipeline.jl` | Planted clique + F6 pool + 512 CSV re-check |
| `scripts/python/degeneracy_scan.py` | Emit Mobius/A-B/circulant degeneracy candidates |
| `scripts/julia/search_special_loci.jl` | Certified targeted search → `special_loci_search.csv` |
| `scripts/julia/investigate_candidate.jl` | BigFloat re-certify + PSLQ on found_fourth hits |
| `scripts/julia/retag_512_completeness.jl` | Re-tag 512 sweep with corrected pool_complete_flag |
| `scripts/julia/witness_decomposition.jl` | NID / witness-set probe at F6 and Dita |
| `scripts/julia/verify_quick_run_task1.jl` | Re-certify 20-point quick-run anchors |
| `scripts/julia/audit_argmax_task2.jl` | Pre/post argmax fix audit for found_fourth |
| `scripts/julia/third_mub_locus_dimension.jl` | Parameter-space third-MUB locus dimension |
| `scripts/julia/run_task1_mobius.jl` | Task-1 Mobius candidate investigation |
| `scripts/julia/validate_third_mub_candidates.jl` | Steps 1–3: HP audit, per-basis fourth test, clustering |
| `scripts/julia/lambda_periodicity_dita.jl` | Item 1: full-circle λ periodicity at Dita |
| `scripts/julia/locus_dimension_confirm.jl` | Item 2: 1D vs 2D locus confirmation |
| `scripts/julia/locus_phi_sweep_full_circle.jl` | Item A: φ-offset sweep at 9 λ values (full circle) |
| `scripts/julia/f6_boundary_reverify.jl` | Item B: F6 θ=0 λ boundary re-verify (400-bit HP) |
| `scripts/julia/csv_reconciliation.jl` | Item 3: 588 vs 471 CSV gap |
| `scripts/julia/locus_classification.jl` | Track D: HP clique-6 classification |
| `scripts/julia/dita_lambda_fourth_dense.jl` | Track A: dense Dita λ fourth test |
| `scripts/julia/symbolic_elimination.jl` | Real pool + witness M2 export |
| `scripts/docker/run_m2.ps1` | Macaulay2 Docker (host mount) |
| `docs/SEARCH_SETS.md` | S588, S*, region R |
| `REPRODUCE.md` | One-command reproduction |
| `results/locus_classification.json` | Track D deliverable |
| `pool_completeness_report` | In `src/mub_zauner_6d_liang_chen.jl` |

---

## 10. Breakthrough plan artifacts (Phases 0–4, 2026-08-03)

| Phase | Deliverable | Status |
|-------|-------------|--------|
| 0 | `REPRODUCE.md`, `docs/SEARCH_SETS.md`, claim ledger | **Done** |
| 1A | `search_special_loci.jl --resume` toward 588 rows | **Done** (588 rows; `--verify-only` PASS; trimmed from 891 overshoot archive) |
| 1B | `--degen-cap`, `-z4_dev` sort, `--no-refine` | **Done** |
| 1C | `dita_lambda_fourth_dense.jl` | **Done** (628/628 clique≥6, fourth=0; `dita_lambda_fourth_dense.csv`) |
| 2 | `locus_classification.json` | **Done** — F6 arc + Dita circle only |
| 3 | `pool_F6.m2`, `pool_Dita.m2` (real); witness `.m2`; `run_m2.ps1` | **Pool real; Groebner open** |
| 4 | Paper §Track C/D + infrastructure; `FINDINGS.md` | **Done** |

Track D numerical criterion: the only HP-verified third-MUB loci in the current CSV are (up to CHM class) **F6\_theta0 bounded λ-arc** and **Dita full λ-circle** at φ=π/4.
