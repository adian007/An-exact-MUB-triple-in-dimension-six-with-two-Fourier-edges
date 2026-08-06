# Methods: searching for four MUBs in dimension six via the Karlsson CHM family

*Updated 2026-08-03 after rigorous audit closure and certified targeted search.*

## Which question each phase tests

| Phase | Pair under test | Pool definition | What `found_fourth` means |
|-------|-----------------|-----------------|---------------------------|
| **Deprecated (wrong scope)** | Fixed `{I, F6}` | Vectors MU to I and F6, solved **once** | 6-clique also MU to H — tests `{I,F6}` extension (Grassl 2004 caps at 3; **do not use**) |
| **Current (correct)** | `{I, H(θ,φ,λ)}` per grid point | Vectors MU to I **and that H**, pool **re-solved per point** | Two 6-cliques among MU vectors → third and fourth MUB bases |

**Scope rule:** if the polynomial system does not depend on `H(θ,φ,λ)`, the sweep is not testing the Karlsson family.

---

## 1. Problem

Zauner's conjecture: at most three MUBs exist in \(\mathbb{C}^6\). We ask, for each Karlsson CHM \(H(\theta,\phi,\lambda)\):

> Does the pair `{I, H}` extend to **four** mutually unbiased bases?

Equivalently: find unit vectors unbiased to both the computational basis and the columns of \(H/\sqrt{6}\); among them find a **6-clique** (third ONB), then another **6-clique** among vectors also unbiased to that third basis (fourth ONB).

---

## 2. Per-H polynomial system (audited derivation)

Write \(v = z/\sqrt{6}\) with gauge \(z_0 = 1\).

**Computational unbiasedness:** \(|z_i| = 1\) → polynomial \(z_i w_i = 1\) with \(w_i \approx \overline{z_i}\).

**Unbiasedness to column \(k\) of \(H\):** with \(c_k = H_{:,k}/\sqrt{6}\),

\[
\Bigl|\sum_{j=0}^{5} \overline{H_{j,k}}\, z_j\Bigr|^2 = 6
\quad\Longleftrightarrow\quad
\Bigl(\sum_j \overline{H_{j,k}} z_j\Bigr)\Bigl(\sum_j H_{j,k} w_j\Bigr) = 6 .
\]

**Index convention:** column \(k\) uses **`H[j,k]`** (not `H[k,j]`). For \(F_6\) this matches \(\omega^{jk}\) because exponents commute; for general CHMs row/column differ — a row-index bug was fixed 2026-08-02.

Impose \(k = 0..4\) only; \(k=5\) follows on the conjugate locus by unitarity of \(H/\sqrt{6}\) (Parseval).

System: 10 equations, 10 unknowns \((z_{1..5}, w_{1..5})\). Coefficients depend on \(H\); **must be re-solved for each parameter point**.

---

## 3. Solvers

### 3.1 Fresh polyhedral solve (`generate_candidate_pool_fresh`)

Default reliability layer. HomotopyContinuation polyhedral start at fixed \(H\). ~25–40 s per point.

### 3.2 Parameter homotopy (`init_parametric_pool_tracker` + `track_pool_solutions`)

- **Anchor:** \(H_{\text{start}} = F_6\) (156 start paths; consistent root structure).
- **Track** start solutions to `params_from_hadamard(H) = [vec(conj(H)); vec(H)]`.
- **Drop detection:** if tracked count \(<\) start count, or homotopy pool \(<\) fresh pool, flag `drop_suspected`.
- **Fallback:** fresh solve replaces homotopy result when drops detected (`cross_check=true`).

Parameter homotopy alone is **not trusted**; path tracking loses solutions at many Karlsson points (observed: 83/156 paths at generic test point before fallback).

### 3.3 Certification (`certify_pool_at_H`)

HomotopyContinuation interval arithmetic certifies all tracked roots at fixed \(H\).
At F6: 156/156 certified, 156 distinct, 48 deduplicated pool vectors.

### 3.4 Pool completeness gate (`pool_completeness_report`)

Compares `n_tracked`, `n_certified`, `n_distinct_certified`, `n_verified`, `n_dedup_pool`,
and mixed volume. Sets `pool_complete_flag=true` when **`n_certified == n_tracked`**
(and `n_certified > 0`, `n_verified > 0`). Required before trusting negative clique results.

**512-point re-tag (2026-08-03):** all 512 generic grid points satisfy both the corrected
gate and the legacy 90%-of-mixed-volume threshold (512/512 each; 0 flips). See
`scripts/julia/retag_512_completeness.jl` → `results/karlsson_perH_sweep_retagged.csv`.

---

## 4. Extension test (`check_four_mub_extension`)

Input: pool of vectors already MU to `{I, H}`.

1. Build orthogonality graph (`ortho_tol`).
2. `found_third` := some clique has size \(\ge 6\).
3. For each 6-clique \(B_3\), filter pool to vectors MU to \(B_3\); search for another 6-clique → `found_fourth`.

**Julia 1.12 fix:** `argmax(f, cliques)` returns the maximizing clique directly (not an index).

Decoupled tolerances: `mu_tol`, `ortho_tol`. Sweeps report `found_fourth` at \(10^{-6}, 10^{-8}, 10^{-10}\).

---

## 5. Targeted search (`scripts/julia/search_special_loci.jl`)

**Do NOT** repeat the 8×8×8 uniform 512-point grid.

1. **Anchors:** F6, θ=0, Dita, Mobius audit-fail, λ=kπ/6.
2. **Degeneracy candidates:** from `scripts/python/degeneracy_scan.py` (Mobius z4_dev, A/B eigenvalue degeneracy, circulant match).
3. **Certification:** `certify_pool_at_H` + `pool_completeness_report` on **every** searched point.
4. **Adaptive refinement:** 5×5×5 local grid (scale 0.005–0.01) when max_clique≥6, pool count jumps, or cert/tracked mismatch.
5. **Counterexample escalation:** `scripts/julia/investigate_candidate.jl` triggered on `found_fourth=true`.

**Certified search set S:** points with `pool_complete_flag=true` and `n_certified==n_tracked`.

Output: `results/special_loci_search.csv` (UTF-8 direct write; do **not** use PowerShell `Tee-Object`).

**Refinement guard:** adaptive grids spawn only from **primary** queue points (`from_primary`), not from refinement hits — prevents nested `ref_ref_*` explosion.

Set `SPECIAL_LOCI_CSV` env var to write alternate output path.

---

## 5b. Third-MUB validation (`scripts/julia/validate_third_mub_candidates.jl`)

For each CSV row with `max_clique>=6`:

1. Locus distances (Fourier θ=0, Dita, Mobius z4\_dev).
2. BigFloat clique re-verify (`verify_clique_hp`, 256-bit).
3. Per-basis fourth-MUB report (`fourth_mub_per_basis_report`).
4. Parameter-space SVD rank (clustering).

Flags: `--exclude-deep-ref` (drop `ref_ref_*`), `--dedupe-ref` (0.005 cell representatives).

---

## 6. Validation and audit

| Check | Expected | Rigor |
|-------|----------|-------|
| **V1:** `H = F_6` | 156 raw, 48 pool, pool_complete | Certified + literature match |
| **V2:** Planted F6 6-clique | max_clique≥6 recovered | Pipeline audit PASS |
| **V3:** 512-point CSV | 512/512 max_clique=2 | Read-only re-check; grid missed special loci |
| **V4:** Dita / θ=0 | max_clique=6, found_fourth=false | Numerical anchor |

Run audit: `julia --project=. scripts/julia/audit_clique_pipeline.jl`

---

## 7. Symbolic elimination (`scripts/julia/symbolic_elimination.jl`)

Tractable subsystems only:
- θ=0 slice numerical scan (108 points)
- **Real** per-H pool export (`symbolic_export/pool_F6.m2`, `pool_Dita.m2`) — 10 equations with numeric H substituted
- Fourth-MUB witness equation skeleton (`build_fourth_mub_witness_equations`) → `fourth_mub_w_elimination_*.m2`
- Macaulay2 via Docker: `scripts/docker/run_m2.ps1 pool_F6.m2` (host path bind-mount)
- Oscar.jl optional for Groebner probes

Full family elimination is **not attempted**.

### 7.0 Oscar.jl install (optional, self-directed)

Oscar is **not** in `[deps]` by default (large binary). To enable `oscar_groebner_probe()` in Julia:

```powershell
cd "d:\IBM-Quantum\MUBs in 6-dimension"
julia --project=. -e "using Pkg; Pkg.add(\"Oscar\")"
```

First install may take 30–60 minutes (pulls `libpolymake_jll`, `Singular`, etc.). Verify:

```powershell
julia --project=. -e "using Oscar; println(Oscar.versioninfo())"
```

Tutorial entry point: [Oscar Groebner bases](https://docs.oscar-system.org/stable/CommutativeAlgebra/GroebnerBases/groebner_bases/) — load exported ideals from `symbolic_export/pool_F6.m2` via `evalfile` or manual transcription into `QQ`/`GF` rings.

If Oscar install fails, use Docker Macaulay2 only (`run_m2.ps1`); the 10-eq per-H pool ideal is the minimum symbolic content.

### 7.1 Witness decomposition (`scripts/julia/witness_decomposition.jl`)

Before Groebner/M2, probe the per-H pool system at F6 and Dita:

| Anchor | eqs×vars | mixed_volume | tracked/certified | inferred dimension |
|--------|----------|-------------|---------------------|-------------------|
| F6 | 10×10 | 252 | 156/156 | 0 (isolated) |
| Dita | 10×10 | 252 | 240/240 | 0 (isolated) |

Use NID/monodromy for irreducible factor degrees; local Groebner feasible when degree ≤300.
Full-family elimination remains deprioritized.

### 7.2 Third-MUB locus in parameter space (2026-08-03, Phase 1 resolved)

Scripts: `gauge_equivalence_analysis.jl`, `chm_equivalence.py`, `third_mub_locus_sampling.jl`,
`extended_axis_sweep.jl`, `fourth_mub_locus_sweep.jl`, `nid_probe.jl`.

**Gauge-equivalence methodology (Phase 1.1):**
1. Analytic: at θ=0, sin θ=0 ⇒ A independent of φ; H(0,φ,λ) entrywise identical.
2. CHM test: search D₁·H·P·D₂ over 720 column permutations + diagonal phases (+ transpose/conjugate).
3. Branch: φ-at-θ=0 equivalent (artifact); λ-at-Dita inequivalent (genuine family).

**A1 — Isotropic Gaussian** (100 samples/radius): 0/100 clique-6 at all radii (both anchors).

**A2 — Axis boundaries** (400-bit HP, binary search → `phase1_boundaries.txt`):

| Anchor | Axis | Boundary | Gauge? |
|--------|------|----------|--------|
| Dita | λ | **Full circle** [0,2π]: 126/126 (step 0.05), 7/7 spot checks | No (CHM-inequivalent) |
| Dita | 20×20 (φ,λ) grid | 0/400 clique≥6 (nearest Δφ=2.6×10⁻³ from π/4) | — |
| Dita | Δφ at λ∈{0.3,0.4,0.5} | drops at **10⁻³** | No |
| F6_theta0 | λ | ±0.055 / ±0.063 | No |
| F6_theta0 | φ | redundant | Yes |
| Both | θ | drops at 1e-6 | No |

**2D grid:** 20×20 at Dita θ fixed → `locus_2d_grid.csv` (400 points). **0/400** clique≥6 — grid omits exact φ=π/4. Targeted φ sweep (`locus_dimension_confirm.jl`) confirms 1D curve: clique≥6 only at φ=π/4; |Δφ|≥10⁻³ → clique=2.

**Item 1 — λ periodicity:** `lambda_periodicity_dita.jl` → `lambda_periodicity_dita.txt`. **VERDICT: full circle** (no arc boundary).

**Item 3 — CSV reconciliation:** `csv_reconciliation.jl` → `csv_reconciliation.txt`. 588 target vs 471 baseline: gap **117** = incomplete ref run + degen cap 80/1845 + dedup.

**Item 4 — Macaulay2:** Docker `unlhcc/macaulay2:1.24.05` ran stub export; no Groebner (`symbolic_elimination_attempt.txt`).

**A3 — NID (HC 2.22.1):** fixed-H pool at 7 locus λ values: dim=0 degree=1 each.
`witness_set(dim=0)` fails on square 10-var systems; use NID → witness_sets(N).

**Result:** Per-H MU-pool is 0-dimensional at each H on the Dita λ-curve.
Global third-MUB locus in (θ,φ,λ): 1D λ-family at Dita (inequivalent CHMs); φ redundant at θ=0.

Legacy `pool_complete_flag` (pre-2026-08-03) also required `n_tracked ≥ 0.9×mixed_volume`,
incorrectly excluding F6_theta0 (204/252 tracked but fully certified). Corrected gate:
`n_certified == n_tracked && n_verified > 0` only.

---

## 8. Software

| File | Purpose |
|------|---------|
| `src/mub_zauner_6d_liang_chen.jl` | Core pipeline + pool completeness |
| `scripts/julia/audit_clique_pipeline.jl` | Planted clique audit + 512 CSV parse |
| `scripts/python/degeneracy_scan.py` | Degeneracy candidate generator |
| `scripts/julia/search_special_loci.jl` | Certified targeted search |
| `scripts/julia/investigate_candidate.jl` | Counterexample protocol |
| `scripts/julia/retag_512_completeness.jl` | 512-point pool completeness re-tag |
| `scripts/julia/witness_decomposition.jl` | Witness/NID at F6 and Dita |
| `scripts/julia/third_mub_locus_sampling.jl` | Multi-radius/directional locus probe |
| `scripts/julia/extended_axis_sweep.jl` | Extended axis sweeps for drop points |
| `scripts/julia/nid_probe.jl` | HC NID/witness_set with error capture |
| `scripts/julia/third_mub_locus_dimension.jl` | Initial locus probe (superseded by above) |
| `scripts/julia/verify_quick_run_task1.jl` | Quick-run anchor re-certification |
| `scripts/julia/audit_argmax_task2.jl` | Argmax fix audit |
| `scripts/julia/lambda_periodicity_dita.jl` | Item 1: full-circle λ test at Dita |
| `scripts/julia/locus_dimension_confirm.jl` | Item 2: 1D vs 2D confirmation |
| `scripts/julia/locus_classification.jl` | Track D: HP clique-6 cluster catalog |
| `scripts/julia/dita_lambda_fourth_dense.jl` | Track A: dense Dita λ fourth test |
| `scripts/docker/run_m2.ps1` | Macaulay2 Docker with host mount |
| `docs/SEARCH_SETS.md` | S588, S*, region R definitions |
| `REPRODUCE.md` | One-command reproduction |
| `results/locus_classification.json` | Track D deliverable |

Julia 1.12.6, HomotopyContinuation 2.22.1, Graphs 1.14.0.

---

## 9. Rigor summary

| Claim | Status |
|-------|--------|
| F6 pool = 48 (156 certified distinct) | **Certified** at F6 anchor |
| Per-H pool at generic Karlsson | Heuristic fresh solve |
| `found_fourth = false` on certified set S | **Numerical** — not family-wide |
| "Karlsson never extends to four MUBs" | **Not proved** |

---

## 10. References

Karlsson (2011); McNulty–Weigert (2024) Secs. 7.1, 8.1; Grassl (2004); Brierley–Weigert (2009); Bengtsson et al. (2007); Haagerup / Björck–Fröberg (cyclic 6-roots).
