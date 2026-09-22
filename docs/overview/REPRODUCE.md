# Reproduction guide (Phase 0 baseline)

One-command instructions to regenerate key tables and audit artifacts for the Karlsson \(K_6^{(3)}\) MUB search. Search-set definitions: [`REPRODUCE.md`](REPRODUCE.md). Claim ledger: [`../results/final_honest_status.md`](../results/final_honest_status.md). Pipeline bugs: [`docs/ENGINEERING.md`](docs/ENGINEERING.md).

---

## Environment (Windows)

From the project root (`d:\MUBs in 6-dimension`):

```powershell
. .\setup_julia_env.ps1
```

This sets:

| Variable | Value |
|----------|-------|
| Julia binary | `C:\Users\adian\AppData\Local\Programs\Julia-1.12.6\bin\julia.exe` |
| `JULIA_DEPOT_PATH` | `<project-root>\.julia-depot` |
| `PATH` | Julia `bin` prepended |

First-time package install (once per depot):

```powershell
julia --project=. -e "using Pkg; Pkg.instantiate()"
```

Dependencies are pinned in [`Project.toml`](Project.toml) (HomotopyContinuation 2.22, Graphs, CSV, DataFrames). Oscar.jl is optional and not required for Phase 0 scripts.

All commands below assume you have sourced `setup_julia_env.ps1` in the same PowerShell session.

---

## Windows Smart App Control blocker (2026-09) — use WSL

Since ~September 2026, **Windows Smart App Control (state: On) blocks Julia's
unsigned core DLLs** (`libblastrampoline-5.dll`), so Julia cannot initialize on
the Windows host at all (crash at `ijl_init_`). SAC has no exclusion mechanism
and was not modified. All Julia computation now runs in **WSL2 Ubuntu**:

```bash
# one-time setup inside WSL (Ubuntu)
curl -sL -o ~/julia.tar.gz https://julialang-s3.julialang.org/bin/linux/x64/1.12/julia-1.12.6-linux-x86_64.tar.gz
tar xzf ~/julia.tar.gz -C ~/ && rm ~/julia.tar.gz
mkdir -p ~/mub6 && cp "/mnt/d/MUBs in 6-dimension/Project.toml" ~/mub6/
cp "/mnt/d/MUBs in 6-dimension/Manifest.toml" ~/mub6/
cp -r "/mnt/d/MUBs in 6-dimension/src" ~/mub6/
cp -r "/mnt/d/MUBs in 6-dimension/scripts" ~/mub6/
cd ~/mub6 && JULIA_DEPOT_PATH=~/mub-depot ~/julia-1.12.6/bin/julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Then run scripts with `MUB_GIT_DIR='/mnt/d/MUBs in 6-dimension'` so provenance
records the canonical commit. Copy outputs back to `/mnt/d/.../results/`.

---

## Regression tests and benchmarks (Part XVII / VIII)

```bash
cd ~/mub6 && export JULIA_DEPOT_PATH=~/mub-depot
MUB_LONG_TESTS=0 ~/julia-1.12.6/bin/julia --project=. test/runtests.jl   # fast (~12 s)
~/julia-1.12.6/bin/julia --project=. test/runtests.jl                    # full (~3 min)
~/julia-1.12.6/bin/julia --project=. scripts/julia/run_benchmarks.jl     # ~10 min
```

Benchmark pins (do not adjust results to pass): F6 pool = **48** vectors;
D₀-equivalent point λ=π/2 = **120** vectors / **10** third bases / N_p = 0 /
W₁ empty on all 10 cliques; Tao S₆ = no third MUB. Artifacts:
`results/benchmarks/benchmarks.{json,md}` with full PoolAudit bookkeeping and
provenance. The audit and bug ledger: [`docs/METHODOLOGY_AUDIT.md`](docs/METHODOLOGY_AUDIT.md).

---

## One-command regeneration table

Run from **project root**. Approximate wall times are for a typical laptop; homotopy solves are ~25–40 s per fresh pool solve.

| Artifact | Command | Primary outputs | Approx. time |
|----------|---------|-----------------|--------------|
| Dita \(\lambda\)-circle periodicity | `julia --project=. scripts/julia/lambda_periodicity_dita.jl` | `results/lambda_periodicity_dita.txt` | 1–3 h (126 circle samples + spot checks) |
| Dita \(\phi\)-offset sweep (9 \(\lambda\) values) | `julia --project=. scripts/julia/locus_phi_sweep_full_circle.jl` | `results/locus_phi_sweep_full_circle.txt` | 30–90 min |
| F6 \(\theta{=}0\) \(\lambda\)-arc boundaries | `julia --project=. scripts/julia/f6_boundary_reverify.jl` | `results/f6_boundary_reverify.txt` | 15–45 min |
| Third-MUB candidate audit + fourth test | `julia --project=. scripts/julia/validate_third_mub_candidates.jl` | `results/third_mub_audit.csv`, `results/fourth_mub_per_basis.csv`, `results/validation_summary.txt` | 30–120 min (scales with CSV row count) |
| CSV row-count reconciliation | `julia --project=. scripts/julia/csv_reconciliation.jl` | `results/csv_reconciliation.txt` | \< 1 min |
| Planted-clique pipeline audit | `julia --project=. scripts/julia/audit_clique_pipeline.jl` | console + logs under `results/` | ~5 min |
| Per-\(H\) validation gate (V1–V4) | `julia --project=. scripts/julia/validate_perH_pipeline.jl` | console | ~10 min |
| Dense Dita \(\lambda\) fourth test | `julia --project=. scripts/julia/dita_lambda_fourth_dense.jl --test10` | `results/dita_lambda_fourth_dense.txt` | ~5–30 min |
| Track D locus classification | `julia --project=. scripts/julia/locus_classification.jl` | `results/locus_classification.json`, `.txt` | 30–120 min |
| Per-H pool + witness M2 export | `julia --project=. scripts/julia/symbolic_elimination.jl` | `symbolic_export/pool_*.m2`, `fourth_mub_w_elimination_*.m2`, `fourth_mub_reduced_*.m2` | ~10 min |
| **Theorem T1** gauge lemmas | `julia scripts/julia/formalize_gauge_lemmas.jl` | `results/formalize_gauge_lemmas.txt` | \< 1 min (no HC) |
| **Finding T2** Dita construction (HP, not a theorem) | `julia --project=. scripts/julia/verify_dita_construction.jl` | `results/verify_dita_construction.txt` | ~30–60 min |
| **Theorem T3** \(W_1\) all-clique certify (seven \(\lambda\)) | `julia --project=. --compiled-modules=no scripts/julia/certify_nwit1_all_cliques_four_classes.jl` then `certify_nwit1_all_cliques_lambdapi.jl` | `results/certify_nwit1_all_cliques_four_classes.txt`, `certify_nwit1_all_cliques_lambdapi.txt` | hours |
| Locus geometry L5–L6 | `julia --project=. scripts/julia/locus_geometry_probes.jl` | `results/locus_geometry_probes.txt` | ~1–3 h |
| Finish audit (Python) | `python scripts/python/finish_project_audit.py` | `results/project_finish_audit.txt` | \< 1 min |
| S\* batch 1 (degen200, no-refine) | `julia --project=. scripts/julia/search_special_loci.jl --degen-cap 200 --no-refine` | `results/special_loci_degen200.csv` | 30–90 min |
| Extended HP validation (928-row CSV) | `julia --project=. scripts/julia/validate_third_mub_candidates.jl --csv results/special_loci_search.csv` | `results/third_mub_audit.csv`, `validation_summary.txt` | 30–60 min |
| Oscar install probe (optional) | `julia --project=. scripts/julia/install_oscar_probe.jl` | `results/track_c_elimination/oscar_install.log` | 30–60 min |
| Macaulay2 Groebner (Docker) | `.\scripts\docker\run_m2.ps1 pool_F6.m2` | `results/track_c_elimination/*.log` | 1–10 min |
| **Theorem T4** exact W1 Groebner at D0-equivalent point | `julia --project=. scripts/julia/verify_field_degree_w1.jl` then `.\scripts\docker\run_m2.ps1 w1_D0_groebner.m2` | `results/track_c_elimination/w1_D0_field_degree.txt`, `m2_w1_D0_groebner.log` | < 2 min total |

### Optional: regenerate the certified search CSV (long)

**Not part of the quick audit bundle.** Completing search set \(\mathcal{S}_{588}\) to 588 rows requires an overnight run:

```powershell
julia --project=. scripts/julia/search_special_loci.jl --resume
julia --project=. scripts/julia/search_special_loci.jl --verify-only
# Full degen load (Track D / S*):
julia --project=. scripts/julia/search_special_loci.jl --degen-cap 500 --no-refine --resume
```

Current certified search CSV: **928** rows (899 pool-complete) in `results/special_loci_search.csv`. Design target \(\mathcal{S}_{588}\) is a subset; see `results/csv_reconciliation.txt` for the older 471/588 gap. \(\mathcal{S}^*\) primary: `special_loci_degen1845.csv` (**1863/1865**, regenerated 2026-08-13).

Degeneracy candidates (input to search queue):

```powershell
python scripts/python/degeneracy_scan.py
```

---

## Regenerate theorem-path artifacts (Phase 1–4)

```powershell
$J = "C:\Users\adian\AppData\Local\Programs\Julia-1.12.6\bin\julia.exe"
$env:JULIA_DEPOT_PATH = "d:\MUBs in 6-dimension\.julia-depot"

# T1 — no HomotopyContinuation required
& $J "d:\MUBs in 6-dimension\scripts\julia\formalize_gauge_lemmas.jl"

# Finding T2, T3 certify, symbolic (requires HC if App Control allows)
& $J --project="d:\MUBs in 6-dimension" "d:\MUBs in 6-dimension\scripts\julia\verify_dita_construction.jl"
& $J --project="d:\MUBs in 6-dimension" "d:\MUBs in 6-dimension\scripts\julia\symbolic_elimination.jl"

# Finish audit (Python, no HC)
python scripts/python/finish_project_audit.py
```

Proofs: `../paper/proofs/gauge_and_locus.tex`, `dita_third_mub.tex`, `fourth_mub_obstruction.tex` (included in `../paper/main_theorems.tex` appendix).

---

## Regenerate all Phase 0 audit artifacts (sequential)

```powershell
. .\setup_julia_env.ps1
julia --project=. scripts/julia/csv_reconciliation.jl
julia --project=. scripts/julia/f6_boundary_reverify.jl
julia --project=. scripts/julia/locus_phi_sweep_full_circle.jl
julia --project=. scripts/julia/lambda_periodicity_dita.jl
julia --project=. scripts/julia/validate_third_mub_candidates.jl
```

For a faster smoke check (anchors only, no full circle):

```powershell
julia --project=. scripts/julia/search_special_loci.jl --quick
julia --project=. scripts/julia/validate_third_mub_candidates.jl --csv results/special_loci_search.csv
```

---

## Key result files cross-reference

| Result | Script | Used in paper / ledger |
|--------|--------|------------------------|
| `lambda_periodicity_dita.txt` | `lambda_periodicity_dita.jl` | Dita full \(\lambda\) circle (126/126) |
| `locus_phi_sweep_full_circle.txt` | `locus_phi_sweep_full_circle.jl` | Table: \(\phi\)-offset at 9 \(\lambda\) values |
| `f6_boundary_reverify.txt` | `f6_boundary_reverify.jl` | F6 \(\theta{=}0\) arc width \(\approx 0.118\) |
| `third_mub_audit.csv` | `validate_third_mub_candidates.jl` | HP audit of clique-6 CSV rows |
| `fourth_mub_per_basis.csv` | same | Per-basis fourth-MUB test |
| `validation_summary.txt` | same | Cluster rank, HP survival count |
| `csv_reconciliation.txt` | `csv_reconciliation.jl` | \(\mathcal{S}_{588}\) gap analysis |
| `locus_classification.json` | `locus_classification.jl` | Track D component catalog |
| `dita_lambda_fourth_dense.txt` | `dita_lambda_fourth_dense.jl` | Track A dense λ circle |
| `REPRODUCE.md` | manual | \(\mathcal{S}_{588}\), \(\mathcal{S}^*\), \(\mathcal{R}\) |

---

## arXiv v1 readiness checklist

Use before submitting [`../paper/main_theorems.tex`](../paper/main_theorems.tex) to arXiv (`quant-ph`).

- [ ] Replace placeholder `\author{...}` block with real name, affiliation, email
- [ ] Confirm abstract scope: **negative heuristic evidence only**; no global \(N(6)=3\) or full-family impossibility claim
- [ ] Cite reproducibility: point readers to this file and `../results/final_honest_status.md`
- [ ] Verify key numbers against regenerated artifacts (928-row CSV, 36/45 HP clique-6, 219/219 degen200, 126/126 + 628/628 Dita)
- [ ] Run `pdflatex` + `bibtex` cycle (see header comment in `.tex`)
- [ ] Optional: upload `results/` CSV subset or Zenodo archive with commit hash
- [ ] **Deferred to Phase 4:** move Tables III–IV to appendix; Track C theorem on region \(\mathcal{R}\); Track D classification JSON

Paper already uses first-person singular and documents honest scope in the abstract and Discussion. No structural appendix migration is required for arXiv v1.

---

## Troubleshooting

| Issue | Fix |
|-------|-----|
| `julia` not found | Run `. .\setup_julia_env.ps1` or edit Julia path in that script |
| `libsymengine` / HC InitError after moving project | Set `JULIA_DEPOT_PATH` to `<project-root>\.julia-depot`, delete `.julia-depot\compiled`, run `julia --project=. -e "using Pkg; Pkg.precompile()"` |
| Precompile / depot errors | Ensure `JULIA_DEPOT_PATH` points at project `.julia-depot`, not a removed LearningHub tree |
| Rounded Dita coords give `max_clique=2` | Use full-precision `acos(1/sqrt(3))`, `pi/4`, exact \(\lambda\); see `FINDINGS.md` §3b |
| `--verify-only` fails at 588 | CSV incomplete; run `--resume` or accept current 491-row baseline |
| Docker Macaulay2 | Start Docker Desktop, then run `.\scripts\docker\run_m2.ps1 pool_F6.m2`; see `results/track_c_elimination/` |

---

## Phase 3–4 quick bundle

```powershell
. .\setup_julia_env.ps1
julia --project=. scripts/julia/symbolic_elimination.jl
.\scripts\docker\run_m2.ps1 pool_F6.m2
julia --project=. scripts/julia/locus_classification.jl
julia --project=. scripts/julia/dita_lambda_fourth_dense.jl --test10
```
