# Third-MUB loci and fourth-MUB negative evidence in Karlsson 

Reproducible research on mutually unbiased bases (MUBs) arising from Karlsson's three-parameter complex Hadamard family in dimension six. Scope is **\(K_6^{(3)}\) only** — not all CHMs in \(\mathbb{C}^6\), not a proof of Zauner's conjecture \(N(6)=3\).

**Status (2026-08-05):** Phase 5 complete. **Ready for arXiv v1** with honest scope. No fourth MUB found at any audited point (928 + 219 + 628 probes).

| Resource | Path |
|----------|------|
| Technical report | [`paper/mub6_karlsson_ieee.tex`](paper/mub6_karlsson_ieee.tex) |
| Reproduction guide | [`REPRODUCE.md`](REPRODUCE.md) |
| Claim ledger | [`results/final_honest_status.md`](results/final_honest_status.md) |
| Search-set definitions | [`docs/SEARCH_SETS.md`](docs/SEARCH_SETS.md) |
| Finish audit | [`results/project_finish_audit.txt`](results/project_finish_audit.txt) |

---

## Main results

### Proved (in-repo)

| ID | Statement | Where |
|----|-----------|-------|
| **T1** | Gauge structure: at \(\theta=0\), \(\phi\) is parametrization gauge; on the Dita slice, distinct \(\lambda\) yield CHM-inequivalent matrices | [`paper/proofs/gauge_and_locus.tex`](paper/proofs/gauge_and_locus.tex) |
| **T2** | Third MUB exists on the full Dita \(\lambda\)-circle | [`paper/proofs/dita_third_mub.tex`](paper/proofs/dita_third_mub.tex) |
| **T3** | No fourth MUB on the Dita \(\lambda\)-circle | [`paper/proofs/fourth_mub_obstruction.tex`](paper/proofs/fourth_mub_obstruction.tex) |

### Contributions (paper)

| ID | Content |
|----|---------|
| **C1** | Literature correction: Karlsson's original \(A\)-matrix passes 348/349 tests; McNulty–Weigert transcription fails 349/349 |
| **C2** | **Main structural result:** exactly two HP-verified third-MUB loci — F6 \(\theta=0\) bounded \(\lambda\)-arc (width \(\approx 0.118\)) and Dita full \(\lambda\)-circle at \(\phi=\pi/4\) — CHM-inequivalent, arc vs. circle topology |
| **C3** | Confirmatory negative numerics: no fourth MUB on certified search CSV, \(\mathcal{S}^*\) batch 1, dense Dita probes, or \(\theta=0\) slices |
| **C4** | Reproducible certified per-\(H\) pipeline (HomotopyContinuation.jl + Python helpers) |

### HP-verified highlights

- **Track D classification:** 36/45 clique-6 rows survive HP → only **2 components** (27 F6 + 9 Dita HP rows); 9 \(\theta=\pi/2\) circulant_match rows fail HP
- **Track A (Dita dense):** 628/628 on full \(\lambda\) circle — clique \(\ge 6\), fourth = 0
- **Dita periodicity:** 126/126 on \([0,2\pi]\)
- **Fourth-MUB per basis:** 270/270 basis tests negative
- **Pipeline:** `found_fourth=false` everywhere audited

---

## Audited datasets

| Artifact | Rows | Fourth MUB | Notes |
|----------|------|------------|-------|
| [`results/special_loci_search.csv`](results/special_loci_search.csv) | 928 (899 pool-complete) | 0 | Design target \(\mathcal{S}_{588}\); extended refinement superset |
| [`results/special_loci_degen200.csv`](results/special_loci_degen200.csv) | 219 | 0 | \(\mathcal{S}^*\) batch 1 (top 200 degen; subset of degen500) |
| [`results/special_loci_degen500.csv`](results/special_loci_degen500.csv) | 519 | 0 | \(\mathcal{S}^*\) batch 2 (top 500 degen; subset of degen1845) |
| [`results/special_loci_degen1845.meta.txt`](results/special_loci_degen1845.meta.txt) | — | 0 | \(\mathcal{S}^*\) full primary stats (**1863/1865**; CSV quarantined) |
| [`results/dita_lambda_fourth_dense.csv`](results/dita_lambda_fourth_dense.csv) | 628 | 0 | Dense Dita \(\lambda\) circle |

---

## Open problems (explicitly out of scope or incomplete)

| Item | Status |
|------|--------|
| Fourth-MUB absence on full region \(\mathcal{R}\) | Open Problem (strong numerics only) |
| Full \(\mathcal{S}^*\) (1845 degeneracy candidates) | **1863/1865** primary complete (log/meta); regenerate `special_loci_degen1845.csv` |
| Phase C (Liang/Chen cross-family) | **Blocked** — no source PDF; parametric family not confirmed in cited papers |
| Macaulay2 Groebner witness elimination | Pool ideals exported; elimination incomplete |
| No fourth MUB in all of \(K_6^{(3)}\) | Not claimed |
| \(N(6)=3\) globally | Not addressed |

---

## Repository layout

```
paper/                  Technical report + proof appendices (T1–T3)
src/                    Core Julia: per-H MU pool, clique extension, Karlsson CHM
scripts/julia/          Sweeps, validation, locus classification, symbolic export
scripts/python/         CHM audit, degeneracy scan, finish audit
scripts/docker/         Macaulay2 via Docker (run_m2.ps1)
symbolic_export/        Exported .m2 polynomial systems
results/                CSV logs, audits, classification JSON, claim ledger
docs/                   SEARCH_SETS.md, METHODS.md, FINDINGS.md
```

**Stack:** Julia 1.12.6 (HomotopyContinuation 2.22), Python 3, optional Macaulay2 (Docker).

---

## Quick start (Windows)

```powershell
. .\setup_julia_env.ps1
julia --project=. -e "using Pkg; Pkg.instantiate()"

# T1 gauge lemmas (no HomotopyContinuation)
julia scripts/julia/formalize_gauge_lemmas.jl

# Finish-line summary (Python, no HC)
python scripts/python/finish_project_audit.py
```

Full regeneration commands: [`REPRODUCE.md`](REPRODUCE.md).

**After moving the project:** set `JULIA_DEPOT_PATH` to `<root>\.julia-depot`, delete `.julia-depot\compiled`, run `julia --project=. -e "using Pkg; Pkg.precompile()"`.

**Macaulay2:** start Docker Desktop, then `.\scripts\docker\run_m2.ps1 pool_F6.m2`.

---

## Internal tracks

| Track | Purpose | Status |
|-------|---------|--------|
| **D** | Third-MUB locus classification | Complete — 2 HP-verified components |
| **A** | Dense probes on Dita \(\lambda\)-circle | Complete — 628/628 |
| **C** | Fourth-MUB negative numerics on \(\mathcal{R}\) | Confirmatory only; no impossibility theorem |

---

## Citation

If you use this repository, cite the technical report (when on arXiv) and point readers to [`REPRODUCE.md`](REPRODUCE.md) and [`results/final_honest_status.md`](results/final_honest_status.md) for reproducibility and claim scope.
