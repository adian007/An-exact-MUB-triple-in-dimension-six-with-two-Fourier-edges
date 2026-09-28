# Third-MUB loci and fourth-MUB negative evidence in Karlsson \(K_6^{(3)}\)

Reproducible search for mutually unbiased bases arising from Karlsson's three-parameter complex Hadamard family in dimension six. Scope is **\(K_6^{(3)}\) only** — not all CHMs in \(\mathbb{C}^6\), not a proof of Zauner's conjecture \(N(6)=3\).

**Labeling:** **Theorem** is reserved for T1 (algebra), T3 (certified \(W_1\) emptiness at seven Dita-\(\lambda\) points), and T4 (exact Groebner \(W_1\) unit ideal at the \(D_0\)-equivalent pair \(\lambda\in\{\pi/2,3\pi/2\}\), single clique). T2 and the Track D locus classification are **Findings** (HP computation). See [`findings.md`](findings.md) and [`scientific_status.md`](../scientific_status.md).

**Engineering log** (bugs diagnosed and fixed): [`docs/ENGINEERING.md`](docs/ENGINEERING.md).

| Resource | Path |
|----------|------|
| Technical report | [`../paper/main_theorems.tex`](../paper/main_theorems.tex) |
| Reproduction guide | [`reproduce.md`](reproduce.md) |
| Claim ledger | [`../results/final_honest_status.md`](../results/final_honest_status.md) |
| Search-set definitions | [`reproduce.md`](reproduce.md) and [`../scientific_status.md`](../scientific_status.md) |
| Audit findings | [`findings.md`](findings.md) |
| Algebraic T3 map (D0 already a theorem; E0–E5) | [`docs/ALGEBRAIC_ATTACK.md`](docs/ALGEBRAIC_ATTACK.md) |

---

## Results (status 2026-08-14)

### Theorems

| ID | Statement | Where |
|----|-----------|-------|
| **T1** | At \(\theta=0\), \(\phi\) is a parametrization gauge; on the Dita slice, distinct \(\lambda\) yield CHM-inequivalent matrices | [`../paper/proofs/gauge_structure.tex`](../paper/proofs/gauge_structure.tex) |
| **T3** | No fourth MUB at **seven** Dita-\(\lambda\) points \(\{0,0.4,\pi/3,2\pi/3,\pi,4\pi/3,5\pi/3\}\): \(W_1\) empty at every recovered-pool 6-clique (**40/40**) | [`../paper/proofs/fourth_mub_theorems.tex`](../paper/proofs/fourth_mub_theorems.tex) |
| **T4** | Independent exact re-verification of a BW 2009 fact: \(W_1(D_{\mathrm{bc}},F_D)\) is the **unit ideal** over \(\mathbb{Q}(\zeta_{24},\sqrt5)\) (GB \(\{1\}\)). BW already had the stronger complete-pool result (\(N_v{=}120\), \(N_t{=}10\), \(N_p{=}0\)); T4 covers **one** third basis via a different method (fixed-\((H,B_3)\) exact GB). Scope: \(\lambda\in\{\pi/2,3\pi/2\}\) at \(D_{\mathrm{bc}}\) only — not an extension of BW | [`../paper/proofs/fourth_mub_theorems.tex`](../paper/proofs/fourth_mub_theorems.tex), [`findings.md`](findings.md) §17 |

T3 is **not** a certificate on the full \(\lambda\)-circle or on covering region \(\mathcal{R}\). T4 is an independent exact re-check of one BW instance (not all-clique / not an extension of BW) and is not a theorem on the Dita circle.

### Findings (HP / certified homotopy, not theorems)

| ID | Statement | Where |
|----|-----------|-------|
| **T2** | Third MUB at \(126/126\) Dita periodicity samples and \(628/628\) dense probes; all-\(\lambda\) extension is informal | [`../paper/proofs/dita_third_mub_methods.tex`](../paper/proofs/dita_third_mub_methods.tex) |
| **C2 / Track D** | Two reproducing HP-verified third-MUB loci in the audited search: F6 \(\theta=0\) bounded \(\lambda\)-arc (\(\approx 0.118\) rad) and Dita full \(\lambda\)-circle | Arc extent: [`results/f6_boundary_reverify.txt`](../../results/f6_boundary_reverify.txt) (400-bit bracket search, the authoritative number). Locus catalogue: [`results/locus_classification.json`](results/locus_classification.json) — note its `topology` interval is the sampled-point range, **not** the arc boundary (see §"Defect resolution" in `research/reports/audit_followup_2026-09-26.md`). |

> **Path B / C2 status (2026-09-26) — read before explaining the arc/circle contrast.**
> The two loci are not in question; their *explanation* is. Two candidate mechanisms have
> now been disposed of, in this order:
> 1. **Retired (false).** "F6 A-block non-Hermitian vs Dita Hermitian." Both blocks are
>    **non-Hermitian**: `results/arc_circle_asymmetry.txt` records `A Hermitian=false` for
>    both loci (lines 7, 13) while its prose asserts the contrast (lines 35, 50–51). The
>    Dita block `[[i,-1],[-1,i]]` is symmetric and satisfies `AA† = 2I`, but `i ≠ -i`, so it
>    is not Hermitian. Lemma L3 is correct and never asserted Hermiticity; the false wording
>    was added downstream of it.
> 2. **Refuted as a mechanism (coincidental two-point match).** The replacement invariant
>    `D = |α|² − |β|²` is a function of `(θ,φ)` only, and both loci hold `(θ,φ)` fixed along
>    their `λ`-loop. So `D` is **constant along each locus** and takes the *same* value at
>    surviving (in-arc) and non-surviving (out-of-arc) points of the F6 locus. It cannot
>    generate the boundary. `D` remains a correct stratum/degeneracy classifier (Karlsson
>    arXiv:1003.4177 Sec. 5 defines it as the Möbius liveness condition), but the link from
>    the stratum to `κ(H(λ)) ≥ 6` is **open**.
> 3. **The phenomenon itself is prior art.** Bounded parameter regions for third-MUB
>    existence in non-affine order-six CHM families are published: Goyeneche, *J. Phys. A*
>    **46**, 105301 (2013), **arXiv:1209.4126**, Table 1 (`K6^(3)` black/white regions;
>    `M6^(1)` bounded intervals; `D6^(1)`, `F6^(2)` full-range). **Do not present the
>    arc/circle contrast as a new structural discovery**, and do not cite this repo for the
>    phenomenon. What would be new is a *mechanism*, which is not published.
>
> Details: [`research/reports/audit_followup_2026-09-26.md`](../../research/reports/audit_followup_2026-09-26.md),
> section "D-invariant dichotomy – resolved".


### Other contributions

| ID | Content |
|----|---------|
| **C1** | Karlsson's original \(A\)-matrix passes 348/349 tests; McNulty–Weigert transcription fails 349/349 |
| **C3** | Confirmatory negatives on larger sets (CSV, \(\mathcal{S}^*\), dense Dita). Same class of method as the 512-point grid that missed every known third-MUB locus — weaker than T3 |
| **C4** | Per-\(H\) HomotopyContinuation pipeline; bug history in [`docs/ENGINEERING.md`](docs/ENGINEERING.md) |

---

## Audited datasets

| Artifact | Rows | Fourth MUB | Notes |
|----------|------|------------|-------|
| [`results/special_loci_search.csv`](results/special_loci_search.csv) | 928 (899 pool-complete) | 0 | Design target \(\mathcal{S}_{588}\); extended-refinement superset |
| [`results/special_loci_degen200.csv`](results/special_loci_degen200.csv) | 219 | 0 | \(\mathcal{S}^*\) batch 1 |
| [`results/special_loci_degen500.csv`](results/special_loci_degen500.csv) | 519 | 0 | \(\mathcal{S}^*\) batch 2 |
| [`results/special_loci_degen1845.csv`](results/special_loci_degen1845.csv) | **1863** | 0 | \(\mathcal{S}^*\) full primary (**regenerated 2026-08-13**; 1863/1865, 21 clique-6, 1852 pool-complete) |
| [`results/dita_lambda_fourth_dense.csv`](results/dita_lambda_fourth_dense.csv) | 628 | 0 | Dense Dita \(\lambda\) circle (numerical; not T3) |

---

## Open problems

| Item | Status |
|------|--------|
| Fourth-MUB absence on full region \(\mathcal{R}\) | **Open** — T3 settles seven points; exact Gröbner witness is the central obstacle |
| \(\mathcal{S}^*\) primary (1845 degeneracy candidates) | **Done** — `special_loci_degen1845.csv` regenerated (1863/1865) |
| Liang/Chen cross-family | **Closed gap** — [`docs/LIANG_CHEN_FAMILY_ASSESSMENT.md`](docs/LIANG_CHEN_FAMILY_ASSESSMENT.md) |
| Macaulay2 Groebner of the two-vector witness | **Inconclusive** (\(\dim I=-1\) over inexact `CC`); single-vector witness now **exactly eliminated at the \(D_0\)-equivalent point** (T4) — non-\(D_0\) λ still open |
| No fourth MUB in all of \(K_6^{(3)}\); \(N(6)=3\) | **Not claimed** |

---

## Repository layout

```
docs/paper/              Technical report + proofs (T1, T3) and findings (T2, L5/L6)
src/MubSearch.jl        Audited module: Karlsson variants + validator, PoolAudit,
                        projective dedup, exhaustive clique enumeration, W1 witness,
                        claim tiers, provenance (src/Karlsson|Pool|Cliques|Certification|
                        Benchmarks|Provenance.jl; legacy core in mub_zauner_6d_liang_chen.jl)
src/                    Per-H MU pool, clique extension, Karlsson CHM (legacy scripts)
scripts/julia/          Sweeps, validation, T3 certify, symbolic export, run_benchmarks.jl
scripts/python/         CHM audit, degeneracy scan, finish audit
scripts/docker/         Macaulay2 via Docker
symbolic_export/        Exported .m2 systems
test/runtests.jl        Regression suite (Karlsson/pool/graph/fourth-MUB + benchmarks)
results/                CSV logs, audits, claim ledger, benchmarks/
results/benchmarks/     Regression benchmark artifacts (A: F6=48, B: D0=120/10/Np=0, C: Tao)
docs/                   SEARCH_SETS, METHODS, FINDINGS, ENGINEERING, ALGEBRAIC_ATTACK,
                        METHODOLOGY_AUDIT, LITERATURE_2026
```

**Stack:** Julia 1.12.6 (HomotopyContinuation 2.22), Python 3, optional Macaulay2 (Docker). On Windows hosts with Smart App Control, Julia runs in WSL2 (see `reproduce.md`).

---

## Quick start (Windows)

```powershell
. .\setup_julia_env.ps1
julia --project=. -e "using Pkg; Pkg.instantiate()"

julia scripts/julia/formalize_gauge_lemmas.jl
python scripts/python/finish_project_audit.py
```

If Windows Application Control blocks precompiled DLLs under `.julia-depot`, add `--compiled-modules=no`. Full commands: [`reproduce.md`](reproduce.md).

---

## Citation

Cite the technical report (when on arXiv). For claim scope use [`../results/final_honest_status.md`](../results/final_honest_status.md); for the scientific interpretation use [`../scientific_status.md`](../scientific_status.md).
