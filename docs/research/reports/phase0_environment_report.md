# Phase 0 — Environment and Repository Audit Report

Date: 2026-09-17
Repository: `D:\MUBs in 6-dimension`
Commit: `7bf60d4aaf2881030b5d376c1757706f30107ff1` (branch `main`)
Remote: `https://github.com/adian007/MUB.git`

Scope: **Audit only.** No symbolic elimination, no mathematical computation,
no modification of any existing source file.

---

## 1. Host environment (Windows)

| Item | Value |
|------|-------|
| OS | Microsoft Windows 11 Pro, build 10.0.26200 |
| CPU | Intel(R) Core(TM) Ultra 7 258V, 8 cores / 8 logical processors |
| RAM | 32,073 MB total physical memory |
| Disk C: | 28.8 GB free |
| Disk D: | 548.5 GB free |

## 2. Software versions

### Windows side

| Software | Version | Status |
|----------|---------|--------|
| Python | 3.13.14 (`C:\Users\adian\AppData\Local\Programs\Python\Python313\python.exe`) | Available |
| SymPy | 1.14.0 | Available |
| Julia | **not installed on Windows** | Not available (see WSL) |
| Sage / Mathematica / Magma | not installed | Not available |
| Docker CLI | 29.7.2 | Installed; daemon not running |
| WSL | Ubuntu 26.04 LTS distro present | Available |

### WSL (Ubuntu 26.04 LTS)

| Software | Version | Status |
|----------|---------|--------|
| Julia | **1.12.6** at `/home/adian/julia-1.12.6/bin/julia` | Available (not on default PATH) |
| Macaulay2 | **1.26.06** (`/usr/bin/M2`) | Available |
| Singular | **4.4.1** (x86_64-Linux, 64-bit, build 44105) | Available |
| Sage / Mathematica / Magma | not installed | Not available |
| CPUs / RAM in WSL | 8 cores, ~15 GB visible | Configured |

**Conclusion:** all three critical systems — Julia (HomotopyContinuation),
Macaulay2 (exact Gröbner), Singular (backup elimination) — are available in
WSL. **No software blocker; Phase 1 may proceed.** Mathematica/Magma/Sage are
unavailable; independent reproduction will use Macaulay2 ↔ Singular ↔
Julia/Nemo cross-checks, which satisfies the "different system" requirement
for exact Gröbner results.

### Julia depot status (caveat)

- Windows depot `.julia-depot` has only `compiled/` and `logs/` — **no
  `packages/`** (artifacts deleted, matching git status).
- WSL `~/.julia`: **no packages installed**.
- Network from WSL to GitHub / pkg.julialang.org: **confirmed working**.

Setup required before any Julia run: `Pkg.instantiate()` from WSL against the
pinned Manifest (julia_version 1.12.6). Not performed in this audit. Windows
Smart App Control blocks unsigned Julia DLLs (repo doc F4); all Julia
computation therefore runs in WSL, matching repo policy. The Docker/M2 path is
unnecessary since M2 1.26.06 is native in WSL.

---

## 3. Repository audit (inspection only)

### 3.1 Working tree vs commit

`git status`: `.julia-depot/artifacts` staged-deleted (benign binary cleanup);
large untracked set = the effective current codebase: `src/MubSearch.jl`,
`src/Karlsson.jl`, `src/Pool.jl`, `src/Cliques.jl`, `src/Certification.jl`,
`src/Benchmarks.jl`, `src/Provenance.jl`, `src/FourthMUBHPC.jl`,
`src/hpc/*.jl`, ~40 scripts in `scripts/julia/`, `test/runtests.jl`, 16
`symbolic_export/*.m2` files, `paper/proofs/*.tex`, many `results/` files.

**Implication:** the effective codebase is ahead of the last commit; all
future experiments must hash-pin their inputs (already repo policy).

### 3.2 Existing claim structure

From readme.md, docs/findings.md, results/final_honest_status.md:

- **T1 (Proved):** θ=0 gauge structure; Dita-slice λ CHM-inequivalence.
- **T3 (Certified):** no fourth MUB at **seven** Dita-λ points
  {0, 0.4, π/3, 2π/3, π, 4π/3, 5π/3}; W₁ empty at all 40/40 pool cliques.
  Scope: those seven points only — not the circle, not region 𝓡.
- **T4 (Proved, exact Gröbner):** W₁(D_bc, F_D) = unit ideal over ℚ(ζ₂₄, √5)
  at λ∈{π/2, 3π/2} — re-verification of Brierley–Weigert 2009 for **one**
  third basis, not an extension of BW's complete-pool result.
- **T2 (HP finding):** third MUB at 126/126 Dita samples, 628/628 dense probes.
- **Open per repo:** full region 𝓡; non-D₀ λ exact W₁; all of K₆⁽³⁾; N(6)=3 —
  **explicitly not claimed**.

This maps onto claim categories C1–C11. **C10, C11 remain OPEN.**

### 3.3 Literature context in-repo

`docs/LITERATURE_2026.md` catalogs: Wuttig–Tindall 2026
(arXiv:2608.18053, claimed complete order-6 CHM classification — repo treats
as **not established** until formalization is public); "Triplets of MUBs"
(2025, Hadamard-cube program — conjectural); Brierley–Weigert 2009 (the D₀
anchor result that T4 re-verifies).

### 3.4 Known certificates, numerical results, and flagged risks

**Exact certificates present:**
- `symbolic_export/w1_D0_exact.m2`, `w1_D0_groebner.m2`,
  `w1_D0_dim_degree.m2` (T4 exact GB exports)
- `results/track_c_elimination/w1_D0_groebner.log` (raw log)
- `results/track_c_elimination/w1_D0_field_degree.txt`
- `paper/proofs/fourth_mub_obstruction.tex`, `gauge_and_locus.tex`,
  `dita_third_mub.tex`

**Numerical result sets present:**
- `results/special_loci_search.csv` (928), `special_loci_degen200.csv` (219),
  `special_loci_degen500.csv` (519), `special_loci_degen1845.csv` (1863),
  `dita_lambda_fourth_dense.csv` (628) — all 0 fourth-MUB hits
  (HEURISTIC_NUMERICAL class).
- `results/benchmarks/benchmarks.json` (CERTIFIED_NUMERICAL: F6 pool 48,
  D₀ pool 120 / 10 third bases / N_p=0, Tao S₆ pool 90).

**Risks / items flagged (prompt §4.6):**

1. **Unexecuted templates:** Hadamard-cube prototype "NOT started"
   (LITERATURE_2026.md §2); G₆⁽⁴⁾ (Szöllősi) family identified but not
   implemented.
2. **Corrupted/duplicate artifacts:** `special_loci_degen1845_CORRUPTED592.csv`
   marked "do not cite"; `special_loci_degen500_run2.csv` is a mislabeled
   duplicate (cap 500). degen1845 CSV regenerated 2026-08-13, log-backed.
3. **Refuted sub-claim preserved:** `circulant_match_anomaly.txt` documents a
   suspected π/2 third-MUB match as a **tolerance artifact** — correctly
   refuted, not silently dropped.
4. **Circularity check:** anchors built from source formulas
   (`exact_karlsson_dita_H.jl`, `verify_closed_form_d0.jl`); repo rejects the
   McNulty–Weigert transcription vs Karlsson's original A-matrix via
   regression test. No anchor currently defined equal to a candidate.
5. **Conjugation-convention caveat (prompt rule 1.7):** the W₁ witness systems
   in `symbolic_export/*.m2` must be re-inspected in Phase 3 to confirm
   whether `xbar` means complex conjugation, torus complexification
   (`zbar = z^{-1}`), or real coordinates. T4 docs say field ℚ(ζ₂₄, √5) with
   an involution — must be re-verified exactly, not assumed.
6. **Working tree ≠ commit:** many load-bearing files untracked; hash-pinning
   per experiment is mandatory.
7. **Depot not instantiated:** Julia packages must be installed before any
   numerical run (network confirmed).

### 3.5 Source-file hashes

All 129 relevant source files (`src/**`, `scripts/**`, `symbolic_export/**`,
`test/**`, `paper/**` — .jl, .py, .m2, .tex, .md) were SHA-256 hashed. Full
list written to `results/phase0_source_hashes.txt` (129 lines).

SHA-256 of that hash file:
`F5E097F66FDF7E84730C9BEB11B322C1619C295DB1C2BE8FD988D4901338C68A`

Copy placed at `research/hashes/phase0_source_hashes.txt`.

---

## 4. Artifact policy compliance

Required `research/` tree created: `literature/ definitions/ constructions/
anchors/ algebraic/ numerical/ intervals/ independent_reproduction/ referee/
reports/ logs/ certificates/ hashes/ claim_ledger/`.

Audit artifacts under `research/`:
- `research/reports/phase0_environment_report.md` (this file)
- `research/hashes/phase0_source_hashes.txt`
- `research/logs/phase0_wsl_probe.sh` (probe script; raw output in §2)

**Overwrite policy:** nothing pre-existing was overwritten. New files created
by this audit: `results/phase0_source_hashes.txt`,
`results/phase0_wsl_probe.sh`, and the new `research/` tree. No collisions
with existing artifacts.

---

## 5. Verdict and next-step gate

| Gate | Status |
|------|--------|
| Required software available? | **YES** — Julia 1.12.6 + M2 1.26.06 + Singular 4.4.1 in WSL Ubuntu 26.04; Python 3.13 + SymPy 1.14 on Windows |
| Any blocking failure? | **NO** |
| Repo internally consistent? | **YES** — ledger/engineering/audit docs consistent; one known corrupted CSV already flagged |
| Ready for Phase 1? | **YES**, with one setup step |

**Phase 1 gate conditions:**
1. Instantiate the Julia environment in WSL (`Pkg.instantiate()` against the
   pinned Manifest).
2. Begin Phase 1 (definitions/conventions exact test suite). Existing T1/T3/T4
   claims are treated as prior results to be **independently re-checked**, not
   trusted.

**Stopping-criteria check:** criterion C (implementation error blocking work)
**not** triggered; criterion D (budget exhausted) **not** triggered. Research
proceeds.

---

*End of Phase 0 report. No mathematical computation was performed.*

