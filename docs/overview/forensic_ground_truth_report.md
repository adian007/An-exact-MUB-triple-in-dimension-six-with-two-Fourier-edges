# Forensic ground-truth audit

**Repository:** `D:\MUBs in 6-dimension`
**Audit date:** 2026-09-23
**Scope:** reconstruction only; no research code or paper claims were fixed.

## Executive finding

The repository is version-controlled, but its history is dominated by automated
checkpoint commits: 241 commits are reachable, including many `cline
checkpoint`, `untracked files on ...`, `index on main`, and `Agent host session`
messages. The current branch is `main`, one commit ahead of `origin/main`, and
the repository was clean before this audit apart from the auditer's temporary
inventory files.

The defensible ground truth is narrower than the paper/status wording:

* The fast Julia regression suite passes **101/101** in WSL Julia 1.12.6.
* The T1 gauge script runs successfully and writes its result.
* The stored T3 artifacts support local certified-numerical W1 emptiness at
  seven Dita lambda points, covering **40/40 recovered-pool six-cliques**
  (22 cliques at `{0, 0.4, pi/3, 2pi/3}` and 18 native partner cliques at
  `{pi, 4pi/3, 5pi/3}`). This is not a full-circle theorem.
* The stored T4 Macaulay2 output reports an exact unit ideal for one
  `(D_bc, F_D)` witness over `Q(zeta_24, sqrt(5))`, degree 16. This is a
  local re-verification of the stronger published D0 result, not a new
  full-family theorem.
* The Python variant audit reproduces **348/349** for the original Karlsson
  transcription and **0/349** for the literal review transcription.
* A fresh Python table-validation run reports **848 rows / 15 candidates**,
  whereas the documented Julia artifacts report **928 rows / 45 candidates**.
  The two pipelines are not checking the same data scope and must not be
  conflated.
* Julia is not installed at the documented Windows path. WSL Julia and its
  project/depot are available and were used for the tests above. Docker is
  installed; the stored Macaulay2 log is available.

## 1. Inventory and classification

There are **542 files** under the repository tree at audit time, including
generated dependency/cache files. The machine-generated appendix at the end of
this document lists every file, its last-write timestamp, and byte length.

| Class | Contents |
|---|---|
| Source code | `src/` Julia modules; `scripts/julia/`, `scripts/python/`, `scripts/docker/`, `scripts/wolfram/`; tests; phase research prototypes |
| Generated results/data | `results/` CSV/JSON/log/TXT artifacts; `runs/`; symbolic exports; `.julia-depot/`; Python bytecode; cached PDF/source material |
| Paper/LaTeX | `paper/` `.tex` and paper workflow/status files; PDFs and extracted literature text under `pdf/` |
| Status/notes | root `readme.md`, `reproduce.md`, `findings.md`, `claim_audit.md`, `known_results_table.md`, `open_problems.md`, `literature_audit.md`, `matrix_catalogue.md`, `research_gap_report.md`, `source_quotes.md`, plus `research/` reports and ledgers |
| Unclear/provenance risk | `_spawn_cache/`, `pdf/` helper scripts and downloaded material, duplicate run outputs, `.env`, compiled depot artifacts, archived/corrupted CSVs, and checkpoint-generated research files |

The tree contains 93 `.jl`, 40 `.py`, 17 `.m2`, 9 `.tex`, 35 `.csv`, 113
`.txt`, 69 `.log`, 45 `.md`, 20 `.pdf`, and many environment/cache files.
Several generated or downloaded files are committed, so timestamps alone do not
identify the last computational producer.

## 2. Claims checked against executable evidence

### Supported, with scope

**T1 gauge lemmas.** Command:

```text
wsl.exe -d Ubuntu -- bash -lc 'export JULIA_DEPOT_PATH=$HOME/mub-depot; cd "/mnt/d/MUBs in 6-dimension"; /home/adian/julia-1.12.6/bin/julia --project=. scripts/julia/formalize_gauge_lemmas.jl'
```

Output: `Wrote /mnt/d/MUBs in 6-dimension/results/formalize_gauge_lemmas.txt`.
The stored result says `ALL PASS`; it is an algebraic/no-HomotopyContinuation
script.

**T3 local certified-numerical result.** Script:
`scripts/julia/certify_nwit1_all_cliques_four_classes.jl`, followed by
`scripts/julia/certify_nwit1_all_cliques_lambdapi.jl`. Stored outputs say:

```text
total_cliques=22 total_empty=22 all_pass=true
PASS: W1 empty for every size-6 clique at every CHM-class representative.
total_cliques=18 total_empty=18 all_pass=true
PASS: W1 empty for every native size-6 clique at pi, 4pi/3, 5pi/3.
```

The underlying implementation uses HomotopyContinuation square systems,
certification, and full residual checks, but the matrices, pools, and
cliques are numerically represented. This is certified local numerical
evidence, not an exact symbolic proof and not a statement for every lambda.

**T4 exact local algebra.** Script/export:
`scripts/julia/verify_field_degree_w1.jl`,
`symbolic_export/w1_D0_groebner.m2`, and
`scripts/docker/run_m2.ps1`. Stored output:

```text
numgens witness = 17
numvars R       = 10
field: Q(zeta_24, sqrt5), [K:Q] = 16
Stage B: dim (exact) = -1
Stage B: witness is the UNIT IDEAL over Q(zeta_24, sqrt5).
=> W1 empty at D_bc, lambda in {pi/2, 3pi/2} (exact-coefficient certificate).
```

The field-degree output independently records irreducibility and degree 16.
The separate expected `results/track_c_elimination/export_w1_exact.txt` is
absent even though `results/export_w1_exact.txt` exists; this is a reproducibility
path inconsistency. The exact GB log itself is present.

**Stored finite searches.** `python scripts/python/finish_project_audit.py`
ran successfully and reported:

```text
special_loci_search.csv: rows=928 ok=899 clique6=45 fourth=0 pool_complete=899
special_loci_degen1845.csv: rows=1863 ok=1852 clique6=21 fourth=0 pool_complete=1852
dita_lambda_fourth_dense.csv: rows=628 ok=628 clique6=628 fourth=0 pool_complete=0
VERDICT: Paper-scope claims reconciled; T3 numerical only; ... exact Groebner witness remains future work.
```

These are finite sampled outputs. They do not prove absence on the full Dita
circle, all of `K_6^(3)`, or global `N(6)=3`.

### Contradictions, missing evidence, or overclaims

* The documented T3 headline “40/40” is supported by the two stored all-clique
  logs, but the six-point multi-lambda script separately reports one clique per
  lambda and is not the 40/40 all-clique run. These are different artifacts.
* The fresh Python validator reports `n_csv_rows=848`, `n_candidates=15`,
  `found_fourth_any=False`; the stored Julia/finish audit reports 928 rows,
  45 candidates. The Python script writes a reduced/a different audit table and
  cannot be used to reproduce the Julia count.
* `results/csv_reconciliation.txt` explicitly reports 928 actual rows against a
  588 target, 925 unique coordinate keys, 899 `ok`, and only 45 refinement
  candidates. The README's “current certified CSV” is therefore a superseding
  overshoot artifact, not a 588-row result.
* `special_loci_degen500_run2.meta.txt` labels a requested 1845 run whose actual
  cap was 500. The finish audit correctly calls it a mislabeled duplicate.
* The dense Dita result is 628/628 numerical samples, not a universal theorem.
* The repository itself warns that `dim I=-1` over `CC` is inexact and
  inconclusive. Only the named-number-field T4 output supports an exact claim.
* `scripts/python/audit_karlsson_variants.py` uses floating-point numerical
  residuals. It supports a transcription audit, not exact certification.

## 3. Git history and drift

Git is present. The earliest reachable commit is
`b477404c1974703651eb398f4a9b1938f499cdf4` (`first commit`, 2026-08-06).
The first source-code commit is
`688396fe420dbe9d047bb60413239bf424782062` (`Add project source code and
initial scripts`, 2026-08-06). The current branch includes a large
checkpoint/import burst on 2026-09-17--22. The complete 241-entry log and the
complete earliest-to-HEAD name-status diff are appended below.

Material drift visible in the actual diffs:

1. The later `src/Certification.jl` introduced automatic normalized versus
   unnormalized B3 detection. Its own documentation calls this a
   **2026-09-17 normalization fix**: the prior code used RHS 6 for unit-normalized
   vectors. This is a substantive result-changing correction.
2. The later source added HP verification, tolerance diagnostics, cross-run
   checks, and explicit `karlsson2011` versus literal `mw_review` variants.
   These are methodological changes, not mere formatting.
3. The current repository adds a large research/paper/results corpus and
   committed environment/cache material after the first source commit. The
   checkpoint history does not identify a single authoritative run for every
   result file.
4. The current paper explicitly narrows T3/T4 scope compared with earlier
   broad wording; the honest-status files are newer than many result artifacts.
   The paper/status layer is therefore newer than and partly editorially
   downstream of the computation.

No exact silent constant change was found that can be safely called a theorem
change without rerunning the full Julia history. The normalization correction
is the high-risk change and must be treated as a provenance boundary.

## 4. Duplicate scripts, regenerated artifacts, and loose ends

Duplicate/near-duplicate families include:

* `locus_classification.jl` with `.json` and `.txt` outputs;
* `dita_lambda_fourth_dense.jl` with `.csv` and multiple `.log` outputs;
* `special_loci_degen200`, `degen500`, `degen500_run2`, `degen1845`, and
  archived/backup/overshoot search CSVs;
* `identify_d0_in_karlsson.py` and `.jl`;
* repeated `phase1_*`, `investigate_*`, `validate_*`, `probe*`, and
  `*_run`, `*_run2`, `*_quick` artifacts;
* exact M2 exports paired with multiple stale/partial logs.

The latest producer identifiable from timestamps is generally the 2026-08-13
T3/T4 run for theorem artifacts, 2026-09-17 for newer phase-1/source-hash
research, and 2026-09-22 for the project-finish audit/status files. The
repository does not consistently embed a producing commit hash in result
files, so “latest timestamp” is not proof of which script produced an artifact.

Open or half-finished branches not fully reflected in the top-level status:

* exact W1 elimination away from the D0-equivalent point;
* an exact algebraic third basis at a non-D0 Dita/Karlsson point;
* full-circle and full-region fourth-MUB exclusion;
* exact transition parameters for 120/72/48-vector pool changes;
* reconciliation of the Python reduced audit with the Julia 928-row audit;
* the incomplete 588-target queue and the 1863/1865 S* coverage gap;
* literature/classification verification for the newer Wuttig--Tindall claim;
* the Liang/Chen cross-family route, which the repository labels closed but
  does not implement as a verified parametric family;
* paper build/readiness and stale generated PDF/TeX synchronization.

## 5. Recommended next actions, in priority order

1. Freeze this state with a clean commit and a manifest containing the current
   Git commit, Julia/Macaulay2 versions, exact command lines, and hashes of
   every result used by the paper.
2. Reconcile the two audit pipelines: explain why the Python validator sees
   848/15 while Julia sees 928/45, and designate one authoritative dataset.
3. Re-run T3 from the current source in WSL, after recording the normalization
   fix boundary, and archive stdout/stderr plus commit hash beside each result.
4. Reproduce T4 from the documented path, including the missing export log, and
   verify the Docker/Macaulay2 command from a clean environment.
5. Only after those steps, decide whether to pursue exact non-D0 elimination or
   revise the paper to the strictly local claims already supported.

## Appendix A — commands actually run

```text
python scripts/python/audit_karlsson_variants.py
python scripts/python/validate_third_mub_table.py
python scripts/python/finish_project_audit.py
wsl Julia 1.12.6: test/runtests.jl with MUB_LONG_TESTS=0
wsl Julia 1.12.6: scripts/julia/formalize_gauge_lemmas.jl
```

## Appendix B — complete timestamped inventory


| Path | Last modified (local) | Bytes |
|---|---:|---:|
| `__pycache__/hadamard6.cpython-313.pyc` | 2026-08-02T22:05:53.0031732+02:00 | 3187 |
| `__pycache__/karlsson_k6_3.cpython-313.pyc` | 2026-08-03T10:44:44.0498274+02:00 | 7806 |
| `_all_jl.txt` | 2026-09-18T13:05:02.3483028+02:00 | 5561 |
| `_spawn_cache/fetch.ps1` | 2026-09-18T13:07:38.4022199+02:00 | 427 |
| `_spawn_cache/matrefs.txt` | 2026-09-18T13:08:29.5612705+02:00 | 411 |
| `_spawn_cache/paper_2110.12206.html` | 2026-09-18T13:07:39.8645882+02:00 | 457014 |
| `_spawn_cache/tagcontext.txt` | 2026-09-18T13:08:27.1733520+02:00 | 430 |
| `.gitignore` | 2026-09-14T21:28:36.8090394+02:00 | 408 |
| `.julia-depot/compiled/v1.12/ArgTools/aGHFV_okjxN.ji` | 2026-09-15T20:01:01.7598508+02:00 | 11803 |
| `.julia-depot/compiled/v1.12/ArgTools/aGHFV_okjxN.so` | 2026-09-15T20:00:50.1885182+02:00 | 1582064 |
| `.julia-depot/compiled/v1.12/Base64/D7K0n_okjxN.ji` | 2026-09-15T20:03:11.8654515+02:00 | 16043 |
| `.julia-depot/compiled/v1.12/Base64/D7K0n_okjxN.so` | 2026-09-15T20:01:10.8777257+02:00 | 140776 |
| `.julia-depot/compiled/v1.12/Dates/p8See_okjxN.ji` | 2026-09-15T20:00:38.6280138+02:00 | 133389 |
| `.julia-depot/compiled/v1.12/Dates/p8See_okjxN.so` | 2026-09-15T20:00:37.7965235+02:00 | 1795576 |
| `.julia-depot/compiled/v1.12/Downloads/eiA4B_okjxN.ji` | 2026-09-15T20:01:00.2008505+02:00 | 61760 |
| `.julia-depot/compiled/v1.12/Downloads/eiA4B_okjxN.so` | 2026-09-15T20:01:00.0697513+02:00 | 2157776 |
| `.julia-depot/compiled/v1.12/InteractiveUtils/0TrXF_okjxN.ji` | 2026-09-15T20:03:16.5014414+02:00 | 77380 |
| `.julia-depot/compiled/v1.12/InteractiveUtils/0TrXF_okjxN.so` | 2026-09-15T20:03:14.7610902+02:00 | 692864 |
| `.julia-depot/compiled/v1.12/JuliaSyntaxHighlighting/8OCEv_okjxN.ji` | 2026-09-15T20:03:11.7921203+02:00 | 19201 |
| `.julia-depot/compiled/v1.12/JuliaSyntaxHighlighting/8OCEv_okjxN.so` | 2026-09-15T20:01:24.8025022+02:00 | 1377312 |
| `.julia-depot/compiled/v1.12/LibCURL_jll/9JWaY_okjxN.ji` | 2026-09-15T20:01:00.7518338+02:00 | 2897 |
| `.julia-depot/compiled/v1.12/LibCURL_jll/9JWaY_okjxN.so` | 2026-09-15T20:00:53.7440351+02:00 | 34760 |
| `.julia-depot/compiled/v1.12/LibCURL/ht49g_okjxN.ji` | 2026-09-15T20:01:00.9210509+02:00 | 125412 |
| `.julia-depot/compiled/v1.12/LibCURL/ht49g_okjxN.so` | 2026-09-15T20:00:55.5696244+02:00 | 594120 |
| `.julia-depot/compiled/v1.12/LibGit2_jll/nfCpg_okjxN.ji` | 2026-09-15T20:00:47.4149315+02:00 | 2749 |
| `.julia-depot/compiled/v1.12/LibGit2_jll/nfCpg_okjxN.so` | 2026-09-15T20:00:43.1673652+02:00 | 34688 |
| `.julia-depot/compiled/v1.12/LibGit2/xrYJZ_okjxN.ji` | 2026-09-15T20:00:46.9532296+02:00 | 294239 |
| `.julia-depot/compiled/v1.12/LibGit2/xrYJZ_okjxN.so` | 2026-09-15T20:00:46.8599626+02:00 | 2210664 |
| `.julia-depot/compiled/v1.12/LibSSH2_jll/K6mup_okjxN.ji` | 2026-09-15T20:00:55.8769979+02:00 | 2597 |
| `.julia-depot/compiled/v1.12/LibSSH2_jll/K6mup_okjxN.so` | 2026-09-15T20:00:42.5449806+02:00 | 34624 |
| `.julia-depot/compiled/v1.12/Logging/PWFjL_okjxN.ji` | 2026-09-15T20:03:16.9087208+02:00 | 3446 |
| `.julia-depot/compiled/v1.12/Logging/PWFjL_okjxN.so` | 2026-09-15T20:01:08.6965921+02:00 | 28248 |
| `.julia-depot/compiled/v1.12/Markdown/AREjX_okjxN.ji` | 2026-09-15T20:03:13.8382966+02:00 | 76440 |
| `.julia-depot/compiled/v1.12/Markdown/AREjX_okjxN.so` | 2026-09-15T20:01:31.0702808+02:00 | 3280832 |
| `.julia-depot/compiled/v1.12/MozillaCACerts_jll/XKIUi_okjxN.ji` | 2026-09-15T20:01:00.8163619+02:00 | 2038 |
| `.julia-depot/compiled/v1.12/MozillaCACerts_jll/XKIUi_okjxN.so` | 2026-09-15T20:00:54.8242900+02:00 | 25176 |
| `.julia-depot/compiled/v1.12/NetworkOptions/J8H6s_okjxN.ji` | 2026-09-15T20:00:48.1999040+02:00 | 19734 |
| `.julia-depot/compiled/v1.12/NetworkOptions/J8H6s_okjxN.so` | 2026-09-15T20:00:40.2021791+02:00 | 100328 |
| `.julia-depot/compiled/v1.12/nghttp2_jll/KTGSA_okjxN.ji` | 2026-09-15T20:01:00.6063305+02:00 | 2443 |
| `.julia-depot/compiled/v1.12/nghttp2_jll/KTGSA_okjxN.so` | 2026-09-15T20:00:52.0076642+02:00 | 34496 |
| `.julia-depot/compiled/v1.12/OpenSSL_jll/M3X35_okjxN.ji` | 2026-09-15T20:00:55.8345005+02:00 | 2898 |
| `.julia-depot/compiled/v1.12/OpenSSL_jll/M3X35_okjxN.so` | 2026-09-15T20:00:42.0914446+02:00 | 39680 |
| `.julia-depot/compiled/v1.12/p7zip_jll/dfuGM_okjxN.ji` | 2026-09-15T20:01:03.4379593+02:00 | 3075 |
| `.julia-depot/compiled/v1.12/p7zip_jll/dfuGM_okjxN.so` | 2026-09-15T20:01:03.3897241+02:00 | 44352 |
| `.julia-depot/compiled/v1.12/Pkg/tUTdb_okjxN.ji` | 2026-09-15T20:03:01.8487658+02:00 | 688300 |
| `.julia-depot/compiled/v1.12/Pkg/tUTdb_okjxN.so` | 2026-09-15T20:02:49.8579681+02:00 | 47372264 |
| `.julia-depot/compiled/v1.12/Printf/3FQLY_okjxN.ji` | 2026-09-15T20:00:40.3056784+02:00 | 34964 |
| `.julia-depot/compiled/v1.12/Printf/3FQLY_okjxN.so` | 2026-09-15T20:00:34.2784689+02:00 | 163664 |
| `.julia-depot/compiled/v1.12/Serialization/zGad9_okjxN.ji` | 2026-09-15T20:03:16.6193364+02:00 | 56583 |
| `.julia-depot/compiled/v1.12/Serialization/zGad9_okjxN.so` | 2026-09-15T20:03:13.2828202+02:00 | 293328 |
| `.julia-depot/compiled/v1.12/StyledStrings/UcVoM_okjxN.ji` | 2026-09-15T20:03:11.5736355+02:00 | 98928 |
| `.julia-depot/compiled/v1.12/StyledStrings/UcVoM_okjxN.so` | 2026-09-15T20:01:20.9524428+02:00 | 6105880 |
| `.julia-depot/compiled/v1.12/Tar/G9ZYP_okjxN.ji` | 2026-09-15T20:01:02.5072031+02:00 | 64112 |
| `.julia-depot/compiled/v1.12/Tar/G9ZYP_okjxN.so` | 2026-09-15T20:01:02.4607940+02:00 | 448024 |
| `.julia-depot/compiled/v1.12/Test/JfdTE_okjxN.ji` | 2026-09-15T20:03:25.0982991+02:00 | 100438 |
| `.julia-depot/compiled/v1.12/Test/JfdTE_okjxN.so` | 2026-09-15T20:03:21.1701640+02:00 | 2765760 |
| `.julia-depot/compiled/v1.12/TOML/mjrwE_okjxN.ji` | 2026-09-15T20:00:38.4967385+02:00 | 15744 |
| `.julia-depot/compiled/v1.12/TOML/mjrwE_okjxN.so` | 2026-09-15T20:00:38.4725057+02:00 | 149512 |
| `.julia-depot/compiled/v1.12/UUIDs/SIw1t_okjxN.ji` | 2026-09-15T20:01:04.8756550+02:00 | 7591 |
| `.julia-depot/compiled/v1.12/UUIDs/SIw1t_okjxN.so` | 2026-09-15T20:01:04.8196500+02:00 | 45272 |
| `.julia-depot/compiled/v1.12/Zlib_jll/xjq3Q_okjxN.ji` | 2026-09-15T20:01:00.6761488+02:00 | 2318 |
| `.julia-depot/compiled/v1.12/Zlib_jll/xjq3Q_okjxN.so` | 2026-09-15T20:00:53.1397356+02:00 | 34416 |
| `.julia-depot/logs/manifest_usage.toml` | 2026-09-15T20:03:07.3947556+02:00 | 79 |
| `.vscode/settings.json` | 2026-08-16T18:20:13.7388180+02:00 | 1224 |
| `bibliography_verified.csv` | 2026-09-17T21:22:26.3625289+02:00 | 5246 |
| `claim_audit.md` | 2026-09-17T21:22:27.6514861+02:00 | 4573 |
| `debug-388a52.log` | 2026-08-13T01:15:03.2682890+02:00 | 248 |
| `findings.md` | 2026-08-14T17:31:05.9375204+02:00 | 319 |
| `forensic_ground_truth_report.md` | 2026-09-23T00:11:21.9687799+02:00 | 12168 |
| `known_results_table.md` | 2026-09-17T21:22:27.3142065+02:00 | 2695 |
| `literature_audit.md` | 2026-09-17T21:22:25.1014172+02:00 | 13895 |
| `Manifest.toml` | 2026-08-13T20:20:36.7445791+02:00 | 45116 |
| `matrix_catalogue.md` | 2026-09-17T20:49:31.4057325+02:00 | 3932 |
| `open_problems.md` | 2026-09-17T21:22:28.2863412+02:00 | 1307 |
| `paper/main_theorems_multifile.tex` | 2026-09-15T19:07:33.2078455+02:00 | 22563 |
| `paper/main_theorems.tex` | 2026-08-20T16:49:30.4240976+02:00 | 41137 |
| `paper/methods_audit_multifile.tex` | 2026-09-15T19:08:28.1485152+02:00 | 26522 |
| `paper/methods_audit.tex` | 2026-08-20T16:51:58.4185250+02:00 | 32680 |
| `paper/overleaf.md` | 2026-08-16T18:06:09.6603034+02:00 | 802 |
| `paper/preamble_common.tex` | 2026-08-16T15:27:56.0072389+02:00 | 956 |
| `paper/proofs/dita_third_mub_methods.tex` | 2026-08-16T18:06:17.6180645+02:00 | 2802 |
| `paper/proofs/fourth_mub_theorems.tex` | 2026-09-15T19:08:12.3413941+02:00 | 14717 |
| `paper/proofs/gauge_structure.tex` | 2026-08-16T18:20:23.4199461+02:00 | 5966 |
| `paper/proofs/locus_geometry.tex` | 2026-08-16T18:06:16.8490822+02:00 | 2304 |
| `paper/restructure_proposal.md` | 2026-08-16T15:24:49.1452520+02:00 | 8274 |
| `paper/restructure_split_status.md` | 2026-08-16T18:25:46.1909830+02:00 | 7298 |
| `parse_tindall.ps1` | 2026-09-18T12:38:58.3363913+02:00 | 4650 |
| `pdf/__pycache__/probe.cpython-313.pyc` | 2026-09-18T14:01:13.1145652+02:00 | 2706 |
| `pdf/_all_model_ids.txt` | 2026-09-18T13:56:56.6736449+02:00 | 2446 |
| `pdf/_anal.out` | 2026-09-18T14:21:41.7304774+02:00 | 1746 |
| `pdf/_anal.py` | 2026-09-18T14:21:41.5693281+02:00 | 1673 |
| `pdf/_anal2.py` | 2026-09-18T14:22:25.0838115+02:00 | 989 |
| `pdf/_analm.py` | 2026-09-18T14:22:27.9079444+02:00 | 675 |
| `pdf/_analyze_tmp.py` | 2026-09-18T13:32:19.9089963+02:00 | 760 |
| `pdf/_build_page_nvidia__ai-classification-ocr-2.html` | 2026-09-18T14:43:22.2639124+02:00 | 83735 |
| `pdf/_build_page_nvidia__htc-ocr-2.html` | 2026-09-18T14:43:19.5867094+02:00 | 83679 |
| `pdf/_build_page_nvidia__llava-3.2-8b-vit.html` | 2026-09-18T14:43:20.3470872+02:00 | 83707 |
| `pdf/_build_page_nvidia__mumu-0.5-8b.html` | 2026-09-18T14:43:21.2232959+02:00 | 83687 |
| `pdf/_build_page_nvidia__nemotron-parse-2.0.html` | 2026-09-18T14:43:18.9005056+02:00 | 225190 |
| `pdf/_build_page_nvidia__phi-3.5-vision-instruct.html` | 2026-09-18T14:43:23.2905744+02:00 | 83735 |
| `pdf/_fetchbuild_dump.py` | 2026-09-18T14:43:17.6487096+02:00 | 1172 |
| `pdf/_fetchbuild_out.txt` | 2026-09-18T14:43:23.3183587+02:00 | 8133 |
| `pdf/_probe_results.txt` | 2026-09-18T14:01:52.3360759+02:00 | 791 |
| `pdf/_showids.py` | 2026-09-18T13:58:36.8752455+02:00 | 62 |
| `pdf/_top_files.txt` | 2026-09-18T14:07:28.6865524+02:00 | 1254 |
| `pdf/.env` | 2026-09-18T13:17:04.0605997+02:00 | 83 |
| `pdf/1809.07442.pdf` | 2026-07-28T16:07:41.0214384+02:00 | 389147 |
| `pdf/2110.12206_full.txt` | 2026-09-18T11:59:35.0304500+02:00 | 37103 |
| `pdf/2110.12206.pdf` | 2026-09-18T11:50:18.3862249+02:00 | 215811 |
| `pdf/2110.12206.txt` | 2026-09-18T11:52:09.8000168+02:00 | 37103 |
| `pdf/2503.14752` | 2026-09-18T13:13:29.3173565+02:00 | 488806 |
| `pdf/3488559.pdf` | 2026-07-28T16:07:24.7797579+02:00 | 187410 |
| `pdf/78-files-data-ec5decca5ed3d6b8079e2e7e7bacc9f2-127.pdf` | 2026-07-28T16:12:12.3750510+02:00 | 320178 |
| `pdf/A concise guide to complex Hadamard matrices.pdf` | 2026-07-27T22:55:23.0208648+02:00 | 417708 |
| `pdf/beauchamp_nicoara_2006.pdf` | 2026-07-28T14:06:15.2758750+02:00 | 190237 |
| `pdf/bell_inequalities_mub.pdf` | 2026-09-18T13:14:03.4292327+02:00 | 1898935 |
| `pdf/brierley_weigert_grobner.pdf` | 2026-09-18T13:12:20.3959099+02:00 | 1498508 |
| `pdf/butterley_hall_numerical.pdf` | 2026-09-18T13:13:58.5777773+02:00 | 129800 |
| `pdf/Extension of the Set of Complex Hadamard Matrices of Size 8.pdf` | 2026-07-27T22:58:04.1773240+02:00 | 419881 |
| `pdf/fetchbuild.py` | 2026-09-18T14:03:47.8571958+02:00 | 1170 |
| `pdf/four_most_distant_bases.pdf` | 2026-09-18T13:13:58.3031740+02:00 | 348126 |
| `pdf/full_list.py` | 2026-09-18T14:48:04.9467982+02:00 | 592 |
| `pdf/g3.py` | 2026-09-18T14:05:40.7604995+02:00 | 531 |
| `pdf/grepnode.py` | 2026-09-18T14:04:25.3500144+02:00 | 380 |
| `pdf/grepparse.py` | 2026-09-18T14:05:07.3696685+02:00 | 874 |
| `pdf/grepparse2.py` | 2026-09-18T14:05:23.2560172+02:00 | 431 |
| `pdf/jaming_matolcsi_mora_mub6.pdf` | 2026-09-18T13:13:57.7706328+02:00 | 189041 |
| `pdf/karlsson_hadamards_mub6.pdf` | 2026-09-18T13:13:58.9049505+02:00 | 140864 |
| `pdf/list-models.py` | 2026-09-18T13:36:06.6488563+02:00 | 848 |
| `pdf/matolcsi_triplets_mub.pdf` | 2026-09-18T13:13:57.4159886+02:00 | 488806 |
| `pdf/mconnell_zauner_evidence.pdf` | 2026-09-18T13:13:59.4805907+02:00 | 309292 |
| `pdf/model_report.txt` | 2026-09-18T13:46:06.3999521+02:00 | 1914 |
| `pdf/mub_composite_dims_review.pdf` | 2026-09-18T13:14:01.2448027+02:00 | 2219086 |
| `pdf/MUTUALLY UNBIASED BASES AND.pdf` | 2026-07-27T22:55:53.3890184+02:00 | 251304 |
| `pdf/Mutually Unbiased Bases in Composite Dimensions -- A Review.pdf` | 2026-07-27T22:56:26.6759779+02:00 | 2213987 |
| `pdf/parse_dump.py` | 2026-09-18T14:04:45.0350746+02:00 | 622 |
| `pdf/probe_models.py` | 2026-09-18T13:45:25.9226125+02:00 | 703 |
| `pdf/probe_writer.py` | 2026-09-18T13:46:06.0992126+02:00 | 1451 |
| `pdf/probe.py` | 2026-09-18T13:59:50.4152600+02:00 | 1255 |
| `pdf/probe2.py` | 2026-09-18T14:01:49.6571718+02:00 | 1034 |
| `pdf/probe3.py` | 2026-09-18T14:03:55.6147375+02:00 | 1136 |
| `pdf/PRXQuantum.3.010101.pdf` | 2026-07-28T16:06:03.4845409+02:00 | 1623998 |
| `pdf/search_models.py` | 2026-09-18T14:47:32.7401983+02:00 | 798 |
| `pdf/THREE-PARAMETER COMPLEX HADAMARD.pdf` | 2026-07-28T13:44:38.3048031+02:00 | 187285 |
| `Project.toml` | 2026-09-15T22:01:44.6283613+02:00 | 936 |
| `readme.md` | 2026-09-15T18:33:55.2713049+02:00 | 7485 |
| `reproduce.md` | 2026-09-15T18:33:08.5453035+02:00 | 12382 |
| `research_gap_report.md` | 2026-09-17T20:49:31.6795201+02:00 | 1863 |
| `research/algebraic/algebraic_obstruction_program.md` | 2026-09-17T20:08:37.3343923+02:00 | 2202 |
| `research/anchors/anchor_selection.md` | 2026-09-17T20:08:37.3405443+02:00 | 1510 |
| `research/campaigns/campaign1/campaign1_report.md` | 2026-09-18T01:09:28.3196964+02:00 | 3900 |
| `research/certificates/certificate_design.md` | 2026-09-17T20:08:37.3398229+02:00 | 1957 |
| `research/claim_ledger/research_claim_ledger.md` | 2026-09-17T20:08:37.3405443+02:00 | 1408 |
| `research/claims/claim_ledger_campaigns.csv` | 2026-09-17T22:14:39.9152090+02:00 | 1437 |
| `research/constructions/dimension6_constructions.md` | 2026-09-17T20:08:37.3343923+02:00 | 1920 |
| `research/definitions/__pycache__/phase1_definitions_suite.cpython-313.pyc` | 2026-09-17T21:18:38.9167637+02:00 | 48429 |
| `research/definitions/mub_dimension6_research.md` | 2026-09-17T20:08:37.3343923+02:00 | 1949 |
| `research/definitions/phase1_definitions_suite.out.txt` | 2026-09-18T00:03:49.5128233+02:00 | 8070 |
| `research/definitions/phase1_definitions_suite.py` | 2026-09-17T23:52:40.0430457+02:00 | 36919 |
| `research/exact_systems/agent4/exact_systems_report.md` | 2026-09-17T22:19:35.4063576+02:00 | 4730 |
| `research/exact_systems/agent4/systems.json` | 2026-09-17T22:19:11.3079035+02:00 | 448 |
| `research/exact_systems/agent4/variable_quantifier_table.csv` | 2026-09-17T22:19:09.1900667+02:00 | 1009 |
| `research/exact_systems/agent4/w1_unnormalized_template.m2` | 2026-09-17T22:19:13.9986546+02:00 | 430 |
| `research/four_parameter_model/agent3/report.md` | 2026-09-18T00:47:56.1719312+02:00 | 3822 |
| `research/hashes/phase0_source_hashes.txt` | 2026-09-17T18:42:30.1208132+02:00 | 16259 |
| `research/hashes/research_provenance.md` | 2026-09-17T20:08:37.3416265+02:00 | 1238 |
| `research/independent_reproduction/reproduction_protocol.md` | 2026-09-17T20:08:37.3416265+02:00 | 1377 |
| `research/intervals/interval_arithmetic_research.md` | 2026-09-17T20:08:37.3398229+02:00 | 1908 |
| `research/literature/agent9/applicability_assessment.md` | 2026-09-17T22:20:09.5098847+02:00 | 9957 |
| `research/literature/agent9/method_bibliography.csv` | 2026-09-17T22:19:23.8436890+02:00 | 5843 |
| `research/literature/agent9/method_bibliography.md` | 2026-09-17T22:19:25.7797648+02:00 | 4820 |
| `research/literature/agent9/raw_search_notes.md` | 2026-09-17T22:19:20.6150877+02:00 | 3018 |
| `research/literature/matrix_catalogue_source_audit.md` | 2026-09-17T20:47:12.4097659+02:00 | 8204 |
| `research/literature/order_six_mub_literature.md` | 2026-09-17T20:08:37.3343923+02:00 | 2322 |
| `research/logs/campaign0/campaign0_environment_audit.md` | 2026-09-17T22:14:37.6171083+02:00 | 3801 |
| `research/logs/phase0_wsl_probe.sh` | 2026-09-17T18:44:03.2283586+02:00 | 492 |
| `research/logs/research_log.md` | 2026-09-17T20:08:37.3416265+02:00 | 807 |
| `research/methods/agent2/method_comparison.csv` | 2026-09-17T22:16:48.0271550+02:00 | 4039 |
| `research/methods/agent2/method_comparison.md` | 2026-09-17T22:16:50.4136661+02:00 | 8217 |
| `research/methods/agent2/prototype/command.txt` | 2026-09-17T22:16:55.3274038+02:00 | 65 |
| `research/methods/agent2/prototype/mub_h2_w1.jl` | 2026-09-17T22:16:52.5990589+02:00 | 573 |
| `research/methods/agent2/prototype/NOT_RUN.log` | 2026-09-17T22:16:58.0575107+02:00 | 204 |
| `research/numerical/numerical_strategy_research.md` | 2026-09-17T20:08:37.3398229+02:00 | 1925 |
| `research/optimization/agent6/check_env.py` | 2026-09-17T23:38:01.1826309+02:00 | 1229 |
| `research/parameter_spaces/agent1/logs/literature_and_repository_audit.md` | 2026-09-17T22:16:59.4097793+02:00 | 905 |
| `research/parameter_spaces/agent1/parameter_space_report.md` | 2026-09-17T22:16:52.6126013+02:00 | 11140 |
| `research/parameter_spaces/agent1/parameterization_catalogue.csv` | 2026-09-17T22:16:55.3115874+02:00 | 2518 |
| `research/readme.md` | 2026-09-17T20:08:37.3343923+02:00 | 2637 |
| `research/referee/referee_risk_register.md` | 2026-09-17T20:08:37.3443090+02:00 | 1267 |
| `research/reports/next_research_report.md` | 2026-09-17T20:08:37.3491969+02:00 | 1583 |
| `research/reports/phase0_environment_report.md` | 2026-09-17T18:48:37.8687436+02:00 | 9147 |
| `research/reports/physics_interpretation.md` | 2026-09-17T20:08:37.3495387+02:00 | 1462 |
| `research/symmetry/agent7/symmetry_gauge_report.md` | 2026-09-17T22:19:52.5918606+02:00 | 6892 |
| `research/symmetry/agent7/transformation_catalogue.csv` | 2026-09-17T22:19:50.5232022+02:00 | 5020 |
| `results/arc_circle_asymmetry.txt` | 2026-08-14T15:45:08.9401934+02:00 | 2700 |
| `results/audit_clique_log.txt` | 2026-08-03T10:48:16.0873465+02:00 | 35398 |
| `results/audit_clique_log2.txt` | 2026-08-03T10:51:40.8034668+02:00 | 59574 |
| `results/benchmarks/benchmarks_run.log` | 2026-09-15T18:31:47.8415693+02:00 | 2769 |
| `results/benchmarks/benchmarks.json` | 2026-09-17T20:44:00.4357124+02:00 | 10094 |
| `results/benchmarks/benchmarks.md` | 2026-09-17T20:44:00.5595933+02:00 | 487 |
| `results/certify_fourth_mub_witness_dita.txt` | 2026-08-13T20:17:05.1956484+02:00 | 918 |
| `results/certify_fourth_mub_witness_nwit1_dita.txt` | 2026-08-13T20:48:29.2746587+02:00 | 1423 |
| `results/certify_fourth_mub_witness_nwit1_multi_lambda.txt` | 2026-08-13T21:09:52.9413310+02:00 | 2362 |
| `results/certify_fourth_mub_witness_run.log` | 2026-08-13T15:32:46.3538666+02:00 | 101602 |
| `results/certify_nwit1_all_cliques_four_classes.txt` | 2026-08-13T21:40:29.1637929+02:00 | 3059 |
| `results/certify_nwit1_all_cliques_lambdapi.txt` | 2026-08-14T16:16:51.1539837+02:00 | 2517 |
| `results/chm_equivalence.txt` | 2026-08-03T13:27:00.8538327+02:00 | 1095 |
| `results/circulant_match_anomaly.txt` | 2026-08-14T15:45:08.9401934+02:00 | 1915 |
| `results/csv_reconciliation.txt` | 2026-08-05T15:52:33.4464871+02:00 | 1459 |
| `results/degeneracy_candidates.csv` | 2026-08-03T11:50:26.9764493+02:00 | 266096 |
| `results/degeneracy_candidates.json` | 2026-08-03T10:45:01.8546198+02:00 | 348 |
| `results/dita_lambda_fourth_dense_test.log` | 2026-08-03T15:35:56.8066581+02:00 | 50478 |
| `results/dita_lambda_fourth_dense.csv` | 2026-08-03T16:02:38.9522899+02:00 | 33574 |
| `results/dita_lambda_fourth_dense.log` | 2026-08-03T16:02:38.9551260+02:00 | 30269 |
| `results/drop_detection_stress_test.txt` | 2026-08-13T16:38:42.7918388+02:00 | 289 |
| `results/enumerate_d0_cliques.txt` | 2026-08-15T14:33:38.0812330+02:00 | 1320 |
| `results/exact_karlsson_dita_H.txt` | 2026-09-17T22:02:29.9305720+02:00 | 2060 |
| `results/export_w1_exact.txt` | 2026-08-15T14:48:49.2425430+02:00 | 678 |
| `results/extended_axis_sweep_full.log` | 2026-08-03T14:34:46.6520726+02:00 | 9398324 |
| `results/extended_axis_sweep_item2.log` | 2026-08-03T14:53:42.0934840+02:00 | 8977904 |
| `results/extended_axis_sweep_quick.log` | 2026-08-03T13:42:25.5145387+02:00 | 1797512 |
| `results/extended_axis_sweep_run.log` | 2026-08-03T14:05:27.2464343+02:00 | 6135142 |
| `results/extended_axis_sweep.txt` | 2026-08-03T13:20:39.7698975+02:00 | 125853 |
| `results/f6_arc_boundary_tight.stderr.log` | 2026-09-17T18:59:48.1323401+02:00 | 794 |
| `results/f6_arc_boundary_tight.stdout.log` | 2026-09-17T18:59:42.7752431+02:00 | 2034 |
| `results/f6_boundary_reverify_run.log` | 2026-08-03T14:57:23.2654999+02:00 | 97520 |
| `results/f6_boundary_reverify_run2.log` | 2026-08-03T14:59:48.7505772+02:00 | 59424 |
| `results/f6_boundary_reverify.txt` | 2026-08-03T14:59:48.7245652+02:00 | 419 |
| `results/final_honest_status.md` | 2026-09-15T18:32:41.7181673+02:00 | 10363 |
| `results/formalize_gauge_lemmas.txt` | 2026-09-23T00:11:30.6861205+02:00 | 377 |
| `results/fourth_mub_locus_sweep_quick.log` | 2026-08-03T13:34:38.8590042+02:00 | 9302 |
| `results/fourth_mub_locus_sweep_run.log` | 2026-08-03T13:33:30.0101816+02:00 | 11802 |
| `results/fourth_mub_locus_sweep_run2.log` | 2026-08-03T13:50:18.2399412+02:00 | 81390 |
| `results/fourth_mub_locus_sweep.txt` | 2026-08-03T13:50:18.1979695+02:00 | 2546 |
| `results/fourth_mub_per_basis.csv` | 2026-08-05T15:57:36.6060145+02:00 | 16788 |
| `results/fourth_test_anchors.log` | 2026-08-03T12:13:13.1007390+02:00 | 15904 |
| `results/gauge_analysis.txt` | 2026-08-03T13:27:59.2357623+02:00 | 2496 |
| `results/grid20_run.log` | 2026-08-03T14:18:47.7233895+02:00 | 1099589 |
| `results/identify_d0_in_karlsson.txt` | 2026-08-15T14:33:59.2536273+02:00 | 2430 |
| `results/investigate_advanced_summary.txt` | 2026-08-03T12:37:21.3779677+02:00 | 1194 |
| `results/investigate_degen_circulant.txt` | 2026-08-03T12:44:33.4380117+02:00 | 14132 |
| `results/investigate_Dita.log` | 2026-08-03T14:58:06.6580318+02:00 | 19426 |
| `results/investigate_Dita.txt` | 2026-08-03T12:41:33.2077845+02:00 | 17812 |
| `results/investigate_F6_theta0.log` | 2026-08-03T14:54:11.3727426+02:00 | 15340 |
| `results/investigate_F6_theta0.txt` | 2026-08-03T12:27:36.2341255+02:00 | 22684 |
| `results/karlsson_perH_sweep_retagged.csv` | 2026-08-03T11:40:49.1291169+02:00 | 66690 |
| `results/karlsson_perH_sweep.csv` | 2026-08-03T10:01:07.0460760+02:00 | 58451 |
| `results/karlsson_sweep.csv` | 2026-08-02T22:24:10.6688839+02:00 | 184661 |
| `results/lambda_periodicity_dita_run.log` | 2026-08-03T14:02:11.3385534+02:00 | 904214 |
| `results/lambda_periodicity_dita.txt` | 2026-08-03T14:02:11.3003555+02:00 | 1175 |
| `results/locus_2d_grid.csv` | 2026-08-03T14:56:04.1420846+02:00 | 26971 |
| `results/locus_classification_run.log` | 2026-08-03T15:44:31.5889720+02:00 | 4336 |
| `results/locus_classification.json` | 2026-08-05T15:58:31.5662686+02:00 | 503 |
| `results/locus_classification.txt` | 2026-08-05T15:58:31.7378639+02:00 | 875 |
| `results/locus_dimension_confirm_run.log` | 2026-08-03T14:45:09.8691162+02:00 | 2532 |
| `results/locus_dimension_confirm.txt` | 2026-08-03T14:49:29.8768874+02:00 | 3282 |
| `results/locus_geometry_probes.txt` | 2026-08-06T14:23:47.3394898+02:00 | 1261 |
| `results/locus_phi_sweep_full_circle_run.log` | 2026-08-03T14:58:29.2110095+02:00 | 301718 |
| `results/locus_phi_sweep_full_circle.txt` | 2026-08-03T14:58:29.1638075+02:00 | 4348 |
| `results/nid_probe_run.log` | 2026-08-03T13:43:19.2008085+02:00 | 5744 |
| `results/nid_probe.txt` | 2026-08-03T14:08:46.3545990+02:00 | 77625 |
| `results/phase0_source_hashes.txt` | 2026-09-17T18:42:30.1208132+02:00 | 16259 |
| `results/phase0_wsl_probe.sh` | 2026-09-17T18:44:03.2283586+02:00 | 492 |
| `results/phase1_boundaries.txt` | 2026-08-03T15:00:10.3856009+02:00 | 952 |
| `results/phase1_definitions_suite.stdout.log` | 2026-09-17T19:05:56.6179596+02:00 | 3892 |
| `results/phase1_final_status.md` | 2026-09-18T00:03:37.1683817+02:00 | 2897 |
| `results/phase1_final.txt` | 2026-09-17T20:05:54.9624721+02:00 | 0 |
| `results/phase1_full_run.py` | 2026-09-17T19:15:47.9030658+02:00 | 657 |
| `results/phase1_full_run.txt` | 2026-09-17T19:16:04.5877666+02:00 | 0 |
| `results/phase1_full_run2.txt` | 2026-09-17T19:38:05.5774659+02:00 | 19306 |
| `results/phase1_isolated.py` | 2026-09-17T19:34:58.9158742+02:00 | 467 |
| `results/phase1_pkg_instantiate.sh` | 2026-09-17T18:56:15.6034372+02:00 | 271 |
| `results/phase1_resolution.txt` | 2026-08-03T14:57:00.0421979+02:00 | 4148 |
| `results/phase1_run_dir.txt` | 2026-09-17T19:02:33.0370528+02:00 | 67 |
| `results/phase1_run3.txt` | 2026-09-17T20:45:10.4288481+02:00 | 21650 |
| `results/phase1_simplify_probe.py` | 2026-09-17T19:06:58.6265013+02:00 | 1358 |
| `results/phase1_simplify_probe.txt` | 2026-09-17T19:07:05.4540375+02:00 | 640 |
| `results/phase1_timed.txt` | 2026-09-17T20:07:50.5013157+02:00 | 4732 |
| `results/phase1_timed2.txt` | 2026-09-17T20:09:08.3176867+02:00 | 15720 |
| `results/phase1_timing_probe.py` | 2026-09-17T19:11:36.7798657+02:00 | 873 |
| `results/phase1_timing_probe.txt` | 2026-09-17T19:11:44.5802934+02:00 | 67 |
| `results/phase2_technical.txt` | 2026-08-03T14:11:04.3903274+02:00 | 3440 |
| `results/pi3_float64_smoke.stderr.log` | 2026-09-17T18:30:45.0327561+02:00 | 0 |
| `results/pi3_float64_smoke.stdout.log` | 2026-09-17T18:32:51.8557556+02:00 | 2533 |
| `results/pi3_precision_256.txt` | 2026-09-17T18:54:31.1668149+02:00 | 2470 |
| `results/precision_f6_bounded_20260917/anchor.stderr.log` | 2026-09-17T17:32:51.0375587+02:00 | 1412 |
| `results/precision_f6_bounded_20260917/anchor.stdout.log` | 2026-09-17T17:32:38.8371213+02:00 | 31 |
| `results/precision_f6_bounded_20260917/anchor/metadata.txt` | 2026-09-17T17:32:38.7499173+02:00 | 503 |
| `results/precision_pi3_run.stderr.log` | 2026-09-17T17:38:17.6301218+02:00 | 0 |
| `results/precision_pi3_run.stdout.log` | 2026-09-17T17:38:17.6301218+02:00 | 0 |
| `results/priority_hits_verify.txt` | 2026-08-03T13:17:23.2287051+02:00 | 49814 |
| `results/probe.err.txt` | 2026-09-17T19:06:43.1005826+02:00 | 0 |
| `results/probe.out.txt` | 2026-09-17T19:06:43.4064940+02:00 | 40 |
| `results/project_finish_audit.txt` | 2026-09-23T00:11:30.6861205+02:00 | 1755 |
| `results/reconstruct_b3_algebraic_bound32.txt` | 2026-08-13T20:28:40.3035989+02:00 | 3706 |
| `results/reconstruct_b3_algebraic_lsround_archive.meta.txt` | 2026-08-13T12:17:31.1462369+02:00 | 131 |
| `results/reconstruct_b3_algebraic_lsround_archive.txt` | 2026-08-13T12:17:31.1032603+02:00 | 4787 |
| `results/reconstruct_b3_algebraic_sqrt6.txt` | 2026-08-13T20:28:09.1853826+02:00 | 3539 |
| `results/reconstruct_b3_algebraic.meta.txt` | 2026-08-13T20:28:40.4493456+02:00 | 175 |
| `results/reconstruct_b3_algebraic.txt` | 2026-08-13T20:27:52.9301422+02:00 | 3704 |
| `results/reconstruct_b3_B3_lambda0.4.csv` | 2026-08-13T12:15:25.7235189+02:00 | 1403 |
| `results/reconstruct_b3_run_deg3.log` | 2026-08-13T12:15:30.7557103+02:00 | 14040 |
| `results/reconstruct_b3_run_deg4_sqrt6.log` | 2026-08-13T12:16:20.3330258+02:00 | 7334 |
| `results/reconstruct_b3_run_deg4.log` | 2026-08-13T12:16:20.3800851+02:00 | 7424 |
| `results/refine_pi3_256.stderr.log` | 2026-09-17T18:52:26.8978099+02:00 | 0 |
| `results/refine_pi3_256.stdout.log` | 2026-09-17T18:54:31.1802832+02:00 | 5241 |
| `results/retag_512_run.log` | 2026-08-03T11:40:49.7465610+02:00 | 3378 |
| `results/retag_512_summary.txt` | 2026-08-03T11:40:49.2293843+02:00 | 243 |
| `results/search_full_run_err.log` | 2026-08-03T13:03:34.7807194+02:00 | 1998 |
| `results/search_full_run.log` | 2026-08-03T13:03:34.7807194+02:00 | 966425 |
| `results/special_loci_clean_run.txt` | 2026-08-03T11:50:09.5014975+02:00 | 2150 |
| `results/special_loci_degen1845_CORRUPTED592.csv` | 2026-08-13T01:15:03.3098187+02:00 | 80423 |
| `results/special_loci_degen1845_regen.csv` | 2026-08-13T20:06:07.6239201+02:00 | 277795 |
| `results/special_loci_degen1845.csv` | 2026-08-13T20:09:01.7326953+02:00 | 277795 |
| `results/special_loci_degen1845.log` | 2026-08-13T01:15:03.3053008+02:00 | 4429454 |
| `results/special_loci_degen1845.meta.txt` | 2026-08-14T15:45:08.9381781+02:00 | 942 |
| `results/special_loci_degen1845.timing.txt` | 2026-08-13T01:15:03.3083111+02:00 | 36 |
| `results/special_loci_degen200.csv` | 2026-08-03T16:16:40.2592244+02:00 | 32607 |
| `results/special_loci_degen200.log` | 2026-08-03T16:16:40.2792015+02:00 | 731920 |
| `results/special_loci_degen500_run2.csv` | 2026-08-13T01:15:03.3198603+02:00 | 78190 |
| `results/special_loci_degen500_run2.log` | 2026-08-13T01:15:03.3258573+02:00 | 1700482 |
| `results/special_loci_degen500_run2.meta.txt` | 2026-08-13T01:15:03.3278603+02:00 | 380 |
| `results/special_loci_degen500.csv` | 2026-08-13T01:15:03.3103238+02:00 | 78190 |
| `results/special_loci_degen500.log` | 2026-08-13T01:15:03.3178545+02:00 | 1705486 |
| `results/special_loci_log.txt` | 2026-08-03T10:48:36.8379758+02:00 | 1625638 |
| `results/special_loci_resume_run.log` | 2026-08-03T16:04:19.4213298+02:00 | 1712904 |
| `results/special_loci_run_log.txt` | 2026-08-03T10:52:09.2895725+02:00 | 63604 |
| `results/special_loci_run_phase1a.log` | 2026-08-03T16:02:07.8399206+02:00 | 1909156 |
| `results/special_loci_run.log` | 2026-08-03T14:46:26.3648312+02:00 | 2018106 |
| `results/special_loci_run2_err.txt` | 2026-08-03T11:31:40.0675793+02:00 | 4189 |
| `results/special_loci_run2.txt` | 2026-08-03T11:49:54.9545586+02:00 | 1767643 |
| `results/special_loci_search_891_overshoot.csv` | 2026-08-03T16:02:46.1764167+02:00 | 186690 |
| `results/special_loci_search_896.csv` | 2026-08-03T14:17:18.6972885+02:00 | 128430 |
| `results/special_loci_search_backup848.csv` | 2026-08-03T11:49:54.7226640+02:00 | 112654 |
| `results/special_loci_search_interim.csv` | 2026-08-03T12:23:42.2301353+02:00 | 47034 |
| `results/special_loci_search_partial210.csv` | 2026-08-03T12:38:46.5447816+02:00 | 27448 |
| `results/special_loci_search_v2.csv` | 2026-08-03T12:23:17.9815926+02:00 | 36886 |
| `results/special_loci_search.csv` | 2026-08-03T16:32:11.1153217+02:00 | 192242 |
| `results/sweep_log.txt` | 2026-08-03T10:01:07.0651067+02:00 | 2018258 |
| `results/symbolic_elimination_attempt.txt` | 2026-08-13T01:15:03.3288553+02:00 | 1686 |
| `results/symbolic_elimination_log.txt` | 2026-08-03T10:50:01.2742738+02:00 | 15228 |
| `results/symbolic_elimination_log2.txt` | 2026-08-03T10:53:56.0641116+02:00 | 27172 |
| `results/symbolic_elimination_run.txt` | 2026-08-03T10:58:27.7143896+02:00 | 23105 |
| `results/symbolic_elimination_summary.txt` | 2026-08-13T01:15:03.3298578+02:00 | 968 |
| `results/task1_investigate_candidate.log` | 2026-08-03T11:08:18.8928010+02:00 | 4512 |
| `results/task1_investigate.log` | 2026-08-03T11:16:35.6797021+02:00 | 16328 |
| `results/task1_mobius_investigation.txt` | 2026-08-03T11:07:59.4751219+02:00 | 969 |
| `results/task1_run.log` | 2026-08-03T11:07:59.5031069+02:00 | 10666 |
| `results/task2_argmax_audit.log` | 2026-08-03T12:03:53.8083839+02:00 | 17364 |
| `results/task2_argmax_audit.txt` | 2026-08-03T12:11:19.5072895+02:00 | 1027 |
| `results/task3_locus_run.log` | 2026-08-03T12:14:17.0715685+02:00 | 8992 |
| `results/third_mub_audit.csv` | 2026-08-05T15:56:09.4604129+02:00 | 5518 |
| `results/third_mub_candidates_audit_py.csv` | 2026-09-23T00:09:58.3271883+02:00 | 1407 |
| `results/third_mub_candidates_audit.csv` | 2026-08-03T14:49:18.7297257+02:00 | 2544 |
| `results/third_mub_locus_dimension.txt` | 2026-08-03T12:31:18.7803610+02:00 | 847 |
| `results/third_mub_locus_sampling_run.log` | 2026-08-03T13:14:21.1504921+02:00 | 4562478 |
| `results/third_mub_locus_sampling.txt` | 2026-08-03T13:14:21.1297003+02:00 | 1820 |
| `results/top50_degen_stress.txt` | 2026-08-06T14:23:47.3258084+02:00 | 410 |
| `results/track_c_elimination/arxiv_format_pass/after_dita_claims.txt` | 2026-08-15T15:46:13.7058561+02:00 | 680 |
| `results/track_c_elimination/arxiv_format_pass/after_fourth_claims.txt` | 2026-08-15T15:46:13.7134779+02:00 | 6336 |
| `results/track_c_elimination/arxiv_format_pass/after_gauge_claims.txt` | 2026-08-15T15:46:13.6982906+02:00 | 3243 |
| `results/track_c_elimination/arxiv_format_pass/after_main_source_claims.txt` | 2026-08-15T15:46:13.6912174+02:00 | 4339 |
| `results/track_c_elimination/arxiv_format_pass/before_dita_claims.txt` | 2026-08-15T15:45:14.0967884+02:00 | 680 |
| `results/track_c_elimination/arxiv_format_pass/before_fourth_claims.txt` | 2026-08-15T15:45:14.1008038+02:00 | 6420 |
| `results/track_c_elimination/arxiv_format_pass/before_gauge_claims.txt` | 2026-08-15T15:45:14.0967884+02:00 | 3243 |
| `results/track_c_elimination/arxiv_format_pass/before_main_source_claims.txt` | 2026-08-15T15:45:14.0948970+02:00 | 4339 |
| `results/track_c_elimination/m2_fourth_mub_reduced_Dita_elim.log` | 2026-08-13T01:15:03.3308572+02:00 | 426 |
| `results/track_c_elimination/m2_fourth_mub_reduced_Dita_exact_smoke.log` | 2026-08-13T01:15:03.3308572+02:00 | 66 |
| `results/track_c_elimination/m2_fourth_mub_reduced_Dita.log` | 2026-08-13T01:15:03.3298578+02:00 | 340 |
| `results/track_c_elimination/m2_pool_Dita.log` | 2026-08-03T15:43:09.5701071+02:00 | 364 |
| `results/track_c_elimination/m2_pool_F6.log` | 2026-08-03T15:43:56.2744911+02:00 | 348 |
| `results/track_c_elimination/pool_Dita_exact.log` | 2026-08-06T14:23:47.3394898+02:00 | 284 |
| `results/track_c_elimination/w1_D0_cost_estimate.log` | 2026-08-15T14:49:36.4222815+02:00 | 11613936 |
| `results/track_c_elimination/w1_D0_field_degree.txt` | 2026-08-15T14:56:01.2991586+02:00 | 284 |
| `results/track_c_elimination/w1_D0_groebner.log` | 2026-08-15T14:50:37.7241967+02:00 | 1714 |
| `results/validate_finish_run.log` | 2026-08-05T15:57:50.5642481+02:00 | 599886 |
| `results/validate_interim_err.log` | 2026-08-03T12:35:05.3423164+02:00 | 4193 |
| `results/validate_interim.log` | 2026-08-03T12:35:05.3423164+02:00 | 113780 |
| `results/validate_log.txt` | 2026-08-02T23:11:18.3183823+02:00 | 46466 |
| `results/validate_run.log` | 2026-08-03T14:50:04.5617756+02:00 | 262292 |
| `results/validation_summary_py.txt` | 2026-09-23T00:11:30.6861205+02:00 | 146 |
| `results/validation_summary.txt` | 2026-08-05T15:57:50.5642481+02:00 | 211 |
| `results/verify_chm_b3_transfer.txt` | 2026-08-13T22:51:00.3945136+02:00 | 2844 |
| `results/verify_closed_form_d0.txt` | 2026-08-15T14:34:56.0142753+02:00 | 2476 |
| `results/verify_dita_construction.txt` | 2026-08-06T14:23:47.3394898+02:00 | 961 |
| `results/verify_extra_third_hit_run.log` | 2026-08-03T13:34:51.7009239+02:00 | 2700 |
| `results/verify_nwit1_followup.txt` | 2026-08-13T21:20:03.7368448+02:00 | 1684 |
| `results/verify_nwit1_referee_checks.txt` | 2026-08-13T21:16:08.3906572+02:00 | 3606 |
| `results/w1_D0_mathematica.txt` | 2026-09-17T20:23:32.0413837+02:00 | 237 |
| `results/w1_D0_wolfram_rerun_raw.txt` | 2026-09-17T21:06:08.1270740+02:00 | 148 |
| `results/witness_decomposition_run.log` | 2026-08-03T11:42:52.0066751+02:00 | 12132 |
| `results/witness_decomposition_run2.log` | 2026-08-03T11:43:31.8055591+02:00 | 3054 |
| `results/witness_decomposition.txt` | 2026-08-03T11:50:30.8771644+02:00 | 1284 |
| `results/witness_mixed_volume_profile.txt` | 2026-08-13T20:25:50.4101986+02:00 | 820 |
| `runs/2026-09-17T170232Z_phase1_definitions/phase1_definitions_suite.out.txt` | 2026-09-18T00:03:49.5128233+02:00 | 8070 |
| `runs/2026-09-17T170232Z_phase1_definitions/phase1_definitions_suite.py` | 2026-09-18T00:03:49.5128233+02:00 | 36919 |
| `scan.ps1` | 2026-09-18T13:00:32.0411565+02:00 | 466 |
| `scripts/__pycache__/audit_phase1_anchor.cpython-313.pyc` | 2026-09-16T16:57:01.1237842+02:00 | 8241 |
| `scripts/_run_all.sh` | 2026-09-17T22:21:38.0424699+02:00 | 1641 |
| `scripts/_runner.sh` | 2026-09-17T22:11:17.4259749+02:00 | 1695 |
| `scripts/audit_phase1_anchor.py` | 2026-09-16T16:56:09.2993666+02:00 | 4200 |
| `scripts/docker/run_m2.ps1` | 2026-08-03T15:16:07.5249631+02:00 | 1127 |
| `scripts/julia/_check_env.jl` | 2026-09-17T22:09:10.9805688+02:00 | 499 |
| `scripts/julia/_install_nemo.jl` | 2026-08-14T15:45:08.9381781+02:00 | 45 |
| `scripts/julia/_paths.jl` | 2026-08-03T11:41:23.8679292+02:00 | 277 |
| `scripts/julia/_test_lll_fit.jl` | 2026-08-14T15:45:08.9401934+02:00 | 268 |
| `scripts/julia/analyze_arc_circle_asymmetry.jl` | 2026-08-14T15:45:08.9609987+02:00 | 6749 |
| `scripts/julia/analyze_locus_grid.jl` | 2026-08-03T14:44:21.9677677+02:00 | 2024 |
| `scripts/julia/audit_argmax_task2.jl` | 2026-08-03T12:08:15.6218649+02:00 | 5615 |
| `scripts/julia/audit_clique_pipeline.jl` | 2026-08-14T15:45:08.9401934+02:00 | 11117 |
| `scripts/julia/certify_fourth_mub_witness_nwit1_multi_lambda.jl` | 2026-08-14T15:45:08.9381781+02:00 | 7202 |
| `scripts/julia/certify_fourth_mub_witness_nwit1.jl` | 2026-08-14T15:45:08.9401934+02:00 | 8787 |
| `scripts/julia/certify_fourth_mub_witness.jl` | 2026-08-14T15:45:08.9381781+02:00 | 6562 |
| `scripts/julia/certify_nwit1_all_cliques_four_classes.jl` | 2026-08-14T15:45:08.9381781+02:00 | 3154 |
| `scripts/julia/certify_nwit1_all_cliques_lambdapi.jl` | 2026-08-14T16:18:13.8725495+02:00 | 3206 |
| `scripts/julia/certify_smoke.jl` | 2026-08-03T11:41:54.5218132+02:00 | 128 |
| `scripts/julia/check_anchors.jl` | 2026-08-03T11:41:55.2847924+02:00 | 689 |
| `scripts/julia/compare_dita_coords.jl` | 2026-08-03T11:41:53.6376704+02:00 | 518 |
| `scripts/julia/csv_reconciliation.jl` | 2026-08-03T15:56:57.4136364+02:00 | 10399 |
| `scripts/julia/diagnose_mobius_grid.jl` | 2026-09-14T21:03:46.6250141+02:00 | 974 |
| `scripts/julia/diagnose_theta0_limit.jl` | 2026-09-14T21:11:08.2614900+02:00 | 1419 |
| `scripts/julia/dita_lambda_fourth_dense.jl` | 2026-08-03T15:33:38.4522397+02:00 | 4624 |
| `scripts/julia/drop_detection_stress_test.jl` | 2026-08-14T15:45:08.9609987+02:00 | 1486 |
| `scripts/julia/enumerate_d0_cliques.jl` | 2026-08-15T14:32:39.3009513+02:00 | 8640 |
| `scripts/julia/exact_karlsson_dita_H.jl` | 2026-08-15T14:33:13.3862427+02:00 | 8356 |
| `scripts/julia/export_w1_exact.jl` | 2026-08-15T14:48:05.6756617+02:00 | 7264 |
| `scripts/julia/extended_axis_sweep.jl` | 2026-08-03T14:36:42.2451398+02:00 | 7092 |
| `scripts/julia/f6_arc_boundary_tight.jl` | 2026-09-17T18:58:31.6815801+02:00 | 4244 |
| `scripts/julia/f6_boundary_reverify.jl` | 2026-08-03T14:57:33.1180084+02:00 | 3335 |
| `scripts/julia/finalize_degen1845_csv.jl` | 2026-08-14T15:45:08.9609987+02:00 | 1720 |
| `scripts/julia/formalize_gauge_lemmas.jl` | 2026-08-06T14:23:47.3258084+02:00 | 3507 |
| `scripts/julia/fourth_mub_locus_sweep.jl` | 2026-08-03T13:33:41.2326132+02:00 | 4095 |
| `scripts/julia/gauge_equivalence_analysis.jl` | 2026-08-03T13:30:26.5602355+02:00 | 4170 |
| `scripts/julia/identify_d0_in_karlsson.jl` | 2026-08-14T17:57:58.4127961+02:00 | 6126 |
| `scripts/julia/inspect_witness_coeffs.jl` | 2026-09-17T20:54:15.0987045+02:00 | 1909 |
| `scripts/julia/install_oscar_probe.jl` | 2026-08-06T14:23:47.3394898+02:00 | 1158 |
| `scripts/julia/investigate_candidate.jl` | 2026-08-06T14:23:47.3394898+02:00 | 8189 |
| `scripts/julia/investigate_circulant_match_anomaly.jl` | 2026-08-14T15:45:08.9381781+02:00 | 6379 |
| `scripts/julia/lambda_periodicity_dita.jl` | 2026-08-03T13:55:12.6611368+02:00 | 7671 |
| `scripts/julia/locus_classification.jl` | 2026-08-03T15:34:59.6800346+02:00 | 10494 |
| `scripts/julia/locus_dimension_confirm.jl` | 2026-08-03T14:45:17.0380147+02:00 | 4131 |
| `scripts/julia/locus_geometry_probes.jl` | 2026-08-06T14:23:47.3394898+02:00 | 3427 |
| `scripts/julia/locus_phi_sweep_full_circle.jl` | 2026-08-03T14:54:15.0543874+02:00 | 4458 |
| `scripts/julia/nid_probe.jl` | 2026-08-03T13:30:34.8966814+02:00 | 5502 |
| `scripts/julia/pi3_float64_smoke.jl` | 2026-09-17T17:39:01.5228558+02:00 | 768 |
| `scripts/julia/precision_f6_bounded.jl` | 2026-09-17T17:31:25.8116715+02:00 | 7667 |
| `scripts/julia/precision_f6_launch.sh` | 2026-09-17T17:31:32.6318081+02:00 | 689 |
| `scripts/julia/profile_witness_mixed_volume.jl` | 2026-08-14T15:45:08.9647700+02:00 | 2759 |
| `scripts/julia/reconstruct_b3_algebraic.jl` | 2026-08-14T15:45:08.9226262+02:00 | 17110 |
| `scripts/julia/refine_pi3_256.jl` | 2026-09-17T18:51:57.7958737+02:00 | 9109 |
| `scripts/julia/retag_512_completeness.jl` | 2026-08-03T11:42:02.3211414+02:00 | 6142 |
| `scripts/julia/run_benchmarks.jl` | 2026-09-14T21:24:36.2274958+02:00 | 5786 |
| `scripts/julia/run_fourth_mub_hpc.jl` | 2026-09-15T22:02:01.4896298+02:00 | 1206 |
| `scripts/julia/run_fourth_test_anchors.jl` | 2026-08-03T12:05:20.9284519+02:00 | 1103 |
| `scripts/julia/run_task1_investigate.jl` | 2026-08-03T11:41:53.4144075+02:00 | 2646 |
| `scripts/julia/run_task1_mobius.jl` | 2026-08-03T11:41:54.3670197+02:00 | 5158 |
| `scripts/julia/run_top50_degen.jl` | 2026-08-06T14:23:47.3394898+02:00 | 2472 |
| `scripts/julia/search_special_loci.jl` | 2026-08-13T01:15:03.3338599+02:00 | 18512 |
| `scripts/julia/smoke_test.jl` | 2026-08-03T11:41:55.8599582+02:00 | 863 |
| `scripts/julia/sweep_karlsson.jl` | 2026-08-03T11:42:08.4033960+02:00 | 3231 |
| `scripts/julia/symbolic_elimination.jl` | 2026-08-14T17:55:47.1040469+02:00 | 15016 |
| `scripts/julia/symbolic_probe.jl` | 2026-08-03T11:42:08.1436398+02:00 | 2525 |
| `scripts/julia/third_mub_locus_dimension.jl` | 2026-08-03T12:23:31.2010030+02:00 | 3361 |
| `scripts/julia/third_mub_locus_sampling.jl` | 2026-08-03T12:42:26.3162569+02:00 | 9121 |
| `scripts/julia/validate_perH_pipeline.jl` | 2026-08-03T11:41:55.8928117+02:00 | 2835 |
| `scripts/julia/validate_third_mub_candidates.jl` | 2026-08-03T14:50:19.5269408+02:00 | 11292 |
| `scripts/julia/verify_chm_b3_transfer.jl` | 2026-08-14T15:45:08.9262690+02:00 | 6763 |
| `scripts/julia/verify_closed_form_d0.jl` | 2026-08-15T14:35:36.3367010+02:00 | 7624 |
| `scripts/julia/verify_dita_construction.jl` | 2026-08-06T14:23:47.3394898+02:00 | 1592 |
| `scripts/julia/verify_extra_third_hit.jl` | 2026-08-03T12:33:30.9941577+02:00 | 804 |
| `scripts/julia/verify_field_degree_w1.jl` | 2026-08-15T14:55:38.7176635+02:00 | 1622 |
| `scripts/julia/verify_nwit1_followup.jl` | 2026-08-14T15:45:08.9381781+02:00 | 4701 |
| `scripts/julia/verify_nwit1_referee_checks.jl` | 2026-08-14T15:45:08.9381781+02:00 | 5494 |
| `scripts/julia/verify_priority_hits.jl` | 2026-08-03T13:14:31.7116720+02:00 | 2339 |
| `scripts/julia/verify_quick_run_task1.jl` | 2026-08-03T11:56:18.0929034+02:00 | 3649 |
| `scripts/julia/witness_decomposition.jl` | 2026-08-03T11:43:58.9875796+02:00 | 5275 |
| `scripts/python/__pycache__/chm_equivalence.cpython-313.pyc` | 2026-08-15T14:33:42.9508363+02:00 | 10345 |
| `scripts/python/__pycache__/degeneracy_scan.cpython-313.pyc` | 2026-08-03T11:45:39.9383477+02:00 | 10087 |
| `scripts/python/__pycache__/hadamard6.cpython-313.pyc` | 2026-08-03T11:45:40.3481591+02:00 | 3202 |
| `scripts/python/__pycache__/karlsson_k6_3.cpython-313.pyc` | 2026-08-03T11:45:40.3529975+02:00 | 7821 |
| `scripts/python/audit_karlsson_variants.py` | 2026-08-02T22:05:46.0815231+02:00 | 6214 |
| `scripts/python/audit_latex_warnings.py` | 2026-08-13T01:15:03.3358683+02:00 | 2690 |
| `scripts/python/bundle_main_theorems.py` | 2026-08-16T18:06:01.9895109+02:00 | 2504 |
| `scripts/python/bundle_methods_audit.py` | 2026-08-16T18:06:05.0188219+02:00 | 2499 |
| `scripts/python/check_paper_latex.py` | 2026-08-16T18:06:06.9991076+02:00 | 1518 |
| `scripts/python/chm_equivalence.py` | 2026-08-15T14:32:53.0836149+02:00 | 7010 |
| `scripts/python/degeneracy_scan.py` | 2026-08-03T11:51:21.2262418+02:00 | 6381 |
| `scripts/python/finish_project_audit.py` | 2026-08-13T01:15:03.3388662+02:00 | 7594 |
| `scripts/python/hadamard6.py` | 2026-08-02T22:05:08.4322568+02:00 | 1546 |
| `scripts/python/identify_d0_in_karlsson.py` | 2026-08-15T14:32:55.8486706+02:00 | 6448 |
| `scripts/python/karlsson_k6_3.py` | 2026-07-28T13:45:37.8444322+02:00 | 5671 |
| `scripts/python/validate_third_mub_table.py` | 2026-08-03T12:05:07.5391143+02:00 | 3748 |
| `scripts/wolfram/verify_dita_exact.wl` | 2026-09-17T20:52:38.6270785+02:00 | 2584 |
| `scripts/wolfram/verify_karlsson_audit.wl` | 2026-09-17T20:45:29.9698654+02:00 | 1831 |
| `scripts/wolfram/verify_repository_mathematics.wl` | 2026-09-17T21:00:16.3636072+02:00 | 8009 |
| `scripts/wolfram/w1_D0_exact.wl` | 2026-09-17T20:43:42.4774877+02:00 | 2414 |
| `setup_julia_env.ps1` | 2026-08-03T11:42:09.6984719+02:00 | 1080 |
| `source_quotes.md` | 2026-09-17T21:22:28.7073647+02:00 | 5652 |
| `src/Benchmarks.jl` | 2026-09-14T20:42:03.9319056+02:00 | 9129 |
| `src/brierley_weigert_notes.jl` | 2026-08-14T17:55:41.2677487+02:00 | 5757 |
| `src/Certification.jl` | 2026-09-17T23:42:45.0008077+02:00 | 7394 |
| `src/Cliques.jl` | 2026-09-14T19:11:19.2770376+02:00 | 4734 |
| `src/dita_third_mub_construction.jl` | 2026-08-14T17:55:44.6967776+02:00 | 2980 |
| `src/FourthMUBHPC.jl` | 2026-09-15T22:01:29.7200365+02:00 | 25493 |
| `src/hpc/M1Hardware.jl` | 2026-09-15T19:38:41.8929559+02:00 | 3996 |
| `src/hpc/M2Mubness.jl` | 2026-09-15T19:38:21.3367179+02:00 | 8087 |
| `src/hpc/M3Search.jl` | 2026-09-15T19:47:19.0733628+02:00 | 9008 |
| `src/karlsson_gauge_only.jl` | 2026-08-06T14:23:47.3394898+02:00 | 1845 |
| `src/Karlsson.jl` | 2026-09-18T12:18:29.0283291+02:00 | 16589 |
| `src/LiangChenLongQiu_CHM_Transcription.jl` | 2026-09-18T12:31:12.1983736+02:00 | 6934 |
| `src/mub_zauner_6d_liang_chen.jl` | 2026-09-17T23:42:53.4457970+02:00 | 34109 |
| `src/MubSearch.jl` | 2026-09-14T21:17:04.8452090+02:00 | 2796 |
| `src/Pool.jl` | 2026-09-14T21:23:50.0228534+02:00 | 13182 |
| `src/Provenance.jl` | 2026-09-14T22:40:38.1240728+02:00 | 4558 |
| `src/StaticMUBKernels.jl` | 2026-09-15T21:59:15.3867661+02:00 | 9605 |
| `symbolic_export/fourth_mub_reduced_Dita_elim.m2` | 2026-08-13T01:15:03.3428704+02:00 | 10260 |
| `symbolic_export/fourth_mub_reduced_Dita_exact_smoke.m2` | 2026-08-13T01:15:03.3438695+02:00 | 370 |
| `symbolic_export/fourth_mub_reduced_Dita.m2` | 2026-08-13T01:15:03.3418709+02:00 | 10046 |
| `symbolic_export/fourth_mub_reduced_F6_theta0.m2` | 2026-08-13T01:15:03.3448691+02:00 | 11241 |
| `symbolic_export/fourth_mub_w_elimination_Dita.m2` | 2026-08-13T01:15:03.3458687+02:00 | 26183 |
| `symbolic_export/fourth_mub_w_elimination_F6.m2` | 2026-08-13T01:15:03.3478717+02:00 | 30201 |
| `symbolic_export/fourth_mub_w_elimination_theta0.m2` | 2026-08-13T01:15:03.3493823+02:00 | 28264 |
| `symbolic_export/mobius_degeneracy_locus.m2` | 2026-08-06T14:42:18.3198143+02:00 | 171 |
| `symbolic_export/pool_Dita_exact.m2` | 2026-08-13T01:15:03.3513947+02:00 | 1603 |
| `symbolic_export/pool_Dita.m2` | 2026-08-13T01:15:03.3503938+02:00 | 1590 |
| `symbolic_export/pool_F6_theta0.m2` | 2026-08-13T01:15:03.3533971+02:00 | 1818 |
| `symbolic_export/pool_F6.m2` | 2026-08-13T01:15:03.3523970+02:00 | 2336 |
| `symbolic_export/theta0_fourth_mub_witness.m2` | 2026-08-13T01:15:03.3543995+02:00 | 267 |
| `symbolic_export/w1_D0_dim_degree.m2` | 2026-08-15T14:48:55.1183882+02:00 | 2215 |
| `symbolic_export/w1_D0_exact.m2` | 2026-08-15T14:48:49.2425430+02:00 | 3283 |
| `symbolic_export/w1_D0_groebner.m2` | 2026-08-15T14:50:11.7434747+02:00 | 1307 |
| `test/norm_detection.jl` | 2026-09-18T01:12:40.8937661+02:00 | 11638 |
| `test/runtests.jl` | 2026-09-14T21:16:53.3523141+02:00 | 12277 |
| `tindall_arxiv_records.csv` | 2026-09-18T12:39:04.3931681+02:00 | 9915 |
| `tindall_arxiv.xml` | 2026-09-18T12:37:21.7080133+02:00 | 90493 |
| `tindall_crossref_records.csv` | 2026-09-18T12:39:04.3215956+02:00 | 4179 |
| `tindall_crossref.json` | 2026-09-18T12:37:19.8155082+02:00 | 535828 |
| `tindall_openalex_records.csv` | 2026-09-18T12:39:04.1559542+02:00 | 3 |
| `tindall_openalex.json` | 2026-09-18T12:37:12.8853722+02:00 | 937690 |
| `tindall_orcid_records.csv` | 2026-09-18T12:39:04.3613307+02:00 | 3 |
| `tindall_orcid.json` | 2026-09-18T12:37:23.7888305+02:00 | 40336 |
| `transcribe_C2CHM_source.html` | 2026-09-18T13:15:46.9786175+02:00 | 457009 |


## Appendix C - complete Git commit log (all refs)

```text
2b9a9e3d96c444f2fbc5baded87f6e5895383a33 | 2026-09-22T23:53:31+02:00 | adian007 | Agent host session 6837c324-0387-4763-9000-2785bc773ca3 - baseline checkpoint
54fd2c963e460fcc7b2e1dc233335a610aa02cd6 | 2026-09-22T23:51:04+02:00 | adian007 | Agent host session 19499e7a-a4f3-4ff3-8ce3-7743f931d8b8 - turn 2
dac35adb9217dfc86d8761b5c9cef5764340b08f | 2026-09-22T23:51:04+02:00 | adian007 | Agent host session 19499e7a-a4f3-4ff3-8ce3-7743f931d8b8 - turn 2 start
32d889c968feef75db95715b8304bcb9363743b0 | 2026-09-22T23:50:09+02:00 | adian007 | Agent host session 19499e7a-a4f3-4ff3-8ce3-7743f931d8b8 - turn 1
84a257e8f9a643e538581ef3998f63870a2f7219 | 2026-09-22T23:50:09+02:00 | adian007 | Agent host session 19499e7a-a4f3-4ff3-8ce3-7743f931d8b8 - turn 1 start
4abde4f8253bb951fb075a56cede23084ad0b13a | 2026-09-22T23:46:33+02:00 | adian007 | Agent host session ed97ce84-06ae-414f-bfe4-f7758823eaf8 - baseline checkpoint
e946b198c1bd0a9c7c8848c8dbe5b5746007ec19 | 2026-09-22T23:35:23+02:00 | adian007 | Agent host session 19499e7a-a4f3-4ff3-8ce3-7743f931d8b8 - baseline checkpoint
63ef02f3388ede1064d9f062b4ea7045fd9a6367 | 2026-09-20T22:53:43+02:00 | adian007 | Agent host session 089eb194-b5be-4b57-a45b-48f9eef2f6b9 - baseline checkpoint
7ea7a8032a70d84aca36880df2afff4ee865a332 | 2026-09-20T22:49:14+02:00 | adian007 | third commit
19786e5f924bb027834f2085ec5d0a7c416b53fe | 2026-09-20T21:59:37+02:00 | adian007 | cline checkpoint session=1789724580483_a97wo run=6
1a976f1b064cc7357a8d9431b3de62592f74e26e | 2026-09-20T21:59:37+02:00 | adian007 | untracked files on cline checkpoint
96eb4140e9e24cabe5d40d5af06507928d18489f | 2026-09-20T21:59:36+02:00 | adian007 | index on main: 7bf60d4 idk commit
38fa6b1102babfb589af44e3f5ba541835a38163 | 2026-09-20T21:59:15+02:00 | adian007 | cline checkpoint session=1789724580483_a97wo run=5
23ee064f5a9ec88d8ff4d83db6271939b1417433 | 2026-09-20T21:59:15+02:00 | adian007 | untracked files on cline checkpoint
b9abf47a20367ddcfca1edadb264f6f7c2c3b391 | 2026-09-20T21:59:14+02:00 | adian007 | index on main: 7bf60d4 idk commit
95d22cbd159f42db3f7c2a855e7a660cb8094afe | 2026-09-18T14:53:13+02:00 | adian007 | cline checkpoint session=1789730083230_h8bb1 run=5
8b00e2efe8d972d52d3cc2550f77ce471431f28a | 2026-09-18T14:53:12+02:00 | adian007 | index on main: 7bf60d4 idk commit
2c5a52cddc1de176d4111df98609319e23c8bb84 | 2026-09-18T14:53:12+02:00 | adian007 | untracked files on cline checkpoint
3a0371f40b7939b727b54ba25b1bfb2988541b17 | 2026-09-18T13:56:08+02:00 | adian007 | cline checkpoint session=1789730083230_h8bb1 run=4
f1f40457129d420e89454b93ae745a738929b836 | 2026-09-18T13:56:08+02:00 | adian007 | untracked files on cline checkpoint
4fa7a793c6af1058fd74a7d367fc02d958c213b3 | 2026-09-18T13:56:07+02:00 | adian007 | index on main: 7bf60d4 idk commit
ee02384fd129cf41d654312e3edee0b5705e9806 | 2026-09-18T13:52:59+02:00 | adian007 | cline checkpoint session=1789730083230_h8bb1 run=3
dd56b634b89c4e56032801b19d096bf5e6768ad8 | 2026-09-18T13:52:59+02:00 | adian007 | untracked files on cline checkpoint
44fabc3418dae0a44c615c9fe5d11d175741d2a2 | 2026-09-18T13:52:58+02:00 | adian007 | index on main: 7bf60d4 idk commit
a86317e31c343e0baf50fc11252f3724425e7475 | 2026-09-18T13:38:22+02:00 | adian007 | cline checkpoint session=1789730083230_h8bb1 run=2
f641cff453bfcf3eab99cbda0dc6afc9e961da56 | 2026-09-18T13:38:22+02:00 | adian007 | untracked files on cline checkpoint
42fa04b95edcdeb7dba13cabc2c70c7e3cc86036 | 2026-09-18T13:38:21+02:00 | adian007 | index on main: 7bf60d4 idk commit
693dace46c8f2437add44780fc56c6d0f162a964 | 2026-09-18T13:17:59+02:00 | adian007 | cline checkpoint session=1789730083230_h8bb1 run=1
b1c41b1dc95c521395d782b500ee1a5a23f43e7d | 2026-09-18T13:17:58+02:00 | adian007 | untracked files on cline checkpoint
de8174b10afe02f01fc09b881594be93b088ab92 | 2026-09-18T13:17:57+02:00 | adian007 | index on main: 7bf60d4 idk commit
06e306673e6147cd27e94ed0e4caa593992c1703 | 2026-09-18T12:59:53+02:00 | adian007 | cline checkpoint session=1789724583775_u8hhb run=4
a73278daaf08a7070ac31d4a82830947dc7184cc | 2026-09-18T12:59:53+02:00 | adian007 | untracked files on cline checkpoint
35ff2606232735baf4dda9427c711ec255a94a42 | 2026-09-18T12:59:52+02:00 | adian007 | index on main: 7bf60d4 idk commit
a270340ca9a707b57467283b43701ba10dae1c2e | 2026-09-18T12:35:45+02:00 | adian007 | cline checkpoint session=1789724580483_a97wo run=4
ae5f727b778a13c612ebf3facded69d14fd61e2b | 2026-09-18T12:35:44+02:00 | adian007 | index on main: 7bf60d4 idk commit
e5beae7d55b5762fea300947f8b9320feb7357c3 | 2026-09-18T12:35:44+02:00 | adian007 | untracked files on cline checkpoint
00cf96da0d37f38072c2b6db09479ec04d67c818 | 2026-09-18T12:34:57+02:00 | adian007 | cline checkpoint session=1789724580483_a97wo run=3
5c2bb402fd278dc41be6b3d1a38cfdb56460de77 | 2026-09-18T12:34:57+02:00 | adian007 | index on main: 7bf60d4 idk commit
5b6a49ff690f7a23715bcd22382076108a1f20e8 | 2026-09-18T12:34:57+02:00 | adian007 | untracked files on cline checkpoint
cd529d5fe6c7157be058b8308d6c16431da97f28 | 2026-09-18T12:34:45+02:00 | adian007 | cline checkpoint session=1789724577275_9zhho run=3
d1855c13056911a8159a2602e12d72a4bdd6184b | 2026-09-18T12:34:45+02:00 | adian007 | untracked files on cline checkpoint
d7f50191344a75b3df0a4c127e9795357912e9db | 2026-09-18T12:34:44+02:00 | adian007 | index on main: 7bf60d4 idk commit
6e6b71cd9248a4a14e0527d49573624ecc04386f | 2026-09-18T12:34:37+02:00 | adian007 | cline checkpoint session=1789724583775_u8hhb run=3
6f7839dccb443b8c06f3747a39ce9fb9606e9609 | 2026-09-18T12:34:37+02:00 | adian007 | index on main: 7bf60d4 idk commit
168a551014f3e5846f6a30ecf2495d9fd834c62b | 2026-09-18T12:34:37+02:00 | adian007 | untracked files on cline checkpoint
4ee043d1037554f22d26e01705affa29c4e9e777 | 2026-09-18T12:34:33+02:00 | adian007 | cline checkpoint session=1789724580483_a97wo run=2
70b2e7a1a2ff36b49fdc5f179c3481c7f8f0c128 | 2026-09-18T12:34:33+02:00 | adian007 | index on main: 7bf60d4 idk commit
4db284964d22b139a2538b07df0aad6a657c080c | 2026-09-18T12:34:33+02:00 | adian007 | untracked files on cline checkpoint
ded0979eb17bf8b892263bef75a5b66fe8818124 | 2026-09-18T12:19:44+02:00 | adian007 | cline checkpoint session=1789724577275_9zhho run=2
741eacff144362f1231c46e5d77696eeabcf6fc1 | 2026-09-18T12:19:44+02:00 | adian007 | untracked files on cline checkpoint
c8b6f4e1362646f0664001e19ce6dc8fd4f8a7ca | 2026-09-18T12:19:43+02:00 | adian007 | index on main: 7bf60d4 idk commit
4f07612c803e0ee2fef175d19a4aa7165af8de95 | 2026-09-18T12:04:02+02:00 | adian007 | cline checkpoint session=1789724583775_u8hhb run=2
191a211a69eb630bfb823fd4a194998a3c9e2fe1 | 2026-09-18T12:04:02+02:00 | adian007 | untracked files on cline checkpoint
3710adbc60822f146aaf811a3432490d816fd006 | 2026-09-18T12:04:01+02:00 | adian007 | index on main: 7bf60d4 idk commit
6a8ca9f2ef0a30b38446fbba1972dd9afbc11d7e | 2026-09-18T12:03:33+02:00 | adian007 | cline checkpoint session=1789724576355_tyayb run=5
23d8f36dfd35bc5290602900d293a1b33cc660c8 | 2026-09-18T12:03:33+02:00 | adian007 | untracked files on cline checkpoint
5211fd9ecf474d400ba828f1e5c5b6d51bc0da2f | 2026-09-18T12:03:32+02:00 | adian007 | index on main: 7bf60d4 idk commit
428902cce02077924ae5e2767003c688464ddb47 | 2026-09-18T12:03:14+02:00 | adian007 | cline checkpoint session=1789724576355_tyayb run=4
113fdb11f05256b15be291cceee0319f84045186 | 2026-09-18T12:03:14+02:00 | adian007 | untracked files on cline checkpoint
f10cee1fc8642675f3829715a4314214f2e3b225 | 2026-09-18T12:03:13+02:00 | adian007 | index on main: 7bf60d4 idk commit
ea2e262e391bc6b927896de73958970b65f37689 | 2026-09-18T12:03:01+02:00 | adian007 | cline checkpoint session=1789724576355_tyayb run=3
4ea3e2d3611d5f2b64872cda82cb57351200dc6d | 2026-09-18T12:03:01+02:00 | adian007 | index on main: 7bf60d4 idk commit
687b3fa5b4abe59ee42c10c65d7241f7ec74cb58 | 2026-09-18T12:03:01+02:00 | adian007 | untracked files on cline checkpoint
29518ade22b786abb249c9403b5c9d6246aa0408 | 2026-09-18T12:02:38+02:00 | adian007 | cline checkpoint session=1789724576355_tyayb run=2
6103c524a5f84c25cd261d013479ad80513b5e0a | 2026-09-18T12:02:37+02:00 | adian007 | untracked files on cline checkpoint
c3f0ccdde85bd554c681800448e3f3ad90e6f2d1 | 2026-09-18T12:02:36+02:00 | adian007 | index on main: 7bf60d4 idk commit
b9316150f0ae685624d1c8cde79640682a6fe33c | 2026-09-18T11:46:12+02:00 | adian007 | cline checkpoint session=1789724583775_u8hhb run=1
d058cabb7c8cd02a19a3779fabe5ec2cd0530520 | 2026-09-18T11:46:12+02:00 | adian007 | untracked files on cline checkpoint
acb6dac0001ac09bf51fe6bfc3f0a3f6d0c71043 | 2026-09-18T11:46:11+02:00 | adian007 | index on main: 7bf60d4 idk commit
1aab3f6074856b31a2aab6c6d7b9c664ca8af4c5 | 2026-09-18T11:46:04+02:00 | adian007 | cline checkpoint session=1789724580483_a97wo run=1
8931b3047faed011262cdaa60df94ca75c276b9a | 2026-09-18T11:46:04+02:00 | adian007 | untracked files on cline checkpoint
bcf8be0c3ef0a2b989257108c27cf220379fc221 | 2026-09-18T11:46:03+02:00 | adian007 | index on main: 7bf60d4 idk commit
7d423e0e3053dfea4f0014edd319e0bb452427fc | 2026-09-18T11:45:56+02:00 | adian007 | cline checkpoint session=1789724577275_9zhho run=1
a80d4cdf5964194978b2806059935397a1d271b9 | 2026-09-18T11:45:55+02:00 | adian007 | index on main: 7bf60d4 idk commit
4e5d1d3013a8080d9e77ecf8b70bb4fff6d37542 | 2026-09-18T11:45:55+02:00 | adian007 | untracked files on cline checkpoint
c4b1804a034b54f19701e02d3dea72782f211aa6 | 2026-09-18T11:45:47+02:00 | adian007 | cline checkpoint session=1789724576355_tyayb run=1
1ebf65440d05f6ecbf586f21771dd5d0058c6cf0 | 2026-09-18T11:45:47+02:00 | adian007 | untracked files on cline checkpoint
2fee9ff1f7e2d19ddea3e47c4e03cf0ac59de2eb | 2026-09-18T11:45:46+02:00 | adian007 | index on main: 7bf60d4 idk commit
dfb66f95ff3ed9eb18a9892bfaaa24df2c773568 | 2026-09-18T01:09:13+02:00 | adian007 | cline checkpoint session=1789685138465_h9xpf run=15
106b4650b58f6cc255ca46c6c2383f22dd1395dc | 2026-09-18T01:09:12+02:00 | adian007 | index on main: 7bf60d4 idk commit
d98c32f6de354df0bbe2c60af14228922bf7598e | 2026-09-18T01:09:12+02:00 | adian007 | untracked files on cline checkpoint
d6d315d5b8809ce3a32bd48c07df6aa6d639a8da | 2026-09-18T01:06:02+02:00 | adian007 | cline checkpoint session=1789685138465_h9xpf run=14
e43980676d4621ba1e9e3ea661e939a6f741a8c4 | 2026-09-18T01:06:01+02:00 | adian007 | index on main: 7bf60d4 idk commit
2b4b5235bfc67c095f1f4080bbf83aae69464585 | 2026-09-18T01:06:01+02:00 | adian007 | untracked files on cline checkpoint
8bd84d1dcde5dfc89bf0e2f753f499d77b299cc9 | 2026-09-18T00:59:29+02:00 | adian007 | cline checkpoint session=1789685138465_h9xpf run=13
07ac746bd07dd5a2e5dc885873b46331e860f0d8 | 2026-09-18T00:59:29+02:00 | adian007 | untracked files on cline checkpoint
1ebc31178e161e99596fe015bbd0ae393fbed914 | 2026-09-18T00:59:28+02:00 | adian007 | index on main: 7bf60d4 idk commit
62f8000926ed3555dd5f60b30a421b6ef13bf595 | 2026-09-18T00:58:51+02:00 | adian007 | cline checkpoint session=1789685138465_h9xpf run=12
dd8d64cede32cc23a5e31b06e68dc363692e70ba | 2026-09-18T00:58:51+02:00 | adian007 | index on main: 7bf60d4 idk commit
c9efabdf5cf45fb81852bba172595af59eb47c25 | 2026-09-18T00:58:51+02:00 | adian007 | untracked files on cline checkpoint
6cc80730901b37111335ff9227f98d57308af8fe | 2026-09-18T00:54:51+02:00 | adian007 | cline checkpoint session=1789684899609_guz7c run=5
d68da4482dcb756aad5bf62f8966f1d3347057ef | 2026-09-18T00:54:50+02:00 | adian007 | index on main: 7bf60d4 idk commit
0ea0553945a9ebef1a9fbadc08468d569693c742 | 2026-09-18T00:54:50+02:00 | adian007 | untracked files on cline checkpoint
6fafdb549a5e397bd04cf148549eb9af600f756c | 2026-09-18T00:53:17+02:00 | adian007 | cline checkpoint session=1789684899609_guz7c run=4
914a2b28c3b2445da6cf4b326dcbc39d1fcc7765 | 2026-09-18T00:53:16+02:00 | adian007 | untracked files on cline checkpoint
448248c3e7fbde4553efff7cb8c8bbc07272ecf4 | 2026-09-18T00:53:15+02:00 | adian007 | index on main: 7bf60d4 idk commit
7dbd0bffd8703fb22477c876c48c4ec0c576024b | 2026-09-18T00:53:10+02:00 | adian007 | cline checkpoint session=1789685493217_equra run=1
efaac526eb19a3adcc951940c5932501647b3e48 | 2026-09-18T00:53:10+02:00 | adian007 | untracked files on cline checkpoint
f8783c21fdfe23cd211d19cf1b4749006c22ee97 | 2026-09-18T00:53:09+02:00 | adian007 | index on main: 7bf60d4 idk commit
0515b4afa0216daea8210180f91d6090e4db3ace | 2026-09-18T00:49:48+02:00 | adian007 | cline checkpoint session=1789684899609_guz7c run=3
4b43e04f12c7efc25021e460c05ec92c125f04ae | 2026-09-18T00:49:48+02:00 | adian007 | untracked files on cline checkpoint
b0ed58b3fadb2081bb91c6a15537e38d77217947 | 2026-09-18T00:49:47+02:00 | adian007 | index on main: 7bf60d4 idk commit
96143a33345bf758187da597430307bf09440b8c | 2026-09-18T00:46:01+02:00 | adian007 | cline checkpoint session=1789685138465_h9xpf run=1
d90d036cd6090bcd0005622e8d33571bbcd0e979 | 2026-09-18T00:46:00+02:00 | adian007 | untracked files on cline checkpoint
e064c0854abb2ad46743d2c6639bae84438348e2 | 2026-09-18T00:45:59+02:00 | adian007 | index on main: 7bf60d4 idk commit
680aaf739c08ce62b11510a108e262228079a3a6 | 2026-09-18T00:45:40+02:00 | adian007 | cline checkpoint session=1789685138465_h9xpf run=11
de7a03746818ff86261cf419b050a33e1f55f54a | 2026-09-18T00:45:40+02:00 | adian007 | untracked files on cline checkpoint
037ebcbc55a3b918ad09b65a0090e58897fb8cae | 2026-09-18T00:45:39+02:00 | adian007 | index on main: 7bf60d4 idk commit
ff3d4fa002f20600f3b77d835cd7862cca243bc4 | 2026-09-18T00:43:43+02:00 | adian007 | cline checkpoint session=1789684899609_guz7c run=2
a8fe0e1c7dc682e81eef6c1aa2e3cdde833e7847 | 2026-09-18T00:43:42+02:00 | adian007 | untracked files on cline checkpoint
494e3277a8083ed75995d5c96130afbb7d4f47ae | 2026-09-18T00:43:41+02:00 | adian007 | index on main: 7bf60d4 idk commit
4701098c699ac0cc375bc1d16d3ee41d49587b53 | 2026-09-18T00:42:32+02:00 | adian007 | cline checkpoint session=1789684899609_guz7c run=1
53eb834ae261f6db02ea918bd5262685860e3d27 | 2026-09-18T00:42:32+02:00 | adian007 | untracked files on cline checkpoint
8f50a9edd959a19d62475e61cb6f12d32dac0ea1 | 2026-09-18T00:42:31+02:00 | adian007 | index on main: 7bf60d4 idk commit
3f5daaf57ce0ac7fa78c06316741b4247c83291b | 2026-09-18T00:39:24+02:00 | adian007 | cline checkpoint session=1789683931419_uch1w run=3
2124bf503d6c2f58fa53d744f6b7332fdd25ba7a | 2026-09-18T00:39:24+02:00 | adian007 | untracked files on cline checkpoint
5be6e9e658411141fc00f6424244748e21fe3994 | 2026-09-18T00:39:23+02:00 | adian007 | index on main: 7bf60d4 idk commit
dd080e45e92b853bbea89dec548263be0999ecf2 | 2026-09-18T00:25:58+02:00 | adian007 | cline checkpoint session=1789683931419_uch1w run=2
bfaa2cefdcb807158bcb0f529f92c14242ce0a41 | 2026-09-18T00:25:57+02:00 | adian007 | index on main: 7bf60d4 idk commit
bc389b32cffa41fbddd0df9df00d57b01766f23f | 2026-09-18T00:25:57+02:00 | adian007 | untracked files on cline checkpoint
aebbe44f217288e8acfd27f96965c6a8ac2de63c | 2026-09-18T00:25:34+02:00 | adian007 | cline checkpoint session=1789683931419_uch1w run=1
765e02b22377905eb51c9131f0a042fe45e7e8fe | 2026-09-18T00:25:33+02:00 | adian007 | index on main: 7bf60d4 idk commit
259ca3c4a0cc0779e19bbfddcda18f8480df7513 | 2026-09-18T00:25:33+02:00 | adian007 | untracked files on cline checkpoint
babd489e96913cd6969da53453e2fcf23a93e348 | 2026-09-18T00:24:57+02:00 | adian007 | cline checkpoint session=1789683842457_ziksu run=2
be75c874c2997a7a5bca6f4e05b6d2ffd5f0b39e | 2026-09-18T00:24:57+02:00 | adian007 | untracked files on cline checkpoint
3b1aea7e49df6df24667d6aa33dee3df94d23351 | 2026-09-18T00:24:56+02:00 | adian007 | index on main: 7bf60d4 idk commit
67aa4e4cdfc1d332a7f14b0a850038059b80c187 | 2026-09-18T00:24:17+02:00 | adian007 | cline checkpoint session=1789683842457_ziksu run=1
7e87302a4d5712bccd76df760af2872a29e8f5ba | 2026-09-18T00:24:17+02:00 | adian007 | untracked files on cline checkpoint
316366ec497c2be929b0c3c47ca10d862bc7f99d | 2026-09-18T00:24:16+02:00 | adian007 | index on main: 7bf60d4 idk commit
8ced12d080d222e4abc32bc0146950edae5c1171 | 2026-09-18T00:23:51+02:00 | adian007 | cline checkpoint session=1789683790648_cbkaz run=1
487807aee5cbae09024b8c54babed053f9354b84 | 2026-09-18T00:23:51+02:00 | adian007 | untracked files on cline checkpoint
22e2541e88e5c2374f8f75120f0a5ec26efd7e14 | 2026-09-18T00:23:50+02:00 | adian007 | index on main: 7bf60d4 idk commit
c798f87865baae0f2e9371c6d04706d417fd53a5 | 2026-09-18T00:06:21+02:00 | adian007 | cline checkpoint session=1789680383019_laqi9 run=13
fcbeb608cb0106acb07efba4189d1e2d22fd414c | 2026-09-18T00:06:21+02:00 | adian007 | untracked files on cline checkpoint
a1ae6486d88724184170714ec8419947cec53617 | 2026-09-18T00:06:20+02:00 | adian007 | index on main: 7bf60d4 idk commit
075dc96736b9882de7d0e0b296bf013249c415f1 | 2026-09-18T00:01:28+02:00 | adian007 | cline checkpoint session=1789680383019_laqi9 run=12
8ed050bad6e28b85ff3f5fe1da96743966f21214 | 2026-09-18T00:01:28+02:00 | adian007 | untracked files on cline checkpoint
013a46f0b336f272451d1cd1df3a08ad98e84b48 | 2026-09-18T00:01:27+02:00 | adian007 | index on main: 7bf60d4 idk commit
b0906307e80e9f2d5d215602a1b25a269f9b6e9f | 2026-09-18T00:00:24+02:00 | adian007 | cline checkpoint session=1789680671082_28vew run=10
8fb076622ec62238c24236a96b658ac137d7349f | 2026-09-18T00:00:24+02:00 | adian007 | untracked files on cline checkpoint
15e651560b2ecb704841e56dd112b250a4b5126b | 2026-09-18T00:00:23+02:00 | adian007 | index on main: 7bf60d4 idk commit
13e53ffcf8b473bdff9cd29486ba7458b65a5cb1 | 2026-09-17T23:52:07+02:00 | adian007 | cline checkpoint session=1789680383019_laqi9 run=1
b76b6a88059565bb38989e17c57b5064b0e1c9d9 | 2026-09-17T23:52:07+02:00 | adian007 | untracked files on cline checkpoint
fae613c872160ac1a28bb674407dc80b0cf4de75 | 2026-09-17T23:52:06+02:00 | adian007 | index on main: 7bf60d4 idk commit
409c7390f1673e124e3ade13acb1a5c5c55b7aed | 2026-09-17T23:51:48+02:00 | adian007 | cline checkpoint session=1789680383019_laqi9 run=11
0b070820d41d279f074159dd72e52953985ee818 | 2026-09-17T23:51:47+02:00 | adian007 | index on main: 7bf60d4 idk commit
3edfb6e9b759cdd083fc8d25f8c6dafb99b30bf7 | 2026-09-17T23:51:47+02:00 | adian007 | untracked files on cline checkpoint
b0bac75cbd89c846eba8f252c7569bd79b124a91 | 2026-09-17T23:32:37+02:00 | adian007 | cline checkpoint session=1789680671082_28vew run=1
eb0a6736a2945c62bdb08ad3ce0980449383c8fa | 2026-09-17T23:32:36+02:00 | adian007 | index on main: 7bf60d4 idk commit
14fbb8773e0c8656c90a1498570b67946db79c7d | 2026-09-17T23:32:36+02:00 | adian007 | untracked files on cline checkpoint
eab71ecf4c6a6d43928c1879fad29bb50f2530e3 | 2026-09-17T23:31:13+02:00 | adian007 | cline checkpoint session=1789680671082_28vew run=9
13921c57ebaadb2aa3fc83248aaad0c78ed2e39a | 2026-09-17T23:31:13+02:00 | adian007 | untracked files on cline checkpoint
78729c61ba0a4f17561469251cd7653371d709e2 | 2026-09-17T23:31:12+02:00 | adian007 | index on main: 7bf60d4 idk commit
5ebf2d0d95672d49d56b3f190d727a499048e090 | 2026-09-17T23:30:30+02:00 | adian007 | cline checkpoint session=1789676473943_xslik run=8
e11ad7f0fc2bfe9c7165fcd17659a6ad4cb39914 | 2026-09-17T23:30:30+02:00 | adian007 | untracked files on cline checkpoint
f7168a59c32b6de77aa22a11e19ffd5bcfc022e4 | 2026-09-17T23:30:29+02:00 | adian007 | index on main: 7bf60d4 idk commit
70d91d8c2ce1546b1e54742facd1df3576f106ef | 2026-09-17T23:29:06+02:00 | adian007 | cline checkpoint session=1789680383019_laqi9 run=10
8a9c633a1b7db9d3c3da63a4b6da3a9ef547415a | 2026-09-17T23:29:05+02:00 | adian007 | untracked files on cline checkpoint
203b241b5e88cb5cad94de3f972342c667691a50 | 2026-09-17T23:29:04+02:00 | adian007 | index on main: 7bf60d4 idk commit
e4ccc9369f841ca7f714e1f4df7fb0d77092e334 | 2026-09-17T23:28:28+02:00 | adian007 | cline checkpoint session=1789676473943_xslik run=7
3d7f60fdb107620ef50bcb9fca4aca7423005692 | 2026-09-17T23:28:28+02:00 | adian007 | untracked files on cline checkpoint
23da9b665ef6ff0e1acd3bb845ce54e99c3bbf2f | 2026-09-17T23:28:27+02:00 | adian007 | index on main: 7bf60d4 idk commit
b7730ca42b47c3a2fbe37ee018375642d1ead109 | 2026-09-17T23:27:32+02:00 | adian007 | cline checkpoint session=1789676473943_xslik run=6
6e07423180cdc2215c9b2bc781c93758f86d86e0 | 2026-09-17T23:27:32+02:00 | adian007 | untracked files on cline checkpoint
c1507b08f60a5624de869125ddba965ed3375f3b | 2026-09-17T23:27:31+02:00 | adian007 | index on main: 7bf60d4 idk commit
a6bd2cc509c0c8cdebc66cc37ba9b19d1b274358 | 2026-09-17T23:26:43+02:00 | adian007 | cline checkpoint session=1789676473943_xslik run=5
6b91f76fe2ccea0a8bc062656c68a348ad8b71e9 | 2026-09-17T23:26:43+02:00 | adian007 | untracked files on cline checkpoint
66dafe6f1e9a3a8816546cdf2a016ede6e46a727 | 2026-09-17T23:26:42+02:00 | adian007 | index on main: 7bf60d4 idk commit
fa2818bb4a213278cb27aad8177297dca41ab31c | 2026-09-17T23:26:27+02:00 | adian007 | cline checkpoint session=1789680383019_laqi9 run=9
f262d6bc22c8db52d764c17f95d05e151d1d4273 | 2026-09-17T23:26:26+02:00 | adian007 | index on main: 7bf60d4 idk commit
db250db1d296c1f12d9d788374b51402913f2493 | 2026-09-17T23:26:26+02:00 | adian007 | untracked files on cline checkpoint
f4b8ae97645ec2e3369fc5274c4b9cbd44bf4cd3 | 2026-09-17T23:19:02+02:00 | adian007 | cline checkpoint session=1789676473943_xslik run=4
080ff91d6a806b29bd0b5bc28c0f236bc2948d70 | 2026-09-17T23:19:02+02:00 | adian007 | index on main: 7bf60d4 idk commit
430191627eccbcc3d0954875034111778b5eb39d | 2026-09-17T23:19:02+02:00 | adian007 | untracked files on cline checkpoint
2f8ba0590e273190a5d07a4c70d4c6d56a5f7624 | 2026-09-17T22:27:05+02:00 | adian007 | cline checkpoint session=1789670401960_rsrgu run=12
7fb951eb0ccdd82bff7389442e483da29e86d0f8 | 2026-09-17T22:27:05+02:00 | adian007 | untracked files on cline checkpoint
14786bcff73888b346ae93d712040f7848f81797 | 2026-09-17T22:27:04+02:00 | adian007 | index on main: 7bf60d4 idk commit
3ca0460a289e9e1360f7bf846fa483414fc2ebbd | 2026-09-17T22:26:44+02:00 | adian007 | cline checkpoint session=1789676473943_xslik run=3
662993cd2d8f2dce72bacc6033fe6ffaa502571f | 2026-09-17T22:26:44+02:00 | adian007 | index on main: 7bf60d4 idk commit
520e316a7df982d3cefa8f5f7e9cc8184964ca8c | 2026-09-17T22:26:44+02:00 | adian007 | untracked files on cline checkpoint
c7872582c6c0364ec2a1c244bc5e4b1deb579569 | 2026-09-17T22:26:28+02:00 | adian007 | cline checkpoint session=1789676473943_xslik run=2
59c5aa274d09fc888b60a36590cac8d50f8b8f0b | 2026-09-17T22:26:28+02:00 | adian007 | untracked files on cline checkpoint
479514f8e0e42c11a200f75ec53c9b4f9e8a255f | 2026-09-17T22:26:27+02:00 | adian007 | index on main: 7bf60d4 idk commit
9e8e8fcd377b1472dc031bb535e4632526d16922 | 2026-09-17T22:22:32+02:00 | adian007 | cline checkpoint session=1789676473943_xslik run=1
9c0c00009ae9b42cf97503a2175d6158d835c1cd | 2026-09-17T22:22:32+02:00 | adian007 | untracked files on cline checkpoint
ef92ec009ae88a91fb86f88b431d78dba6396b40 | 2026-09-17T22:22:31+02:00 | adian007 | index on main: 7bf60d4 idk commit
e591fdb4beffdd528168806840da0026f04bb3a3 | 2026-09-17T22:07:23+02:00 | adian007 | cline checkpoint session=1789670401960_rsrgu run=11
04d76bdd98a1788f2d2fe79dd0737c4d3fde7ef7 | 2026-09-17T22:07:22+02:00 | adian007 | index on main: 7bf60d4 idk commit
00f205d351dd72aa9aec3a2138980bf989c897fa | 2026-09-17T22:07:22+02:00 | adian007 | untracked files on cline checkpoint
3ad0ea4778c7ff9071df4b2373f37ed96eb5e100 | 2026-09-17T21:47:15+02:00 | adian007 | cline checkpoint session=1789670401960_rsrgu run=10
f4fde0a2ccbbe8aa88bac627fa82cb65a4620c48 | 2026-09-17T21:47:15+02:00 | adian007 | untracked files on cline checkpoint
b9a450b272df5dd216b8ad43271cc88fd584b16e | 2026-09-17T21:47:14+02:00 | adian007 | index on main: 7bf60d4 idk commit
e260d02f835d8e47ee1f9ddb188443666947681b | 2026-09-17T21:45:12+02:00 | adian007 | cline checkpoint session=1789670401960_rsrgu run=9
6a92c09cfffa159cdb5eb7446034fde7d8fcea57 | 2026-09-17T21:45:11+02:00 | adian007 | untracked files on cline checkpoint
dcaa8430e2c21486146350b19773974c5b7147be | 2026-09-17T21:45:10+02:00 | adian007 | index on main: 7bf60d4 idk commit
b38d38e23c7463e2fd5548cda01377afe4bc4524 | 2026-09-17T20:23:20+02:00 | adian007 | cline checkpoint session=1789662989700_vltw2 run=8
7c7974bb38c7a6e3e86d4c77a9fe5a4e2122fdfc | 2026-09-17T20:23:19+02:00 | adian007 | index on main: 7bf60d4 idk commit
36bb0d29c53290393d06e97c05e1a7acd533ba7c | 2026-09-17T20:23:19+02:00 | adian007 | untracked files on cline checkpoint
d9834d193c27b88584ba3cf758a638e1a43fc168 | 2026-09-17T19:55:25+02:00 | adian007 | cline checkpoint session=1789662989700_vltw2 run=7
dcbb06b1271be6c4af67c313d0b204b3af734bea | 2026-09-17T19:55:25+02:00 | adian007 | index on main: 7bf60d4 idk commit
9f67016b299ecb904ac52ad89d53f03cbe4abfc8 | 2026-09-17T19:55:25+02:00 | adian007 | untracked files on cline checkpoint
d2e1e70df84697bffac40f469cb59e25d1e87e40 | 2026-09-17T19:33:28+02:00 | adian007 | cline checkpoint session=1789662989700_vltw2 run=6
eca2c5d61094eb1c1d12160917eee4b355a0b151 | 2026-09-17T19:33:27+02:00 | adian007 | index on main: 7bf60d4 idk commit
fdd4ee702346b12ea380313a2f76424952dc78af | 2026-09-17T19:33:27+02:00 | adian007 | untracked files on cline checkpoint
0ef3e61542b8836a405a7d8fc26ecd7f32e8795c | 2026-09-17T19:12:47+02:00 | adian007 | cline checkpoint session=1789662989700_vltw2 run=5
dc467134695451f3c7707620ae13ceb4d6e27a8e | 2026-09-17T19:12:47+02:00 | adian007 | untracked files on cline checkpoint
55015845bcdf940e50f0d392334b941d9bbfda33 | 2026-09-17T19:12:46+02:00 | adian007 | index on main: 7bf60d4 idk commit
32cdfd1476deddaace9325e481ffbbf895b13877 | 2026-09-17T19:01:27+02:00 | adian007 | cline checkpoint session=1789662989700_vltw2 run=4
4a8236a9d5ef5bf55798d726b62a657f182cce41 | 2026-09-17T19:01:26+02:00 | adian007 | index on main: 7bf60d4 idk commit
e8a08fd0631e38489f3f34b0f8c561ae17d6775b | 2026-09-17T19:01:26+02:00 | adian007 | untracked files on cline checkpoint
c06f1ee319a377f76162818d6fbe7a47ddaec4b8 | 2026-09-17T18:55:48+02:00 | adian007 | cline checkpoint session=1789662989700_vltw2 run=3
22be255545ff20a72cb19709ce96668626bdc336 | 2026-09-17T18:55:48+02:00 | adian007 | untracked files on cline checkpoint
16e074f14b6132386ec18a5ed0ee3f5c5494bb5f | 2026-09-17T18:55:47+02:00 | adian007 | index on main: 7bf60d4 idk commit
c4fc7595cf178f2c78b31f0184fdad2f672dcd23 | 2026-09-17T18:52:40+02:00 | adian007 | cline checkpoint session=1789662989700_vltw2 run=2
4913032911ba749be5953cc31e0a00c9df7bd648 | 2026-09-17T18:52:40+02:00 | adian007 | untracked files on cline checkpoint
fef4360e88cf7a22477d545123ea427508df28e8 | 2026-09-17T18:52:39+02:00 | adian007 | index on main: 7bf60d4 idk commit
79f32f7f2898669007d413359166945a8f880bf5 | 2026-09-17T18:40:45+02:00 | adian007 | cline checkpoint session=1789662989700_vltw2 run=1
557a68b5a05eacba56658d2f7d8ed0b3f8c584f3 | 2026-09-17T18:40:45+02:00 | adian007 | untracked files on cline checkpoint
b66a411b7957cd70e0feac7874c5b9f456133161 | 2026-09-17T18:40:44+02:00 | adian007 | index on main: 7bf60d4 idk commit
3b18eacae4937e3c4e156437673c706eba32f94e | 2026-09-17T18:30:13+02:00 | adian007 | cline checkpoint session=session_1789657722170_ylnwh run=4
eddfd07ef9ceaf1178d34cffc420ecc149435c6e | 2026-09-17T18:30:13+02:00 | adian007 | untracked files on cline checkpoint
b9906f5155b9341a7683d2cbd8ca52dc89580b1a | 2026-09-17T18:30:12+02:00 | adian007 | index on main: 7bf60d4 idk commit
139006df184cdab667d9d9331b71d1bd99d85b93 | 2026-09-17T17:17:59+02:00 | adian007 | cline checkpoint session=session_1789657722170_ylnwh run=3
b3145b703e51ffe13854aeb1c48d94d7a7ccfa6f | 2026-09-17T17:17:58+02:00 | adian007 | index on main: 7bf60d4 idk commit
d5456ff2088c8cd90ffc2e8cf61b05a41454ddc9 | 2026-09-17T17:17:58+02:00 | adian007 | untracked files on cline checkpoint
65c4b4aca715c74f4f586d5e9c2af06cbe3f63f5 | 2026-09-17T17:15:11+02:00 | adian007 | cline checkpoint session=session_1789657722170_ylnwh run=2
3cbf69fb604e4a0a8c06e35ec62f7a1d7d487839 | 2026-09-17T17:15:10+02:00 | adian007 | index on main: 7bf60d4 idk commit
d23e0831a70c33d4d3e471259765bcb90154aecf | 2026-09-17T17:15:10+02:00 | adian007 | untracked files on cline checkpoint
4e70c546e810c4b525456b80a73d46ed16968ef4 | 2026-09-17T17:08:44+02:00 | adian007 | cline checkpoint session=session_1789657722170_ylnwh run=1
32092f062b0658e7ad13ffde2b37c20d3c0b7053 | 2026-09-17T17:08:43+02:00 | adian007 | index on main: 7bf60d4 idk commit
71ae4b1193d1f91050b6a7a49b0d2f3c99e6403d | 2026-09-17T17:08:43+02:00 | adian007 | untracked files on cline checkpoint
7bf60d4aaf2881030b5d376c1757706f30107ff1 | 2026-08-13T01:11:51+02:00 | adian007 | idk commit
2c60c76b0d153cdf5bfd66394f6e727270879d52 | 2026-08-13T01:15:01+02:00 | adian007 | On main: temp: deleted vscode settings
0f73361180bb14a91cd12ce714b55ce38dfc1089 | 2026-08-13T01:15:01+02:00 | adian007 | index on main: ef9100f idk commit
9dc51847f8f236e22c57d1e347602d9e98e49f79 | 2026-08-13T01:15:01+02:00 | adian007 | untracked files on main: ef9100f idk commit
ef9100fc9b8a06ef5957637e112c132f696633cd | 2026-08-13T01:11:51+02:00 | adian007 | idk commit
48761753bf7f1292778e18a2b31f16c21eabe4fb | 2026-08-06T14:30:48+02:00 | adian007 | Fix heading formatting in readme.md
6a85692acb2e13fa80143bf43ab88b6b6d16067f | 2026-08-06T14:30:28+02:00 | adian007 | Fix title formatting in README
d759c5934b2116c0ad0d3e94e7c144397e143509 | 2026-08-06T14:29:53+02:00 | adian007 | updated readme
688396fe420dbe9d047bb60413239bf424782062 | 2026-08-06T14:26:55+02:00 | adian007 | Add project source code and initial scripts
b477404c1974703651eb398f4a9b1938f499cdf4 | 2026-08-06T14:25:13+02:00 | adian007 | first commit
```

## Appendix D - earliest-to-current Git diff name-status

Earliest commit: b477404c1974703651eb398f4a9b1938f499cdf4

```text
A	.gitignore
A	.vscode/settings.json
A	findings.md
A	Manifest.toml
A	Project.toml
M	readme.md
A	reproduce.md
A	__pycache__/hadamard6.cpython-313.pyc
A	__pycache__/karlsson_k6_3.cpython-313.pyc
A	_all_jl.txt
A	_spawn_cache/fetch.ps1
A	_spawn_cache/matrefs.txt
A	_spawn_cache/paper_2110.12206.html
A	_spawn_cache/tagcontext.txt
A	bibliography_verified.csv
A	claim_audit.md
A	debug-388a52.log
A	known_results_table.md
A	literature_audit.md
A	matrix_catalogue.md
A	open_problems.md
A	paper/overleaf.md
A	paper/restructure_proposal.md
A	paper/restructure_split_status.md
A	paper/main_theorems.tex
A	paper/main_theorems_multifile.tex
A	paper/methods_audit.tex
A	paper/methods_audit_multifile.tex
A	paper/preamble_common.tex
A	paper/proofs/dita_third_mub_methods.tex
A	paper/proofs/fourth_mub_theorems.tex
A	paper/proofs/gauge_structure.tex
A	paper/proofs/locus_geometry.tex
A	parse_tindall.ps1
A	pdf/.env
A	pdf/1809.07442.pdf
A	pdf/2110.12206.pdf
A	pdf/2110.12206.txt
A	pdf/2110.12206_full.txt
A	pdf/2503.14752
A	pdf/3488559.pdf
A	pdf/78-files-data-ec5decca5ed3d6b8079e2e7e7bacc9f2-127.pdf
A	pdf/A concise guide to complex Hadamard matrices.pdf
A	pdf/Extension of the Set of Complex Hadamard Matrices of Size 8.pdf
A	pdf/MUTUALLY UNBIASED BASES AND.pdf
A	pdf/Mutually Unbiased Bases in Composite Dimensions -- A Review.pdf
A	pdf/PRXQuantum.3.010101.pdf
A	pdf/THREE-PARAMETER COMPLEX HADAMARD.pdf
A	pdf/_all_model_ids.txt
A	pdf/_anal.out
A	pdf/_anal.py
A	pdf/_anal2.py
A	pdf/_analm.py
A	pdf/_analyze_tmp.py
A	pdf/_build_page_nvidia__ai-classification-ocr-2.html
A	pdf/_build_page_nvidia__htc-ocr-2.html
A	pdf/_build_page_nvidia__llava-3.2-8b-vit.html
A	pdf/_build_page_nvidia__mumu-0.5-8b.html
A	pdf/_build_page_nvidia__nemotron-parse-2.0.html
A	pdf/_build_page_nvidia__phi-3.5-vision-instruct.html
A	pdf/_fetchbuild_dump.py
A	pdf/_fetchbuild_out.txt
A	pdf/_probe_results.txt
A	pdf/_showids.py
A	pdf/_top_files.txt
A	pdf/beauchamp_nicoara_2006.pdf
A	pdf/bell_inequalities_mub.pdf
A	pdf/brierley_weigert_grobner.pdf
A	pdf/butterley_hall_numerical.pdf
A	pdf/fetchbuild.py
A	pdf/four_most_distant_bases.pdf
A	pdf/full_list.py
A	pdf/g3.py
A	pdf/grepnode.py
A	pdf/grepparse.py
A	pdf/grepparse2.py
A	pdf/jaming_matolcsi_mora_mub6.pdf
A	pdf/karlsson_hadamards_mub6.pdf
A	pdf/list-models.py
A	pdf/matolcsi_triplets_mub.pdf
A	pdf/mconnell_zauner_evidence.pdf
A	pdf/model_report.txt
A	pdf/mub_composite_dims_review.pdf
A	pdf/parse_dump.py
A	pdf/probe.py
A	pdf/probe2.py
A	pdf/probe3.py
A	pdf/probe_models.py
A	pdf/probe_writer.py
A	pdf/search_models.py
A	research/readme.md
A	research/algebraic/algebraic_obstruction_program.md
A	research/anchors/anchor_selection.md
A	research/campaigns/campaign1/campaign1_report.md
A	research/certificates/certificate_design.md
A	research/claim_ledger/research_claim_ledger.md
A	research/claims/claim_ledger_campaigns.csv
A	research/constructions/dimension6_constructions.md
A	research/definitions/mub_dimension6_research.md
A	research/definitions/phase1_definitions_suite.out.txt
A	research/definitions/phase1_definitions_suite.py
A	research/exact_systems/agent4/exact_systems_report.md
A	research/exact_systems/agent4/systems.json
A	research/exact_systems/agent4/variable_quantifier_table.csv
A	research/exact_systems/agent4/w1_unnormalized_template.m2
A	research/four_parameter_model/agent3/report.md
A	research/hashes/phase0_source_hashes.txt
A	research/hashes/research_provenance.md
A	research/independent_reproduction/reproduction_protocol.md
A	research/intervals/interval_arithmetic_research.md
A	research/literature/agent9/applicability_assessment.md
A	research/literature/agent9/method_bibliography.csv
A	research/literature/agent9/method_bibliography.md
A	research/literature/agent9/raw_search_notes.md
A	research/literature/matrix_catalogue_source_audit.md
A	research/literature/order_six_mub_literature.md
A	research/logs/campaign0/campaign0_environment_audit.md
A	research/logs/phase0_wsl_probe.sh
A	research/logs/research_log.md
A	research/methods/agent2/method_comparison.csv
A	research/methods/agent2/method_comparison.md
A	research/methods/agent2/prototype/NOT_RUN.log
A	research/methods/agent2/prototype/command.txt
A	research/methods/agent2/prototype/mub_h2_w1.jl
A	research/numerical/numerical_strategy_research.md
A	research/optimization/agent6/check_env.py
A	research/parameter_spaces/agent1/logs/literature_and_repository_audit.md
A	research/parameter_spaces/agent1/parameter_space_report.md
A	research/parameter_spaces/agent1/parameterization_catalogue.csv
A	research/referee/referee_risk_register.md
A	research/reports/next_research_report.md
A	research/reports/phase0_environment_report.md
A	research/reports/physics_interpretation.md
A	research/symmetry/agent7/symmetry_gauge_report.md
A	research/symmetry/agent7/transformation_catalogue.csv
A	research_gap_report.md
A	results/arc_circle_asymmetry.txt
A	results/audit_clique_log.txt
A	results/audit_clique_log2.txt
A	results/benchmarks/benchmarks.json
A	results/benchmarks/benchmarks.md
A	results/benchmarks/benchmarks_run.log
A	results/certify_fourth_mub_witness_dita.txt
A	results/certify_fourth_mub_witness_nwit1_dita.txt
A	results/certify_fourth_mub_witness_nwit1_multi_lambda.txt
A	results/certify_fourth_mub_witness_run.log
A	results/certify_nwit1_all_cliques_four_classes.txt
A	results/certify_nwit1_all_cliques_lambdapi.txt
A	results/chm_equivalence.txt
A	results/circulant_match_anomaly.txt
A	results/csv_reconciliation.txt
A	results/degeneracy_candidates.csv
A	results/degeneracy_candidates.json
A	results/dita_lambda_fourth_dense.csv
A	results/dita_lambda_fourth_dense.log
A	results/dita_lambda_fourth_dense_test.log
A	results/drop_detection_stress_test.txt
A	results/enumerate_d0_cliques.txt
A	results/exact_karlsson_dita_H.txt
A	results/export_w1_exact.txt
A	results/extended_axis_sweep.txt
A	results/extended_axis_sweep_full.log
A	results/extended_axis_sweep_item2.log
A	results/extended_axis_sweep_quick.log
A	results/extended_axis_sweep_run.log
A	results/f6_arc_boundary_tight.stderr.log
A	results/f6_arc_boundary_tight.stdout.log
A	results/f6_boundary_reverify.txt
A	results/f6_boundary_reverify_run.log
A	results/f6_boundary_reverify_run2.log
A	results/final_honest_status.md
A	results/formalize_gauge_lemmas.txt
A	results/fourth_mub_locus_sweep.txt
A	results/fourth_mub_locus_sweep_quick.log
A	results/fourth_mub_locus_sweep_run.log
A	results/fourth_mub_locus_sweep_run2.log
A	results/fourth_mub_per_basis.csv
A	results/fourth_test_anchors.log
A	results/gauge_analysis.txt
A	results/grid20_run.log
A	results/identify_d0_in_karlsson.txt
A	results/investigate_Dita.log
A	results/investigate_Dita.txt
A	results/investigate_F6_theta0.log
A	results/investigate_F6_theta0.txt
A	results/investigate_advanced_summary.txt
A	results/investigate_degen_circulant.txt
A	results/karlsson_perH_sweep.csv
A	results/karlsson_perH_sweep_retagged.csv
A	results/karlsson_sweep.csv
A	results/lambda_periodicity_dita.txt
A	results/lambda_periodicity_dita_run.log
A	results/locus_2d_grid.csv
A	results/locus_classification.json
A	results/locus_classification.txt
A	results/locus_classification_run.log
A	results/locus_dimension_confirm.txt
A	results/locus_dimension_confirm_run.log
A	results/locus_geometry_probes.txt
A	results/locus_phi_sweep_full_circle.txt
A	results/locus_phi_sweep_full_circle_run.log
A	results/nid_probe.txt
A	results/nid_probe_run.log
A	results/phase0_source_hashes.txt
A	results/phase0_wsl_probe.sh
A	results/phase1_boundaries.txt
A	results/phase1_definitions_suite.stdout.log
A	results/phase1_final.txt
A	results/phase1_final_status.md
A	results/phase1_full_run.py
A	results/phase1_full_run.txt
A	results/phase1_full_run2.txt
A	results/phase1_isolated.py
A	results/phase1_pkg_instantiate.sh
A	results/phase1_resolution.txt
A	results/phase1_run3.txt
A	results/phase1_run_dir.txt
A	results/phase1_simplify_probe.py
A	results/phase1_simplify_probe.txt
A	results/phase1_timed.txt
A	results/phase1_timed2.txt
A	results/phase1_timing_probe.py
A	results/phase1_timing_probe.txt
A	results/phase2_technical.txt
A	results/pi3_float64_smoke.stderr.log
A	results/pi3_float64_smoke.stdout.log
A	results/pi3_precision_256.txt
A	results/precision_f6_bounded_20260917/anchor.stderr.log
A	results/precision_f6_bounded_20260917/anchor.stdout.log
A	results/precision_f6_bounded_20260917/anchor/metadata.txt
A	results/precision_pi3_run.stderr.log
A	results/precision_pi3_run.stdout.log
A	results/priority_hits_verify.txt
A	results/probe.err.txt
A	results/probe.out.txt
A	results/project_finish_audit.txt
A	results/reconstruct_b3_B3_lambda0.4.csv
A	results/reconstruct_b3_algebraic.meta.txt
A	results/reconstruct_b3_algebraic.txt
A	results/reconstruct_b3_algebraic_bound32.txt
A	results/reconstruct_b3_algebraic_lsround_archive.meta.txt
A	results/reconstruct_b3_algebraic_lsround_archive.txt
A	results/reconstruct_b3_algebraic_sqrt6.txt
A	results/reconstruct_b3_run_deg3.log
A	results/reconstruct_b3_run_deg4.log
A	results/reconstruct_b3_run_deg4_sqrt6.log
A	results/refine_pi3_256.stderr.log
A	results/refine_pi3_256.stdout.log
A	results/retag_512_run.log
A	results/retag_512_summary.txt
A	results/search_full_run.log
A	results/search_full_run_err.log
A	results/special_loci_clean_run.txt
A	results/special_loci_degen1845.csv
A	results/special_loci_degen1845.log
A	results/special_loci_degen1845.meta.txt
A	results/special_loci_degen1845.timing.txt
A	results/special_loci_degen1845_CORRUPTED592.csv
A	results/special_loci_degen1845_regen.csv
A	results/special_loci_degen200.csv
A	results/special_loci_degen200.log
A	results/special_loci_degen500.csv
A	results/special_loci_degen500.log
A	results/special_loci_degen500_run2.csv
A	results/special_loci_degen500_run2.log
A	results/special_loci_degen500_run2.meta.txt
A	results/special_loci_log.txt
A	results/special_loci_resume_run.log
A	results/special_loci_run.log
A	results/special_loci_run2.txt
A	results/special_loci_run2_err.txt
A	results/special_loci_run_log.txt
A	results/special_loci_run_phase1a.log
A	results/special_loci_search.csv
A	results/special_loci_search_891_overshoot.csv
A	results/special_loci_search_896.csv
A	results/special_loci_search_backup848.csv
A	results/special_loci_search_interim.csv
A	results/special_loci_search_partial210.csv
A	results/special_loci_search_v2.csv
A	results/sweep_log.txt
A	results/symbolic_elimination_attempt.txt
A	results/symbolic_elimination_log.txt
A	results/symbolic_elimination_log2.txt
A	results/symbolic_elimination_run.txt
A	results/symbolic_elimination_summary.txt
A	results/task1_investigate.log
A	results/task1_investigate_candidate.log
A	results/task1_mobius_investigation.txt
A	results/task1_run.log
A	results/task2_argmax_audit.log
A	results/task2_argmax_audit.txt
A	results/task3_locus_run.log
A	results/third_mub_audit.csv
A	results/third_mub_candidates_audit.csv
A	results/third_mub_candidates_audit_py.csv
A	results/third_mub_locus_dimension.txt
A	results/third_mub_locus_sampling.txt
A	results/third_mub_locus_sampling_run.log
A	results/top50_degen_stress.txt
A	results/track_c_elimination/arxiv_format_pass/after_dita_claims.txt
A	results/track_c_elimination/arxiv_format_pass/after_fourth_claims.txt
A	results/track_c_elimination/arxiv_format_pass/after_gauge_claims.txt
A	results/track_c_elimination/arxiv_format_pass/after_main_source_claims.txt
A	results/track_c_elimination/arxiv_format_pass/before_dita_claims.txt
A	results/track_c_elimination/arxiv_format_pass/before_fourth_claims.txt
A	results/track_c_elimination/arxiv_format_pass/before_gauge_claims.txt
A	results/track_c_elimination/arxiv_format_pass/before_main_source_claims.txt
A	results/track_c_elimination/m2_fourth_mub_reduced_Dita.log
A	results/track_c_elimination/m2_fourth_mub_reduced_Dita_elim.log
A	results/track_c_elimination/m2_fourth_mub_reduced_Dita_exact_smoke.log
A	results/track_c_elimination/m2_pool_Dita.log
A	results/track_c_elimination/m2_pool_F6.log
A	results/track_c_elimination/pool_Dita_exact.log
A	results/track_c_elimination/w1_D0_cost_estimate.log
A	results/track_c_elimination/w1_D0_field_degree.txt
A	results/track_c_elimination/w1_D0_groebner.log
A	results/validate_finish_run.log
A	results/validate_interim.log
A	results/validate_interim_err.log
A	results/validate_log.txt
A	results/validate_run.log
A	results/validation_summary.txt
A	results/validation_summary_py.txt
A	results/verify_chm_b3_transfer.txt
A	results/verify_closed_form_d0.txt
A	results/verify_dita_construction.txt
A	results/verify_extra_third_hit_run.log
A	results/verify_nwit1_followup.txt
A	results/verify_nwit1_referee_checks.txt
A	results/w1_D0_mathematica.txt
A	results/w1_D0_wolfram_rerun_raw.txt
A	results/witness_decomposition.txt
A	results/witness_decomposition_run.log
A	results/witness_decomposition_run2.log
A	results/witness_mixed_volume_profile.txt
A	runs/2026-09-17T170232Z_phase1_definitions/phase1_definitions_suite.out.txt
A	runs/2026-09-17T170232Z_phase1_definitions/phase1_definitions_suite.py
A	scan.ps1
A	scripts/_run_all.sh
A	scripts/_runner.sh
A	scripts/audit_phase1_anchor.py
A	scripts/docker/run_m2.ps1
A	scripts/julia/_check_env.jl
A	scripts/julia/_install_nemo.jl
A	scripts/julia/_paths.jl
A	scripts/julia/_test_lll_fit.jl
A	scripts/julia/analyze_arc_circle_asymmetry.jl
A	scripts/julia/analyze_locus_grid.jl
A	scripts/julia/audit_argmax_task2.jl
A	scripts/julia/audit_clique_pipeline.jl
A	scripts/julia/certify_fourth_mub_witness.jl
A	scripts/julia/certify_fourth_mub_witness_nwit1.jl
A	scripts/julia/certify_fourth_mub_witness_nwit1_multi_lambda.jl
A	scripts/julia/certify_nwit1_all_cliques_four_classes.jl
A	scripts/julia/certify_nwit1_all_cliques_lambdapi.jl
A	scripts/julia/certify_smoke.jl
A	scripts/julia/check_anchors.jl
A	scripts/julia/compare_dita_coords.jl
A	scripts/julia/csv_reconciliation.jl
A	scripts/julia/diagnose_mobius_grid.jl
A	scripts/julia/diagnose_theta0_limit.jl
A	scripts/julia/dita_lambda_fourth_dense.jl
A	scripts/julia/drop_detection_stress_test.jl
A	scripts/julia/enumerate_d0_cliques.jl
A	scripts/julia/exact_karlsson_dita_H.jl
A	scripts/julia/export_w1_exact.jl
A	scripts/julia/extended_axis_sweep.jl
A	scripts/julia/f6_arc_boundary_tight.jl
A	scripts/julia/f6_boundary_reverify.jl
A	scripts/julia/finalize_degen1845_csv.jl
A	scripts/julia/formalize_gauge_lemmas.jl
A	scripts/julia/fourth_mub_locus_sweep.jl
A	scripts/julia/gauge_equivalence_analysis.jl
A	scripts/julia/identify_d0_in_karlsson.jl
A	scripts/julia/inspect_witness_coeffs.jl
A	scripts/julia/install_oscar_probe.jl
A	scripts/julia/investigate_candidate.jl
A	scripts/julia/investigate_circulant_match_anomaly.jl
A	scripts/julia/lambda_periodicity_dita.jl
A	scripts/julia/locus_classification.jl
A	scripts/julia/locus_dimension_confirm.jl
A	scripts/julia/locus_geometry_probes.jl
A	scripts/julia/locus_phi_sweep_full_circle.jl
A	scripts/julia/nid_probe.jl
A	scripts/julia/pi3_float64_smoke.jl
A	scripts/julia/precision_f6_bounded.jl
A	scripts/julia/precision_f6_launch.sh
A	scripts/julia/profile_witness_mixed_volume.jl
A	scripts/julia/reconstruct_b3_algebraic.jl
A	scripts/julia/refine_pi3_256.jl
A	scripts/julia/retag_512_completeness.jl
A	scripts/julia/run_benchmarks.jl
A	scripts/julia/run_fourth_mub_hpc.jl
A	scripts/julia/run_fourth_test_anchors.jl
A	scripts/julia/run_task1_investigate.jl
A	scripts/julia/run_task1_mobius.jl
A	scripts/julia/run_top50_degen.jl
A	scripts/julia/search_special_loci.jl
A	scripts/julia/smoke_test.jl
A	scripts/julia/sweep_karlsson.jl
A	scripts/julia/symbolic_elimination.jl
A	scripts/julia/symbolic_probe.jl
A	scripts/julia/third_mub_locus_dimension.jl
A	scripts/julia/third_mub_locus_sampling.jl
A	scripts/julia/validate_perH_pipeline.jl
A	scripts/julia/validate_third_mub_candidates.jl
A	scripts/julia/verify_chm_b3_transfer.jl
A	scripts/julia/verify_closed_form_d0.jl
A	scripts/julia/verify_dita_construction.jl
A	scripts/julia/verify_extra_third_hit.jl
A	scripts/julia/verify_field_degree_w1.jl
A	scripts/julia/verify_nwit1_followup.jl
A	scripts/julia/verify_nwit1_referee_checks.jl
A	scripts/julia/verify_priority_hits.jl
A	scripts/julia/verify_quick_run_task1.jl
A	scripts/julia/witness_decomposition.jl
A	scripts/python/__pycache__/degeneracy_scan.cpython-313.pyc
A	scripts/python/__pycache__/hadamard6.cpython-313.pyc
A	scripts/python/__pycache__/karlsson_k6_3.cpython-313.pyc
A	scripts/python/audit_karlsson_variants.py
A	scripts/python/audit_latex_warnings.py
A	scripts/python/bundle_main_theorems.py
A	scripts/python/bundle_methods_audit.py
A	scripts/python/check_paper_latex.py
A	scripts/python/chm_equivalence.py
A	scripts/python/degeneracy_scan.py
A	scripts/python/finish_project_audit.py
A	scripts/python/hadamard6.py
A	scripts/python/identify_d0_in_karlsson.py
A	scripts/python/karlsson_k6_3.py
A	scripts/python/validate_third_mub_table.py
A	scripts/wolfram/verify_dita_exact.wl
A	scripts/wolfram/verify_karlsson_audit.wl
A	scripts/wolfram/verify_repository_mathematics.wl
A	scripts/wolfram/w1_D0_exact.wl
A	setup_julia_env.ps1
A	source_quotes.md
A	src/Benchmarks.jl
A	src/Certification.jl
A	src/Cliques.jl
A	src/FourthMUBHPC.jl
A	src/Karlsson.jl
A	src/LiangChenLongQiu_CHM_Transcription.jl
A	src/MubSearch.jl
A	src/Pool.jl
A	src/Provenance.jl
A	src/StaticMUBKernels.jl
A	src/brierley_weigert_notes.jl
A	src/dita_third_mub_construction.jl
A	src/hpc/M1Hardware.jl
A	src/hpc/M2Mubness.jl
A	src/hpc/M3Search.jl
A	src/karlsson_gauge_only.jl
A	src/mub_zauner_6d_liang_chen.jl
A	symbolic_export/fourth_mub_reduced_Dita.m2
A	symbolic_export/fourth_mub_reduced_Dita_elim.m2
A	symbolic_export/fourth_mub_reduced_Dita_exact_smoke.m2
A	symbolic_export/fourth_mub_reduced_F6_theta0.m2
A	symbolic_export/fourth_mub_w_elimination_Dita.m2
A	symbolic_export/fourth_mub_w_elimination_F6.m2
A	symbolic_export/fourth_mub_w_elimination_theta0.m2
A	symbolic_export/mobius_degeneracy_locus.m2
A	symbolic_export/pool_Dita.m2
A	symbolic_export/pool_Dita_exact.m2
A	symbolic_export/pool_F6.m2
A	symbolic_export/pool_F6_theta0.m2
A	symbolic_export/theta0_fourth_mub_witness.m2
A	symbolic_export/w1_D0_dim_degree.m2
A	symbolic_export/w1_D0_exact.m2
A	symbolic_export/w1_D0_groebner.m2
A	test/norm_detection.jl
A	test/runtests.jl
A	tindall_arxiv.xml
A	tindall_arxiv_records.csv
A	tindall_crossref.json
A	tindall_crossref_records.csv
A	tindall_openalex.json
A	tindall_openalex_records.csv
A	tindall_orcid.json
A	tindall_orcid_records.csv
A	transcribe_C2CHM_source.html
```
