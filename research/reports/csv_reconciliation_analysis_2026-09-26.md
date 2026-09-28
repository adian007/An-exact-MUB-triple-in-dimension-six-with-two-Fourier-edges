# CSV reconciliation analysis — is `results/csv_reconciliation.txt` a real partial reconciliation?

**Agent:** csv-scout · **Date:** 2026-09-26 · **Branch:** `main` · **Commit under audit:** `b1c81e1`
**Scope:** `results/csv_reconciliation.txt`, `scripts/julia/csv_reconciliation.jl`,
`scripts/python/validate_third_mub_table.py`, `scripts/julia/validate_third_mub_candidates.jl`.
**Paper files touched:** none.

---

## 0. Verdict up front

| Question | Verdict |
|---|---|
| Is `588` a hard-coded **design target** (`S_588`), not a measurement? | **YES — confirmed.** |
| Is `588 vs 928` therefore a *pipeline discrepancy*? | **NO.** It is target-vs-achieved **overshoot**. |
| Is `csv_reconciliation.txt` a partial reconciliation of the **848 vs 928** (15 vs 45) issue? | **NO — it does not address that issue at all.** It is a **RED HERRING** w.r.t. it. |
| Is the `588 vs 928` relationship itself *explained*? | **NO — the "GAP EXPLANATION" section is sign-incoherent and self-refuting.** |
| Is `848` explained? | **YES — a mundane cause: a different, older input CSV + a different `ref_ref` filter.** |

**Bottom line for the docs agent:** the open discrepancy in `docs/scientific_status.md` §6 is
**substantially de-risked but NOT yet closed**, because the one-line fix that reconciles
848/15 with 928/45 is a *code* change (`CSV_PATH` in the Python validator) that has **not been
made and not been re-run**. Until that patch is made and re-run, the discrepancy is *explained*,
not *resolved*.

---

## 1. What does `scripts/julia/csv_reconciliation.jl` actually compute?

Read in full (254 lines). It is a **text report generator**, not a search. It reads three files:

| Constant | File | Role |
|---|---|---|
| `DEGEN_CSV` | `results/degeneracy_candidates.csv` | candidate pool size (1845 rows) |
| `SEARCH_CSV` | `results/special_loci_search.csv` | the artifact under analysis (928 rows) |
| `OUT` | `results/csv_reconciliation.txt` | **hard-coded output — the script overwrites its own artifact** |

Evidence tags for this section: **`exact symbolic`** (source reading; the literals below are
visible in the file text and are not derived from any run).

### 1.1 Is `588` a hard-coded target? — CONFIRMED, three ways

```julia
# line 137
push!(lines, @sprintf("  expected (verify_csv): 588"))
# line 138
push!(lines, @sprintf("  gap vs 588: %d (negative = more rows than target)", 588 - n_csv))
# line 139
push!(lines, @sprintf("  historical baseline (phase2): 471 rows, gap=117"))
```

* `588` is a **bare integer literal** in a `%d` format string. It is not computed, not read
  from a config file, not derived from `degeneracy_candidates.csv`, and not passed as an
  argument. `parse_degen_cap()` (lines 9–18) parameterises only `degen_cap` (default 80) —
  **there is no `--expected-rows` flag and no code path that can change 588.**
* `471` and `gap=117` on line 139 are **likewise hard-coded literals**, not measured. They
  are a transcription of `results/phase2_technical.txt` §2.1.
* The script itself labels 588 as a *target*: line 138 prints the parenthetical
  **"(negative = more rows than target)"**. The script's own author therefore intended
  `588 = target`, `n_csv = achieved`.
* Independent corroboration that 588 is a design target, not a measurement:
  * `docs/overview/reproduce.md:104` — "Completing search set `\\mathcal{S}_{588}` to 588 rows
    requires an overnight run"; line 113 — "Design target `\\mathcal{S}_{588}` is a subset".
  * `docs/overview/readme.md:53` — "928 … Design target `\\mathcal{S}_{588}`; extended-refinement superset".
  * `docs/results/final_honest_status.md:43` — "**928 rows** … design target `S_588`".
  * `results/phase2_technical.txt:8` — `search_special_loci.jl --verify-only, expected_rows=588`.
    So 588 originates as an `expected_rows` assertion in the *search* script, i.e. a target.
  * `results/project_finish_audit.txt:11` — "--- S588 primary clique-6 (exclude ref_ref) ---".

**Definitive answer to Q1: `588` is a hard-coded design target `S_588`, asserted by
`search_special_loci.jl --verify-only` as `expected_rows=588` and re-asserted as a string
literal in `csv_reconciliation.jl`. It is not a measurement. `928` is a measurement of
`results/special_loci_search.csv`. Therefore `588 vs 928` is a target-vs-achieved OVERSHOOT of
+340 rows, and is categorically NOT the Python-vs-Julia pipeline discrepancy.**

### 1.2 Broken cross-reference

Line 3 of the artifact points to `docs/SEARCH_SETS.md`. That file **does not exist**
(verified: `Test-Path docs/SEARCH_SETS.md` -> `False`; no `docs/SEARCH*` file exists at all).
Evidence tag: **`sampled/numerical only`** (filesystem check). `docs/overview/readme.md:89`
lists `SEARCH_SETS` in a directory sketch, so the document appears to have been lost or
never committed. The S588 definition is therefore only recoverable from
`reproduce.md`, `readme.md` and `phase2_technical.txt` — not from the file the artifact cites.


---

## 2. Does `csv_reconciliation.txt` address the 848 vs 928 discrepancy? — **NO. CONFIRMED ORTHOGONAL.**

Evidence tag: **`exact symbolic`** (full-text reading of both the 45-line artifact and the
254-line script).

The artifact contains the string `588` four times and the strings `848` / `15` **zero times**.
The script contains `588` three times (lines 137, 138, 201) and `848`/`backup848` **zero times**.
Neither file ever opens `results/special_loci_search_backup848.csv`, never opens
`scripts/python/validate_third_mub_table.py`, and never references `validation_summary_py.txt`.
The two subjects share only the number 928, which is a coincidence of both being analyses of
`special_loci_search.csv`.

**Therefore: `csv_reconciliation.txt` explains only the 588/928 target question. It is entirely
orthogonal to the 848-vs-928 Python/Julia discrepancy, and is a RED HERRING with respect to it.**
The audit's phrasing that it "partially explains 588→928" is literally true but describes a
*different* discrepancy from the one §6 is about.

---

## 3. Tracing BOTH pipelines' row-selection logic

### 3.1 Julia pipeline — `scripts/julia/validate_third_mub_candidates.jl`

Evidence tag: **`exact symbolic`** (source) + **`sampled/numerical only`** (recount).

```julia
const DEFAULT_CSV = joinpath(RESULTS_DIR, "special_loci_search.csv")   # line 15
cand = [r for r in rows if r.status == "ok" && r.max_clique >= 6]      # line 43
# exclude_deep_ref = false  and  dedupe_ref = false  by default  (lines 247-248)
println(io, "n_csv_rows=$(length(rows))")                              # line 233
println(io, "n_candidates=$(length(cand))")                            # line 234
```

* Input: **`results/special_loci_search.csv`** (mtime 2026-08-03 16:32:11).
* Row selection: every row, no dedup, no `ref_ref` exclusion, no `pool_complete_flag` filter.
* Row acceptance: `status == "ok" && max_clique >= 6`.
* Malformed-row handling (line 26): `length(parts) < length(header) && continue` — **lenient**
  (keeps short rows, unlike the Python `DictReader`).
* Output: `results/validation_summary.txt` (928 / 45), `results/third_mub_audit.csv`,
  `results/fourth_mub_per_basis.csv`.

### 3.2 Python pipeline — `scripts/python/validate_third_mub_table.py`

Evidence tag: **`exact symbolic`** (source) + **`sampled/numerical only`** (recount).

```python
CSV_PATH = RESULTS / "special_loci_search_backup848.csv"               # line 10  <-- THE BUG
def third_candidates(rows, exclude_deep_ref=True):                    # line 20
    if row.get("status") != "ok":            continue                 # line 23
    if mc < 6:                               continue                 # line 29
    if exclude_deep_ref and "ref_ref" in row.get("name",""): continue # lines 31-32
```

* Input: **`results/special_loci_search_backup848.csv`** (mtime 2026-08-03 **11:49:54**).
* Row selection: same `status=="ok" && max_clique>=6` core, **but** `exclude_deep_ref=True`
  is the **default and is hard-coded** (no CLI flag, no env var, no argparse).
* Malformed-row handling: `csv.DictReader` — will raise/skip differently from the Julia
  lenient parser.
* Output: `results/validation_summary_py.txt` (848 / 15),
  `results/third_mub_candidates_audit_py.csv` (10 columns vs the Julia table's 12).

### 3.3 The decisive experiment: cross-applying both predicates to both files

Recount performed (read-only, deterministic; see
`results/NOT_RUN_csvscout_julia_reconciliation_rerun.log`, Entry 3). Evidence tag:
**`sampled/numerical only`** (exact integer counting over stored artifacts — no
interval-arithmetic certification is claimed, and none is needed, since these are counts).

| Input CSV | Julia predicate<br>(ok ∧ mc≥6, ref_ref **kept**) | Python predicate<br>(ok ∧ mc≥6 ∧ **no** ref_ref) |
|---|---:|---:|
| `special_loci_search.csv` (928 rows) | **45** | **20** |
| `special_loci_search_backup848.csv` (848 rows) | 64 | **15** |

* **45** on the 928-row file reproduces `results/validation_summary.txt` `n_candidates=45` exactly. ✔
* **15** on the 848-row file reproduces `results/validation_summary_py.txt` `n_candidates=15` exactly. ✔
* `results/third_mub_audit.csv` has **45** data rows; `results/third_mub_candidates_audit_py.csv`
  has **15** data rows. ✔ (independent corroboration of the candidate tables, not just the summaries)

**This is a complete, closed decomposition. The discrepancy is fully attributable:**

```
15  = Python predicate  on backup848      (stored Python result)
20  = Python predicate  on search.csv     (+5 primary-level candidates from the later run)
45  = Julia predicate   on search.csv     (stored Julia result; +25 ref_ref candidates)
64  = Julia predicate   on backup848      (the 848 file has MORE clique-6 rows, not fewer!)
```


### 3.4 Complete enumeration of every divergence point

| # | Divergence | Julia | Python | Impact |
|---|---|---|---|---|
| **D1** | **Input CSV path** | `special_loci_search.csv` (hard-coded const, line 15) | `special_loci_search_backup848.csv` (hard-coded, line 10) | **Dominant.** +80 rows. *This is the whole story of the 848.* |
| **D2** | `ref_ref` (2nd-level refinement) rows | **INCLUDED** (`exclude_deep_ref=false` default) | **EXCLUDED** (hard-coded default `True`) | On `search.csv`: 45 vs 20 → **+25 candidates** |
| **D3** | Deduplication | none | none | *No divergence today* — but Julia exposes `--dedupe-ref` (0.005-cell rep, line 51) and Python has no equivalent. A future run *with* the flag would silently change 45. Latent trap. |
| **D4** | Malformed-row policy | lenient: skips only rows **shorter** than header, keeps longer ones | strict `DictReader` | Would diverge on ragged rows; **no divergence today** (`n_bad` line 217 was never emitted) |
| **D5** | Column-name→index binding | `Dict` from header, order-independent | `DictReader` | No divergence (both 22 identical columns, verified) |
| **D6** | `max_clique` parse | `parse(Int, ...)`; empty → `0`; **throws** on non-integer | `int(... or 0)`, `ValueError` → **silently `continue`** | Diverges only on a corrupt field: Julia crashes, Python drops the row. Latent. |
| **D7** | `found_fourth` | recomputes from scratch (`fourth_mub_per_point`, HP) | reads the stored `found_fourth` column | Different *meaning*, same *value* today (both `False`). Not a count divergence. |
| **D8** | HP verification | runs `hp_verify_third` (256-bit) on every candidate | **no HP step at all** | Julia reports `hp_survived=36`; Python has no such field. This is a *scientific* asymmetry, not a row-count one. |
| **D9** | Output table | `third_mub_audit.csv`, 12 cols, `hp_survival` | `third_mub_candidates_audit_py.csv`, 10 cols | This is the audit's "reduced/different audit table" remark — accurate but **not the cause** of 15 vs 45. |
| **D10** | Clustering | `svd` on centred matrix, `MersenneTwister(20260803)` | `np.linalg.svd` on centred matrix | `cluster_rank=3` in both — coincidentally equal despite different inputs. Cosmetic. |

**Which pipelines could NOT diverge (verified equal):** the `status=="ok"` criterion, the
`max_clique >= 6` threshold, the absence of any dedup, and the 22-column schema. So the audit's
suggestion that "the Python script writes a reduced/different audit table [and] cannot be used to
reproduce the Julia count" is **misleading as a causal explanation** — the reduced table is a
*symptom*; the *cause* is D1 (different input file) and D2 (different `ref_ref` filter).

### 3.5 An unexpected and important finding

`special_loci_search.csv` has **more rows (928 > 848) but FEWER clique-6 candidates
under Julia's own predicate (45 < 64)**. So `search.csv` is **not** a superset of `backup848.csv`
in any scientifically meaningful sense — later is *not* more here. Composition explains it:

| | `search.csv` (928) | `backup848.csv` (848) |
|---|---:|---:|
| anchor / non-`ref_` names | 21 | 20 |
| `degen_*` rows | **108** | **0** |
| `ref_*` level-1 | 372 | 248 |
| `ref_ref_*` level-2+ | **428** | **580** |
| `ok` rows | 899 | 794 |
| unique coord keys | **925** (3 dupes) | 848 (0 dupes) |
| `pool_complete_flag=true` | 899 | — |

`backup848.csv` is **degenerate-candidate-free** (zero `degen_*` rows) and is 68% second-level
refinements. The later `search.csv` introduced 108 `degen_*` primaries and *fewer* `ref_ref_`
rows. Whatever produced them, they are **structurally different runs**, not the same run at two
lengths. This is a second, independent reason not to treat the two as one dataset.

---

## 4. Is the 848/45 issue better explained by a different input file? — **YES, decisively.**

Evidence tag: **`sampled/numerical only`** (filesystem mtimes + deterministic recount).

| Artifact | mtime (local) | bytes | rows |
|---|---|---:|---:|
| `special_loci_search_backup848.csv` | **2026-08-03 11:49:54** | 112,654 | 848 |
| `special_loci_search_v2.csv` | 2026-08-03 12:23:17 | 36,886 | — |
| `special_loci_search_interim.csv` | 2026-08-03 12:23:42 | 47,034 | — |
| `special_loci_search_partial210.csv` | 2026-08-03 12:38:46 | 27,448 | — |
| `special_loci_search_896.csv` | 2026-08-03 14:17:18 | 128,430 | 896 |
| `special_loci_search_891_overshoot.csv` | 2026-08-03 16:02:46 | 186,690 | 891 |
| **`special_loci_search.csv`** | **2026-08-03 16:32:11** | 192,242 | **928** |

`backup848.csv` is the **OLDEST artifact in the family — 4 h 42 min older** than
`special_loci_search.csv`, and it is also the earliest file in the entire `results/special_loci*.csv`
set. Its very name ("backup848") advertises that it is a snapshot.

**The mundane explanation is confirmed:** the 848 is *not* a different selection algorithm applied
to the same data. It is **a stale, superseded input file** hard-coded into the Python validator,
combined with a **hard-coded `ref_ref` exclusion** that the Julia validator does not apply. The
Python validator has not been pointed at the authoritative CSV.

**Caveat (must not be overstated):** this explains the *arithmetic* completely, but it does **not**
establish that `search.csv` is scientifically *correct* — only that it is the designated
authoritative artifact (`docs/scientific_status.md:128`, `reproduce.md:113`, `readme.md:53`).
`results/phase2_technical.txt` §2.2 already records that at least 4 `ref_degen_*` rows carry
**stale `max_clique=6` labels from a pre-argmax-fix run**, and one of them,
`ref_ref_degen_circulant_match`, is present in `search.csv` with `max_clique=2`. So the 45
candidates themselves are not automatically sound.


---

## 5. The sign-coherence problem — **CONFIRMED, and worse than suspected**

Evidence tag: **`exact symbolic`** (arithmetic on the literals in the source) +
**`sampled/numerical only`** (the counts they refer to).

The artifact's section header is
`=== GAP EXPLANATION (588 - 928 = -340) ===` and the reported gap is **negative**, i.e. the CSV
contains **340 MORE rows than the target**. Every one of the four listed reasons, however, is a
statement about things that were **NOT** produced. Each one, if true, pushes the count **DOWN**,
i.e. in the direction of a **positive** gap. **Not one of the four can produce a negative gap.**

| # | Listed reason | Sign of its effect on row count | Can it explain gap = −340? |
|---|---|---|---|
| 1 | "degen cap: only 80/1845 unique degen candidates loaded (−1765 never queued)" | **negative** (fewer rows) | **No** |
| 2 | "refinement shortfall: ~1701 ref points not written (crash/incomplete run/dedup skips)" | **negative** (fewer rows) | **No** — and see contradiction C2 |
| 3 | "second-level `ref_ref` blocked by `can_spawn_refinement`" | **negative** (fewer rows) | **No** — and see contradiction C2 |
| 4 | "coordinate dedup: refinements overlapping anchors/degen skipped" | **negative** (fewer rows) | **No** — and see contradiction C3 |

### Verified arithmetic behind each claim

* **Reason 1**: `1845 − 80 = 1765`. ✔ arithmetic correct. But the premise is already dead:
  the script's own model says `primary total = 20 anchors + 80 degen = 100`, whereas the file
  actually contains **21 + 108 = 129** primary rows. **The CSV has 28 MORE `degen_*` rows than the
  cap of 80 permits.** The `degen_cap = 80` model simply does not describe the run that produced
  `search.csv`. (`degeneracy_candidates.csv` does contain exactly 1845 data rows — verified — so
  the *pool* count is right; the *cap* is what fails to apply.)
* **Reason 2**: `est_ref_typ = 20 × 5^3 = 2500`, then `2500 − 799 = 1701`. ✔ arithmetic correct.
  (20 = the 45 spawn-eligible rows minus the 25 `ref_ref` ones, per line 186; all 20 have
  `max_clique ≥ 6` so all contribute `125`.) But this compares the CSV's `ref_*` count against a
  **theoretical upper bound with no dedup** — the script labels it as much
  (`"theoretical max, no dedup"`, line 195). A shortfall against a *maximum* is not evidence of a
  crash. Note also the estimator uses the *Julia* spawn predicate (which **includes** `ref_ref`)
  and then *excludes* `ref_ref` when summing — internally inconsistent bookkeeping.
* **Reason 3**: "second-level `ref_ref` blocked." **Directly contradicted by the artifact's own
  data**: it reports `ref_ref_* (level 2+): 428` in the very same file. `ref_ref` is not blocked;
  it is the single largest category after `ref_*` level-1. And `can_spawn_refinement` **does not
  exist anywhere in the repository** — a whole-codebase grep for it returns exactly one hit: the
  text of `results/csv_reconciliation.txt` itself. It is a **phantom function name**, invented by
  this script and never cross-checked against `search_special_loci.jl`.
* **Reason 4**: "coordinate dedup skipped refinements." The artifact's own line
  `unique coordinate keys: 925 (duplicates: 3)` means dedup removed **at most 3 rows** (verified
  independently: 925 unique keys from 928 rows). It cannot account for hundreds of rows.

### Additional internal inconsistencies found

* **C1 — the row-category block does not sum to the row count.** The artifact lists
  `21 + 108 + 372 + 428 = 929`, but reports `928` data rows. Verified. The cause is a predicate
  mismatch in the source: `n_anchor`/`n_degen`/`n_ref1` use `startswith`, while `n_ref2` uses
  `occursin("ref_ref", …)` (line 147, *substring anywhere*). The clean prefix partition
  `21 + 108 + 799 = 928` ✔ is exact, so exactly **one** row is double-counted — the row whose
  name contains `ref_ref` without carrying the `ref_` prefix (rendered as
  `ref_ref_degen_circulant_match`). *Mechanism note: the raw field appears to carry leading
  padding, so the prefix test fails while the substring test passes; I could not fully confirm
  the padding with the tools available. Immaterial to every other conclusion here.*
  Tag for this micro-detail: **`sampled/numerical only`**.
* **C2 — reasons 3 and 4 contradict the file's own numbers**, as shown above.
* **C3 — the "degen rows missing" and "primary shortfall" guards never fired.** Lines 207–216 emit
  those reasons only when `exp.n_degen_loaded − n_degen > 0` (here `80 − 108 = −28`) and
  `n_primary_csv < exp.n_primary` (here `129 < 100` is false). Their absence is the script's own
  evidence that the 100-primary model does not fit — yet the file still prints the model's
  "primary total: 100" as if it were the plan of record.
* **C4 — the two conditional reasons that DID fire were computed against a different baseline than
  the reported gap.** `missing_degen`, `missing_ref` and `n_primary_csv` are all measured against
  `exp` (the S588/80-degen model), while the header gap is `588 − n_csv`. Mixing two baselines
  under one heading is the structural cause of the incoherence.

**Conclusion on the 588-vs-928 relationship: the relationship itself is REAL and correctly signed
(−340 = a genuine 340-row overshoot of the S588 target), but its "GAP EXPLANATION" is
INCOHERENT — it is a sign-flipped list of under-production diagnostics pasted under an
over-production heading, and two of its four items are falsified by the very file it analyses.**


---

## 6. Answers to the five questions, consolidated

**Q1 — Is 588 a hard-coded design target?** YES. Bare literals at `csv_reconciliation.jl:137-139`
and `:201`; no flag, no config, no computation. It originates as `expected_rows=588` in
`search_special_loci.jl --verify-only` (`phase2_technical.txt:8`) and is documented as
`\mathcal{S}_{588}` in `reproduce.md:104,113`, `readme.md:53`, `final_honest_status.md:43`.
**So 588-vs-928 is target-vs-achieved overshoot and is NOT the Python/Julia discrepancy.**
`exact symbolic`.

**Q2 — Does the file address 848 vs 928?** NO. It never mentions 848, 15, `backup848`, or the
Python validator. **Confirmed orthogonal; it is a red herring w.r.t. the §6 discrepancy.**
`exact symbolic`.

**Q3 — Do the pipelines agree on row selection?** NO, on two counts: **(D1)** different input
CSV, and **(D2)** different `ref_ref` treatment. Their `status=="ok"` and `max_clique>=6` criteria
DO agree, and neither dedups. Ten divergence points are enumerated in §3.4; only D1 and D2 are
active today.

**Q4 — Is 848 from a different/older input file rather than different logic?** YES, and *both*:
it is an older file (4 h 42 min) *and* a different filter. The mundane explanation is the
correct one, and it is now **quantitatively complete** (15 → 20 → 45 → 64, §3.3).

**Q5 — Verdict.** `csv_reconciliation.txt` is a **red herring** with respect to the
848/15-vs-928/45 discrepancy, and a **sign-incoherent, partly self-refuting** artifact with
respect to the 588/928 question it does address. It should **not** be cited as partial progress
on §6. Its correct use is narrow: it is a 2026-08-05 *target-vs-achieved* report proving the S588
queue was overshot by 340 rows, and it also usefully documents the CSV's composition
(925 unique keys, 899 `ok`, 3 duplicates).

---

## 7. Residual open items (do NOT overclaim closure)

1. **The Python validator still points at the stale file.** `scripts/python/validate_third_mub_table.py:10`
   must be changed to `RESULTS / "special_loci_search.csv"` (or, better, take a `--csv` argument
   mirroring the Julia script). **Until this is patched and re-run, §6 is *explained*, not *resolved*.**
   *(Not done here: modifying the script is outside this agent's read-only brief, and the
   re-run would overwrite two `results/` artifacts, which is forbidden without backup+campaign gate.)*
2. **The `ref_ref` semantics must be fixed in writing, not just in code.** Which is correct —
   counting or excluding 2nd-level refinements — is a **scientific** decision, not a
   reproducibility one. On the current authoritative file the two readings give **45 vs 20**.
   The paper and `results/project_finish_audit.txt:11` ("S588 primary clique-6 (exclude ref_ref) —
   count=45") are **internally inconsistent on this exact point** and must be reconciled by a human.
3. **`ref_ref` staleness.** `phase2_technical.txt` §2.2 records ≥4 `ref_degen_*` rows retaining
   pre-fix `max_clique=6` labels. The 45 should be re-derived from a post-fix run.
4. **`docs/SEARCH_SETS.md` is missing** though cited by the artifact and sketched in `readme.md:89`.
5. **Julia was unavailable**, so no Julia code path was executed. All Julia claims are
   **source-reading (`exact symbolic`)** or **recount (`sampled/numerical only`)**. No claim here
   is tagged `certified numerical (interval arithmetic)`.

---

## 8. RECOMMENDED WORDING for `docs/scientific_status.md` §6

> Ready-to-paste. The docs agent should replace the current §6 sentence *"First reconcile the
> Python audit's 848-row/15-candidate view with the Julia audit's 928-row/45-candidate view.
> Until that discrepancy is explained, aggregate counts should not be used as a single
> scientific dataset."* with the following. I have **not** edited that file.

```markdown
**Reconciliation status (2026-09-26).** The 848/15-vs-928/45 difference is **explained but
not yet closed**, and it is *not* a disagreement about the mathematics. It decomposes exactly
into two independent, non-mathematical causes:

1. **Different input file.** `scripts/python/validate_third_mub_table.py:10` hard-codes
   `results/special_loci_search_backup848.csv` (2026-08-03 11:49:54, 848 rows), whereas
   `scripts/julia/validate_third_mub_candidates.jl:15` defaults to
   `results/special_loci_search.csv` (2026-08-03 16:32:11, 928 rows). The Python validator is
   pointed at a superseded snapshot roughly 4 h 42 min older than the authoritative artifact.
2. **Different `ref_ref` handling.** The Julia predicate is
   `status == "ok" && max_clique >= 6` and **keeps** second-level `ref_ref_*` rows. The Python
   predicate is the same core test but **hard-codes** `exclude_deep_ref = True`. Neither
   pipeline deduplicates, and both use an identical 22-column schema and identical
   `status`/`max_clique` thresholds.

Cross-applying both predicates to both files reproduces both stored numbers exactly:

| input CSV | Julia predicate (keeps `ref_ref`) | Python predicate (drops `ref_ref`) |
|---|---:|---:|
| `special_loci_search.csv` (928) | **45** = stored `n_candidates=45` | 20 |
| `special_loci_search_backup848.csv` (848) | 64 | **15** = stored `n_candidates=15` |

So `15 -> 20 -> 45` is fully accounted for, and no computational disagreement remains between
the two scripts. **The two CSV snapshots are nevertheless not interchangeable:** the 848-row
file contains **zero** `degen_*` rows and 580 `ref_ref_*` rows, whereas the 928-row file has
108 `degen_*` and 428 `ref_ref_*` rows; and the 848-row file yields **64** clique-6 candidates
under Julia's own predicate versus 45 for the authoritative file. The later artifact is not a
superset of the earlier one.

**Two items remain open and are *not* resolved by this decomposition:**

* `scripts/python/validate_third_mub_table.py` still points at the stale file. Until it is
  repointed at `special_loci_search.csv` and re-run (writing to new output files), the §6
  discrepancy is **explained, not resolved**, and the stale
  `results/validation_summary_py.txt` (848/15) and `results/third_mub_candidates_audit_py.csv`
  (15 rows) must not be cited.
* The **scientific** choice of whether second-level `ref_ref_*` refinements count as
  third-MUB candidate loci is undecided: on the authoritative file the two readings give
  **45 (inclusive)** versus **20 (exclusive)**. `results/project_finish_audit.txt:11` labels the
  inclusive count as "S588 primary clique-6 (**exclude** ref_ref) — count=45", which is
  self-contradictory and must be resolved by hand.
```

```markdown
**A separate, unrelated discrepancy.** `results/csv_reconciliation.txt` does **not** bear on the
848/928 question. It compares the *hard-coded design target* `S_588 = 588`
(`scripts/julia/csv_reconciliation.jl:137-139` and `:201`, bare integer literals, not a
measurement; originally an `expected_rows=588` assertion in `search_special_loci.jl --verify-only`)
against the *measured* 928 rows. That relationship is a **+340-row target overshoot**, not a
pipeline mismatch. Its `GAP EXPLANATION` section is additionally **sign-incoherent**: it lists
only under-production causes ("never queued", "not written", "skipped") under a header whose
gap is negative, and two of its four items are falsified by the file itself — the CSV contains
428 `ref_ref_*` rows despite the claim that second-level refinement is "blocked", and coordinate
dedup removed at most 3 rows (925 unique keys from 928 rows). The named helper
`can_spawn_refinement` does not exist anywhere in the repository outside that text. The artifact
should be cited only as evidence that the S588 queue was overshot, never as a reconciliation of
the Python/Julia count difference. Its citation of `docs/SEARCH_SETS.md` is also a dead link —
that file is absent.

Until the validator is repointed and re-run, **aggregate counts should still not be used as a
single scientific dataset**; however the reason is now a stale input path and an undecided
`ref_ref` convention, not an unexplained numerical disagreement.
```

---

## 9. Evidence tag ledger

| Claim | Tag |
|---|---|
| `588`, `471`, `gap=117` are hard-coded literals; no code path changes them | `exact symbolic` |
| `CSV_PATH` = `backup848.csv`; Julia `DEFAULT_CSV` = `special_loci_search.csv` | `exact symbolic` |
| Python excludes `ref_ref`; Julia includes it; neither dedups | `exact symbolic` |
| 15 / 20 / 45 / 64 cross-predicates; 928 vs 848 rows; 899 vs 794 `ok`; 925 vs 848 keys | `sampled/numerical only` (deterministic integer recount of stored CSVs) |
| `third_mub_audit.csv` = 45 rows; `third_mub_candidates_audit_py.csv` = 15 rows | `sampled/numerical only` |
| `can_spawn_refinement` absent from the codebase (1 hit, in the artifact's own text) | `exact symbolic` (whole-repo grep) |
| `docs/SEARCH_SETS.md` does not exist | `sampled/numerical only` (filesystem) |
| File mtimes / byte sizes | `sampled/numerical only` (filesystem) |
| Row-category block sums to 929, not 928; offending row `ref_ref_degen_circulant_match` | `sampled/numerical only` |
| `1845−80=1765`, `20×125=2500`, `2500−799=1701` | `exact symbolic` (integer arithmetic on verified counts) |
| All Julia runtime behaviour; `estimate=1701` internal bookkeeping | `not run` (no Julia on host) |
| `phase2_technical.txt` §2.2 stale `ref_degen_*` clique labels | `sampled/numerical only` (read the file; did **not** re-verify the 4 rows) |
| `cluster_sv` differences; `hp_survived=36` | `not run` (Julia unavailable) |

**No claim in this report is tagged `certified numerical (interval arithmetic)`.** Nothing here
required interval certification — the disputed quantities are row counts, not real numbers.

**External citations:** none. No arXiv/DOI sources are relied upon; the entire analysis is
internal to the repository, so no external search was performed or needed.






