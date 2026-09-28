# Independent verification of four headline numbers

**Verifier:** `verifier` (independent re-derivation from stored artifacts + source)
**Date:** 2026-09-26
**Branch/commit:** `main` @ `b1c81e1` (working tree clean except untracked `IDEA.md`)
**Method:** static re-derivation only. **No numerical experiment was executed.** No file in
`results/` was modified or overwritten; one new file `results/NOT_RUN_verifier_round1.log`
was added to record what was deliberately *not* run.

**Evidence tags used below (exactly one per claim):** `exact symbolic`,
`certified numerical (interval arithmetic)`, `sampled/numerical only`, `not run`.

---

## Tooling warning that invalidates negative greps (read this first)

`search_codebase` in this environment **does not index `.jl` files** (neither `src/*.jl` nor
`scripts/julia/*.jl`); it indexes `.py`, `.md`, `.txt`, `.csv`, `.json`. Proof: a query for
`WITNESS_MV_SOLVE_CAP` returned only `results/witness_mixed_volume_profile.txt`, yet the symbol
is defined at `src/mub_zauner_6d_liang_chen.jl:877` and used at
`scripts/julia/profile_witness_mixed_volume.jl:49`. Likewise `build_numeric_fourth_mub_witness_system`
returned zero `.jl` hits although it is defined at `src/mub_zauner_6d_liang_chen.jl:695`.

**Consequence for the team: any negative result from `search_codebase` about Julia source is
unreliable and must not be reported as "not present in the code".** Use the file reader instead.

---

## (4) Cross-component CHM residual — **VERDICT: QUALIFIED** (both numbers CONFIRMED; the attribution error is in the *brief*, not in the repo docs)

### 4.1 What each file actually contains

| File | Contains | Does **not** contain |
|---|---|---|
| `results/locus_classification.txt` | `5.241449e+00`, twice: line 20 (`CHM F6_theta0 vs Dita: residual=5.241449e+00 equiv=false`) and line 22 (`F6 vs Dita cross-locus: residual=5.241449e+00 equiv=false (distinct components)`) | any within-locus pair scan |
| `results/chm_equivalence.txt` | line 17 `Dita lambda pairs: min_residual=2.309398e-02 max_residual=3.163720e+00 all_equiv=False`; line 18 `F6 phi pairs: max_residual=0.000000e+00 all_equiv=True` | any `5.24…` value, any cross-locus comparison |

- `5.241449` occurs in **exactly one file in the whole repository** (`results/locus_classification.txt`,
  lines 20 and 22) — verified by full-tree regex search.
- Producer: `scripts/julia/locus_classification.jl` — the pair loop at lines 244–253 writes the
  line-20 string, and lines 278–280 write the line-22 "cross-locus" string. Tolerance `CHM_TOL = 1e-10`
  (line 19).
- Producer of the other file: `scripts/python/chm_equivalence.py` (`out = …"chm_equivalence.txt"`, line 153;
  header string line 156). It scans **within** a locus: Dita λ-pairs against the fixed anchor λ=0.4
  (lines 4–9) and F6 φ-pairs at λ=0.3 (lines 12–14).

**Scopes genuinely differ and the files are not interchangeable:**
`chm_equivalence.txt` = *within-locus pair sampling* (Dita λ-variation, F6 φ-gauge);
`locus_classification.txt` = *cross-locus component comparison* between the two HP-verified
cluster representatives, F6 arc `(0, 0.5, 0.3)` and Diţă circle `(0.95531661812450919, 0.78539816339744828, 0.4)`.

- Evidence tag for both numbers: `sampled/numerical only` (float64 CHM-residual minimisation over
  720 column permutations × 4 matrix variants × diagonal phases).

### 4.2 The attribution error — corrected

The brief states: *"the status doc cites `chm_equivalence.txt` for a number that file does not contain."*
**This is FALSE as stated, and I decline to repeat it.** Checked directly:

- `docs/scientific_status.md` **never mentions `chm_equivalence.txt` at all** (full-file read; the
  file is 157 lines and cites only the T1/T3/T4 scripts, `w1_D0_groebner.m2`, and the `results/` logs
  at lines 103–122). It makes no 5.24 claim.
- `docs/results/final_honest_status.md` cites `chm_equivalence.txt` twice, at **line 30** ("φ is
  gauge at θ=0", evidence `gauge_analysis.txt`, `chm_equivalence.txt`) and **line 31** ("λ at Diţă:
  CHM-inequivalent 1D family", evidence `chm_equivalence.txt`, Item 1). **Both citations are correct**:
  line 18 of the artifact (`F6 phi pairs: max_residual=0.000000e+00 all_equiv=True`) supports line 30,
  and lines 4–9 + 17 (all six Diţă λ-pairs inequivalent) support line 31.

**So the mis-attribution lives in the incoming task brief, not in any repository document.** The
correct citation for 5.241449e+00 is `results/locus_classification.txt:20,22`. No repo file needs
editing on this account. (I am instructed not to edit docs without a mandate, and none is needed.)

### 4.3 Additional qualification the repo itself supports

`(distinct components)` on line 22 is a **label, not a proof**. The residual only excludes the tested
subgroup (column permutations × {H, Hᵀ, conj H, conj Hᵀ} × diagonal phases). The repo is aware that
this test is narrow: `scripts/python/chm_equivalence.py:100` defines a *second* routine
`chm_equivalence_residual_full` ("Full monomial CHM equivalence: row perm + column perm + diagonal
phases + variants") which is **not** used for the F6-vs-Diţă comparison. So the strongest supportable
statement is: *no element of the column-only tested subgroup maps H_F6 to H_Diţă.* Do not upgrade to
"the two loci are different MUB structures". Tag: `sampled/numerical only`.

Also note the two files use **different implementations and different tolerances** (Python
`TOL_EQUIV` vs Julia `CHM_TOL = 1e-10`) and are never cross-validated against each other.

---
## (1) n_wit=2 mixed volume 119210 vs cap 5000 — **VERDICT: QUALIFIED**

**Headline answer to the priority question: `mixed_volume()` is genuinely COMPUTED in both
artifacts. 119210 is NOT a hardcoded constant and NOT an assumed Bézout-style bound.** However
three qualifications materially weaken Blocker A, one of them severe — see 1.4.

### 1.1 Artifact fields — all CONFIRMED

`results/certify_fourth_mub_witness_dita.txt`:
line 8 `n_eqs = 35`; line 9 `n_vars = 20`; line 10 `mixed_volume = 119210`; line 11
`skipped_solve = true`; line 12 `skip_reason = mixed_volume=119210 exceeds cap 5000`; line 13
`n_paths = 0`; line 21 `verdict = SKIPPED_MV_119210`; line 25 `certify() was NEVER CALLED — the
solve was skipped before tracking.`
`results/witness_mixed_volume_profile.txt`: line 10 `n_wit=2 reduced witness: n_eqs=35 n_vars=20
mixed_volume=119210`; line 11 `WITNESS_MV_SOLVE_CAP = 5000`; line 12 `tractable …? NO`; line 15
`The ALREADY-REDUCED 2-vector witness has mv=119210.`
Every one of these strings is reproduced by the corresponding `println` in the generating scripts,
so artifacts and code are consistent. Tag for the string content: `exact symbolic`
(trivial transcription); the *value* is discussed in 1.2–1.3.

### 1.2 Is 119210 computed? YES — two independent call sites

- `scripts/julia/profile_witness_mixed_volume.jl:45–48`
  `sys2 = build_numeric_fourth_mub_witness_system(H, B3; n_wit=2)` then
  `mv2 = mixed_volume(sys2)`, printed as `mixed_volume=%d`. `sys2` is the **raw overdetermined
  witness system (35 eqs / 20 vars)**. This is the number quoted in the profile file and cited in
  `docs/results/final_honest_status.md:52`.
- `src/mub_zauner_6d_liang_chen.jl:793` `mv = mixed_volume(solve_system)`, inside
  `certify_fourth_mub_witness_at_H`; line 776 additionally computes `mv_full = mixed_volume(system)`.
  The artifact's `mixed_volume` field is this `mv` (written at `certify_fourth_mub_witness.jl:110`).

No hardcoded `119210` and no Bézout-product bound exists anywhere in the Julia source; the number is
a HomotopyContinuation mixed-volume (BKK) computation. Tag: `exact symbolic` for the fact that the
computation is invoked; the resulting integer is a numerical invariant of one formulation, see 1.3.

### 1.3 QUALIFICATION A — which system the artifact's 119210 belongs to (seeded random square system)

In `certify_fourth_mub_witness_at_H` (`src/mub_zauner_6d_liang_chen.jl:773–793`):

```julia
system = build_numeric_fourth_mub_witness_system(H, B3; n_wit=n_wit)   # 35 eq, 20 var
n_eqs = length(system); n_vars = nvariables(system)
mv_full = mixed_volume(system)                                          # line 776
squared_for_certify = n_eqs != n_vars                                  # line 783  -> true
Random.seed!(20260813 + n_wit)                                         # line 788  -> seed 20260815
R = randn(ComplexF64, n_vars, n_eqs)                                   # line 789  -> 20x35 Gaussian
squared_eqs = [sum(R[i,j]*eqs[j] for j in 1:n_eqs) for i in 1:n_vars]  # line 790
solve_system = System(squared_eqs; variables=vars)                     # line 791  -> 20x20
mv = mixed_volume(solve_system)                                        # line 793  <- the artifact value
```

So the `mixed_volume = 119210` printed in `certify_fourth_mub_witness_dita.txt` is the BKK count of a
**seeded random linear-combination (Bertini) square system**, not of the 35-equation witness system
as such. `mv_full` (the true 35-eq count) is computed but **never written to the artifact**, so a
reader cannot distinguish the two from the artifact alone.

**Important precision (correcting an over-strong reading I was asked to consider):** the value is
*not* merely a seed artifact, because the profile script computes the mixed volume of the **raw
35-eq system** and obtains the **same integer 119210**. Equality of the BKK count of a generic linear
section with that of the original overdetermined system is the expected behaviour here. The correct
statement is:

> 119210 is a genuine computed mixed volume of the n_wit=2 witness system *and* of its canonical
> Bertini square section (seed 20260815); the two coincide. It is nonetheless a bound for **one
> choice of equations** — an upper bound on the isolated-root count of that formulation, not an
> invariant of the underlying non-extension question.

**The cap is the only genuinely arbitrary number:** `src/mub_zauner_6d_liang_chen.jl:877`
`const WITNESS_MV_SOLVE_CAP = 5000` is a hardcoded, undocumented tractability threshold with no
derivation, no sensitivity study, and no record of what 5000 paths actually costs. "Not tractable" is
a *policy* statement, not a measurement.

---
### 1.4 QUALIFICATION B (SEVERE) — the n_wit=2 system is RETIRED by repo policy, so Blocker A targets a dead end

`src/Certification.jl:28–40`, verbatim:

```
# Part X — the one-vector witness W1(B3). The logical statement:
#   a fourth MUB extending {I, H, B3} ⇒ W1(H, B3) != empty
#   contrapositive:  W1(H, B3) = empty ⇒ the triple {I, H, B3} does NOT extend
#
# Emptiness of W1 is therefore a COMPLETE certificate of non-extension for
# that triple — strictly stronger than needed. The n_wit=2 system is
# logically redundant for this claim structure (two orthogonal W1-solutions
# are neither necessary nor sufficient for a fourth ONB) and is retired:
# build_numeric_fourth_mub_witness_system(n_wit=2) remains in the legacy
# file for archival reproduction only. It must not be used in new claims.
```

Scope/date: an undated comment block in `src/Certification.jl`; the only dated marker in the same
file is `NORMALIZATION FIX (2026-09-17)` (line 79). The file was last committed in `7ea7a80`
("third commit", 2026-09-20). The retirement therefore **post-dates both artifacts**
(2026-08-13 / 2026-08-14) and applies to them retroactively.

**Implication for Blocker A, stated loudly:** making the n_wit=2 system tractable would mean solving a
system the repository has already declared logically redundant and forbidden for new claims. Blocker A
as written ("get mv below 5000 for n_wit=2") is aimed at the wrong target. The live path is the
n_wit=1 witness (`Certification.jl:42–72`, `w1_witness_certificate`), which returns
`EMPTY_CERTIFIED / CERTIFIED_NUMERICAL` or `INCOMPLETE / OPEN` and enforces
`complete = (n_tracked == mv && n_cert == n_tracked)` at line 148 — a formulation already tractable at
mv=252 (item 2). Note also that `docs/scientific_status.md:136–140` records the normalization hazard
that forced this rewrite; the item-(2) artifacts predate it, which bears on how much weight they
carry (see 2.4).

### 1.5 "certify() was NEVER CALLED" — CONFIRMED, and the wording is honest

The cap test is at line 795 and returns at lines 796–806; `solve()` is line 808 and `certify()` is
line 812 — both **after** the return. So the artifact's line 25 and lines 26–28 ("This is NOT an
interval-arithmetic non-existence certificate… Interval exclusion remains open") are literally true
of the code. The generator branches on exactly this: `certify_fourth_mub_witness.jl:126` selects
that interpretation block when the verdict starts with `SKIPPED_MV`.
**`n_certified = 0` here means "not attempted", not "attempted and found nothing."** The two must
never be conflated in any downstream table. Tag: `exact symbolic` (a control-flow fact).

---
## (2) T3 all-clique tally 22 + 18 = 40/40 — **VERDICT: CONFIRMED (numbers) / QUALIFIED (independence)**

### 2.1 The numbers, pool sizes and clique counts — all CONFIRMED

`results/certify_nwit1_all_cliques_four_classes.txt` (line 3 `reps = {0, 0.4, π/3, 2π/3}`):

| λ | pool | n_6cliques | per-clique lines |
|---|---|---|---|
| 0.0000000000 | 120 | 10 | 6–15 |
| 0.4000000000 | 72 | 4 | 19–22 |
| 1.0471975512 (π/3) | 72 | 4 | 26–29 |
| 2.0943951024 (2π/3) | 72 | 4 | 33–36 |

line 40 `total_cliques=22 total_empty=22 all_pass=true`; line 41 `PASS: …`

`results/certify_nwit1_all_cliques_lambdapi.txt` (line 3 `lambdas = {π, 4π/3, 5π/3}`):

| λ | pool | n_6cliques | per-clique lines |
|---|---|---|---|
| 3.1415926536 (π) | 120 | 10 | 6–15 |
| 4.1887902048 (4π/3) | 72 | 4 | 19–22 |
| 5.2359877560 (5π/3) | 72 | 4 | 26–29 |

line 33 `total_cliques=18 total_empty=18 all_pass=true`; line 34 `PASS: …`

22 = 10+4+4+4 and 18 = 10+4+4, so **40 = 22 + 18** arithmetically. The λ sets are disjoint
({0, 0.4, 1.047, 2.094} vs {3.142, 4.189, 5.236}) → seven distinct λ values, none reused.
Tag: `certified numerical (interval arithmetic)` for the emptiness statements; the *pool
enumeration* is `sampled/numerical only` (see 2.4).

### 2.2 THE DOUBLE-COUNTING QUESTION — unambiguous answers

**(a) Is the `lambdapi` run an INDEPENDENT solve at λ+π, or the SAME cliques transported?**
→ **It is an INDEPENDENT solve, with natively regenerated cliques. NOT a transport.**
Proof from `scripts/julia/certify_nwit1_all_cliques_lambdapi.jl`:
- line 52 `H = build_karlsson_family(DITA_THETA, DITA_PHI, λ)` — fresh H at λ+π;
- line 53 `pool, cliques = all_six_cliques(H)` → lines 18–29: `generate_candidate_pool_fresh(H)`,
  `deduplicate_pool`, `_orthogonality_graph`, `maximal_cliques` — a **fresh pool and fresh cliques
  computed at λ+π**, not a permuted copy of the λ-cliques;
- lines 59–60 `certify_fourth_mub_witness_at_H(H; n_wit=1, clique_indices=idx, B3=B3)` — a **fresh
  HC solve + certify** per clique.

The clique index lists in the two files are indeed disjoint (λ=0 clique 1 = `[32,52,51,59,31,17]`;
λ=π clique 1 = `[60,56,79,43,36,46]`).

The **transported** experiment is a *different script*: `scripts/julia/verify_chm_b3_transfer.jl`
line 149 `pool, cliques = all_six_cliques(H_src)` (cliques enumerated at λ_src), line 155
`B3p = transform_b3(B3, mp.D_r)`, lines 162–163 `certify_fourth_mub_witness_at_H(H_tgt; n_wit=1,
clique_indices=idx, B3=B3p)` — transporting the λ_src cliques to λ+π and re-solving there.

**(b) Is 40/40 really 40 independent certificates?**
→ **YES for the tally as published: 40/40 = 40 independent HC solves, with NO double counting in
any current document.** Every place the tally appears — `docs/scientific_status.md:121–122`,
`docs/overview/forensic_ground_truth_report.md:21–23` and `76–80`,
`docs/overview/claim_audit.md:13`, `docs/overview/known_results_table.md:15`, and the T3 statement
in `docs/paper/proofs/fourth_mub_theorems.tex` — derives 40 from the two all-clique logs only,
22 + 18. `verify_chm_b3_transfer.txt` is nowhere added to the 40.

**Two real qualifications, and one live double-counting hazard:**

**Q1 — independence holds at the level of *computations*, not of *CHM classes*.** The repo already
says this in its own voice: `results/verify_nwit1_referee_checks.txt:22` reads verbatim
`FAIL: some certified λ are CHM-equivalent — certificates are not independent.`, with lines 7 and 11
giving `λ=1.047198 vs λ=4.188790: residual=1.430483e-15 equiv=true` and
`λ=2.094395 vs λ=5.235988: residual=1.596843e-15 equiv=true`; `results/verify_nwit1_followup.txt:5–9`
gives all five `λ vs λ+π` residuals at 4.6e-16 … 2.2e-15, `equiv=true`, and line 10 names the
representatives `{0, 0.4, π/3, 2π/3}`. The generator admits it too: the header of
`certify_nwit1_all_cliques_four_classes.jl:3` reads "(λ+π duplicates are CHM-equivalent and
skipped.)". **So: 40 independent solves, but only 4 independent CHM classes.** The paper's own
scope remark already words it correctly — "seven explicit Dita-λ values (**four CHM classes**)"
(`results/track_c_elimination/arxiv_format_pass/before_fourth_claims.txt:72`). T3 stays "local
certified numerical, 7 λ points / 4 CHM classes"; the 18 λ+π points add numerical confirmation of
an already-covered class, **not** new mathematical content. The mirror-image confirmation — pool
sizes match class-by-class (120 at λ=0 and λ=π; 72 at the other five) and clique counts match
(10, 4, 4) — is consistent with the equivalence but is not itself a proof of it.

**Q2 — LIVE DOUBLE-COUNTING HAZARD: the transfer file's 18 are the SAME cliques as the 22.**
`results/verify_chm_b3_transfer.txt` covers 3 pairs with 10+4+4 = 18 cliques, and its per-pair
headers (`src pool=120 n_6cliques=10`; `src pool=72 n_6cliques=4`; `src pool=72 n_6cliques=4`;
lines 9, 25, 35) are *identical* to the four-classes pool/clique counts — because they are literally
the same recovered cliques in transported coordinates. **If anyone adds the transfer file's 18 to
the 40, that is double counting.** Rule for the team: *40 = four_classes(22) + lambdapi(18) only;
the transfer file corroborates the 22 and is never an additional 18.* The transfer file is itself a
genuine independent *solve* (line 162), so the precise description is "18 further solves of the
same 18 bases in new coordinates" — corroborating, not additive.

---
### 2.3 What `cert=252 fullpass=0` means — precise, without inflation

Per clique the artifacts record `mv=252 tracked=252 cert=252 fullpass=0 found4=false
verdict=EMPTY_SQUARE_252`. Decoded against the code:
- n_wit=1 ⇒ 10 vars, 17 eqs ⇒ `squared_for_certify = true`
  (`src/mub_zauner_6d_liang_chen.jl:783`); a 10×10 Bertini square system is formed from the 17
  equations with `Random.seed!(20260813 + 1) = 20260814` (lines 788–791).
- `mv=252` = mixed volume of that square system; `tracked=252` ⇒ **all 252 paths tracked, none lost**;
  `cert=252` ⇒ **all 252 endpoints certified by interval arithmetic**. Complete tracking and complete
  certification — this is the strong part.
- `fullpass=0` = `n_certified_pass_full_residual = 0`: **none** of the 252 certified square-system
  roots satisfies all 17 original equations (residual gate `rmax < 1e-8`, lines 817–830).
- Inference chain: because `R` is a 10×17 Gaussian (generically rank 10), every solution of the square
  system is a solution of the 17-eq system; if no certified square root passes the full residual, the
  square system has no solutions, hence neither does the 17-eq system, hence `W1 = ∅`.

**Strength, stated exactly:** the interval-arithmetic certification is applied to the **square**
system, not to the 17-eq overdetermined witness ideal. The final step rests on the rank of `R`, a
genericity property that is **not recorded or checked in any artifact** (no rank or conditioning
output exists anywhere). So the honest tag is `certified numerical (interval arithmetic)` **plus an
unverified rank assumption** — not an exact proof of ideal emptiness, and not an interval
*exclusion* certificate. The artifact name `EMPTY_SQUARE_252` is therefore accurate, and these
artifacts do not overstate themselves. `found4=false` is combinatorial (clique graph), not certified.

### 2.4 Further caveats that bound this claim (all pre-existing, all material)

- **Pool completeness is not certified.** `generate_candidate_pool_fresh` + `deduplicate_pool` is a
  float64 numerical pool; "every size-6 clique" means every clique *in the recovered pool graph*
  (`ortho_tol=1e-8`). The paper's own scope remark says so: "it does not by itself certify against
  hypothetical third MUBs missed by an incomplete root track"
  (`results/track_c_elimination/arxiv_format_pass/before_fourth_claims.txt:74`).
- **The clique set is not deterministic.** `results/verify_nwit1_referee_checks.txt:26–55` re-runs the
  selection 3× per λ and reports `NONDETERMINISTIC: clique set varies` at all six λ tested (counts
  4, 4, 4, 10, 4, 4 are stable; the index sets are not), and line 57 warns "clique indices vary
  across runs". "40/40 recovered-pool six-cliques" therefore denotes a *per-run recovered* set, not
  40 canonically pinned objects. Taking *all* cliques, as these two scripts do, mitigates but does not
  remove this.
- **Pre-normalization-fix.** These artifacts are dated 2026-08-13/14; `src/Certification.jl:79–82`
  records a 2026-09-17 correction of the B3 RHS (RHS 1 vs 6) that "fixes the prior defect where
  unit-normalized pool vectors were incorrectly tested with RHS=6". The legacy path used here does
  auto-detect (`_witness_mu_to_basis`, `src/mub_zauner_6d_liang_chen.jl:669–684`), but comparisons
  with post-fix results must be checked for convention consistency, as
  `docs/scientific_status.md:136–140` warns.
- **Scope unchanged:** not a full-circle result, not `K_6^(3)`, not `N(6)=3`. No upgrade.

---
## (3) Arc width 0.117643222422 rad — **VERDICT: CONFIRMED** (one caveat is *not* reproducible)

`results/f6_boundary_reverify.txt`:
- line 2 `Date: 2026-08-03T14:58:19.856  HP_bits=400`;
- line 3 `theta=0, phi=0.5 [gauge], lambda0=0.3`;
- line 5 `anchor (0,0.5,0.3): clique=6 hp=true`;
- line 7 `sign=+1: boundary |Delta lambda| in [0.0550448030233, 0.0550448030466]  mid=0.055044803035`;
- line 8 `sign=-1: boundary |Delta lambda| in [0.0625984193757, 0.062598419399]  mid=0.0625984193874`;
- line 10 `arc width (|Delta+| + |Delta-|) = 0.117643222422`.

**Arithmetic re-derivation (done by hand; no code was run):**
0.055044803035 + 0.0625984193874 = 0.1176432224224 → printed as `0.117643222422` under `%.12g`.
Midpoints also check: (0.0550448030233+0.0550448030466)/2 = 0.05504480303495 ≈ 0.055044803035;
(0.0625984193757+0.062598419399)/2 = 0.06259841938735 ≈ 0.0625984193874. Both bracket widths are
2.33e-11, consistent with the script's stopping rule `rel_tol = 1e-10` against `lm0 = 0.3` (3e-11).
**CONFIRMED: the number is exactly |Δ+| + |Δ−| from the two bracketed boundaries.**

`scripts/julia/f6_boundary_reverify.jl` confirms provenance: line 6 `const HP_BITS = 400` (printed at
line 70); line 91 `arc_width = sum(abs.(boundaries))`, where `boundaries` collects only the midpoints
(line 87); line 93 prints with `%.12g`; the anchor probe is lines 75–76. Tag: `sampled/numerical only`
for the width — the point probes use 400-bit HP interval arithmetic, but the pool and clique search
that define "has a 6-clique" are float64, so this is not an exact symbolic value of the true locus.

**The `phase1_boundaries.txt` caveat — reproducible as a code choice, NOT as an empirical result.**
`results/phase1_boundaries.txt:12` states: `NOTE: coarse bracket step (0.05) required; fine step
(0.01) finds spurious early dips.` Status:
- **Reproducible from the script as a parameter choice:** `f6_boundary_reverify.jl:33`
  `bracket_axis(lm0, sign; max_scan = 0.15, step = 0.05)` hardcodes `step = 0.05`, and lines 31–32
  carry the comment "Coarse step (0.05) avoids spurious early clique dips seen with step=0.01
  (matches extended_axis_sweep.jl --quick bracket that yields +0.055 / -0.063)". The stored brackets
  are exactly what a 0.05-step scan produces (`inner = 0.05`, `outer = 0.10` for both signs, given
  true boundaries 0.0550 and 0.0626), so the code path is self-consistent.
- **NOT reproducible as an empirical finding:** no stored artifact demonstrates "spurious early dips"
  at step = 0.01 for this scan. The only stored step-0.01 fine sweeps are in
  `results/lambda_periodicity_dita.txt:16–20`, a *different* scan, and they show **no** dips
  (`points=14 clique>=6=14`; `points=11 clique>=6=11`). Tag: `not run` by me; unsupported by any
  stored run.
- **A concrete mechanism does exist in the code**, which supports the concern without proving the
  claim: `bracket_axis` initialises `prev6 = true` unconditionally (line 37), so if the very first
  probe already returned `clique < 6` the routine would report `outer = 0.05` and silently emit a
  boundary that is too small. The probe is tolerance-fragile: `probe_point` uses `MU_TOL = 1e-12` and
  `ORTHO_TOL = 1e-12` (lines 7–8) — two orders of magnitude tighter than the 1e-8 used in the
  item-(2) clique scripts — on a float64 pool, so an isolated dip at a fine scan point is a
  plausible failure mode. Follow-up: seed `prev6` from an actual probe at d → 0, and log the
  per-d `clique` values so any dip is visible in the artifact.

---
## Summary table

| # | Number | Verdict | Primary evidence file:line | Evidence tag |
|---|---|---|---|---|
| 4 | `5.241449e+00` cross-locus CHM residual | **CONFIRMED**; attribution in the brief is wrong (it lives in `locus_classification.txt`, not `chm_equivalence.txt`); "distinct components" is a label from a column-only test | `results/locus_classification.txt:20,22`; producer `scripts/julia/locus_classification.jl:244–253,278–280` | `sampled/numerical only` |
| 4 | `chm_equivalence.txt` summary (Dita max 3.163720e+00; F6 max 0.0) | **CONFIRMED** as a *different, within-locus* scan | `results/chm_equivalence.txt:9,17,18`; producer `scripts/python/chm_equivalence.py:153` | `sampled/numerical only` |
| 1 | `mixed_volume = 119210` | **CONFIRMED as genuinely COMPUTED** (not hardcoded, not a Bézout bound) — holds for the raw 35-eq system *and* its seeded Bertini square section; **QUALIFIED**: cap 5000 is arbitrary, and n_wit=2 is retired by policy | `results/certify_fourth_mub_witness_dita.txt:10`; `results/witness_mixed_volume_profile.txt:10,15`; `scripts/julia/profile_witness_mixed_volume.jl:46`; `src/mub_zauner_6d_liang_chen.jl:776,793` | `exact symbolic` (computation invoked) / `sampled/numerical only` (the BKK integer) |
| 1 | n_eqs=35, n_vars=20, skipped_solve=true, n_paths=0, verdict=SKIPPED_MV_119210, certify() never called | **CONFIRMED** | `results/certify_fourth_mub_witness_dita.txt:8–13,21,25`; `src/mub_zauner_6d_liang_chen.jl:795–806` (return) vs `808,812` (solve/certify) | `exact symbolic` (control flow) |
| 2 | 22 + 18 = 40/40 over 7 λ | **CONFIRMED** (pools 120/72/72/72 and 120/72/72; counts 10+4+4+4 and 10+4+4) | `results/certify_nwit1_all_cliques_four_classes.txt:5,18,25,32,40`; `…_lambdapi.txt:5,18,25,33` | `certified numerical (interval arithmetic)` + `sampled/numerical only` (pool) |
| 2 | Double counting? | **NO double counting in the published 40/40.** lambdapi is an independent solve with fresh cliques; the *transfer* file's 18 ARE the 22 transported — adding them to 40 **would** double count. Only **4 independent CHM classes**, per the repo's own `verify_nwit1_referee_checks.txt:22` | `scripts/julia/certify_nwit1_all_cliques_lambdapi.jl:52,53,59–60` vs `scripts/julia/verify_chm_b3_transfer.jl:149,155,162`; `results/verify_chm_b3_transfer.txt:9,25,35` | `exact symbolic` (code path) |
| 2 | `fullpass=0` | **CONFIRMED**; = no certified *square*-system root satisfies all 17 eqs; the inference to `W1 = ∅` additionally needs rank(R)=10, which no artifact records | `src/mub_zauner_6d_liang_chen.jl:783–793,817–830` | `certified numerical (interval arithmetic)` + unverified rank assumption |
| 3 | arc width `0.117643222422` = 0.055044803035 + 0.0625984193874 | **CONFIRMED**; HP_bits=400 and anchor `(0,0.5,0.3) clique=6 hp=true` confirmed | `results/f6_boundary_reverify.txt:2,5,7,8,10`; `scripts/julia/f6_boundary_reverify.jl:6,87,91,93` | `sampled/numerical only` (400-bit HP probes, float64 pool) |
| 3 | "coarse step 0.05 required; fine step finds spurious dips" | **QUALIFIED**: confirmed as a hardcoded parameter + comment; **not** reproducible from any stored artifact | `scripts/julia/f6_boundary_reverify.jl:31–33,37`; `results/phase1_boundaries.txt:12`; counter-evidence `results/lambda_periodicity_dita.txt:16–20` | `not run` |

**No claim was upgraded. T3 remains "local certified numerical, 7 λ points / 4 CHM classes"; no
full-circle or full-`K_6^(3)` certificate exists, and none was produced.**

---
## Open questions for the team

1. **Blocker A target (highest priority).** Given `src/Certification.jl:34–39` retires n_wit=2 as
   "logically redundant … must not be used in new claims", is Blocker A (lower the n_wit=2 mixed
   volume below 5000) still the intended target, or should it be re-aimed at the live n_wit=1 path
   (`w1_witness_certificate`, mv=252)? My reading: re-aim it. Who owns that decision?
2. **Where is `mv_full` (the 35-eq system's own mixed volume) recorded?** It is computed at
   `src/mub_zauner_6d_liang_chen.jl:776` and then dropped. `certify_fourth_mub_witness_dita.txt`
   exposes only the randomized square section's count, and neither `mv_full`, nor
   `squared_for_certify`, nor the seed, is written to any artifact. Recommend adding them to the
   artifact header (as a new file, or with a timestamped `.bak_2026-09-26` archive, per the
   no-silent-overwrite rule).
3. **What justifies `WITNESS_MV_SOLVE_CAP = 5000`?** No derivation, no measured cost per path, no
   sensitivity analysis. Either document the basis or make it a CLI parameter, so that "not
   tractable" stops being an unfalsifiable policy statement.
4. **Should the rank/conditioning of the Bertini matrix `R` be recorded?** The `W1 = ∅` inference in
   the 40/40 chain silently needs rank(R) = n_vars. One line of output would close the gap between
   "certified numerical" and a defensible exclusion argument.
5. **Confirm the 40/40 accounting rule in writing** somewhere durable: 40 = four_classes(22) +
   lambdapi(18); `verify_chm_b3_transfer.txt` corroborates the 22 and is never additive. This is the
   one place where a future edit could silently double count.
6. **Should the 40/40 be restated as "40 solves over 4 CHM classes"?** The repo's own referee-check
   artifact says "certificates are not independent" (`verify_nwit1_referee_checks.txt:22`) while the
   status doc and paper say "seven points". Both are true but invite opposite readings; the paper's
   "four CHM classes" parenthetical is the right fix and should propagate to
   `docs/scientific_status.md:38,121–122`.
7. **F6-vs-Diţă "distinct components":** the comparison is column-only. Should we run
   `chm_equivalence_residual_full` (`scripts/python/chm_equivalence.py:100`, row+column) before that
   label is used in prose? Currently it excludes only the tested subgroup.
8. **The arc-width caveat needs evidence or deletion.** Either store a step=0.01 run that shows the
   dips, or drop the claim; and fix the unconditional `prev6 = true` initialisation at
   `f6_boundary_reverify.jl:37`, which can silently under-report a boundary.
9. **Tooling:** `search_codebase` does not index `.jl` files, so negative grep results about Julia
   source from earlier rounds are unreliable and should be re-checked with the file reader before
   being repeated in any report.
10. **Pool-completeness and clique nondeterminism** (§2.4), not the count 40, are the real limiting
    factors on T3. Is a certified pool-completeness argument — or a fixed, hashed clique list per λ —
    the right next increment?

---

## Method / reproducibility notes for this report

- No file in `results/` was modified or overwritten. The only new files written by this agent are
  `results/NOT_RUN_verifier_round1.log` (a record of experiments deliberately **not** run) and this
  report at `research/reports/verification_2026-09-26.md`.
- No Julia/Python computation was executed. Every number above was re-derived by reading the artifact
  text and the generating source, plus hand-arithmetic for item (3).
- No external theorem is asserted from memory. Where the repo itself cites literature
  (Brierley–Weigert, Phys. Rev. A 79, 052316 (2009), arXiv:0901.4051, as referenced in
  `src/dita_third_mub_construction.jl:4`), I reproduce the repo's own attribution and make no new
  claim about it. **No new literature claim should be added to any document on the basis of this
  report without an actual search-and-cite pass.**







