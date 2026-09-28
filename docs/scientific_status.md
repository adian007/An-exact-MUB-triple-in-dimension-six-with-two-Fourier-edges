# Scientific calculation status

## 0. Changelog

**2026-09-26 — documentation reconciliation. No new computation was executed
for this revision.** Every number added below was re-derived by reading the
stored artifact together with the source that produced it, and each carries an
explicit evidence tag. No claim strength was upgraded. T3 remains a local
certified numerical result at seven \(\lambda\) points; no full-circle or full
\(K_6^{(3)}\) certificate exists, and none was produced here. The theorem
statements in `paper/` were not modified.

Reconciled against the forensic ground-truth audit
`docs/overview/forensic_ground_truth_report.md` and two independent round-one
reports, `research/reports/verification_2026-09-26.md` and
`research/reports/csv_reconciliation_analysis_2026-09-26.md`. The previous
revision is archived as `docs/scientific_status.md.bak_2026-09-26`.

* §2 and §4: recorded the exact T3 arithmetic (pools 120/72/72/72 and
  120/72/72; clique counts 10+4+4+4 and 10+4+4; seven disjoint \(\lambda\)),
  and the limitation that those 40 certificates span only four independent
  CHM classes.
* §4 and §5: added the negative findings at equal visibility — the
  `n_wit=2` system is retired by repository policy; the mixed-volume cap of
  5000 is an undocumented policy threshold rather than a measurement; the
  `fullpass=0` inference rests on a \(\mathrm{rank}\,R=10\) assumption that no
  artifact records; and the F6 coarse-step caveat is unsupported by any stored
  run.
* §4: added the F6 arc width and the cross-locus CHM residual, both tagged
  `sampled/numerical only`.
* §6: rewrote the 848/15-versus-928/45 note. The discrepancy is now **explained,
  not resolved**; `results/csv_reconciliation.txt` turns out to be irrelevant
  to it. The guard on aggregate counts is retained unchanged.
* `research/claims/claim_ledger_campaigns.csv`: appended thirteen rows covering
  the verified items and the negative findings. The eleven Campaign 0 rows are
  untouched (verified: all eleven original lines still present verbatim).

* **Part 2 cross-reference (flagged for your review, no tier added).** The Blocker C
  investigation in `research/reports/audit_followup_2026-09-26.md` proved the
  arc-versus-circle mechanism currently attributed to an Hermitian/non-Hermitian
  difference in the A-block to be **false**: both A-blocks are non-Hermitian, and the
  repository's own artifact `results/arc_circle_asymmetry.txt` records
  `A Hermitian=false` for both loci (lines 7, 13) while asserting the contrast in its
  prose (lines 35, 50–51). Tag: `exact symbolic` for the refutation. This file was
  **not** changed to carry that correction, because `docs/results/final_honest_status.md`
  is outside the set of files cleared for direct editing; the affected row is
  `docs/results/final_honest_status.md:35` ("Dita vs F6_theta0 topology contrast |
  HP-supported structural | A-block geometry"). That row should be downgraded and
  re-worded in a separate reviewable diff. Lemma L3 itself is correct and needs no
  change; the false wording was added downstream of it.

## 1. Problem and mathematical objects

The project studies mutually unbiased bases (MUBs) in \(\mathbb C^6\), mainly
after fixing the computational basis \(I\). A second basis is represented by an
order-six complex Hadamard matrix \(H\), with

\[
|H_{jk}|=1,\qquad H H^\ast = 6I.
\]

The principal parameter space is Karlsson's three-parameter family
\(K_6^{(3)}\), with a Diţă slice obtained by fixing
\[
\theta=\arccos(1/\sqrt 3),\qquad \phi=\pi/4,
\]
and varying \(\lambda\). The computational question is whether a third basis
\(B_3\) exists that is unbiased to both \(I\) and \(H\), and whether a fourth
basis can extend \(I,H,B_3\).

For a flat vector \(v\), the projective constraints are:

\[
|v_j|^2=1/6,\qquad |\langle h_k,v\rangle|^2=1/6
\]

for every column \(h_k\) of \(H\), with analogous equations for every column
of \(B_3\). The repository's `W_1` witness system encodes one such vector;
emptiness of `W_1` is a necessary-condition obstruction to a fourth MUB.

## 2. What has actually been calculated

| ID | Calculation | Method | Evidence strength | Current scope |
|---|---|---|---|---|
| T1 | Gauge structure at \(\theta=0\), and separation of \(\lambda\) values on the Diţă slice | Exact symbolic/algebraic identities in the T1 script | Exact, assuming the stated formulas | Karlsson family and specified slice |
| T2 | Third-MUB construction and periodicity on the Diţă circle | HomotopyContinuation/high-precision numerical computation plus residual checks | Certified/numerical findings, not an analytic family proof | 126 periodicity samples and 628 dense probes |
| T3 | No fourth-MUB witness at seven selected Diţă points | Square-system homotopy solving, interval certification of a seeded \(10\times10\) Bertini section, and full residual checks over recovered six-cliques | Local certified numerical result; emptiness inferred via a square section, with an unrecorded \(\mathrm{rank}\,R=10\) assumption (see §4) | 40/40 recovered-pool six-cliques at 7 disjoint \(\lambda\) values, spanning only 4 independent CHM classes (see §4) |
| T4 | Exact elimination of the one-vector fourth-MUB witness at a \(D_0\)-equivalent point | Macaulay2 Gröbner basis over \(\mathbb Q(\zeta_{24},\sqrt5)\) | Exact local algebra | One \((D_{\rm bc},F_D)\) pair; \(\lambda\in\{\pi/2,3\pi/2\}\) |
| C-search | Finite scans of special loci and the Diţă circle | Floating/high-precision numerical pool and clique computations | Sampled evidence only | 928-row search, 1863-row S\* run, 628-point dense circle; row counts from the two 848/928 views are not yet reconciled (§6) |
| V1 | Formula/transcription audit of Karlsson variants | Double-precision residual tests | Numerical implementation audit | 349 parameter points per variant |

## 3. Calculation pipeline

### 3.1 Constructing \(H\)

`src/Karlsson.jl` and the legacy family code construct \(H(\theta,\phi,\lambda)\)
from the \(2\times2\) blocks \(A,B\), Möbius-derived phases, and the
four-block Hadamard assembly. The current code distinguishes the original
Karlsson transcription from a literal review transcription because the latter
fails the unitary/Hadamard checks in the repository's 349-point audit.

### 3.2 Finding third-basis vectors

The pool pipeline solves the MU equations for vectors unbiased to \(I\) and
\(H\), rejects candidates failing conjugacy or residual checks, then performs
projective deduplication and constructs an orthogonality graph. A six-clique
is interpreted as a candidate third basis. Pool sizes depend on the point:
the stored T3 logs show 120-vector pools at \(\lambda=0,\pi\) and 72-vector
pools at the other selected points.

### 3.3 Testing a fourth basis

For a selected third basis \(B_3\), the `W_1` system searches for one flat
vector unbiased to \(I,H,B_3\). A fourth MUB would contain six such vectors,
so:

\[
W_1=\varnothing\quad\Longrightarrow\quad
\text{no fourth MUB extends } \{I,H,B_3\}.
\]

This implication is one-way. A nonempty one-vector witness system would not
prove that six mutually orthogonal witnesses exist.

The numerical T3 workflow squares a witness system for a homotopy solve,
tracks all roots, certifies roots, and then checks the unsquared/full residual.
The exact T4 workflow instead uses named algebraic coefficients and a
Macaulay2 Gröbner basis. These are different evidentiary categories and must
not be merged under the word “proof”.

### 3.4 Exact T4 field

The T4 coefficients lie in
\[
K=\mathbb Q(\zeta_{24},\sqrt5).
\]

The stored field-degree calculation gives
\[
[\mathbb Q(\zeta_{24}):\mathbb Q]=8,\qquad [K:\mathbb Q]=16.
\]

The Macaulay2 log reports 17 generators in 10 variables and Gröbner basis
\(\{1\}\), hence the ideal is the unit ideal and the exact \(W_1\) variety is
empty for that fixed pair. This verifies a local obstruction; it does not
eliminate all third bases or all \(\lambda\).

## 4. Evidence ledger

### Exact or algebraic

* T1 gauge identities: `scripts/julia/formalize_gauge_lemmas.jl`,
  `results/formalize_gauge_lemmas.txt`.
* T4 field and Gröbner calculation:
  `scripts/julia/verify_field_degree_w1.jl`,
  `symbolic_export/w1_D0_groebner.m2`,
  `results/track_c_elimination/w1_D0_field_degree.txt`,
  `results/track_c_elimination/w1_D0_groebner.log`.
* The exact T4 result is a re-verification of a stronger published
  \(D_0\) nonextendability result, not a new global theorem.

### Certified numerical

* T3 scripts:
  `scripts/julia/certify_nwit1_all_cliques_four_classes.jl`,
  `scripts/julia/certify_nwit1_all_cliques_lambdapi.jl`.
* Outputs:
  `results/certify_nwit1_all_cliques_four_classes.txt`,
  `results/certify_nwit1_all_cliques_lambdapi.txt`.
* The stored totals are 22/22 and 18/18 empty cliques, respectively.
  Together they cover the seven selected points and 40/40 cliques. The
  arithmetic is explicit in the artifacts: pools of 120/72/72/72 and
  120/72/72 vectors, clique counts \(10+4+4+4\) and \(10+4+4\), over the
  disjoint sets \(\{0,0.4,\pi/3,2\pi/3\}\) and \(\{\pi,4\pi/3,5\pi/3\}\).
  Tag for the emptiness statements: `certified numerical (interval
  arithmetic)`.
* Four limitations bound this entry, and all four are material.
  **(i) Four CHM classes, not seven independent points.** The 40 certificates
  are 40 independent *computations* but only four independent CHM classes:
  the repository's own referee check states `FAIL: some certified lambda are
  CHM-equivalent — certificates are not independent`
  (`results/verify_nwit1_referee_checks.txt:22`), with residuals \(1.43\cdot
  10^{-15}\) and \(1.60\cdot10^{-15}\) for the pairs \(\lambda=1.047198\) vs
  \(\lambda=4.188790\) and \(\lambda=2.094395\) vs \(\lambda=5.235988\) (lines
  7, 11), and all five \(\lambda\) vs \(\lambda+\pi\) residuals at
  \(4.6\cdot10^{-16}\) to \(2.2\cdot10^{-15}\) in
  `results/verify_nwit1_followup.txt:5–9`. The generating script says the same
  in its header (`certify_nwit1_all_cliques_four_classes.jl:3`). The \(\lambda+\pi\)
  run therefore adds numerical confirmation of an already-covered class, not
  new mathematical content. Tag: `sampled/numerical only`.
  **(ii) Emptiness is inferred, not directly certified.** Per clique the
  artifacts record `mv=252 tracked=252 cert=252 fullpass=0`; the interval
  certification is applied to a seeded \(10\times10\) square section
  (`src/mub_zauner_6d_liang_chen.jl:783–793`), and `fullpass=0` means no
  certified square-system root satisfies all 17 original equations (residual
  gate `rmax < 1e-8`, lines 817–830). The step from there to \(W_1=\varnothing\)
  additionally requires the \(10\times17\) Gaussian \(R\) to have rank 10. That
  rank is a genericity property which **no artifact records**; it is an
  unverified assumption, and these artifacts are not an interval *exclusion*
  certificate. The artifact name `EMPTY_SQUARE_252` is accurate.
  **(iii) Pool completeness is not certified.** "Every six-clique" means every
  clique in the float64 recovered pool graph at `ortho_tol=1e-8`; nothing here
  certifies against a hypothetical third MUB missed by an incomplete root
  track. **(iv) The clique set is not deterministic.** Re-running the selection
  three times per \(\lambda\) varies the index sets at all six \(\lambda\) tested
  (`results/verify_nwit1_referee_checks.txt:26–55`), so "40/40 recovered-pool
  six-cliques" denotes a per-run recovered set, not 40 canonically pinned
  objects. Taking all cliques mitigates this without removing it.
* No double counting occurs in the published 40/40, and one future edit could
  introduce it. The \(\lambda+\pi\) run is an *independent* solve with freshly
  generated cliques (`scripts/julia/certify_nwit1_all_cliques_lambdapi.jl:52,
  53, 59–60`), not a transport. The separate transport experiment
  `results/verify_chm_b3_transfer.txt`, by contrast, re-solves the *same* 18
  cliques in transported coordinates, so its counts must never be added to the
  40. Standing rule: 40 = 22 (four classes) + 18 (\(\lambda+\pi\)) only.
  Tag: `exact symbolic` (control-flow fact).

### Numerical or sampled

* `results/dita_lambda_fourth_dense.csv`: 628 samples, no fourth result;
  not a theorem for the continuous circle.
* `results/special_loci_search.csv`: 928 rows, 899 pool-complete.
* `results/special_loci_degen1845.csv`: 1863 rows, 1852 pool-complete,
  21 six-clique loci, 0 fourth-MUB detections.
* `scripts/python/audit_karlsson_variants.py`: 348/349 for the original
  transcription and 0/349 for the literal review transcription.
* F6 arc width at \(\lambda_0=0.3\), \(\theta=0\): the \(\lambda\)-arc on which
  the \(\phi=0.5\) locus retains a six-clique has width
  \[
  0.117643222422 = 0.055044803035 + 0.0625984193874
  \]
  radians, from the two bracketed boundaries in
  `results/f6_boundary_reverify.txt:7,8,10`; the anchor probe
  `(0,0.5,0.3)` is recorded as `clique=6 hp=true` at line 5 with
  `HP_bits=400` at line 2. The probe is 400-bit high-precision interval
  arithmetic, but the pool and clique search that define "has a six-clique"
  are float64, so this is not an exact value for the true locus. Tag:
  `sampled/numerical only`.
* The caveat attached to that arc width — "coarse bracket step (0.05)
  required; fine step (0.01) finds spurious early dips"
  (`results/phase1_boundaries.txt:12`) — is **not reproducible from any stored
  artifact**. It exists in the repository only as a hardcoded parameter plus a
  code comment (`scripts/julia/f6_boundary_reverify.jl:31–33`), and the stored
  step-0.01 sweeps that do exist belong to a different scan and show no dips
  (`results/lambda_periodicity_dita.txt:16–20`). Tag: `not run`. A concrete
  mechanism for the worry does exist in code — `bracket_axis` initialises
  `prev6 = true` unconditionally at line 37, so a dip at the very first probe
  would be silently under-reported — but that is a code reading, not evidence
  that the dips occur.
* Cross-locus CHM residual between the F6 arc representative and the Diţă
  circle representative: \(5.241449\), reported twice in
  `results/locus_classification.txt:20,22`. Tag: `sampled/numerical only`
  (float64 CHM-residual minimisation). The parenthetical "distinct
  components" on line 22 is a **label from a column-only CHM test**, not an
  independently proven component decomposition: the residual excludes only the
  tested subgroup of column permutations, matrix variants and diagonal
  phases. The broader routine `chm_equivalence_residual_full`
  (`scripts/python/chm_equivalence.py:100`, row *and* column) exists but is
  not used for this comparison. The supportable statement is that no element
  of the tested subgroup maps \(H_{\mathrm{F6}}\) to \(H_{\mathrm{Diță}}\).
  This figure is unrelated to `results/chm_equivalence.txt`, which holds a
  different, within-locus scan (Diţă max \(3.163720\), F6 max \(0\)).
* The following three entries are **negative or non-reproducible results**,
  recorded here with the same visibility as the positive findings above.
* **The witness mixed volume was never attempted, not attempted and empty.**
  The n_wit=2 witness system reports `mixed_volume = 119210` against
  `WITNESS_MV_SOLVE_CAP = 5000`, with `n_eqs = 35`, `n_vars = 20`,
  `skipped_solve = true`, `n_paths = 0` and verdict `SKIPPED_MV_119210`
  (`results/certify_fourth_mub_witness_dita.txt:8–13,21`). The artifact's own
  note "`certify()` was NEVER CALLED — the solve was skipped before tracking"
  (line 25) is literally true of the code: the cap test returns at
  `src/mub_zauner_6d_liang_chen.jl:795–806`, before `solve()` at line 808 and
  `certify()` at line 812. So `n_certified = 0` means **not attempted**, and
  must never be tabulated alongside a genuine zero count. Tag: `exact
  symbolic` (control flow). The value 119210 *is* genuinely computed, not
  hardcoded and not a Bézout bound, at
  `scripts/julia/profile_witness_mixed_volume.jl:46` and
  `src/mub_zauner_6d_liang_chen.jl:776,793`; tag for the call `exact
  symbolic`, for the BKK integer `sampled/numerical only`. It is an upper
  bound for one choice of equations, not an invariant of the underlying
  non-extension question.
* **The cap of 5000 is a policy threshold, not a measurement.**
  `WITNESS_MV_SOLVE_CAP = 5000`
  (`src/mub_zauner_6d_liang_chen.jl:877`) is hardcoded, has no derivation, no
  measured cost per path and no sensitivity analysis. "Not tractable" is
  therefore a policy statement and not an experimental result.
* **The n_wit=2 system is retired by repository policy.** `src/Certification.jl:28–40`
  states that emptiness of the one-vector witness \(W_1\) is already a
  complete certificate of non-extension for a fixed triple, that the n_wit=2
  system is "logically redundant for this claim structure (two orthogonal
  W1-solutions are neither necessary nor sufficient for a fourth ONB)", and
  that it "must not be used in new claims", remaining in the legacy file for
  archival reproduction only. This comment post-dates the 2026-08-13/14
  artifacts and applies to them retroactively, so the skipped 119210-path
  result above should be read as historical rather than as a live blocker.
  The live formulation is the n_wit=1 witness `w1_witness_certificate`
  (`src/Certification.jl:42–72`), which is already tractable at `mv=252`. Tag:
  `exact symbolic` (repository policy text).

## 5. Important methodological boundary

The repository contains a documented normalization correction in
`src/Certification.jl`: unit-normalized third-basis inputs must use RHS 1,
while unnormalized Hadamard-style inputs use RHS 6. Earlier calculations may
therefore not be comparable to later calculations until rerun with the same
normalization convention.

The current evidence does **not** prove:

* no fourth MUB on the full Diţă circle;
* no fourth MUB throughout \(K_6^{(3)}\);
* no fourth MUB in every order-six complex Hadamard family;
* \(N(6)=3\).

One further boundary is structural rather than numerical. The T3 emptiness
statements are certified for a *square section* of the overdetermined witness
system, and the inference to \(W_1=\varnothing\) passes through the rank of the
Gaussian matrix \(R\), which no artifact records. Combined with the four-CHM-
class reduction, the honest description of T3 is therefore *local certified
numerical over four independent CHM classes, with one unverified rank
assumption in the exclusion step* — not a certificate of non-extension over a
family.

## 6. Scientific next step

Before extending the research, establish one authoritative reproducibility
bundle: current commit, environment versions, exact command, stdout/stderr,
result hashes, normalization convention, pool completeness, and whether the
calculation is exact, certified numerical, or sampled.

The row-count discrepancy between the Python audit's 848-row/15-candidate view
and the Julia audit's 928-row/45-candidate view is now **explained but not
resolved**, and two earlier characterizations of it must be corrected.
`results/csv_reconciliation.txt` is a red herring with respect to this
discrepancy: it never mentions 848, 15, `backup848`, or the Python validator,
and it never opens `results/special_loci_search_backup848.csv`. It addresses only
the unrelated question of the \(\mathcal{S}_{588}\) row target, so any
description of it as *partially explaining* the 588-to-928 gap is incorrect on
its own terms and should not be repeated. That gap is not a pipeline
discrepancy either: \(588\) is a hard-coded design target written as a bare
integer literal at `scripts/julia/csv_reconciliation.jl:137–139`, with no flag
or configuration able to change it, and the script's own text — "negative =
more rows than target" — fixes the direction. So \(588\) against \(928\) is a
target-versus-achieved overshoot of 340 rows, and the artifact's "GAP
EXPLANATION" section is sign-incoherent: all four of its listed reasons describe
rows that were *not* produced, each of which would push the count down, while
the stated gap is negative, that is, more rows than target. Two of the four are
further contradicted by the artifact's own data.

The real 848-versus-928 discrepancy does have a complete arithmetic
explanation. `scripts/python/validate_third_mub_table.py:10` reads
`results/special_loci_search_backup848.csv`, a superseded snapshot some four
and three-quarter hours older than the authoritative
`results/special_loci_search.csv`, and its candidate predicate excludes
second-level `ref_ref` rows by a hard-coded default (line 20 and lines 31–32)
where the Julia validator includes them. Cross-applying both predicates to
both files closes the chain exactly: 15 is the Python predicate on the 848-row
file, 20 is the Python predicate on `special_loci_search.csv`, 45 is the Julia
predicate on `special_loci_search.csv`, and 64 is the Julia predicate on the
848-row file, which holds more clique-6 rows, not fewer. Tag: `exact symbolic`
for the source reading, `sampled/numerical only` for the recount. The
explanation is nonetheless **not a resolution**: the Python validator still
points at the stale file, so the one-line change to its input path has not been
made and the script has not been re-run. Until that patch lands, the
discrepancy is explained only. Nothing here certifies that
`special_loci_search.csv` is scientifically correct, merely that it is the
designated authoritative artifact; at least four `ref_degen_*` rows in it carry
stale `max_clique` labels from a pre-argmax-fix run.

The guard that previously stood therefore remains in force: aggregate counts
from these two views should not be used as a single scientific dataset. This is
reinforced rather than lifted, since the two files are structurally different
runs and not one run at two lengths — `search.csv` carries 108 `degen_*` rows
and 428 `ref_ref` rows, where the backup carries none and 580 respectively.

One scientific question inside this material is deliberately left open. Whether
`ref_ref`, the second-level refinements, should be counted or excluded is a
question of search-set semantics and not of reproducibility, and the two
readings give 45 against 20 on the authoritative file. The paper and
`results/project_finish_audit.txt:11`, which reports `S588 primary clique-6
(exclude ref_ref) count=45`, are internally inconsistent on precisely this
point. This document does not decide between them; the choice needs a human
decision, and the 45 should in any case be re-derived from a post-fix run.

Full derivations, the ten-point divergence table between the two validators,
and the phantom-function and double-counting findings are in
`research/reports/csv_reconciliation_analysis_2026-09-26.md`; the numerical
claims behind §2 and §4 are re-derived in
`research/reports/verification_2026-09-26.md`.
