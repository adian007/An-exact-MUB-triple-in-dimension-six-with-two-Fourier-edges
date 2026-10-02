# Plan: Next Steps Toward Determining Whether a 4th MUB Exists in Dimension 6

**Context**: This plan builds directly on the repository at `D:\MUBs in 6-dimension` and the two newly added packages (`b3_a4_i4_results/`, `provenance/fold_a4_final.zip`). It assumes the reader has the full context from the prior report.

**Problem restated**: Does there exist a 4th mutually unbiased basis extending {I, H, B3} for any order-6 complex Hadamard H? Equivalently: does W1(H,B3) ≠ ∅ for any B3? Or can we prove W1 = ∅ for all relevant B3 and all relevant H?

**Current evidence (from prior report)**:
- T1: PROVED (gauge structure)
- T2: NUMERICALLY_SUPPORTED (third-MUB on Diţă circle, 126+628 samples)
- T3: CERTIFIED_NUMERICAL (W1 empty at 7 λ-points, 40/40 cliques, 4 CHM classes, 4 caveats)
- T4: EXACT_ALGEBRAIC (W1 unit ideal at one D0-equivalent pair, re-verification of BW 2009)
- Large sampled searches: 0 fourth-MUB hits in 928 + 1863 + 628 points
- Fold/A4: 120→96→72 transition at λ\*, 24 singular roots = 12+12 A4 orbits (numerical)
- b3_a4_i4_results: 4 B3 cliques collapse to 1 A4 orbit at λ=0.4, π/3, 2π/3; fourth-vector probe found no zero residual (residuals 0.49–0.64)
- k3 audit: all 12 B3 transition matrices numerically fit Fourier family F/F^T; Jaming Theorem 1.4 bridge conditionally applicable but not exact

**The evidence is overwhelmingly negative but not a proof.** The plan below structures the path from "strong numerical evidence" to "exact or certified theorem."

---

## Phase 0: Freeze and certify the current state

**Goal**: Make the current evidence reproducible and hash-bound before extending it. The forensic audit (`docs/overview/forensic_ground_truth_report.md`) already flagged that the repo does not consistently embed producing commit hashes in result files.

### 0.1 Clean commit + manifest
- Commit the two new packages (`b3_a4_i4_results/`, `provenance/`) with a descriptive message.
- Generate a top-level `REPRODUCIBILITY_MANIFEST.json` recording:
  - Git commit hash
  - Julia version (1.12.6), WSL/Windows platform
  - Python version, NumPy/SciPy versions
  - HomotopyContinuation version (2.22.1)
  - SHA-256 of every result file cited in the paper/status docs
  - Exact command lines for each run
- Cross-reference with existing `results/benchmarks/run_manifest.json` and `src/Provenance.jl` patterns.

### 0.2 Resolve the 848/15 vs 928/45 discrepancy
- The Python validator `scripts/python/validate_third_mub_table.py` still reads `results/special_loci_search_backup848.csv` (stale, ~4.75 hours older than the authoritative `special_loci_search.csv`).
- One-line fix: update the input path in `validate_third_mub_table.py:10` to point at `special_loci_search.csv`.
- Re-run. Expect 20 candidates with the Python predicate (already known from the cross-application arithmetic).
- Decide: which predicate is correct — include `ref_ref` (45) or exclude (20)? The paper and `results/project_finish_audit.txt:11` are internally inconsistent on this point. This needs a human decision, then a consistent re-derivation.
- Tag the resolution in the claim ledger.

### 0.3 Re-verify T3/T4 post-normalization-fix
- The 2026-09-17 normalization fix (RHS=1 vs RHS=6 auto-detection in `src/Certification.jl`) was a substantive correction.
- Re-run `certify_nwit1_all_cliques_four_classes.jl` and `certify_nwit1_all_cliques_lambdapi.jl` from clean state, recording stdout/stderr, commit hash, and seed beside each result file.
- Re-run the Macaulay2 T4 Gröbner computation via `scripts/docker/run_m2.ps1` and verify the unit-ideal result, including the missing `results/track_c_elimination/export_w1_exact.txt` path inconsistency noted in the forensic audit.

**Exit criteria**: A frozen, hash-bound state where every number in `docs/scientific_status.md` and `docs/results/final_honest_status.md` traces to a specific artifact with a specific commit.

---

## Phase 1: Resolve the fold discrepancy and certify the singular locus

**Goal**: Understand whether the 24 singular roots at λ\* are genuine singular MU-vector roots, and whether the pool count at λ\* is 72 or 96. This matters because the b3_a4_i4_results work uses the 72-vector pool at non-critical λ, but the fold structure informs the global picture.

### 1.1 Diagnose the 72 vs 96 discrepancy
- The Phase 2 fold-sweep (`docs/results/fold_sweep_phase2_2026-09-29.md`) found Julia's `MubSearch.generate_candidate_pool_fresh` returned 72 at λ\*, while the Python 1200-start search (in the provenance archive) found 96.
- The 24-difference matches the singular candidates exactly. The question: are the 24 singular roots legitimate physical MU vectors that the Julia pool builder misses (because it uses a non-singular solver), or are they artifacts of the phase-coordinate least-squares approach?
- Action: Re-run the Julia pool builder at λ\* with explicit diagnostic: track ALL 252 HomotopyContinuation paths, report how many diverge/fail, and check whether the 24 missing vectors correspond to the singular locus.
- Action: Re-run the Python 1200-start search from the provenance archive scripts (`scripts/refine_singular_augmented.py`, `scripts/a4_orbit_refined.py`) to confirm the 96 count is reproducible.

### 1.2 Certify or refute the 24 singular roots
- The provenance archive's `FINAL_RESULTS.md` states: "These are numerical singular-root results, not interval-certified singularity proofs." This is accurate and must stay accurate.
- The 24 roots satisfy the augmented system F=0, Ju=0, ‖u‖²=1 to numerical precision. The next step is a genuine singular-root certification.
- **Option A — Deflation + interval certification**: Construct a deflated system that separates the singular roots, then use interval/Krawczyk methods to certify. This requires custom implementation; HomotopyContinuation.jl `certify()` is designed for nonsingular isolated solutions.
- **Option B — Algebraic verification at λ\***: If λ\* has an exact algebraic description (it appears to be a specific algebraic number), work over a number field and attempt an exact singularity analysis. This is harder but would upgrade the result from numerical to exact.
- **Option C — Accept as numerical evidence**: If A and B are infeasible, document the 24 roots as high-quality numerical evidence with explicit uncertainty, and do not promote them to certified status.

### 1.3 A4 orbit certification for the 24 singular roots
- The provenance archive reports 24 = 12 + 12 under the G2, G3 generators, with max matching error 8.5×10⁻¹⁰.
- Action: Verify the A4 closure is complete (not just that the 12 generated elements permute the 24 roots, but that there are no additional monomial automorphisms). This requires checking the full monomial automorphism group at λ\*, not just the fixed-z subgroup.
- The fixed-z group of order 12 was verified at a generic numerical sample in the `mub_complete_computational_package`. At the singular λ\*, the stabilizer could be larger.

**Exit criteria**: A documented, honest status for the λ\* singular locus: either (a) certified singular roots with A4 orbit decomposition, or (b) numerically resolved with explicit uncertainty, or (c) identified as a dead end with reasoning.

---

## Phase 2: Exact W1 elimination at one non-D0 algebraic Diţă point (the smallest meaningful theorem)

**Goal**: This is the single highest-value next result identified in `docs/research/research_gap_report.md` and `docs/overview/open_problems.md`. Prove that at one algebraic Diţă/Karlsson point NOT CHM-equivalent to D0, the W1 witness ideal is the unit ideal over a named number field.

### 2.1 Choose the target point
- **Preferred**: λ = 0 or λ = π/3 on the Diţă slice (θ = arccos(1/√3), φ = π/4).
- λ = 0: has a 120-vector pool (the larger pool, 10 cliques), and is CHM-inequivalent to D0. The b3_a4_i4_results work has pool_0_4.npz but not pool at λ=0 — that's a separate run.
- λ = π/3: has 72 vectors, 4 cliques, and the b3_a4_i4_results work has pool_pi_over_3.npz with the cliques already enumerated. The k3 audit already numerically fit these B3 matrices to the Fourier family.
- **λ = π/3 is the more tool-ready target** because the pool and cliques already exist in b3_a4_i4_results/.

### 2.2 Reconstruct an exact B3 at the target point
- The b3_a4_i4_results B3 vectors are numerical (float64 from the pool solver). For an exact certificate, you need exact B3 entries.
- **Route 1 — LLL reconstruction**: The repo already attempted this (`results/reconstruct_b3_algebraic.meta.txt`: Nemo.lll found 6/36 trivial row-1 only, bound ≤ 8). This failed. Try with higher precision, a larger basis, or a different lattice construction.
- **Route 2 — Algebraic Ansatz**: If the B3 vectors at λ=π/3 are exactly Fourier-family matrices (as the k3 audit numerically suggests with residuals ~10⁻¹⁵), then the exact B3 entries are algebraic numbers in a known field. The Fourier family F(a,b) has entries that are roots of unity and simple algebraic numbers. Fit the exact parameters (a,b) and reconstruct B3 exactly.
- **Route 3 — High-precision certified B3**: Use 400-bit (or higher) interval arithmetic to compute B3 to sufficient precision that the W1 emptiness margin is provable. This is a CERTIFIED_NUMERICAL route, not EXACT_ALGEBRAIC, but would be a substantial upgrade from the current float64 B3.
- The `docs/research/certificates/certificate_design.md` explicitly warns about the epsilon-transfer issue: W1 emptiness for a numerical B3 is not W1 emptiness for the exact intended B3 unless the perturbation is bounded and the exclusion has a margin.

### 2.3 Set up the exact W1 system
- The exact W1 system is documented in `docs/research/algebraic/algebraic_obstruction_program.md` and `src/Certification.jl`.
- 10 variables (z₁,...,z₅, w₁,...,w₅ with z₀=w₀=1), 17 equations:
  - zⱼwⱼ = 1 (j=1..5)
  - |Σⱼ conj(H[j,k]) zⱼ|² = 6 (k=1..6) — with exact H entries
  - |Σⱼ conj(B3[j,k]) zⱼ|² = 1 (for unit-normalized B3) — with exact B3 entries
- The field: depends on the exact H and B3 entries. For λ=π/3 on the Diţă slice, H has entries in ℚ(ζ₂₄) or a similar cyclotomic field. B3 may enlarge the field.
- `verify_field_degree_w1.jl` already exists for the D0 case; adapt it for the new point.

### 2.4 Run the Gröbner basis computation
- Export the exact system to Macaulay2 (`symbolic_export/w1_D0_groebner.m2` is the template).
- Run via `scripts/docker/run_m2.ps1` (Docker Desktop required).
- Target: GB = {1}, dim = -1, ideal is unit ideal.
- Record: field, monomial order, number of generators, number of variables, computation time, and the exact certificate hash.
- If GB ≠ {1}: analyze the dimension. A positive-dimensional component would be a major discovery (it would contain a candidate 4th-MUB vector). A finite nonzero number of solutions would be a finite set of candidate vectors to check.

### 2.5 Independent verification
- Repeat the Gröbner computation in a different CAS (Mathematica via `scripts/wolfram/w1_D0_exact.wl` exists as a template; Wolfram Engine is available).
- Cross-check the field degree calculation independently.

**Exit criteria**: Either (a) an exact unit-ideal certificate at λ=π/3 (or λ=0) — a genuine new theorem, or (b) a documented non-certification with the reason ( GB ≠ {1}, field too large, computation infeasible) and the next-best result.

---

## Phase 3: The Fourier family bridge — from numerical fit to exact exclusion

**Goal**: The k3 audit (`docs/results/k3_family_membership_2026-09-29.md`) found that all 12 B3 transition matrices at λ=0.4, π/3, 2π/3 numerically fit the Jaming et al. Fourier family F (or F^T) with residuals ~10⁻¹⁵. Jaming Theorem 1.4 proves that for ANY Fourier family member F(a,b), the pair {I, F(a,b)} cannot be extended to a quartet. If the B3 matrices are EXACTLY Fourier family members, then {I, B3} cannot be extended to a quartet — which would mean {I, H, B3, B4} cannot exist for any B4, because {I, B3} is already a non-extendable pair.

This is potentially a much more powerful result than the pointwise W1 certificate, because it would apply to ALL Fourier-family B3s at those λ values, not just the recovered-pool ones.

### 3.1 Upgrade the numerical fit to exact membership
- The current fit uses `fit_fourier_family(M)` at tolerance 10⁻⁸. Residuals are ~10⁻¹⁵. This is strong numerical evidence but not exact membership.
- Action: For each of the 12 B3 matrices, attempt to prove exact membership in the Fourier family. This requires:
  - Identifying the exact (a,b) parameters as algebraic numbers
  - Verifying that the matrix entries match the Fourier family formula exactly
- If the B3 matrices at λ=π/3 are exactly Fourier, the parameters (a,b) should be in a computable number field. The k3 audit already found numerical (a,b) values like (.521141+.853470i, -.753539-.657403i) for clique 0. These look like they could be roots of unity or simple algebraic numbers.

### 3.2 Apply Jaming Theorem 1.4 as an exact exclusion
- If exact membership is established, Jaming Theorem 1.4 (pair-to-quartet nonextension for all Fourier parameters) applies exactly.
- This would prove: for the B3 matrices that are exactly Fourier, {I, B3} cannot be extended to a quartet. Therefore {I, H, B3} cannot be extended to a 4th MUB (because any such extension would contain {I, B3} as a sub-pair).
- Scope: this is a theorem about the specific B3 matrices that are exactly Fourier, not about all possible B3s at that λ. But combined with the A4 orbit result (all 4 cliques are A4-equivalent), it covers all 4 recovered cliques at that λ.

### 3.3 Check the orientation/placement issue
- The k3 audit found that the Jaming paper uses a transpose placement relative to the repo convention. C0c in the audit proves an extension symmetry lemma: (I,H) extends iff (I,H^†) extends iff (I,H^T) extends. This covers both orientations.
- Verify this lemma is correctly applied and that the exact membership proof respects the orientation.

**Exit criteria**: Either (a) exact Fourier membership proved for at least one B3 clique, enabling an exact Jaming-based exclusion, or (b) documented that the residuals, while small, do not constitute exact membership, and the bridge remains numerical-only.

---

## Phase 4: Full-circle and full-region coverage (the long-term theorem targets)

**Goal**: Extend from pointwise results to continuum results. These are the open problems #4 and #5 from `docs/overview/open_problems.md`.

### 4.1 Interval covering on the Diţă circle
- **Objective**: Prove W1 = ∅ for ALL λ on the Diţă circle (or at least on a substantial interval).
- **Method** (from `docs/research/intervals/interval_arithmetic_research.md`):
  1. Subdivide the λ-circle into intervals
  2. For each interval, bound the variation of H(λ) and B3(λ) coefficients
  3. Use interval Newton/Krawczyk exclusion or certified homotopy continuation from a box interior
  4. Retain a numerical margin separating every candidate endpoint from the full witness equations
  5. Recurse until every box is excluded or promoted to an exact/algebraic exceptional stratum
- **Hazards** (documented in the intervals doc):
  - Interval dependency near Dita degeneracies
  - Moving B3 is not a fixed coefficient matrix
  - Clique identities change when pool vectors collide
  - Parameter denominators create false boxes
- **Realistic first target**: Prove an interval exclusion on a small compact subarc of the Diţă circle AWAY from the four special constellation-change windows (near λ=0, π/2, π, 3π/2 where the pool size changes). This would test whether the local certified point method has a viable uniform upgrade.

### 4.2 Complete-pool certification
- The current T3 result is "every clique in the recovered pool." It does not certify that the recovered pool is complete (i.e., that there are no MU vectors missed by the homotopy solve).
- **Objective**: Certify that the MU pool at a given H is complete — that the homotopy solve found ALL solutions.
- HomotopyContinuation.jl tracks all paths from a start system with known number of solutions (the mixed volume). If all 252 paths are tracked and certified, and the start system is a valid upper bound, this is a completeness certificate for the polynomial system.
- The T3 artifacts report `mv=252 tracked=252 cert=252` — all paths tracked and certified. This IS a completeness certificate for the square system, but the step to the full 17-equation system requires the rank-10 assumption on the Gaussian matrix R (unverified).
- Action: Verify the rank-10 assumption explicitly, or find a way to certify W1 emptiness without it.

### 4.3 The Wuttig–Tindall classification
- Wuttig–Tindall v2 (2026, arXiv:2608.18053) claims a complete classification of order-6 CHMs: H₆ = G₆⁽⁴⁾ ∪ K₆⁽³⁾ ∪ T₆.
- If this classification is correct and independently verified, then K₆⁽³⁾ (the Karlsson family) is only one of three sectors. The repo's Karlsson-only results would still not cover G₆⁽⁴⁾ or T₆.
- Action: Track the Wuttig–Tindall preprint. If it gains acceptance and independent verification, update the repo's scope statements. If it's refuted, note that.
- The repo's `docs/research/literature/agent9/applicability_assessment.md` already assesses this.

**Exit criteria**: For Phase 4, there is no single exit. Each sub-phase produces a documented result (theorem, numerical evidence, or dead-end with reasoning).

---

## Phase 5: The positive case — what would a 4th MUB look like?

**Goal**: If a 4th MUB exists, what would we expect to see? This is the constructive counterpart to the exclusion searches. It's important to have a positive vision, because the negative evidence could be missing something.

### 5.1 What the searches have ruled out
- The 928 + 1863 + 628 point searches found 0 fourth-MUB hits. These are finite samples, not exhaustive.
- The T3 result rules out 4th MUBs at 7 λ-points for the recovered-pool B3s.
- The T4 result rules out 4th MUBs at one D0-equivalent pair exactly.
- The fold/A4 work characterizes the singular structure of the MU pool at λ\*.

### 5.2 Where a 4th MUB could hide
- **Outside K₆⁽³⁾**: If the Wuttig–Tindall classification is correct, a 4th MUB could involve a G₆⁽⁴⁾ or T₆ matrix as H, or a B3 that is not in K₆⁽³⁾.
- **At a non-sampled λ**: The Diţă circle is continuous. The 7 T3 points and the 628 dense samples don't cover every λ.
- **With a non-recovered B3**: The pool completeness caveat means there could be a B3 that the homotopy solve missed.
- **With a non-clique B3**: The repo searches for B3 as 6-cliques in the orthogonality graph. If a B3 exists that is not a 6-clique (e.g., because the vectors aren't all mutually orthogonal), it would be missed.

### 5.3 Constructive search strategies
- **Random CHM generation**: Generate random order-6 CHMs (not just from K₆⁽³⁾) and test for 4th-MUB extendibility. This is a shotgun approach but could hit a positive case.
- **Optimization-based search**: Formulate the 4th-MUB existence as an optimization problem (minimize the MU defect over H, B3, B4) and use global optimization to search for a zero. The b3_a4_i4_results work already does a local nonlinear least-squares probe for the 4th vector; a global optimizer with multiple starts could be more thorough.
- **Symmetry-constrained search**: Use the known symmetries (A4, gauge freedoms) to reduce the search space and make a global search feasible.

### 5.4 The positive evidence threshold
- If a 4th MUB is found: it must be verified exactly (or at certified numerical precision) to rule out a numerical artifact.
- The verification must include: CHM residual for H, ONB verification for B3 and B4, MU verification for all 6 pairs, and independence from gauge choices.
- A single positive case would settle N(6) ≥ 4, refuting Zauner's conjecture.

**Exit criteria**: Either (a) a verified 4th MUB construction (would be a major result), or (b) a documented exhaustive search over a well-defined domain with zero hits, or (c) a clear statement of where a 4th MUB could still hide given current evidence.

---

## Phase 6: Synthesis and publication readiness

**Goal**: Bring the results together into a coherent, honest, publication-ready form.

### 6.1 Claim ledger update
- Update `docs/results/final_honest_status.md` and `docs/scientific_status.md` with all new results, maintaining the six-tier vocabulary.
- Every new claim must carry: tier, evidence file, method, scope, and explicit limitations.
- No promotional language. If a result is numerical, say so. If it's local, say so. If it has an unverified assumption, say so.

### 6.2 Paper draft
- The existing LaTeX drafts in `docs/paper/` (main_theorems.tex, methods_audit.tex, proofs/*) should be updated to reflect the current evidence.
- The paper should clearly state what is proved, what is numerically supported, and what is open.
- The T3 and T4 scopes should be stated as local, not as family-wide.
- The fold/A4 results should be included as numerical evidence for the pool structure, not as a theorem.

### 6.3 Literature reconciliation
- The McNulty–Weigert 2025 correction (arXiv:2504.13067) invalidated a lemma used in later exclusion theorems. Any later exclusion claims based on that lemma must be reconciled.
- The Wuttig–Tindall v2 classification should be cited as a preprint claim, not as established fact.
- The Jaming et al. Theorem 1.4 bridge should be cited as a conditional numerical application, not as an exact exclusion (unless exact membership is proved).

### 6.4 Reproducibility package
- Every result cited in the paper should have a companion reproducibility entry: script, command, commit hash, seed, tolerance, and output hash.
- The `run_manifest.json` pattern from `scripts/julia/run_benchmarks.jl` should be extended to all major runs.

**Exit criteria**: A publication-ready preprint with honest claim tiers, full reproducibility metadata, and a clear open-problems section.

---

## Priority ordering and dependencies

| Priority | Phase | Depends on | Why first |
|---|---|---|---|
| **P0** | 0 (Freeze + certify) | Nothing | Without this, new results can't be trusted to be reproducible |
| **P1** | 2 (Exact W1 at non-D0 point) | 0 | Highest-value single theorem; uses existing b3_a4_i4_results pool at π/3 |
| **P1** | 3 (Fourier bridge) | 0, 2.2 (exact B3) | Potential for a stronger, family-level exclusion via Jaming Theorem 1.4 |
| **P2** | 1 (Fold discrepancy + singular certification) | 0 | Resolves an open discrepancy; informs the global picture |
| **P3** | 4 (Interval covering, complete-pool) | 0, 1, 2 | Long-term theorem targets; require the pointwise results as base cases |
| **P3** | 5 (Positive search) | 0 | Parallel track; a positive result would be transformative |
| **P4** | 6 (Synthesis + publication) | All of the above | Brings everything together |

**The critical path is P0 → P1 (Phase 2) → P3 (Phase 6).** The exact W1 at λ=π/3 (or λ=0) is the single most valuable next result, and it's achievable with existing infrastructure.

---

## Resource and feasibility notes

- **Julia/WSL**: All Julia computation runs in WSL2 Ubuntu with JULIA_DEPOT_PATH pointing to `/home/adian/mub-depot`. Windows Smart App Control blocks unsigned Julia DLLs natively. Documented in `docs/overview/reproduce.md`.
- **Macaulay2**: Requires Docker Desktop running. Invoked via `scripts/docker/run_m2.ps1`.
- **Wolfram**: Wolfram Engine scripts exist in `scripts/wolfram/` for exact verification.
- **Python**: Used for the k3 audit, fold computation, and auxiliary analysis. NumPy/SciPy stack.
- **Time estimates**:
  - Phase 0: 1–2 days (mostly documentation + one-line Python fix + re-runs)
  - Phase 2 (exact W1 at π/3): 1–4 weeks (B3 reconstruction is the uncertain part; Gröbner computation is hours to days once the system is set up)
  - Phase 3 (Fourier bridge): 1–2 weeks (parallel with Phase 2; depends on B3 reconstruction)
  - Phase 1 (fold): 1–2 weeks
  - Phase 4 (interval covering): 1–3 months (substantial implementation work)
  - Phase 5 (positive search): ongoing, parallel

---

## What would constitute success at each tier

| Tier | Success looks like |
|---|---|
| PROVED | Exact Gröbner unit ideal at λ=π/3 (or λ=0), with independent CAS verification |
| EXACT_ALGEBRAIC | Same as above; or exact Fourier membership proof for a B3 clique enabling Jaming exclusion |
| CERTIFIED_NUMERICAL | Interval-certified W1 emptiness on a λ-interval (not just points); certified singular roots at λ\* |
| NUMERICALLY_SUPPORTED | Higher-precision B3 reconstruction with certified margin; denser λ sampling with consistent negative results |
| HEURISTIC | Current state: numerical fits, local probes, sampled searches |
| OPEN | N(6)=3 globally; no-4th-MUB on all of K₆⁽³⁾; full Diţă circle; non-Karlsson families |

---

## The honest bottom line

The repository has done something genuinely rigorous: it has accumulated strong numerical evidence that no 4th MUB exists on the Diţă circle, with certified numerical certificates at 7 points and exact algebraic certificates at one D0-equivalent pair. The evidence is negative and consistent across multiple independent methods (homotopy continuation, interval certification, Gröbner bases, CHM equivalence checks, Fourier family fitting).

But the evidence is not a proof. The gaps are:
1. Pointwise, not continuum (T3 covers 7 points, not the full circle)
2. Recovered-pool, not complete-pool (T3 covers recovered cliques, not all possible B3s)
3. One unverified rank assumption in the W1 emptiness inference
4. Only one exact certificate (T4 at D0-equivalent, which is corollary-grade relative to published BW 2009)
5. Karlsson-only, not all CHM families (K₆⁽³⁾ is only the H₂-reducible sector)

The plan above addresses these gaps in priority order, with the exact W1 at a non-D0 point as the single highest-value next step. If that succeeds, it would be a genuine new theorem. If it fails (GB ≠ {1}), the failure mode itself is informative and should be analyzed.
