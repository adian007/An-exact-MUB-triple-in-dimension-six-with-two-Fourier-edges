# Phase 1: Dita \(X_6\) identity and certified-clique audit

Date: 2026-09-29

Scope: Phase 1.0 and the available portions of Phase 1.1–1.2. The requested
all-40-clique family-fit sweep could not be completed because the repository
does not contain the clique matrices for four of the seven certificate
parameters. No Phase 2 work was performed.

## Files and provenance

**Created:** this report only.

**Modified:** none.

**Untouched result artifacts (SHA-256, checked before and after this report
was created):**

| File | SHA-256 |
|---|---|
| `results/campaigns/i3_singular_locus/third_mub_cliques.json` | `EA8176FCDB853A15A4B1D710B8954E5A59CB59C55C7270165A55A53651C2A2A5` |
| `results/campaigns/i3_singular_locus/root_certification.json` | `0DC30C8ED615BC323B77975F143D7A7482E59CFA6057C3C499AEF871AC0C8E34` |
| `results/certify_nwit1_all_cliques_four_classes.txt` | `04EC7E96BB8D21D9C19334F2D63DCB6C505B3B986A0646430D6BE581B4376C65` |
| `results/certify_nwit1_all_cliques_lambdapi.txt` | `DDEAACAB63630EED96EF020B6F9644D5CA38F9CCBEF243865BF77985AF3B6831` |
| `results/verify_nwit1_referee_checks.txt` | `12CF22B46DA72340EC025C7FCFE77D7650A2DE4652D33E52720BAEFEAAB53A2C` |
| `results/verify_nwit1_followup.txt` | `D4F17678265C4AEA4E107280AF7386751723AFE932128A4A6D4C4493A174E3B6` |

The worktree was clean before this report was added; `HEAD` was
`5330de2 Add K3 family-fit audit and pi/3 addendum`.

## Phase 1.0: closed form for the Dita \(X_6\) parameters

**Method and precision.** Direct inspection of double-precision \(X_6\) fits
at the three required parameters and ten additional values, followed by an
exact SymPy identity check over Laurent polynomials in a formal variable
\(t\). The 13-point numerical scan used \(\lambda=0.4,\pi/3,2\pi/3\) and
\(0.1,0.2,0.6,0.8,1.2,1.5,2.4,2.8,4.0,5.5\). Each numerical fit returned
FIT; matrix and constraint residuals were below \(6\times10^{-16}\).
The sequential 13-point run's peak working set was 83,402,752 bytes.

The fits at \(\lambda=1.2\) and \(1.5\) shared the chart
\[
p=(0,2,1,5,4,3),\qquad q=(0,5,1,2,3,4)
\]
(zero-based row and column orders) and exposed the branch
\[
t=e^{i\lambda/3},\quad z=t^3=e^{i\lambda},\qquad
\beta=t,\quad \gamma=-t^{-1},\quad
\epsilon=i\,t,\quad \phi=i\,t^{-1}.
\]
The other sampled fits sometimes selected different equivalent parameter
charts; the fixed-chart identity below, rather than continuity of those
chart choices, establishes the formula.

**Exact symbolic check.** Let \(H_D(z)\) be the unnormalized Dita matrix in
`dita_hadamard`'s convention, with \(\bar z=z^{-1}\), and let
\(C=H_D(z)[p,q]\). Dephase \(C\) by its first row and column:
\[
C^{\rm dep}_{rc}=\frac{C_{rc}C_{00}}{C_{r0}C_{0c}}.
\]
Substitution of the four parameter functions above into the repository's
\(X_6\) matrix gives
\[
C^{\rm dep}=X_6(\beta,\gamma,\epsilon,\phi)
\]
entry-by-entry. SymPy simplified every one of the 36 Laurent-polynomial
differences to zero in \(\mathbb Q(i,t,t^{-1})\). The source constraint (3),
\[
\beta\gamma\epsilon^2+\beta\gamma\phi+\beta^2\epsilon\phi+
\gamma\epsilon\phi+\beta\gamma^2\epsilon\phi+
\beta\gamma\epsilon\phi^2,
\]
also simplified identically to zero. Since \(t\ne0\), this is an exact
identity for generic \(z=t^3\), and for every real \(\lambda\) with the
specified branch \(t=e^{i\lambda/3}\). The common \(1/\sqrt6\) normalization
cancels in dephasing. Thus the Dita circle has an exact projective
permutation-equivalence to this \(X_6\) branch. Tier: **EXACT** (symbolic
identity under the stated parameterization and chart).

This is a statement about the displayed Dita matrix's \(X_6\) form. It does
not establish anything about completeness of the pools or additional MUBs.

## Phase 1.1: seven certificate parameters and CHM classes

The two native all-clique records list these parameters and recovered-pool
clique counts (10 + 4 + 4 + 4 + 10 + 4 + 4 = 40). Counts are recovered-pool
enumeration results, not completeness certificates.

| \(\lambda\) | Pool size | Recovered 6-cliques | CHM class from recorded equivalence checks | Native witness-log provenance | Tier |
|---:|---:|---:|---|---|---|
| \(0\) | 120 | 10 | class \(\{0,\pi\}\) | 2026-08-13, pre-2026-09-17 | NUMERICAL |
| \(0.4\) | 72 | 4 | representative class \(\{0.4\}\) within this seven-point sample | 2026-08-13, pre-2026-09-17 | NUMERICAL |
| \(\pi/3\) | 72 | 4 | class \(\{\pi/3,4\pi/3\}\) | 2026-08-13, pre-2026-09-17 | NUMERICAL |
| \(2\pi/3\) | 72 | 4 | class \(\{2\pi/3,5\pi/3\}\) | 2026-08-13, pre-2026-09-17 | NUMERICAL |
| \(\pi\) | 120 | 10 | class \(\{0,\pi\}\) | 2026-08-14, pre-2026-09-17 | NUMERICAL |
| \(4\pi/3\) | 72 | 4 | class \(\{\pi/3,4\pi/3\}\) | 2026-08-14, pre-2026-09-17 | NUMERICAL |
| \(5\pi/3\) | 72 | 4 | class \(\{2\pi/3,5\pi/3\}\) | 2026-08-14, pre-2026-09-17 | NUMERICAL |

The equivalence residuals in the follow-up/referee records are approximately
\(4.6\times10^{-16}\)–\(2.2\times10^{-15}\) for the tested
\(\lambda\leftrightarrow\lambda+\pi\) pairs. The four representative classes
are pairwise inequivalent in the recorded checks. Tier: **NUMERICAL** for the
class-equivalence checks and recovered-pool counts.

The seven native witness logs are dated August 13–14 and therefore predate
the 2026-09-17 normalization correction in `src/Certification.jl`. A separate
2026-09-28 `third_mub_cliques.json` export provides normalized clique-vector
matrices only for \(0.4,\pi/3,2\pi/3\); the Phase 1.2 calculations below use
those post-correction exported matrices. This date distinction does not
retroactively repair or upgrade the older witness calculations.

The old native logs report `mv=252 tracked=252 cert=252 fullpass=0` and
`EMPTY_SQUARE_252` for every listed clique. Those figures certify the tracked
square-system paths, but do not by themselves certify emptiness of the full
overdetermined witness system: the logs do not record the needed rank
verification. Pool completeness is also not established.

## Phase 1.2: fits on the 12 stored clique matrices

**Sequential run.** Loaded the four exported cliques at each of
\(0.4,\pi/3,2\pi/3\), and tested \(M=B_3\) and
\(M=H_D^\dagger B_3\), one matrix at a time. For every matrix, the run checked
flatness, unitarity, the three-distinct-columns predicate on \(M\) and \(M^T\),
direct \(F(a,b)\), the separate \(F(a,b)\)-via-transpose call, and \(X_6\) in
both orientations. Peak working set: 85,942,272 bytes.

Each cell below counts the 8 matrices at that parameter (4 cliques times 2
transitions). The same fit residuals were classified at all three requested
tolerances; the fitter's best chart and numerical residuals are independent
of the acceptance threshold.

| \(\lambda\) | Matrices | Flatness/unitarity max | Direct \(F\) fits at \(10^{-6},10^{-8},10^{-10}\) | \(F\)-transpose fits at each tolerance | \(X_6\), either orientation, at each tolerance | Theorem predicate on \(M\) / \(M^T\) at each tolerance |
|---:|---:|---:|---:|---:|---:|---:|
| \(0.4\) | 8 | \(<2\times10^{-15}\) | 8 / 8 / 8 | 0 / 0 / 0 | 0 / 0 / 0 | 0 / 8 |
| \(\pi/3\) | 8 | \(<3\times10^{-15}\) | 8 / 8 / 8 | 0 / 0 / 0 | 0 / 0 / 0 | 0 / 8 |
| \(2\pi/3\) | 8 | \(<5\times10^{-15}\) | 8 / 8 / 8 | 0 / 0 / 0 | 0 / 0 / 0 | 0 / 8 |

Across the 24 matrices, direct-\(F\) fit scores
\(\max(\text{consistency},\text{fixed-entry residual},\text{total residual})\)
were below \(5\times10^{-15}\). The best \(F\)-transpose residuals were at
least \(0.415\); the best \(X_6\) matrix/constraint residuals were at least
\(0.193\). The three-distinct-columns predicate was false on every direct
matrix and true on every transpose, with transpose residuals below
\(9\times10^{-16}\). These numerical results agree with the prior
three-parameter family-fit audit at
`docs/results/k3_family_membership_2026-09-29.md`.

For these 12 stored matrices only, there is no matrix satisfying the
theorem predicate in an orientation while lacking an \(F\), \(F^T\), or
\(X_6\) fit; thus no such anomaly was observed. Tier: **NUMERICAL**.

## Gate 1 blocker and stop

The full requested Phase 1.2 cannot be executed from the checked-in data:

* `third_mub_cliques.json` contains four 6-by-6 clique matrices at each of
  only \(0.4,\pi/3,2\pi/3\) (12 matrices total).
* The August 13–14 native witness logs contain 40 clique index lists and
  certification counts, but no pool vectors, bases, or transition matrices.
* `root_certification.json` contains certified root approximations, not the
  missing third-basis clique matrices.

Consequently, the other four \(\lambda\) values contribute 28 missing clique
bases, and the 12 available matrices cannot be represented as the requested
40-clique sweep. The existing native pool builder could regenerate data, but
that would be a new computation rather than loading the requested artifacts;
it has not been run here. **Stop at Gate 1.** No Phase 2 work was performed.

### What is not established

* The family fits do not cover all 40 certificate cliques or all seven
  \(\lambda\) values.
* The recovered pools are not certified complete.
* The August `EMPTY_SQUARE_252` records do not establish full witness-system
  emptiness without the missing rank step.
* The exact \(X_6\) identity does not establish pool completeness,
  non-extendability, or any conclusion below the fold.

### Commands and exit status

| Command/action | Exit |
|---|---:|
| Sequential 13-point \(X_6\) scan with peak working-set logging | 0 |
| Exact SymPy Laurent-polynomial comparison and constraint simplification | 0 |
| Sequential 24-matrix fit/predicate run on stored cliques | 0 |
| `git status --short` before report creation | 0 (clean) |
| SHA-256 checks of the six untouched artifacts before report creation | 0 |

### Contradictions and proposed next action

The supplied state says there are 40 certified recovered-pool clique records
at seven parameters; that is consistent with the two August text logs'
counts. It is not sufficient to perform the requested fits: only 12 clique
matrices are stored. The state that the three selected parameters have
four cliques each is confirmed by the September 28 matrix export. There is
no evidence here that the other 28 matrices are available in the referenced
artifacts.

**Proposed next action (not performed):** approve or provide fresh,
post-2026-09-17 pool-vector exports for \(\lambda=0,\pi,4\pi/3,5\pi/3\).
Then rerun Phase 1.2 for all seven parameters and all 40 clique matrices,
retaining new dated artifacts. Do not start Phase 2 until that Phase 1 gate
is resolved.
