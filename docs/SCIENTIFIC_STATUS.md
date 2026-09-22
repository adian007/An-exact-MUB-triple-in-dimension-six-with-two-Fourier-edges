# Scientific calculation status

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
| T3 | No fourth-MUB witness at seven selected Diţă points | Square-system homotopy solving, certification, and full residual checks over recovered six-cliques | Local certified numerical result | 40/40 recovered-pool six-cliques at 7 \(\lambda\) values |
| T4 | Exact elimination of the one-vector fourth-MUB witness at a \(D_0\)-equivalent point | Macaulay2 Gröbner basis over \(\mathbb Q(\zeta_{24},\sqrt5)\) | Exact local algebra | One \((D_{\rm bc},F_D)\) pair; \(\lambda\in\{\pi/2,3\pi/2\}\) |
| C-search | Finite scans of special loci and the Diţă circle | Floating/high-precision numerical pool and clique computations | Sampled evidence only | 928-row search, 1863-row S\* run, 628-point dense circle |
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
  Together they cover the seven selected points and 40/40 cliques.

### Numerical or sampled

* `results/dita_lambda_fourth_dense.csv`: 628 samples, no fourth result;
  not a theorem for the continuous circle.
* `results/special_loci_search.csv`: 928 rows, 899 pool-complete.
* `results/special_loci_degen1845.csv`: 1863 rows, 1852 pool-complete,
  21 six-clique loci, 0 fourth-MUB detections.
* `scripts/python/audit_karlsson_variants.py`: 348/349 for the original
  transcription and 0/349 for the literal review transcription.

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

## 6. Scientific next step

Before extending the research, establish one authoritative reproducibility
bundle: current commit, environment versions, exact command, stdout/stderr,
result hashes, normalization convention, pool completeness, and whether the
calculation is exact, certified numerical, or sampled. First reconcile the
Python audit's 848-row/15-candidate view with the Julia audit's
928-row/45-candidate view. Until that discrepancy is explained, aggregate
counts should not be used as a single scientific dataset.
