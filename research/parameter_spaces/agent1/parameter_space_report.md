# Agent 1 — Parameter-space and classification audit

**Campaign:** pre-computation literature/repository audit  
**Date:** 2026-09-17  
**Scope:** order-six complex Hadamard matrices (CHMs) and the meaning of
“extend dimension six into four parameters.” No CHM construction, numerical
sweep, homotopy, or elimination was run by this agent.

## Executive conclusion

There is a mathematically valid interpretation of “four parameters”: the
published, non-affine four-real-parameter CHM family
\(G_6^{(4)}\) constructed by Szöllősi. This is **not** a verified
four-parameter extension of Karlsson’s \(K_6^{(3)}\). The source calls the
former *generic* and the latter *degenerate*, but does not establish a
set-inclusion \(K_6^{(3)}\subset G_6^{(4)}\), and the repository contains no
such derivation. Accordingly, no “Karlsson plus one parameter” builder is
justified.

The recommended Campaign-3 object, subject to an independent transcription
and CHM validation, is Szöllősi’s \(G_6^{(4)}\) as its **own** four-real-
parameter sector. It must be labelled `G6_4`, not `Karlsson_4`.

The recent Wuttig–Tindall classification preprint gives a second, broader
meaning: four initial dephased-corner phases from which it claims every
order-six class is reconstructed algebraically (with stated exceptions). It
is a **coordinate/reconstruction statement**, not a four-parameter CHM
formula and not evidence that every physical moduli component is a smooth
four-manifold. It is a v2 arXiv preprint dated 2026-08-28; this audit did not
independently check its proof or any claimed formalization. Treat its global
classification as `LITERATURE_RESULT (unindependently reviewed preprint)` for
planning, never as an internal theorem.

## Definitions: dimensions that must not be conflated

For an unnormalised CHM \(H\in\mathbb C^{6\times6}\),
\[
  |H_{jk}|=1,\qquad HH^\dagger=6I_6.
\]

* **Raw phase/gauge dimension.** Standard Hadamard equivalence acts by left
  and right monomial unitaries (row/column phase multipliers and
  permutations). Its continuous row/column-phase part has \(2n-1=11\)
  effective real phases at \(n=6\): the simultaneous opposite global row/
  column phase is redundant. These are not inequivalent CHM parameters.
* **Dephased family dimension.** Fix the first row and first column to one;
  then quotient only the continuous phase gauge. A “three-/four-parameter
  family” below means real parameters of a dephased family of inequivalent
  matrices, subject still to discrete identifications and singular overlaps.
* **Defect.** The (dephased) defect comes from the real linearization of the
  phase/orthogonality equations. It bounds the dimension of a smooth local
  deformation orbit; it neither proves integrability nor gives a global
  parametrization. It must not be substituted for family dimension.
* **MUB variables.** Parameters of a third basis \(B_3\), a fourth-basis
  candidate, or phase coordinates for MU vectors live in a different
  incidence variety. Adding them produces a multi-parameter *search space*,
  not a four-parameter CHM family.

## What “four parameters” can validly mean

| Interpretation | Valid? | Precise formulation | Research status / warning |
|---|---|---|---|
| 1. A four-real-parameter CHM family | **Yes** | Szöllősi’s \(G_6^{(4)}\), whose entries are algebraic functions involving roots of sextics | Genuine published family; formulas are not yet transcribed or validated in this repository. |
| 1a. A four-parameter *extension of Karlsson* | **No justification found** | Would require an explicit map/family with a proved relation to \(K_6^{(3)}\) | Do not build or name such an ansatz. The literature describes K and G as separate sectors/families. |
| 2. Add a parameter to a chosen MUB/triple orbit | **Conditionally valid** | A proved map \(t\mapsto (H(t),B_3(t))\) satisfying all pairwise MUB identities for every \(t\) in a stated domain | Current repository has sampled/high-precision loci, not a universal exact orbit proof. Its parameter is not a CHM-classification parameter unless separately shown. |
| 3. Four parameters for a prospective fourth basis | **Valid coordinate choice, not a family claim** | Gauge-fixed variables for a candidate \(B_4\), constrained by its 36 unimodularity and 36 orthogonality/MU equations | Usually vastly underparameterized if only four reals; any such reduced ansatz needs an explicit completeness disclaimer. |
| 4. A four-real-dimensional incidence/search space | **Valid if defined** | e.g. \((p_1,p_2,p_3,t)\in D_K\times I\), with a separately defined candidate-vector/basis branch | Does not imply a four-parameter CHM family or coverage of all third bases. |
| 5. Four dephased corner phases | **Provisional literature coordinate chart** | Wuttig–Tindall’s claimed reconstruction from a dephased \(3\times3\) corner specified by four phases | Need exact transcription, branch accounting, and independent proof audit before using as an exhaustive domain. |

## Families/sectors and repository coverage

Karlsson’s Theorem 11 gives a three-real-parameter description of the
\(H_2\)-reducible sector. Its paper explicitly calls it a three-parameter
family and says it contains previously described one- and two-parameter
families as subfamilies. It does not claim to be all order-six CHMs.

Szöllősi constructs a distinct four-real-parameter family \(G_6^{(4)}\), with
entries described through algebraic roots of sextics. The paper conjectured
that \(G_6^{(4)}\), \(K_6^{(3)}\), and Tao’s isolated spectral matrix exhaust
order six. The 2026 Wuttig–Tindall preprint claims a complete classification
and an algebraic four-corner-phase reconstruction. The preprint is relevant
but has not been independently validated here.

Repository inspection shows:

* `src/Karlsson.jl` implements and validates only Karlsson’s original
  three-coordinate constructor \((\theta,\phi,\lambda)\). It records the
  Fourier seam and Möbius degeneracies; they are strata, not a fourth
  coordinate.
* The repository’s own `docs/LITERATURE_2026.md` and
  `research/literature/matrix_catalogue_source_audit.md` explicitly identify
  \(G_6^{(4)}\) as unimplemented and Karlsson plus Tao as insufficient scope.
* The Liang/Chen/Long–Qiu path is a non-operational stub; the local audit
  found special-matrix results, not a published Karlsson-style parametric
  family. This is not a proof that no other parametrization exists.
* Tao \(S_6\) is an isolated anchor and is not a point of a documented local
  four-parameter implementation.

Thus existing Karlsson computations cover neither \(G_6^{(4)}\) nor all
order-six CHM sectors. They also do not cover all MUB triples even within
Karlsson without a proof that every third basis has been quantified.

## Recommended domains and reductions

### Karlsson (existing, three-coordinate program)

Use Karlsson’s source-reduced parameter ranges when transcribing a theorem,
not merely the exploratory \([0,2\pi)^3\) coordinates. In code and all
reports, retain explicit periodic identifications, discrete row/column
permutations, conjugation/transpose only when their exact action has been
proved, and branch labels for Möbius denominators. The known \(\theta=0\)
seam and Diţă indeterminacy must be separately included; no open-cell result
extends across them without a specialization/limit proof.

### Szöllősi \(G_6^{(4)}\) (recommended next CHM family)

Do **not** infer a rectangular fundamental domain from the phrase
“four-parameter family.” Its formula is algebraic/multibranch. The safe
domain is: the exact parameter domain and branch conditions stated in the
primary paper, augmented by discriminant/denominator exclusions; each
excluded boundary/branch locus is an explicit stratum. Dephase every output,
then quotient only exact standard-equivalence certificates (monomial phases
and permutations). A future builder must emit parameters, branch choices,
minimal-polynomial/root-isolation data, and exact or certified CHM checks.

### Four-corner reconstruction (preprint-driven option)

Until independently audited, use it only as a hypothesis-generating chart:
four torus phases modulo exact dephasing/permutation symmetries, with every
quadratic/cubic root and every vanishing leading coefficient retained as a
separate branch. It is not yet an authorized replacement for \(G_6^{(4)}\)
or a compact fundamental domain for a universal theorem.

## Recommended next action

Agent 3 may proceed only under interpretation 1: formulate
`G6_4_transcription_validation_plan`, **not** an implementation, until the
primary formulas and their branch domain are independently checked. A useful
four-dimensional project meanwhile is the explicitly named incidence domain
\(D_{G}\) (four G-family coordinates) crossed with separately quantified MU
variables. It must never be described as “Karlsson extended to four
parameters.”

## Sources (direct, checked)

1. B. R. Karlsson, *Three-parameter complex Hadamard matrices of order 6*,
   LAA 434 (2011), 247–258; [arXiv:1003.4177](https://arxiv.org/abs/1003.4177).
   Abstract and source provenance directly support the three-parameter
   construction; the repository uses its Theorem 11 for the \(H_2\)-reducible
   sector.
2. F. Szöllősi, *Complex Hadamard matrices of order 6: a four-parameter
   family*, JLMS 85 (2012), 616–632,
   [DOI](https://doi.org/10.1112/jlms/jdr052),
   [arXiv:1008.0632](https://arxiv.org/abs/1008.0632). The abstract explicitly
   says it constructs a previously unknown four-parameter family and describes
   the sextic algebraic dependence.
3. P. Diţă, *Four-parameter families of complex Hadamard matrices of order
   six*, [arXiv:1207.2593](https://arxiv.org/abs/1207.2593). A separate claimed
   construction route. It was not formula-audited here; do not assume it is a
   disjoint sector or an extension of G/K without an equivalence analysis.
4. M. Cárdenes Wuttig and J. Tindall, *A Complete Classification of Complex
   Hadamard Matrices of Order Six*, v2 (2026-08-28),
   [arXiv:2608.18053](https://arxiv.org/abs/2608.18053). Recent preprint;
   its completeness claim and four-corner reconstruction are recorded but
   not independently verified.
5. W. Tadej and K. Życzkowski, *Defect of a unitary matrix*,
   [arXiv:math/0702510](https://arxiv.org/abs/math/0702510). Supports the
   distinction between defect and an actual smooth orbit dimension.
6. D. McNulty and S. Weigert, *Mutually Unbiased Bases in Composite
   Dimensions — A Review*, [arXiv:2410.23997](https://arxiv.org/abs/2410.23997).
   Used only as a secondary catalogue pointer; the repository flags a local
   Karlsson-transcription discrepancy, so Karlsson’s primary paper controls
   formulas.

## Limitations

This is a source/repository audit, not a proof of classification or of an MUB
statement. The local exact seven-point results, sampled third-MUB loci, and
any asserted global CHM classification were not recomputed here. Nothing in
this report establishes a fourth MUB, its nonexistence, nonextendability on a
continuous orbit, or \(N(6)=3\).
