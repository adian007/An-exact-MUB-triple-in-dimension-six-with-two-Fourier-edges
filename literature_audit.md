# Literature Audit: MUBs in Dimension Six

Audit date: 2026-09-17. Primary-source URLs were inspected through the arXiv abstract/HTML versions listed below. This report does not claim a global nonexistence theorem.

## Executive finding

The MUB problem is open. In dimension six, three MUBs exist and the general upper bound is seven. Zauner's conjecture is the stronger statement `N(6)=3`; neither `N(6)=3` nor global fourth-MUB nonexistence is proved. Wuttig--Tindall v2 (2026) now claims and formally audits a complete order-six CHM classification, but that preprint has not been independently reproduced in this repository. The precise local status is mixed:

- Grassl proved maximality for the Heisenberg/Weyl pair in dimension six: 48 candidate vectors, 16 third bases, and no vector extending any of those triples.
- Brierley--Weigert gave an exact D0 computation: 120 vectors MU to `{I,D0}`, 10 third bases, and no pair of those third bases is MU. Hence no four MUBs containing `D0`.
- Brierley--Weigert also surveyed known order-six families. Some cases are exact, while non-affine-family results are explicitly approximate and not universal.
- Karlsson's `K_6^(3)` is a complete parametrization of the `H_2`-reducible sector, not by itself a classification of all order-six CHMs. Wuttig--Tindall v2 claims the exhaustive decomposition `H_6=G_6^(4) union K_6^(3) union T_6` in Corollary 25. The repository has not independently audited that classification, and its computations remain Karlsson-only.
- This repository has independently reproduced the F6 and D0-equivalent benchmark counts, and has local certified-numerical fourth-vector obstructions at seven non-D0-equivalent Karlsson-Dita parameter points. Those are finite/local results, not a theorem on the Dita circle, all Karlsson matrices, or `N(6)`.

## A. Global problem

### A1. What is open?

There are several related questions, not one question:

1. Does a fourth MUB exist in `C^6`?
2. Does a complete set of seven MUBs exist in `C^6`?
3. What is the exact value of `N(6)`?
4. Is the 2026 Wuttig--Tindall classification correct and independently reproducible?
5. Which MUB triples are equivalent and which are extendible?

A fourth MUB would imply at least four MUBs, while seven is the complete-set question. The implication is one-way: proving no fourth MUB would prove `N(6)=3`, but evidence about selected CHM families does not address all MUB triples.

Bengtsson et al. state that for `N=6` the known bounds are at least 3 and at most 7, and that they found no MUB quartets. Their conclusion is not a global theorem. Grassl states that the maximal size `N(6)` was still open and records Zauner's conjecture `N(6)=3`.

### A2. Zauner's conjecture

In the dimension-six MUB context used by Grassl, Zauner's conjecture is `N(6)=3`: no set of four MUBs exists in `C^6`. It is not the SIC-POVM conjecture that also appears in Zauner's thesis and later literature. The two Zauner conjectures must not be conflated.

### A3. Rigorous bounds

The universal bounds verified in the primary sources are:

- lower bound: `N(6) >= 3`, by the general three-basis construction;
- upper bound: `N(6) <= 7`, by the standard dimension bound;
- no cited primary source inspected here proves `N(6)=3`.

## B. CHMs and MUBs

For an unnormalised order-six CHM `H`, `H H^dagger = 6 I` and `|H_ij|=1`. The columns of `H/sqrt(6)` form a basis unbiased to the computational basis. If `H1` and `H2` represent two non-computational bases, they are MU exactly when every entry of `H1^dagger H2` has modulus `sqrt(6)` in the unnormalised convention, equivalently modulus `1/sqrt(6)` after normalization.

The exact CHM equivalence used by Bengtsson et al. and Karlsson is

`H2 = D2 P2 H1 P1 D1`,

where `D1,D2` are diagonal unitary matrices and `P1,P2` are permutation matrices. Thus row phases, column phases, row permutations, and column permutations are allowed. For unordered MUB pairs, the pair `{I,H}` is also equivalent to `{I,H^dagger}`; transpose is not an independent generic equivalence operation, though particular families may relate to their transpose.

## C. Known order-six CHMs

The 2007 primary survey lists the Fourier family `F(x1,x2)`, its transpose, the circulant Björck matrix, the Diță family `D(x)`, the Hermitian family `B(theta)`, and Tao's isolated spectral matrix `S`. It explicitly says the complete classification is open at that time. Wuttig--Tindall v2 is a later preprint claiming that the classification has now been completed.

Karlsson's 2011 Theorem 11 proves that every `H_2`-reducible order-six CHM is equivalent to a member of the three-parameter family `K_6^(3)`. The theorem does not say every order-six CHM is `H_2`-reducible. Karlsson's Section 7 says the full characterization is only partial and that evidence points to an additional parameter. Maxwell--Brierley likewise state that the complete classification is open and that Karlsson excludes the isolated spectral matrix.

The source-aware current statement is:

- `K_6^(3)` covers the `H_2`-reducible subclass.
- It contains the previously known explicit one- and two-parameter families in that subclass.
- Karlsson's theorem alone does not cover all order-six CHMs.
- Wuttig--Tindall v2 claims that the remaining classes form the `G_6^(4)` and Tao sectors.
- The repository's Karlsson-only results still cannot be promoted to a global MUB theorem, even if that new classification is accepted.

## D. Exact source results

### D1. Grassl

Grassl's Theorem 2 gives exactly 48 normalized vectors unbiased to the two Heisenberg eigenbases in dimension six. They form 16 orthonormal third bases, and every resulting triple is maximal: no vector is unbiased to all three bases. The scope is the Weyl--Heisenberg pair and its Jacobi-group equivalents, not arbitrary CHM pairs.

### D2. Bengtsson et al.

The 2007 survey gives:

- Eq. (11): the Diță family `D(x)` and `D0=D(0)`;
- Eqs. (20)--(21): the MUB/CHM conventions and equivalence operations;
- Eqs. (61)--(62): the twisted Fourier form `F_D`;
- Eqs. (68), (74): four explicit third bases for a Fourier-family pair and the associated circulant blocks;
- Eqs. (78)--(80): the explicit triple `{I,F_D,D_bc}`, with `D_bc` equivalent to `D0`.

The paper reports only triplets in its search and explicitly says that complete classification of order-six CHMs is unavailable.

### D3. Brierley--Weigert 2009

Their exact D0 computation is in Section 4.1 and Table 1. They report 120 vectors MU to `{I,D0}`, 60 of which form 10 bases, with no two of those bases MU to one another. They state: sets of four MU bases containing `D0` do not exist. Their phase set is Eq. (20), with `tan(alpha)=2`.

Their general survey is not a universal proof: Section 4 distinguishes rigorous special/affine cases from approximate non-affine computations. The conclusion says the calculations cover nearly 6,000 sampled CHMs, not all CHMs.

### D4. Maxwell--Brierley 2015

Theorem 4.1 proves the additional Karlsson-family Fourier constraint `g_j(rho)=0` for permutations of `(1,1,1,-1,-1,-1)`. Section 5 reports that the resulting linear programs did not yield a contradiction. This is a structural partial result, not a proof that Karlsson matrices cannot occur in a complete MUB set.

### D5. McNulty--Weigert review and 2025 comment

The 2026 review is a secondary source useful for navigation, not a substitute for the original papers. It identifies Karlsson as the known three-parameter family and distinguishes it from the four-parameter `G_6^(4)` family. The 2025 McNulty--Weigert comment says a lemma used in later exclusion theorems contains an error and that three later theorems are invalidated. Those later exclusion claims should not be used as foundations without separate verification.

## E. Repository audit

### E1. Correctly supported

- The per-H pool equations and the CHM/MUB convention are consistent with Brierley--Weigert Section 3, Eqs. (6)--(8).
- `D0`, `D_bc`, and `F_D` are reconstructed as separate source objects, not defined by equality to the repository's candidate.
- F6 benchmark: 48 physical MU vectors and 16 third bases, with no fourth basis, agrees with Grassl/Brierley--Weigert.
- D0-equivalent benchmark: 120 vectors, 10 third bases, and no MU pair among third bases agrees with Brierley--Weigert.
- The exact `W1(D_bc,F_D)` unit-ideal computation is executed over `Q(zeta_24,sqrt(5))`, but its mathematical no-fourth content is corollary-grade relative to the stronger Brierley--Weigert D0 classification.
- The Karlsson builder matches the original Karlsson A-block and the repository's regression tests reject the literal McNulty--Weigert review transcription as a valid CHM block. This is an independently checked transcription issue, not evidence that the review's mathematical family is wrong.

### E2. Local repository findings

- Seven Dita-Karlsson parameter points have all recovered-pool third-basis cliques tested by a certified numerical `W1` procedure: 40/40 empty. This is a finite certified-numerical statement under the documented solver, residual, and completeness conditions.
- The dense Dita-circle and larger special-locus datasets are numerical/certified-numerical evidence, not a universal interval or exact-parametric certificate.
- The 512-point interior grid is only a sampled grid and deliberately misses named special loci; it is not a family theorem.
- `Dita lambda=pi/2,3pi/2` is numerically CHM-equivalent to `D0` under the repository's full row+column check; the seven T3 points are not matched to `D0` by the narrower column-only dephased comparison used in E0. A full monomial-equivalence check for every T3 point is still required before calling them definitively non-D0-equivalent.

## F. Terminology and scope rules

Use `no single fourth vector` for an exact or certified-empty `W1` statement. Use `no fourth basis` only after the candidate-vector system is complete and the logical implication is stated. Use `no fourth MUB globally` only for a proof covering every possible third basis and every possible initial pair, which is not present here.

Do not call a finite sweep exhaustive unless the finite set itself is the theorem's quantified domain. Do not call HomotopyContinuation `certify()` an ideal-emptiness proof unless all roots of a complete polynomial system are covered by a valid exclusion argument. Do not use `dim I=-1` over inexact `CC` as an exact certificate.

## Primary sources inspected

- Bengtsson et al., [Mubs and Hadamards of Order Six](https://arxiv.org/abs/quant-ph/0610161), HTML [v3](https://arxiv.org/html/quant-ph/0610161v3).
- Brierley and Weigert, [Constructing Mutually Unbiased Bases in Dimension Six](https://arxiv.org/abs/0901.4051), HTML [v2](https://arxiv.org/html/0901.4051v2).
- Karlsson, [Three-parameter complex Hadamard matrices of order 6](https://arxiv.org/abs/1003.4177), HTML [v1](https://arxiv.org/html/1003.4177v1).
- Grassl, [On SIC-POVMs and MUBs in Dimension 6](https://arxiv.org/abs/quant-ph/0406175), HTML [v2](https://arxiv.org/html/quant-ph/0406175v2).
- Maxwell and Brierley, [On properties of Karlsson Hadamards and sets of Mutually Unbiased Bases in dimension six](https://arxiv.org/abs/1402.4070), HTML [v2](https://arxiv.org/html/1402.4070v2).
- McNulty and Weigert, [Mutually Unbiased Bases in Composite Dimensions -- A Review](https://arxiv.org/abs/2410.23997).
- McNulty and Weigert, [Comment on Product states and Schmidt rank of mutually unbiased bases in dimension six](https://arxiv.org/abs/2504.13067).
- Mateo Cárdenes Wuttig and Joseph Tindall, [A Complete Classification of Complex Hadamard Matrices of Order Six](https://arxiv.org/abs/2608.18053), v2. This source was inspected as a current preprint claim; its classification and formal certificates were not independently reproduced here.

Page-number note: the arXiv HTML was inspected directly. Where this report does not give a PDF page, the corresponding section/equation/theorem is the controlling locator; a PDF page number is `UNVERIFIED` rather than inferred.

## Mandatory summary

| Question | Answer | Evidence | Scope |
|---|---|---|---|
| Is N(6)=3 proved? | No. | Grassl Sec. 4; Bengtsson Sec. 1; no global theorem found. | Global. |
| Is a fourth MUB known to exist? | No construction is known in the inspected primary literature. | Bengtsson abstract/conclusion; Brierley--Weigert abstract/conclusion. | Evidence, not proof. |
| Is a fourth MUB known not to exist globally? | No. | Primary papers state the problem remains open. | Global. |
| Is nonextendability known for D0? | Yes. | Brierley--Weigert Sec. 4.1, Table 1, Eq. (20). | Exact, fixed D0 pair. |
| Is nonextendability known for selected Diță points? | Yes, in several senses. | BW for D0; repository T3 certified-numerical at seven non-D0-equivalent Karlsson points. | Pointwise. |
| Is nonextendability known for all Diță parameters? | Not by an exact universal proof inspected here. | BW sampled/affine analysis; repository dense circle is numerical. | Full circle open. |
| Is nonextendability known for all Karlsson matrices? | No. | Maxwell--Brierley LP inconclusive; Karlsson is only H2-reducible. | `K_6^(3)` global question open. |
| Does Karlsson's family cover all order-six CHMs? | Not by itself. | Karlsson Theorem 11 covers H2-reducible matrices; Wuttig--Tindall v2 Corollary 25 claims the remaining sectors are `G_6^(4)` and `T_6`. | Full CHM space; new classification claim not independently audited here. |
| What is the strongest defensible new target? | Exact no-fourth at one algebraic, non-D0-equivalent Dita/Karlsson point, preferably with complete-pool or unit-ideal coverage. | Repository `docs/ALGEBRAIC_ATTACK.md`, E0/E3b. | Local exact theorem. |
