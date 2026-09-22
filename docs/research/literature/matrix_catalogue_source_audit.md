# Matrix Catalogue and Claim Audit

**Date:** 2026-09-17  
**Mode:** literature and source verification only  
**Computations:** none run for this audit; no Gröbner basis, sweep, or parameter elimination was used.

## Review boundary

This audit verifies the named matrix catalogue, local source/transcription locations, primary literature, and the quantifier of the associated claims. It does not certify the local numerical outputs. The project remains scoped to Karlsson's \(K_6^{(3)}\) family unless a claim explicitly names another anchor.

## Catalogue

| Object or family | Local source of truth | Primary source checked | Source-aligned status |
|---|---|---|---|
| Fourier \(F_6\) | `src/mub_zauner_6d_liang_chen.jl`; `src/Benchmarks.jl` | Grassl, arXiv:quant-ph/0406175, Eq. (15); Brierley--Weigert, arXiv:0901.4051v2 | Aligned up to Fourier-sign/conjugation convention, which is standard-equivalence compatible. The local pool claim is fixed-pair, not family-wide. |
| Diță \(D_0\) and \(D(x)\) | `src/brierley_weigert_notes.jl`; `src/Benchmarks.jl`; `paper/proofs/fourth_mub_theorems.tex` | Brierley--Weigert, arXiv:0901.4051v2, Eq. (30), Sec. 4.1, Eq. (20), Appendix A.2 | Aligned. The source reports 120 MU vectors, 10 additional Hadamards/bases, and no pairwise-MU pair among those ten. The local claim must retain the distinction between vectors and bases. |
| Karlsson \(K_6^{(3)}\) | `src/Karlsson.jl`; `src/mub_zauner_6d_liang_chen.jl` | Karlsson, arXiv:1003.4177v1, Theorem 11, Eqs. (2.4)--(2.5), (5.1); LAA 434 (2011) | Aligned. The local \(A\)-block, \(B=-F_2-A\), block assembly, and four Möbius relations match the primary formula. The source parameter domain uses reduced ranges; local exploratory coordinates may be redundant. |
| Karlsson Fourier seam / degenerate branches | `src/Karlsson.jl`; `docs/METHODOLOGY_AUDIT.md` | Karlsson, arXiv:1003.4177v1, Sec. 5 and Appendices 1--2 | Conceptually aligned. The source explicitly warns that degenerate Möbius maps require limiting/seam treatment and that limit families can depend on approach direction. Any exact-seam formula must not be presented as the \(\theta\to0^+\) limit without a separate argument. |
| Tao isolated matrix \(S_6(0)\) | `src/Benchmarks.jl`; `research/constructions/dimension6_constructions.md` | Wuttig--Tindall, arXiv:2608.18053v2, Definition 5 / Eq. (10); Tao (2004) | Aligned in role and source. It is an isolated order-six CHM sector and is a negative-control anchor, not a Karlsson-family point. |
| General order-six catalogue | `docs/LITERATURE_2026.md`; `docs/LIANG_CHEN_FAMILY_ASSESSMENT.md` | Wuttig--Tindall, arXiv:2608.18053v2, Theorem 6, Corollary 25; McNulty--Weigert, arXiv:2410.23997v2 | Updated source context. The 2026 preprint states the exact three-sector classification \(\mathcal G_6^{(4)}\cup\mathcal K_6^{(3)}\cup\mathcal T_6\). The current repository intentionally implements only \(\mathcal K_6^{(3)}\), with Tao as an anchor; it does not implement \(\mathcal G_6^{(4)}\). |
| Liang/Chen/Long--Qiu special matrices | `docs/LIANG_CHEN_FAMILY_ASSESSMENT.md`; `src/mub_zauner_6d_liang_chen.jl` stub | Liang et al., arXiv:2110.12206; QIP 23 (2024), 278; Wuttig--Tindall v2, Ref. [27] | No evidence in the checked catalogue that these papers provide a Karlsson-style three-parameter family. The local stub is therefore appropriately non-operational, but the negative assessment should remain explicitly source-reviewed rather than called a mathematical impossibility. |
| Szöllősi / \(\mathcal G_6^{(4)}\) sector | Not implemented; discussed in `docs/LIANG_CHEN_FAMILY_ASSESSMENT.md` and `docs/LITERATURE_2026.md` | Szöllősi (2012); Wuttig--Tindall v2, Sec. V and Corollary 25 | Catalogue entry is acknowledged but not locally transcribed. No Karlsson-family claim may be generalized to this sector. |

## Claim audit

| Local claim | Audit result | Required wording discipline |
|---|---|---|
| Grassl gives a Fourier-pair obstruction | Supported for the fixed pair \(\{I,F_6\}\): the source states that no more than one additional Hadamard can complement the pair. | Do not upgrade this to a statement about every order-six CHM or every MUB triple. |
| Brierley--Weigert gives the Dita benchmark | Supported: 120 vectors, 10 bases, no two of those bases mutually unbiased. | “Ten third bases” means ten candidate bases from the 120-vector pool; it does not mean ten mutually compatible bases. |
| Brierley--Weigert sampled known families | Supported, but the paper itself describes a sampled/computational result and distinguishes rigorous special cases from approximate non-affine calculations. | Keep “sampled,” “fixed matrix,” and “upper bound” qualifiers. |
| Karlsson is a three-parameter family | Supported by Karlsson's Theorem 11 for the \(H_2\)-reducible sector. | Do not call it the full order-six CHM space. Wuttig--Tindall v2 explicitly places it alongside \(\mathcal G_6^{(4)}\) and Tao. |
| The local \(A\)-block is the authoritative transcription | Supported against Karlsson v1. The McNulty--Weigert printed equations are a separate review transcription and must be labelled as such until corrected by the authors or an erratum. | State “the local implementation follows Karlsson's primary formula”; avoid implying that the survey authors endorsed the local correction. |
| Seven Dita points or 40/40 W1 cases exclude a fourth MUB | This is a local certified-computation claim, not a literature claim. | Keep its finite-point and fixed-triple quantifiers; it does not establish a Dita-circle, Karlsson-family, or global \(N(6)=3\) theorem. |
| “No fourth MUB on all of \(K_6^{(3)}\)” | Not established by the checked literature or by the local claim ledger. | Keep **Open**. |
| “No fourth MUB in all of dimension six” or \(N(6)=3\) | Not established. The primary sources and the 2026 classification both leave the MUB existence question open. | Keep outside project scope. |

## Independent review findings

1. The catalogue is internally coherent for a Karlsson-focused project: \(F_6\), \(D_0\), the Karlsson family, Tao's isolated matrix, and the unimplemented \(\mathcal G_6^{(4)}\) sector are distinguishable.
2. The strongest literature-supported local anchors are the exact/special Brierley--Weigert counts: \(F_6\) has 48 MU vectors and 16 candidate bases with no pairwise-MU pair; \(D_0\) has 120 MU vectors and 10 candidate bases with no pairwise-MU pair; Tao's matrix has 90 MU vectors and no third basis in that paper's reported analysis.
3. The local documentation is correct to separate the 2008 constellation paper from the 2009 algebraic-construction paper. Citation keys and wording should continue to preserve that separation.
4. Wuttig--Tindall v2 materially changes the background catalogue as of 2026: \(\mathcal G_6^{(4)}\) is no longer merely a conjectural or loosely described sector in this project context. It remains unimplemented here, so no computational result in this repository covers it.
5. Before any later elimination or sweep, a second human review should verify the exact local transcriptions of `src/brierley_weigert_notes.jl`, `src/Benchmarks.jl`, and `src/Karlsson.jl` against the cited equations. This audit identifies the source targets; it is not a substitute for that line-by-line sign and indexing review.

## Sources

- Karlsson, “Three-parameter complex Hadamard matrices of order 6,” arXiv:1003.4177v1, LAA 434 (2011), 247--258.
- Brierley and Weigert, “Constructing Mutually Unbiased Bases in Dimension Six,” arXiv:0901.4051v2, Phys. Rev. A 79, 052316 (2009).
- Grassl, “On SIC-POVMs and MUBs in Dimension 6,” arXiv:quant-ph/0406175v2.
- Bengtsson et al., “MUBs and Hadamards of Order Six,” arXiv:quant-ph/0610161v3, J. Math. Phys. 48, 052106 (2007).
- McNulty and Weigert, “Mutually Unbiased Bases in Composite Dimensions -- A Review,” arXiv:2410.23997v2, Quantum 10, 2051 (2026).
- Wuttig and Tindall, “A Complete Classification of Complex Hadamard Matrices of Order Six,” arXiv:2608.18053v2.
- Liang, Chen, Long, and Qiu, “Some special complex Hadamard matrices of order six,” arXiv:2110.12206.
