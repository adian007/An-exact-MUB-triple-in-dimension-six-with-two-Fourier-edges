# Agent 3 — Four-Parameter Model Report

## Executive Summary

**NO_VALID_FOUR_PARAMETER_EXTENSION_OF_KARLSSON_IDENTIFIED.**

A genuine four-parameter complex Hadamard family exists in the literature — Szöllősi's G_6^(4) — but it is a **separate published family**, not a verified extension of Karlsson's K_6^(3). The repository contains no derivation establishing K_6^(3) ⊂ G_6^(4).

## What "Four Parameters" Can Validly Mean

| Interpretation | Valid? | Source |
|---|---|---|
| Szöllősi G_6^(4) — a genuine 4-parameter CHM family | YES | Szöllősi 2012, arXiv:1008.0632 |
| Karlsson K_6^(3) extended to 4 parameters | NO justification | No source found |
| Diţă 4-parameter families | YES (separate) | Diţă 2012, arXiv:1207.2593 |
| Wuttig-Tindall 4-corner reconstruction | UNVERIFIED | Preprint arXiv:2608.18053 |

## Szöllősi G_6^(4) — Key Facts

- **Published in**: JLMS 85 (2012), 616–632
- **Parameter count**: 4 real parameters
- **Structure**: Entries are algebraic functions involving roots of sextics
- **Relation to Karlsson**: Described as "generic" vs Karlsson's "degenerate" — no verified set inclusion
- **Status**: Not implemented in this repository

## Parameter Domain

The exact parameter domain is defined in the primary paper (arXiv:1008.0632). Key features:
- Four real parameters with algebraic dependence through sextic roots
- Branch cuts from the sextic root choices
- No rectangular fundamental domain may be assumed
- Singular strata occur where the sextic discriminant vanishes

## Gauge Freedom

Standard CHM equivalence (left/right monomial unitaries) applies. The dephased gauge fixes first row and first column to 1, reducing the continuous gauge from 11 to (at most) 4 parameters — consistent with the family dimension.

## MUB Equations

The MUB equations (MU constraints, orthogonality, witness systems) remain polynomial/rational in the parameters once the sextic root branches are explicitly tracked. However:
- The algebraic degree is higher than Karlsson's family
- Branch choices must be recorded for reproducibility
- The mixed volume of associated polynomial systems will be larger

## Recommendation

**Use G_6^(4) as a separate four-parameter sector, NOT as a Karlsson extension.** If this family is to be used for MUB research:
1. First independently transcribe and validate the formulas from arXiv:1008.0632
2. Verify every generated matrix satisfies CHM constraints
3. Identify all singular/boundary strata explicitly
4. Determine whether the MUB pool equations remain tractable
5. Label all results as "G_6^(4) sector" — never as "Karlsson extended"

## What Remains Open

- Whether G_6^(4) contains any MUB triples at all
- Whether any G_6^(4) triple extends to four MUBs
- Whether the G_6^(4) and K_6^(3) sectors overlap or are disjoint
- Whether a complete classification of order-six CHMs exists (Wuttig-Tindall preprint claims this but is unverified)

## Status

| Item | Status |
|---|---|
| G_6^(4) formulas transcribed | NOT_RUN |
| G_6^(4) CHM validation | NOT_RUN |
| Parameter domain mapped | NOT_RUN |
| Singular strata identified | NOT_RUN |
| MUB equations analyzed | NOT_RUN |
| Karlsson extension verified | **REJECTED** — no source supports it |

## Sources

1. F. Szöllősi, "Complex Hadamard matrices of order 6: a four-parameter family," JLMS 85 (2012), 616–632, arXiv:1008.0632
2. B.R. Karlsson, "Three-parameter complex Hadamard matrices of order 6," LAA 434 (2011), 247–258, arXiv:1003.4177
3. P. Diţă, "Four-parameter families of complex Hadamard matrices of order six," arXiv:1207.2593
4. M. Cárdenes Wuttig and J. Tindall, "A Complete Classification of Complex Hadamard Matrices of Order Six," v2 (2026-08-28), arXiv:2608.18053