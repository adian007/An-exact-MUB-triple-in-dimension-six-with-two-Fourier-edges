# An Exact MUB Triple in Dimension Six with Two Fourier Edges

Research code, symbolic computations, numerical experiments, and reproducibility
artifacts for the study of mutually unbiased bases (MUBs) in dimension six.

> **TL;DR.** We construct an explicit MUB triple in $\mathbb{C}^6$,
> $\lbrace I,\ D(\rho)/\sqrt{6},\ B \rbrace$ with $\rho = e^{i\pi/3}$, built from the
> Diţă slice of Karlsson's three-parameter family of complex Hadamard matrices.
> The two transition matrices involving the new basis $B$ lie (up to equivalence
> and transposition) in the two-parameter **Fourier family**. A known
> non-extendability theorem for Fourier-family pairs then shows this particular
> triple **cannot be extended to a fourth MUB**.

---

## Contents

1. [Scope: what is and is not claimed](#1-scope-what-is-and-is-not-claimed)
2. [Background](#2-background)
3. [Main result](#3-main-result)
4. [The Diţă slice](#4-the-diţă-slice)
5. [Why the Fourier family matters](#5-why-the-fourier-family-matters)
6. [Evidence status](#6-evidence-status)
7. [Repository layout](#7-repository-layout)
8. [Reproducing the results](#8-reproducing-the-results)
9. [Provenance and known caveats](#9-provenance-and-known-caveats)
10. [Citation](#10-citation)
11. [References](#11-references)
12. [License](#12-license)

---

## 1. Scope: what is and is not claimed

Read this first.

| Statement | Status |
|---|---|
| The constructed triple $\lbrace I, D(\rho)/\sqrt6, B\rbrace $ is non-extendable to a fourth MUB | **Claimed**, via the Fourier-family obstruction |
| Every MUB triple in dimension six is non-extendable | **Not claimed** |
| Four MUBs do not exist in dimension six ($N(6)=3$) | **Not claimed** |
| Results hold for all complex Hadamard matrices of order 6 | **Not claimed**; the work is confined to the Karlsson $K_6^{(3)}$ family and the Diţă slice |

Only the first row follows from the exact construction plus the Fourier-family
obstruction used here. Please do not cite this work as evidence for the other
three.

---

## 2. Background

A set of orthonormal bases $\mathcal B_1,\dots,\mathcal B_k$ in $\mathbb C^d$ is
**mutually unbiased** if, for every pair of distinct bases,

```math
|\langle u_i, v_j\rangle|^2 = \frac1d .
```

For $d=6$ this reads $|\langle u_i, v_j\rangle| = 1/\sqrt6$.

Whether a complete set of $d+1 = 7$ MUBs exists in dimension six is a
long-standing open problem; the largest set known is three bases.

**Connection to Hadamard matrices.** Fix the computational basis $I$. A basis
unbiased to $I$ is encoded by a complex Hadamard matrix $H$:

```math
H H^\dagger = 6 I, \qquad |H_{jk}| = 1 .
```

Normalizing by $1/\sqrt6$ gives a unitary whose columns form the unbiased basis.
This project works through that correspondence, focusing on Karlsson's
three-parameter family $K_6^{(3)}$ and its one-parameter Diţă slice.

---

## 3. Main result

The central object is the exact triple

```math
\left\lbrace I,\ \frac{D(\rho)}{\sqrt6},\ B \right\rbrace,
\qquad \rho = e^{i\pi/3},
```

where

- $I$ is the computational basis;
- $D(\rho)$ is a complex Hadamard matrix on the Diţă slice;
- $B$ is an explicitly constructed third basis;
- all three bases are pairwise mutually unbiased.

**Key structural observation.** The two transition matrices linking $B$ to the
other two bases belong, up to the relevant equivalences and transposition, to
the two-parameter Fourier family $F_6^{(2)}(a,b)$.

**Consequence.** A known non-extendability result for Fourier-family pairs
applies, so the triple is not contained in any set of four MUBs.

A second exact construction is studied at $z = 1$.

---

## 4. The Diţă slice

The project uses the following one-parameter slice of the Karlsson family:

```math
D(z)=
\begin{pmatrix}
1&1&1&1&1&1\\
1&-1&z&-z&i&-i\\
1&-i&i&i&-i&-1\\
1&i&-z&z&-1&-i\\
1&\bar z&-i&-1&-z&i\\
1&-\bar z&-1&-i&z&i
\end{pmatrix},
\qquad |z|=1,\quad z = e^{i\lambda}.
```

| Point | Role |
|---|---|
| $z = \rho = e^{i\pi/3}$ | Main exact construction |
| $z = 1$ | Second exact construction |
| Other $\lambda$ on the slice | Numerical and certified studies |

---

## 5. Why the Fourier family matters

```text
Exact MUB triple {I, D(ρ)/√6, B}
              │
              ▼
Two transition matrices involving B
              │
              ▼
Membership in the Fourier family F₆⁽²⁾(a,b)
(after equivalences and transposition)
              │
              ▼
Known obstruction: the pair (I, F(a,b)) admits no extension
to a MUB quartet
              │
              ▼
This particular triple is non-extendable
```

The obstruction used is Theorem 1.4 of Jaming, Matolcsi, Móra, Szöllősi and
Weiner (arXiv:0902.0882), whose proof relies on a computer-assisted discretized
search. See [References](#11-references).

**Orientation caveat.** Fourier-family parameter placement differs between
the literature and this repository's fitting code (one is the transpose of the
other). Extension of $(I,H)$ is equivalent to extension of $(I,H^\dagger)$ and,
after complex conjugation, of $(I,H^{T})$, so both orientations are covered, but
every fit in this repository states which convention it uses.

---

## 6. Evidence status

Every claim in this repository carries an explicit tier:

| Tier | Meaning |
|---|---|
| **EXACT** | Symbolic or algebraic proof over an exact field |
| **CERTIFIED** | Interval-arithmetic or otherwise rigorously certified numerics |
| **NUMERICAL** | Double-precision evidence only; residuals reported |
| **OPEN** | Not yet established |

<!-- VERIFY: sync this table with SCIENTIFIC_STATUS.md before publishing. -->

| ID | Claim | Tier |
|---|---|---|
| T1 | Gauge structure of the family | EXACT (symbolic) |
| T2 | Third MUB on the Diţă circle | CERTIFIED / NUMERICAL (sampled; see ledger) |
| T3 | Local fourth-MUB exclusion: 40/40 recovered-pool six-cliques at seven $\Lambda_{\text{cert}}$ points | CERTIFIED (local, not a full-circle theorem) |
| T4 | Exact unit-ideal Gröbner result, one $(D_{bc},F_D)$ witness, degree 16 over $\mathbb Q(\zeta_{24},\sqrt5)$ | EXACT (re-verification of a stronger published result, not a new theorem) |
| F | $B$ and $D^\dagger B$ fit $F_6^{(2)}(a,b)$ directly at tested $\lambda$ (residual $\sim10^{-15}$) | NUMERICAL |
| X | Exact identity placing dephased $D(z)$ in the 2-circulant family $X_6$ | EXACT |
| — | Fourier-family membership for **all** relevant cliques; pool completeness; hypotheses of the quartet theorem verified exactly | OPEN |

If a statement in the paper is stronger than the tier listed here, the ledger
wins and the paper should be corrected.

The authoritative ledger is `SCIENTIFIC_STATUS.md`; claim-level provenance is in
`claim_trace.csv`.

---

## 7. Repository layout

<!-- VERIFY: adjust to the actual tree. -->

```text
.
├── README.md
├── SCIENTIFIC_STATUS.md        # authoritative evidence ledger
├── RETRACTION_NOTICE.md        # retracted claims and why
├── claim_trace.csv             # per-claim evidence-level decision table
├── src/                        # Julia modules (incl. Certification.jl)
├── scripts/
│   ├── julia/                  # certification and search scripts
│   └── python/                 # validation and fitting scripts
├── results/                    # machine-generated outputs (never overwritten)
└── docs/
    ├── overview/               # status, reproduction, forensic report, open problems
    ├── paper/                  # drafts and proofs
    ├── research/               # notes, literature, agent reports
    └── results/                # human-readable result summaries
```

---

## 8. Reproducing the results

### Requirements

| Tool | Used for |
|---|---|
| Julia ≥ 1.12 (with HomotopyContinuation.jl, Arblib.jl) | Numerical solving and certified arithmetic |
| Python 3 (NumPy, SciPy, SymPy) | Fitting and table validation |
| Macaulay2 | Exact Gröbner elimination (T4) |

### Steps

```bash
# 1. Clone
git clone <REPO_URL>
cd <REPO_DIR>

# 2. Julia environment
julia --project=. -e 'using Pkg; Pkg.instantiate()'

# 3. Fast regression suite (expected: 101/101 pass)
julia --project=. <PATH_TO_REGRESSION_ENTRYPOINT>

# 4. Certification and exact checks
#    TODO: list exact commands per claim (T1, T3, T4, family-membership fits)
```

**Windows note.** If Application Control blocks precompiled DLLs, run Julia with
`--compiled-modules=no`. This workaround is also recorded in the relevant script
headers.

Each result file should be traceable to the commit and command that produced it.
Record Julia, Python and Macaulay2 versions alongside any regenerated output.

---

## 9. Provenance and known caveats

- **Normalization boundary (2026-09-17).** A fix in `src/Certification.jl` now
  auto-detects unit-normalized versus $\mathrm{RHS}=6$ third-basis vectors.
  Results produced before and after this change are **not** interchangeable.
- **Retracted diagnostic.** A one-column three-$(-1)$ helper was previously cited
  as Fourier-family evidence. It is only a necessary signature; the full
  criterion (Matszangosz–Szöllősi, Theorem 1) requires three distinct columns each
  containing a $-1$. See `RETRACTION_NOTICE.md`.
- **Retired system.** The $n_{\text{wit}}=2$ system is retired by
  `src/Certification.jl` and must not underpin new claims.
- **Row-count discrepancy (848 vs 928).** Traced to
  `validate_third_mub_table.py` reading a stale backup CSV, not to a pipeline
  logic error. Final reconciliation should be documented.
- **Automated history.** Much of the commit history consists of automated
  checkpoints, and result files do not embed a producing-commit hash. File
  timestamps alone are not proof of which script produced an artifact.
- **Review status.** No human domain-expert review has been completed yet.

---

## 10. Citation

If you use this work, please cite the paper once available:

```bibtex
@misc{TODO_key,
  author = {TODO},
  title  = {An exact {MUB} triple in dimension six with two {F}ourier edges},
  year   = {2026},
  note   = {arXiv:TODO}
}
```

---

## 11. References

1. P. Jaming, M. Matolcsi, P. Móra, F. Szöllősi, M. Weiner, *A generalized
   Pauli problem and an infinite family of MUB-triplets in dimension 6*,
   arXiv:0902.0882 (Theorem 1.4: non-extendability of $(I, F(a,b))$).
2. B. R. Karlsson, *Three-parameter complex Hadamard matrices of order 6*,
   arXiv:1003.4177.
3. G. Matszangosz, F. Szöllősi, Des. Codes Cryptogr. (2024) (Fourier-family
   characterization via the three-$(-1)$ criterion).
4. D. McNulty, S. Weigert, review of MUBs, arXiv:2410.23997.
5. D. Goyeneche, arXiv:1209.4126.
6. W. Tadej, K. Życzkowski, *A concise guide to complex Hadamard matrices*.

---

## 12. License

TODO: add a license (e.g. MIT for code, CC-BY-4.0 for text and data).
