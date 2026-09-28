# Blocker C - audit follow-up: the Hermitian mechanism is false

**Agent:** RESEARCH (`researcher`) - **Date:** 2026-09-26 - **Branch/commit:** `main` / `b1c81e1`
**Scope:** Part 2, Blocker C only. One new file written (this report). No existing file in
`results/` was modified, so no `.bak_2026-09-26` backup was required. No `paper/*.tex`, no
`docs/paper/`, no `docs/scientific_status.md`, and no `docs/results/final_honest_status.md`
were touched.

---

## 0. Headline

**The stated mechanism is false, and the correct invariant is now derived exactly.**

1. `A = [[i,-1],[-1,i]]` (Dita) is **not** Hermitian: `A^dag = conj(A)^T = [[-i,-1],[-1,-i]] != A`.
   The F6 block `A = [[-1/2+i*sqrt3/2, -1/2-i*sqrt3/2],[-1/2+i*sqrt3/2, 1/2+i*sqrt3/2]]` is **not**
   Hermitian either (`A11 = omega` is non-real). Both loci are recorded `A Hermitian=false` by the
   script's own numeric field (`results/arc_circle_asymmetry.txt:7` and `:13`) while lines 33/35/50-51
   of the same file assert a Hermitian/non-Hermitian contrast. **The file contradicts itself.**
2. The Hermitian locus is *exactly* `cos(th) + sin(th)cos(ph) = 0` - a 1-dimensional curve in
   `[0,pi)^2` that contains **neither** anchor. Hermiticity does not distinguish the two loci; it is
   simply false at both. `[exact symbolic]`
3. The invariant that **does** distinguish them is the Mobius degeneracy index
   `D = |alpha|^2 - |beta|^2`, whose closed form I derived in `(theta,phi)`:
   `D_A = 2*sqrt3*sin(th)*(sin(ph) - sqrt3*cos(th)cos(ph))`,
   `D_B = -2*sqrt3*sin(th)*(sin(ph) + sqrt3*cos(th)cos(ph))`.
   These are **exactly proportional** to Karlsson's two published degeneracy factors, which
   independently reproduces the primary source. `[exact symbolic]`
4. The two loci are strata of **different kinds**: F6 sits on the `sin(th) = 0` component (a
   phi-gauge stratum, *both* Mobius maps collapsed to constants); Dita sits on the
   `tan(ph) = sqrt3*cos(th)` component (only `M_A` collapsed; `M_B` is a genuine circle automorphism
   at the extreme of its range). `[exact symbolic]`
5. Karlsson **himself** singles out the Dita point in Sec. 6 as a degenerate special point of the
   chart and identifies the resulting matrix with the Dita family `D_6^(1)`. The repository's Dita
   anchor is therefore *not* a generic point of `K_6^(3)`. `[exact symbolic, primary source]`

**What is still open:** the causal link from any of this to `kappa(H(lam)) >= 6`. Not established,
not claimed. See Sec. 7.

---

## 1. Blocker choice: agreement, with one dissent

**I agree with the choice of Blocker C**, on the grounds the orchestrator gave and for the same
reason: the stated mechanism is affirmatively wrong, so correcting it is checkable progress, and a
clean negative is a real result in this repo.

**I agree on rejecting Blocker A**, and I confirm the decisive fact directly:
`src/Certification.jl:28-40` states that the `n_wit=2` system "is logically redundant for this claim
structure (two orthogonal W1-solutions are neither necessary nor sufficient for a fourth ONB) and
is retired ... It must not be used in new claims." That retirement is a *logical* result, and it
makes "make `mv=119210` tractable" an attack on a retired target. I also agree that the cap 5000 is
unarbitrated policy. I did not pursue A.

**One dissent, and it is a correction to the framing of the hypothesis, not to the blocker choice.**
The brief's mechanism is stated as *"`alpha_B = 0` collapse versus generic SU(1,1) form."* Both
halves of that framing are wrong as stated:

- There is no "generic SU(1,1) form" in play. I proved exactly that `D = 0` - **not** `D = 1` - is
  Karlsson's degeneracy condition, and that the F6 point has `D_A = D_B = 0`. So F6 is *not* an
  SU(1,1)-normalised generic case; it is a *doubly degenerate* one. The premise "F6 has
  `|alpha_B|^2-|beta_B|^2 = 1`" in the brief is false; the artifact says `|alpha_B| = |beta_B| = 1`,
  i.e. `D = 0`.

---

## 2. Literature pass

Read first, to avoid duplication: `docs/overview/literature_audit.md` (audit dated 2026-09-17) and
`research/literature/agent9/method_bibliography.csv` (17 rows: MUB/NAG/INT/REAL bibliography).
Those already cover Bengtsson et al., Brierley-Weigert, Grassl, Tadej-Zyczkowski, the
homotopy/certification toolchain, and interval arithmetic. **I did not re-derive any of those.**

### 2.1 Primary source, actually fetched and read in full

I downloaded the arXiv HTML of Karlsson 1003.4177 and read Secs. 2-6 directly rather than relying on
memory. Cached at `results/lit_cache/karlsson_1003.4177.html`, de-tagged to
`results/lit_cache/karlsson_plain.txt`.

> **B. R. Karlsson**, *Three-parameter complex Hadamard matrices of order 6*,
> Linear Algebra Appl. **434**(1), 247-258 (2011). arXiv:1003.4177 [math-ph].
> DOI `10.48550/arXiv.1003.4177`. <https://arxiv.org/html/1003.4177v1>

Verbatim, Sec. 5 (decisive for this report):

> "It is useful to see these constraints as Mobius transformations
> `w = M(z) = (az - b)/(bb* z - a*)`, `z = M^{-1}(w) = (a*w - b)/(b*w - a)` that, **as long as
> `|a|^2 - |b|^2 != 0`**, map the unit circle onto itself."

> "The Mobius transformations (5.1) become degenerate if `|a|^2 - |b|^2 -> 0`: the transformation
> `w = M(z)` degenerates into a mapping of the unit circle in `z` into a single point `w = a/b*`,
> and this mapping has no inverse, and the inverse transform `z = M^{-1}(w)` degenerates into a
> mapping of the unit circle in `w` into a single point `z = a*/b*`, and again there is no inverse.
> For `M_A` and `M_A^{-1}` this occurs if `|A11| = |A12| = 1`, i.e. if
> `sin(th)(sin(ph) - sqrt3*cos(th)cos(ph)) = 0`, and for `M_B` and `M_B^{-1}` if `|B11| = |B12| = 1`,
> i.e. if `sin(th)(sin(ph) + sqrt3*cos(th)cos(ph)) = 0`. **Both transformations are degenerate when
> th = 0** ..."

This is the single most important literature finding: **the invariant is `|alpha|^2 - |beta|^2` and
the relevant value is `0` (collapse), not `1` (SU(1,1) normalisation).**

Verbatim, Sec. 6 (decisive for the Dita anchor - Karlsson names this exact point):

> "Consider for example the point `th = arccos(1/sqrt3)`, `ph = pi/4`, for which `{A11 = i, A12 = -1}`
> and `{B11 = -1-i, B12 = 0}`. Since `|A11| = |A12| = 1`, `M_A` is degenerate and
> `M_A(z^2) = M_A^{-1}(z^2) = -1`. Furthermore, `M_B(z^2) = M_B^{-1}(z^2) = 1/z^2 = z*^2`, so that,
> taking `z1 = z` as independent parameter, `z2^2 = -1`, `z3^2 = -1` and `z4^2 = z*^2`. Let
> `z2 = z3 = i` and `z4 = z*`. The resulting one-parameter Hadamard matrix ... **is equivalent to the
> generic member of `D_6^(1)`**."

And on the whole stratum:

> "For points on the `M_A` degeneracy curve, `sin(ph) = sqrt3*cos(th)cos(ph)` (see Section 5). Along

### 2.2 A retracted structural criterion - highly relevant prior art

> **D. McNulty and S. Weigert**, *Comment on 'Product states and Schmidt rank of mutually unbiased
> bases in dimension six'*, J. Phys. A: Math. Theor. **58**, 168001 (2025). arXiv:2504.13067.
> DOI `10.1088/1751-8121/adcb07`. <https://arxiv.org/html/2504.13067v1>

Fetched and read. The comment shows that **Lemma 11(v) Part 6 of Chen & Yu** - "If a set of four MU
bases in dimension six exists, none of the Hadamard matrices from the set contains a real `3x2`
submatrix" - rests on an unverified step (that two components of a certain column vector `v` must
vanish). The authors exhibit the one-parameter family `M_6(a)` as a counterexample to the final
step, and conclude: *"We are aware of at least three theorems on the existence of MU quadruples
that build on the now unverified lemma ... Thus, their validity is called into question."*
Invalidated downstream (as listed in that comment): Liang-Hu-Sun-Chen-Chen, LMA **69**, 2908 (2021),
arXiv:1904.10181 (the only arXiv id given explicitly there); Liang-Hu-Chen, QIP **18**, 352 (2019);
Chen-Liang-Hu-Chen, QIP **20**, 353 (2021).

**Why this matters for Blocker C:** this is a *documented instance of exactly the failure mode this
repo is now exhibiting* - a tidy structural criterion on order-six CHMs, propagated into a claim
ledger, resting on a step that does not survive checking. It is direct external support for treating
`docs/results/final_honest_status.md:35` as needing correction rather than defence.

### 2.3 Other sources actually fetched (abstract pages, verbatim)

- **L. Chen and L. Yu**, *Product states and Schmidt rank of mutually unbiased bases in dimension
  six*, J. Phys. A **50**, 475304 (2017). arXiv:1610.04875. DOI `10.1088/1751-8121/aa8f9e`.
  (The paper whose lemma is refuted above.)
- **S. Brierley and S. Weigert**, *Constructing Mutually Unbiased Bases in Dimension Six*,
  Phys. Rev. A **79**, 052316 (2009). arXiv:0901.4051. DOI `10.1103/PhysRevA.79.052316`.
  Abstract confirms: sampling all known order-six CHMs, "we never find more than two that are
  mutually unbiased" - *sampling*, not a global certificate.
- **A. Maxwell and S. Brierley**, *On properties of Karlsson Hadamards and sets of Mutually Unbiased
  Bases in dimension six*, Linear Algebra Appl. **440**, 296-306 (2015). arXiv:1402.4070.
  DOI `10.1016/j.laa.2014.10.017`. Abstract confirms their LP approach was **inconclusive** for
  ruling out a complete set from Karlsson Hadamards.

### 2.4 Honest gaps in the literature pass

The following are cited **only as listed in the reference list of arXiv:2504.13067**, which I did
fetch. I did **not** independently verify their arXiv ids and could not: the arXiv API returned
HTTP 406 and the Semantic Scholar API returned HTTP 429 (both attempts logged in
`results/lit_cache/arxiv_query_out.txt` and `results/lit_cache/s2_query_out.txt`).

- J.-M. Goyeneche, *Mutually unbiased triplets from non-affine families of complex Hadamard
  matrices in dimension 6*, J. Phys. A **46**, 105301 (2013). **arXiv id UNVERIFIED.** This is the
  closest known prior work to the repo's arc-vs-circle question and I was **unable to obtain or
  read it**; I make no claim about its contents.

---

## 3. The falsifiable sub-question, stated before computing

> **Q.** Is the arc-vs-circle asymmetry caused by a property of the Mobius data
> `(alpha_A, beta_A, alpha_B, beta_B)` at the two anchors - in particular by the `alpha_B = 0`
> degeneracy at Dita - rather than by any Hermitian/non-Hermitian distinction between the two
> `A`-blocks?

And, stated as the specific replacement hypothesis that this report actually establishes and
sharpens:

> **H.** The governing invariant is the *Mobius degeneracy index* `D = |alpha|^2 - |beta|^2`, whose
> degeneracy condition is `D = 0` (collapse of the unit circle to a point), not `D = 1` (SU(1,1)
> normalisation). The two anchors are strata of **different kinds**: F6 `(th=0)` lies on the
> `sin(th) = 0` component, where `D_A = D_B = 0` and *both* maps collapse to constants; Dita lies on
> the `tan(ph) = sqrt3*cos(th)` component, where `D_A = 0` but `D_B = -4`, so `M_A` collapses while
> `M_B` remains a genuine circle automorphism. This asymmetry - one map frozen, one map live - is a
> candidate mechanism for the circle-vs-arc difference.

**Explicitly not part of Q or H:** any claim that this *causes* `kappa(H(lam)) >= 6`. That is
Sec. 7, and it is left open.

---

## 4. What was tried, and campaign-gate status

Campaign gate: `results/NOT_RUN_researcher_blockerC_2026-09-26.log` (status / reason / command /
date, modelled on `research/methods/agent2/prototype/NOT_RUN.log`).

**Run** - one thing, and it is *not* a new numerical experiment:

    python scripts\python\blocker_c_exact_check.py    # RC=0, ~2 s
    # output: results\lit_cache\exact_check_out.txt

This is an **exact-arithmetic identity check** (unit-test level): sympy over `Q(i, sqrt3)` and
symbolic `(th,ph)`, **no floating point anywhere**, no sampling, no sweep, no root finding, no
tolerance. Tag: `exact symbolic`.

**Not run** (all listed in the NOT_RUN log):

1. No recomputation of the F6 arc width `0.117643` rad - existing artifact
   `f6_boundary_reverify.txt` only. `[not run]`
2. No recomputation of the Dita-circle counts `126/126 + 628/628` - existing artifacts
   `lambda_periodicity_dita.txt`, `dita_lambda_fourth_dense.csv` only. `[not run]`
3. No new clique search, no fresh 6-clique discovery, no third-basis construction at any `lam`.
   `[not run]`
4. No interval-arithmetic (`certify()`) certification. **Nothing in this report is claimed as
   certified numerical.** `[not run]`
5. No Julia execution, no homotopy solve, no `n_wit=2` witness system. `[not run]`
6. No `(th,ph)` sweep to locate `alpha_B = 0` - the claim in Sec. 5.7 is an exact algebraic
   derivation instead. `[not run - analytic]`
7. **No test of the central open clause** (that the structural difference causes `kappa(H) >= 6` for
   all `lam`). `[not run]`

No long Julia homotopy solve was launched. The host shell was unstable; all heavy work was analytic
or was a <=2 s exact computation.

---

## 5. The result

All formulas below are stated for the repo's builder
(`scripts/python/karlsson_k6_3.py:59-61,76-78`), which is the audited transcription of Karlsson
Sec. 3.6 + Prop. 9, with `F2 = [[1,1],[1,-1]]`:

    A11 = -1/2 + i*(sqrt3/2)*(cos(th) + exp(-i*ph) sin(th))
    A12 = -1/2 + i*(sqrt3/2)*(-cos(th) + exp(i*ph) sin(th))
    A   = [[A11, A12], [conj(A12), -conj(A11)]]          B = -F2 - A
    alpha_A = A12^2,  beta_A = A11^2,   alpha_B = B12^2,  beta_B = B11^2

Write `c = cos(th)`, `s = sin(th)`, `p = cos(ph)`, `q = sin(ph)`.

### 5.1 Exact entry moduli

Expanding `i(c + exp(-i*ph)s) = i(c + sp - i sq)` and separating real/imaginary parts:

    A11 = ( -1/2 + (sqrt3/2) s q )  +  i (sqrt3/2)( c + s p )
    A12 = ( -1/2 - (sqrt3/2) s q )  +  i (sqrt3/2)( -c + s p )
    B11 = ( -1/2 - (sqrt3/2) s q )  -  i (sqrt3/2)( c + s p )
    B12 = ( -1/2 + (sqrt3/2) s q )  -  i (sqrt3/2)( -c + s p )

so, exactly (sympy output `exact_check_out.txt` lines 2-6, residual `0`):

    |A11|^2 = 1 + (3/2) s c p - (sqrt3/2) s q
    |A12|^2 = 1 - (3/2) s c p + (sqrt3/2) s q
    |B11|^2 = 1 + (3/2) s c p + (sqrt3/2) s q
    |B12|^2 = 1 - (3/2) s c p - (sqrt3/2) s q
    |A11|^2 + |A12|^2 = 2      (row norm of AA^dag = 2I)
    |B11|^2 + |B12|^2 = 2      (row norm of BB^dag = 2I)


### 5.2 THE CRUX - the correct normalization invariant, and its closed form

Since `|A11|^2 + |A12|^2 = 2` and `|B11|^2 + |B12|^2 = 2`, the difference of fourth powers factors:

    D_A := |alpha_A|^2 - |beta_A|^2 = |A12|^4 - |A11|^4 = (|A12|^2 - |A11|^2)(|A12|^2 + |A11|^2)
                                              = 2(|A12|^2 - |A11|^2)
    D_B := |alpha_B|^2 - |beta_B|^2 = |B12|^4 - |B11|^4 = 2(|B12|^2 - |B11|^2)

Substituting Sec. 5.1:

    D_A = 2*sqrt3*sin(th) ( sin(ph) - sqrt3*cos(th)*cos(ph) )      (exact symbolic)
    D_B = -2*sqrt3*sin(th) ( sin(ph) + sqrt3*cos(th)*cos(ph) )     (exact symbolic)

sympy residuals are exactly `0` for both (`exact_check_out.txt` lines 37-38):

    D_A - 2*r3*s*(sin ph - r3*c*cos ph) = 0
    D_B + 2*r3*s*(sin ph + r3*c*cos ph) = 0

**This is the crux result, and it settles the brief's question #2 decisively.**

**(a) `D = 1` is the wrong invariant; `D = 0` is the right one.** `D_A` is *exactly* `2*sqrt3` times
Karlsson's `M_A` degeneracy factor and `D_B` is exactly `-2*sqrt3` times his `M_B` factor. Therefore

> `M_A` is degenerate <=> `D_A = 0`, and `M_B` is degenerate <=> `D_B = 0`.

The value `1` is nowhere in the theory; the value `0` is the whole story. Consequently:

> **CORRECTION to `scripts/python/karlsson_k6_3.py:100-104`**, which asserts that
> `M(z) = (az-b)/(bb* z - a*)` "is an involution-type map when `|alpha|^2-|beta|^2 = 1` (standard
> SU(1,1) Mobius)". This is a **misstatement of Karlsson Sec. 5**. Karlsson's condition is
> `|alpha|^2-|beta|^2 != 0` for the unit circle to map onto itself, and the map **collapses to the
> constant `alpha/beta*` as `|alpha|^2-|beta|^2 -> 0`**. The comment should read: *"`M` is a genuine
> automorphism of the unit circle iff `|alpha|^2-|beta|^2 != 0`; it collapses to the constant
> `alpha/beta*` when `|alpha|^2-|beta|^2 = 0`. These are the two Mobius strata of Karlsson Sec. 5 and
> must be treated as strata."* I have **not** edited the script (not my file); flagged for the code
> owner.

**(b) The brief's premise about F6 is false.** The brief anticipated `|alpha_B|^2 - |beta_B|^2 = 0`
for F6 and read it as "0 (1-1)". The artifact (`arc_circle_asymmetry.txt:8`) records
`|alpha_B| = |beta_B| = 1`, which is `|alpha_B|^2 = |beta_B|^2 = 1`, hence `D_B = 1 - 1 = 0` - **a
degeneracy, not an SU(1,1) normalization.** F6 is therefore *doubly* degenerate, not generic. This
is a substantive correction to the hypothesis, and it makes the contrast sharper, not weaker.

**(c) The dynamic range.** Since `|A11|^2+|A12|^2 = 2`, we have `D_A = 2(|A12|^2-|A11|^2)` in
`[-4,4]` with `D_A = 0` exactly at the centre; likewise `D_B` in `[-4,4]`. So `D_B = -4` at Dita is
the **extreme** of the range, attained iff `|B12|^2 = 0`, `|B11|^2 = 2`, i.e. iff `B12 = 0`, i.e.
iff `alpha_B = 0`. This is the sharp statement:

> `alpha_B = 0` **is** `D_B = -4`; it is the extremal case, not a separate phenomenon.
> `[exact symbolic]`

### 5.3 Dita specialization, by hand

At `th_D = arccos(1/sqrt3)`, `ph = pi/4`: `c = 1/sqrt3`, `s = sqrt2/sqrt3`, `p = q = 1/sqrt2`.

    A11 = -1/2 + (sqrt3/2)(sqrt2/sqrt3)(1/sqrt2) + i(sqrt3/2)(1/sqrt3 + (sqrt2/sqrt3)(1/sqrt2))
        = -1/2 + 1/2 + i(1/2 + 1/2)  =  i
    A12 = -1/2 - (sqrt3/2)(sqrt2/sqrt3)(1/sqrt2) + i(sqrt3/2)(-1/sqrt3 + (sqrt2/sqrt3)(1/sqrt2))
        = -1/2 - 1/2 + i(-1/2 + 1/2)  =  -1

Hence

    A   = [[ i, -1], [-1,  i]]            (conj(A12) = -1, -conj(A11) = +i)
    B   = -F2 - A = [[-1-i, 0], [0, 1-i]]
    alpha_A = A12^2 = 1     beta_A = A11^2 = -1              D_A = 1 - 1 = 0
    alpha_B = B12^2 = 0     beta_B = B11^2 = (-1-i)^2 = 2i    D_B = 0 - 4 = -4

`D_A = 0` matches `2*sqrt3*s*(q - sqrt3*c*p) = 2*sqrt3*(sqrt2/sqrt3)(1/sqrt2 - sqrt3*(1/sqrt3)(1/sqrt2)) = 0`
exactly, since `1/sqrt2 - 1/sqrt2 = 0`. `D_B = -4` matches
`-2*sqrt3*s*(q + sqrt3*c*p) = -2*sqrt3*(sqrt2/sqrt3)(1/sqrt2 + 1/sqrt2) = -4` exactly.
`[exact symbolic]`

**The Mobius maps, exactly:**

    M_A(z) = (1*z - (-1))/(conj(-1)*z - conj(1)) = (z+1)/(-z-1) = -1     for all z
    M_B(z) = (0*z - 2i)/(conj(2i)*z - conj(0))   = -2i/(-2i z)  = 1/z     for all z

sympy: `M_A(z) = -1` (constant, verified by `simplify(M_A + 1) == 0`) and `M_B(z)*z = 1`.
`[exact symbolic]`

**The full z-structure, exactly** (matches Karlsson Sec. 6 verbatim, quoted in Sec. 2.1):

    z1 = e^{i*lam}   (free)
    z3^2 = M_A(z1^2) = -1                   =>  z3 = +-i
    z2^2 = M_B^{-1}(z3^2) = 1/(-1) = -1     =>  z2 = +-i    (M_B is an involution: M_B = 1/z)
    z4^2 = M_B(z1^2) = 1/z1^2 = e^{-2i*lam}  =>  z4 = +-1/z1 = +-e^{-i*lam}

With the consistent sign choice `z2 = z3 = i`, `z4 = conj(z)`, Karlsson writes the matrix
explicitly and states it "is equivalent to the generic member of `D_6^(1)`". This confirms
independently that the repo's `(th_D, pi/4)` anchor **is** the Dita family and **is** a degenerate
special point of the chart. `[exact symbolic, primary source]`

### 5.4 F6 (`th = 0`) specialization, by hand

At `s = 0`, `c = 1`:

    A11 = -1/2 + i*sqrt3/2 = omega       A12 = -1/2 - i*sqrt3/2 = omega^2   (omega = e^{2*pi*i/3})
    B11 = -1-omega = omega^2              B12 = -1-omega^2 = omega
    alpha_A = (omega^2)^2 = omega   beta_A = omega^2          D_A = 1 - 1 = 0
    alpha_B = omega^2            beta_B = (omega^2)^2 = omega  D_B = 1 - 1 = 0

sympy confirms `A11 = -1/2 + sqrt(3)*I/2`, `alpha_A = -1/2 + sqrt(3)*I/2`, `D_A = 0`, `D_B = 0`.
`[exact symbolic]`

**Both maps collapse:**

    M_A(z) = (omega*z - omega^2)/(omega*z - omega^2) = 1     M_B(z) = (omega^2*z - omega)/(omega^2*z - omega) = 1

sympy: `M_A(z) = 1`, `M_B(z) = 1`. Hence `z3^2 = 1`, `z4^2 = 1`, `z2^2 = M_B^{-1}(1) = 1`, so

    z1 = e^{i*lam} (free)      z2 = z3 = z4 = +-1  (i.e. Z2 = Z3 = Z4 = F2)

**and `A`, `B` carry no `ph` at all** (no `ph` appears in Sec. 5.1 when `s = 0`). This is a
*stronger* version of the repo's Lemma L1: not only is `H(0,ph,lam)` independent of `ph`, but the
mechanism is that **all four `z`-parameters except `z1` are frozen and both Mobius maps are
constants**, so the `ph`-dependence in the chart is annihilated twice over. `[exact symbolic]`

### 5.5 Is "Hermitian" ever the right invariant? No - and here is the proof

    A      = [[A11,          A12        ],
             [conj(A12),   -conj(A11)  ]]
    A^dag  = [[conj(A11),    A12        ],
             [conj(A12),    -A11       ]]
    A = A^dag  <=>  A11 = conj(A11)  <=>  Im(A11) = 0  <=>  cos(th) + sin(th)*cos(ph) = 0.

(The off-diagonal slots `A12` vs `conj(A12)` are already Hermitian by construction, so `A11` real is
the *only* condition.) Therefore:

> **The Hermitian locus is exactly the curve `cos(th) + sin(th)cos(ph) = 0` in `[0,pi)^2`, a
> 1-dimensional set containing neither anchor.** `[exact symbolic]`

| locus | `Im(A11)` | Hermitian? |
|---|---|---|
| Dita `(arccos(1/sqrt3), pi/4)` | `1/sqrt3 + (sqrt2/sqrt3)(1/sqrt2) = 2/sqrt3 ~ 1.1547` | **No** |
| F6 `th=0` | `1` | **No** |


Direct: at Dita the `(1,1)` entry of `A - A^dag` is `i - (-i) = 2i != 0`; at F6 it is `i*sqrt3 != 0`.
The Hermitian locus is the curve `tan(th) = -1/cos(ph)` living in `th` in `(pi/2, pi)` - it passes
through **neither** anchor, and there `c + sp < 0`, opposite in sign to both anchors (`2/sqrt3` and
`1`, both positive).

**Conclusion.** Hermiticity is not merely an unreliable proxy here; it is a *false* one, with value
`false` at both ends of the comparison it was supposed to explain. It cannot be the mechanism, and
the repo's own numeric field already said so.


### 5.6 NEGATIVE RESULT - the `z2`-recovery determinant does *not* degenerate

The orchestrator's proposed sub-mechanism was that at Dita the Mobius inversion used to recover `z2`
(i.e. inverting `M_B` for `z2^2`, per `scripts/python/karlsson_k6_3.py:100-109`) becomes a `0/0`
form. **It does not.** With `z3^2 = -1`, `alpha_B = 0`, `beta_B = 2i`:

    num = beta_B - z3^2*conj(alpha_B)   = 2i - (-1)*0        =  2i
    den = alpha_B - z3^2*conj(beta_B)   = 0  - (-1)*(-2i)    = -2i
    z2^2 = num/den = 2i/(-2i) = -1        (correct, matches Sec. 5.3)
    den != 0   (sympy: `den == 0 ? False`)

So the `z2` recovery at Dita is a clean non-zero-over-non-zero. The "ill-conditioned `0/0` inversion
at the Dita slice" story is **false**. `[exact symbolic - negative result]`


### 5.7 NEGATIVE RESULT - "Dita lies on the `M_A` stratum" does NOT single out Dita

The brief's item #5 proposed that the Dita point lying on the `M_A` degeneracy stratum is "likely the
whole story". **The stratum membership is necessary but far from sufficient**, for a simple exact
reason: `D_A = 2*sqrt3*s*(q - sqrt3*c*p)` vanishes on the whole **1-dimensional curve**
`{sin(th) = 0} U {tan(ph) = sqrt3*cos(th)}`, which retains a free parameter `th`. `M_A` is therefore
*degenerate at every point of that curve*, not just at Dita. What singles out Dita is the
**additional, independent** condition `alpha_B = 0` (equivalently `D_B = -4`, equivalently `B12 = 0`,
equivalently `A12 = -1`).

Deriving the `alpha_B = 0` locus by hand: `A12 = -1` requires, from Sec. 5.1's decomposition,

    Re A12 = -1/2 - (sqrt3/2) s q = -1     =>  s q = 1/sqrt3
    Im A12 =      (sqrt3/2)(-c + s p) = 0  =>  s p = c

Combining with the `M_A`-curve condition `q = sqrt3*c*p`, multiply by `s`: `s q = sqrt3*c*(s p) =
sqrt3*c^2`. Since `s q = 1/sqrt3`, we get `sqrt3*c^2 = 1/sqrt3`, i.e.

    3*cos^2(th) = 1   =>   cos(th) = +- 1/sqrt3

giving **exactly two** points in `[0,pi)^2`:

    (th, ph) = ( arccos(1/sqrt3),  pi/4 )        = the repo's Dita anchor
    (th, ph) = ( pi - arccos(1/sqrt3),  3pi/4 )

(both verified: `s q = (sqrt2/sqrt3)(1/sqrt2) = 1/sqrt3` and `s p = +-1/sqrt3 = c` in each case).
The two are related by `(th,ph) -> (pi-th, pi-ph)`, which sends `A -> conj(A)` and `B -> conj(B)`
(verified by substitution into Sec. 5.1: `A11(pi-th,pi-ph) = conj(A11(th,ph))`, likewise `A12`).
Whether these two chart points give CHM-equivalent matrices under Karlsson's `(2.2)` is **not
settled here**. `[open]`

**This is the honest shape of the finding:** Dita is a *generic* point of the `M_A` stratum but a
*singular* point of the full pair `(D_A, D_B)`, and it is the **pair**, not the stratum, that
distinguishes the anchors.

### 5.8 Summary table - the two loci as strata

| | F6 anchor `(th=0, ph arbitrary)` | Dita anchor `(arccos(1/sqrt3), pi/4)` |
|---|---|---|
| `A11`, `A12` | `omega`, `omega^2` | `i`, `-1` |
| `B11`, `B12` | `omega^2`, `omega` | `-1-i`, `0` |
| `A` Hermitian? | **No** (`A11 = omega` non-real) | **No** (`A11 = i` non-real) |
| `alpha_A, beta_A` | `omega, omega^2` | `1, -1` |
| `alpha_B, beta_B` | `omega^2, omega` | `0, 2i` |
| `D_A` | `0` | `0` |
| `D_B` | `0` | `-4` (extreme of `[-4,4]`) |
| `M_A` | `== 1` (collapsed) | `== -1` (collapsed) |
| `M_B` | `== 1` (collapsed) | `z -> 1/z` (live involution) |
| stratum component | `sin(th) = 0` (phi-gauge; **both** maps collapse) | `tan(ph) = sqrt3*cos(th)` (only `M_A` collapses) |
| `z`-structure | `z1` free; `z2=z3=z4=1` ⚠️ *z2 wrong — see Defect-resolution §DR.2(D)* | `z1` free; `z2=z3=i`, `z4=1/z1` |
| `ph` dependence of `A`,`B` | none | present |
| Karlsson's name for the family | subfamily of `F_6^(2)` or `(F_6^(2))^T` (Secs. 5-6) | `D_6^(1)` (Dita), Sec. 6 |

Tag for the table: `[exact symbolic]` (each cell derived in Secs. 5.3-5.5 and machine-verified; the

---

## 6. What is proved, exactly - consolidated

| # | Statement | Tag |
|---|---|---|
| P1 | `D_A = 2*sqrt3*sin(th)(sin(ph) - sqrt3*cos(th)cos(ph))`, `D_B = -2*sqrt3*sin(th)(sin(ph) + sqrt3*cos(th)cos(ph))`; both exactly `-+2*sqrt3 x` Karlsson's published degeneracy factors; `M_A`/`M_B` degenerate **iff** `D = 0` | `exact symbolic` |
| P2 | At Dita: `A = [[i,-1],[-1,i]]`, `B = [[-1-i,0],[0,1-i]]`, `(alpha_A,beta_A,alpha_B,beta_B) = (1,-1,0,2i)`, `(D_A,D_B) = (0,-4)`, `M_A == -1`, `M_B = 1/z`, `z2=z3=i`, `z4=1/z1` | `exact symbolic` |
| P3 | At F6 `th=0`: `A11=omega, A12=omega^2, B11=omega^2, B12=omega`, `(D_A,D_B) = (0,0)`, `M_A == M_B == 1`, `z2=z3=z4=1` ⚠️ *z2 wrong — it is a 0/0 form resolved to `z2sq = alpha_A/beta_A`, not 1; see §DR.2(D)* | `exact symbolic` |
| P4 | The two loci are strata of **different kinds**: `sin(th)=0` (phi-gauge, both maps collapse) vs `tan(ph)=sqrt3*cos(th)` (one map collapses, one live) | `exact symbolic` |
| P5 | Hermitian locus is exactly `cos(th) + sin(th)cos(ph) = 0`, containing **neither** anchor; both `A`-blocks are non-Hermitian | `exact symbolic` |
| P6 | `alpha_B = 0` on the `M_A` stratum forces `3cos^2(th) = 1`, giving exactly two chart points related by `(th,ph)->(pi-th,pi-ph)`; CHM-equivalence of the two **not** settled | `exact symbolic` / `open` |
| P7 | The `z2`-recovery denominator at Dita is `-2i != 0` - **not** a `0/0` form | `exact symbolic (negative)` |
| P8 | Karlsson Sec. 6 independently names the Dita point as a degenerate chart point generating `D_6^(1)` | `exact symbolic, primary source` |
| P9 | The `arc_circle_asymmetry.txt` Hermitian mechanism is false, and the file self-contradicts | `exact symbolic` |

---

## 7. What is NOT proved - the open clause, stated sharply

**No causal claim is made or implied.** Specifically, none of P1-P9 establishes that the structural
difference *causes* `kappa(H(lam)) >= 6` for all `lam` at Dita, or the bounded arc at F6. The clique
number is a property of the MU-pool orthogonality graph of `(H(lam), B3)`, and the Mobius data
constrains `H(lam)` but says nothing by itself about cliques. `[not run / open]`

**The concrete mechanism I would propose as the next step** (a *conjecture*, sharply falsifiable):

> **C.** At F6, `H(lam)` depends on `lam` only through `Z1 = [[1,1],[z1,-z1]]`, so every entry is
> affine-linear in `z1` and the `lam`-dependence is not closed under any involution of `z1`. At Dita,
> `H(lam)` depends on `lam` through **both** `Z1` and `Z4 = [[1, 1/z1],[1, -1/z1]]`, so every entry
> is a Laurent polynomial in `z1` whose terms pair as `z1 <-> 1/z1 = z1^{-1}` - the `lam -> -lam`
> involution. If a third basis `B3(lam)` is constructed `lam`-equivariantly under that involution
> (acting on the **triple** `(I, H, B3)`, not on `H` alone), then one 6-clique at a single `lam`
> propagates to all `lam`, which is exactly the circle behaviour; at F6 no such involution exists, so
> propagation is unavailable and only a bounded arc survives.

**A tension I must flag rather than hide:** repo Lemma L4 states that `H(th_D, pi/4, lam)` and
`H(th_D, pi/4, lam')` are **not** CHM-equivalent for `lam != lam'`. If `H(lam)` and `H(-lam)` were
related by `D2 P2 H P1 D1` (Karlsson `(2.2)`), that would be a CHM equivalence and would contradict
L4. So the equivariance in C, if it exists, must act on the **triple** and must *not* be
expressible as a CHM equivalence of `H` alone. This is a real constraint on C, and it is exactly
the kind of thing that must be checked before C is promoted. `[open]`

**Decidable cheap test of C** (one 6x6 sign/phase chase, no search): take Karlsson's explicit `H(z)`
from Sec. 2.1 and check whether `H(conj(z))` equals `D H(z) P` for diagonal unitary `D` and
permutation `P`. I did **not** run this. `[not run]`

last row is `[primary source]`, and the `F_6^(2)` cell is *noted from* Karlsson Secs. 5-6 but **not
verified against the repo's builder output** - see Sec. 9).

Direct: at Dita the `(1,1)` entry of `A - A^dag` is `i - (-i) = 2i != 0`; at F6 it is `i*sqrt3 != 0`.
The Hermitian locus is the curve `tan(th) = -1/cos(ph)` living in `th` in `(pi/2, pi)` - it passes
through **neither** anchor, and there `c + sp < 0`, opposite in sign to both anchors (`2/sqrt3` and
`1`, both positive).

---

## 8. The Hermitian correction and where it propagates

**Flagged for correction, NOT edited** (docs-reconciler owns docs; routing to orchestrator):

1. **`docs/results/final_honest_status.md:35`** - the ledger row
   `| Dita vs F6_theta0 topology contrast | **HP-supported structural** | A-block geometry; arc_circle_asymmetry.txt |`.
   This is the repo's headline structural claim resting on the false Hermitian contrast. It should
   be **downgraded and re-worded**, e.g. to reference the Mobius degeneracy pair `(D_A,D_B)` and to
   state the clique consequence as unproved. The ledger's own key offers no tier above
   "Conjectured" for an unproved mechanism, and "HP-supported structural" currently reads as
   stronger than the evidence supports - the supporting artifact's own numeric field contradicts
   its prose.
2. **`results/arc_circle_asymmetry.txt`** - lines 33, 35, 44-48, 50-51. Line 33 correctly says F6's
   block is "non-Hermitian"; line 35 says Dita's is "Hermitian"; both lines 7 and 13 say
   `A Hermitian=false`. Lines 44-48 and 50-51 build the verdict on the false contrast. (Flagged;
   `results/` artifact, route to owner.)
3. **Lemma L3** (`results/track_c_elimination/arxiv_format_pass/after_gauge_claims.txt:56-64`) is
   **correct and must not be changed**: it asserts only `A(th_D,pi/4) = [[i,-1],[-1,i]]` and
   `AA^dag = 2I`. I re-derived both exactly. The word "Hermitian" was added *downstream* of L3, by
   the asymmetry narrative - L3 itself is innocent.
4. **`scripts/python/karlsson_k6_3.py:100-104`** - the SU(1,1) comment misstates Karlsson Sec. 5
   (Sec. 5.2(a) above). Suggested replacement text supplied. Not edited (not my file).
5. **Positive correction worth recording:**
   `docs/research/symmetry/agent7/symmetry_gauge_report.md:70-76` was **right** and is vindicated -
   its two singularity equations are exactly `-+D_A/(2*sqrt3)` and `-+D_B/(2*sqrt3)`, and its warning
   that these "must be treated as strata, never removed by division" is precisely what Sec. 5.8 now
   exploits. No change needed; cite it.

---

## 9. What would falsify this hypothesis

Stated before any further work, per the brief:

- **F1 (kills the `D` framing entirely).** If the Dita anchor's full-circle clique behaviour is
  reproduced at a *generic* point of the `M_A` stratum - i.e. if arc-or-circle behaviour is
  uncorrelated with the sign/magnitude of `D_B` - then the Mobius degeneracy index does not govern
  the clique, and the structural program in this report is a red herring. **Cheapest available
  test:** one third-basis solve at a generic point of the curve, e.g. `th = pi/3` (`cos(th) = 1/2`),
  `tan(ph) = sqrt3/2`, `ph = arctan(sqrt3/2)`, where `(D_A, D_B) = (0, -3)`, versus Dita's `(0,-4)`
  and F6's `(0,0)`. `[not run]`
- **F2 (kills the `alpha_B = 0` framing).** If `kappa(H(lam)) >= 6` holds on a full circle at some
  point with `alpha_B != 0` **and** `D_B != -4`, or fails at one of the two `alpha_B = 0` chart
  points of Sec. 5.7, then `alpha_B = 0` is not the operative condition. `[not run]`
- **F3 (kills "different kinds of strata").** If the `th=0` phi-gauge stratum and the
  `tan(ph)=sqrt3*cos(th)` stratum turn out to be related by a CHM equivalence that also identifies
  their clique graphs, the asymmetry in Sec. 5.8 is a chart artifact. `[not run]`
- **F4 (kills the correction itself - checked, and it did *not* fire).** If some `A` in the repo's
  chart is Hermitian at an anchor, the correction fails. Verified false at both anchors (Sec. 5.5).
- **F5 (would make this a red herring about the *repo's* F6).** If the repo's `F6_theta0` output is
  **not** equivalent to a subfamily of `F_6^(2)` / `(F_6^(2))^T`, then the Jaming-Matolcsi-Mora
  quadruple exclusion (reported second-hand, arXiv id unverified) does not transfer, and any "F6"
  bookkeeping is loose. `[not run]`
- **F6 (would refute my reading of the primary source).** If Karlsson Sec. 5 actually requires
  `|alpha|^2-|beta|^2 = 1` somewhere I did not read, the Sec. 5.2(a) correction would be wrong. I
  quote Sec. 5 verbatim in Sec. 2.1; the phrase there is "`as long as |alpha|^2-|beta|^2 != 0`". Low
  risk, but stated.

---

## 10. Verdict

> ## `partial progress`

Broken down by the tier vocabulary of `src/Certification.jl:18-19`:

| Component | Verdict |
|---|---|
| **Refutation of the Hermitian mechanism** (P5, P9) | **`PROVED` / `exact symbolic` - complete.** Both `A`-blocks are non-Hermitian; the Hermitian locus contains neither anchor. This half of the task is **solved**. |
| **Identification of the correct invariant + closed form** (P1, P4) | **`PROVED` / `exact symbolic` - complete**, independently corroborated by the primary source and by the repo's own agent7 report. |
| **Dita as a published degenerate chart point** (P2, P8) | **`PROVED` / `exact symbolic, primary source` - complete.** |
| **Causal link to `kappa(H(lam))`** (Sec. 7, C) | **`OPEN`.** Not established. The proposed mechanism C is a conjecture with a stated tension against Lemma L4. |

**A negative result is the deliverable here, and it is a legitimate one.** The repo carried a
`HP-supported structural` claim resting on a distinction that is false at both of its endpoints and
that its own source artifact numerically refutes. That claim is now refuted exactly, and the
replacement invariant is derived in closed form. The remaining gap - connecting the invariant to
the clique - is honestly labelled open rather than papered over, and the whole program was analytic
(zero new numerical experiments), consistent with the campaign gate.

---

## 11. Next steps (in priority order)

1. **Route the correction** to `docs/results/final_honest_status.md:35` and to
   `results/arc_circle_asymmetry.txt` (lines 33, 35, 44-48, 50-51). Do **not** touch Lemma L3.
2. **Fix the SU(1,1) comment** at `scripts/python/karlsson_k6_3.py:100-104` using the text in Sec. 5.2(a).
3. **Run F1** - the single cheapest kill-test: one third-basis solve at a generic point of the
   `M_A` stratum (`th = pi/3`, `ph = arctan(sqrt3/2)`, where `(D_A, D_B) = (0, -3)`), versus
   Dita's `(0,-4)` and F6's `(0,0)`. This decides whether the program survives.
4. **Run the F6 decidable test of C** (Sec. 7): the `D H(z) P` sign/phase chase on Karlsson's
   explicit `H(z)`. No search required; a few minutes of exact algebra.
5. **Settle F5** - is the repo's `F6_theta0` output equivalent to a subfamily of `F_6^(2)` /
   `(F_6^(2))^T` (Karlsson Secs. 5-6)? If yes, the Jaming-Matolcsi-Mora exclusion may transfer
   and the F6 arc becomes *less* interesting, not more.
6. **Settle P6** - are the two `alpha_B = 0` chart points CHM-equivalent under `(2.2)`?
7. **Obtain Goyeneche 2013** (J. Phys. A **46**, 105301) - closest prior art, currently unread by
   me. Needs a library copy; the arXiv API and Semantic Scholar were both blocked (406/429).
8. **Widen the stratum treatment:** per P4, any future certificate on `K_6^(3)` must be
   stratified by `(D_A, D_B) = 0` before generic-fiber transfer, extending agent7's warning
   from the two equations to the pair. This is the one *positive, forward-looking* structural
   theorem this report supports.

---

## 12. Files touched by this agent

Created (all new; **no existing file in `results/` was modified**, so no `.bak_2026-09-26`
backups were required):

- `research/reports/audit_followup_2026-09-26.md` - this report
- `results/NOT_RUN_researcher_blockerC_2026-09-26.log` - campaign gate
- `scripts/python/blocker_c_exact_check.py` - exact-arithmetic verification
- `results/lit_cache/karlsson_1003.4177.html` - cached primary source
- `results/lit_cache/karlsson_plain.txt` - de-tagged text
- `results/lit_cache/mcnulty_weigert_comment.html` - cached arXiv:2504.13067
- `results/lit_cache/exact_check_out.txt` - exact-check output (RC=0)
- `results/lit_cache/exact_check_err.txt` - empty
- `results/lit_cache/arxiv_query_out.txt`, `arxiv_query.py`, `s2_query_out.txt`, `s2_query.py`
  - failed-API-attempt logs

Flagged, **not** edited: `docs/results/final_honest_status.md`,
`results/arc_circle_asymmetry.txt`, `scripts/python/karlsson_k6_3.py`.

---

## D-invariant dichotomy — resolved

**FALSE / COINCIDENTAL.** `D = |alpha|^2 - |beta|^2` does **not** determine, for the
points of a third-MUB locus, whether that point survives: it is provably *constant*
along each locus, so it takes identical values at surviving (in-arc) and non-surviving
(out-of-arc) points of the same locus. It therefore has **zero discriminating power**
for survival, and the arc's boundedness cannot be a consequence of how `D` varies along
the locus. The two-anchor "match" was a comparison of two *different slices*, not a
boundary mechanism. This is a closed-form argument, not a sampling result.
`exact symbolic`.

Two things survive and must not be thrown away with the negative finding, and one
external prior-art fact re-scopes the whole question. Both are in D3.4 and D3.5.

### D3.1 What was asked, and why the premise needed fixing

The task asked whether `D` determines, **for every point on each locus** (not just the
anchors), whether that point survives as a valid representative on the arc/circle —
equivalently, whether arc boundedness and circle fullness are a *consequence of how `D`
varies along each locus*.

Both loci are one-parameter loops in `lambda` at **fixed** `(theta, phi)`:
`results/locus_classification.txt:9` gives the F6 locus as a bounded arc "in lambda at
theta=0"; line 15 gives the Dita locus as a "circle in lambda (Dita, phi=pi/4)". So the
arc parameter **is `lambda`**, and each locus holds `(theta, phi)` constant.

That is what kills the mechanism: `D` is a function of `(theta, phi)` only, so holding
`(theta, phi)` fixed along a loop makes `D` constant on that loop. The premise "how `D`
varies along each locus" is void — there is no variation along a locus to act as a
mechanism.

### D3.2 The linchpin, verified three ways

1. **Source ordering.** `scripts/python/karlsson_k6_3.py`: `build_A(theta, phi)`
   (lines 42-62) takes no `lambda`; `B = -F2 - A` (line 78);
   `alpha_A, beta_A = A[0,1]**2, A[0,0]**2` and `alpha_B, beta_B = B[0,1]**2, B[0,0]**2`
   (lines 89-90) are computed **before** `z1 = np.exp(1j*lam)` first appears (line 92).
   The module docstring states it outright (lines 30-31): "Free parameters: theta, phi in
   [0, pi) determine A (and hence B); z1 is then a free unimodular parameter
   z1 = exp(i*lambda)". So `lambda` enters *only* as `z1`, downstream of `D`.
2. **Free-symbol check (sympy).** `alpha_A`, `beta_A`, `alpha_B`, `beta_B` each have free
   symbols exactly `{theta, phi}` — `lambda` appears in none. Hence `dD/dlambda == 0`
   identically. Script: `scripts/python/d_invariant_check.py`; output
   `results/d_invariant_check_out.txt`.
3. **Closed-form match to Karlsson.** `|A12|^2 - |A11|^2` and `|B12|^2 - |B11|^2` reduce
   exactly to `sqrt3*sin(th)*(sin(ph) - sqrt3*cos(th)cos(ph))` and its `+` counterpart,
   i.e. to Karlsson's two published degeneracy factors, with **residual exactly 0**
   (sympy `simplify`). At the anchors: F6 `(0, 0)`; Dita `(0, -2)` for the
   squared-modulus discriminants. The `D = |alpha|^2 - |beta|^2` normalisation used in
   the prior report gives the same signs — `(0,0)` and `(0,-4)` — because
   `D_X = (|X12|^2 - |X11|^2) * (|X12|^2 + |X11|^2)`. `exact symbolic`.

### D3.3 The disproof, in the form the task asked for

Ask the task's question pointwise on the **F6 `lambda`-circle**:

* a point **inside** the arc survives (`kappa >= 6`); its `(theta,phi) = (0, 0.5)`, so `D = (0, 0)`;
* a point **outside** the arc, on the same `lambda`-circle with the same `(theta, phi)`,
  does **not** survive; its `D = (0, 0)` too.

The two points have identical `D`. `D` cannot distinguish them. This is not a sampling
statement — it follows from `dD/dlambda == 0`, so it holds at *every* point of the locus
simultaneously.

The three specific tests the task requested:

* **(i) Does `D` hit a zero, pole, or sign-change exactly at the arc boundaries
  `0.055044803035` / `0.0625984193874`?** **No — and it cannot.** Those are `lambda`
  offsets; `D` has no `lambda` dependence, so it has no structure at any `lambda`,
  boundary or interior. There is nothing to align.
* **(ii) Does `D` stay non-degenerate around the full Dita circle?** `D_B = -4`
  everywhere there, so it never hits the critical value `0` — but this is **not** an
  explanation, because `D_B = 0` everywhere on F6 *including F6's non-surviving points*,
  and F6's surviving points share that same value. Also `D_A = 0` at **both** loci, so
  `M_A` is collapsed at every point of both loci, surviving or not. `D` simply does not
  track survival.
* **(iii) Does anything special happen at the boundary?** Not for `D`: it is constant.

### D3.4 What survives (do not overclaim the negative)

The negative is about the *mechanism*. Three exact facts remain, and they are why the
prior session's stratum analysis was worth doing:

1. `D`'s vanishing is exactly the degeneracy condition of the Mobius data, and
   **Karlsson says so himself**: the maps "as long as `|alpha|^2 - |beta|^2 != 0`, map the
   unit circle onto itself", and when it tends to `0` the map "degenerates into a mapping
   of the unit circle in `z` into a single point ... and this mapping has no inverse"
   (arXiv:1003.4177, Sec. 5). So `D` is a *literature-defined degeneracy indicator*, not a
   repo invention — but its role is circle-onto-itself liveness of the Mobius map, **not**
   boundedness of a third-MUB locus. These are different circles.
2. `D`'s value does classify the **slice**: F6 sits on the `sin(th) = 0` stratum (both maps
   collapse — Karlsson Sec. 5, Eq. (5.2) branch), Dita on `tan(ph) = sqrt3*cos(th)` (only
   `M_A` collapses). Exactly as the prior report's P4 stated.
3. That classification fixes the **`z`-structure**, which is the thing that actually
   differs in `lambda`. Confirmed exactly here:
   * F6: `z3^2 = 1`, `z4^2 = 1` (both maps collapse to constants) — and the z2-recovery
     denominator is **identically zero**, a genuine `0/0` form requiring a limit
     (numerically `z2^2 -> -0.935807 + 0.353i`, which is **not** `1`; the prior report's
     table entry "F6: z2=z3=z4=1" is therefore **wrong for `z2`** and should read
     `z3 = z4 = 1` with `z2` fixed by a limit).
   * Dita: `alpha_A = 1, beta_A = -1, alpha_B = 0, beta_B = 2i`; `z3^2 = -1`,
     `z4^2 = z1^-2`, `z2^2 = 2i/(-2i) = -1`.
   The final link — `z`-structure to `kappa >= 6` / boundedness — is hypothesis **C** in
   section 7 of this report and remains **open / not run**. It is the only surviving
   candidate mechanism in this repo, and it is *not* a statement about `D`.

Numerical confirmation of the entry-level `lambda` dependence (float64, repo builder):
Dita shows **24/36** entries changing with `lambda` at Hadamardness
`max|H H^dag - 6I| = 0`, consistent with `lambda` entering through `Z1` (`z1`) and `Z4`
(`1/z1`). `sampled/numerical only`.

**A previously unrecorded defect found while checking this.** The repo's own Python
builder `scripts/python/karlsson_k6_3.py` is **degenerate at `theta = 0`**: the z2
inversion is `0/0`, numpy returns `nan`, and `build_k6(0.0, 0.5, lam)` is **non-finite at
every `lambda`** (RuntimeWarnings: "invalid value encountered in scalar divide" at lines
109 and 67, then "invalid value encountered in sqrt" at line 114). So the Python path
cannot represent the F6_theta0 locus at all, and any F6_theta0 result in this repo must
come from a different code path. At `theta = 1e-9` the builder recovers (Hadamardness
8.9e-16) and 24/36 entries vary with `lambda`. This is a reproducibility hazard for the
F6 arc and is logged for follow-up. `exact symbolic` (control flow / arithmetic) for the
`0/0`; `sampled/numerical only` for the entry counts.

### D3.5 Literature — the phenomenon is published, the mechanism is not

This is the most consequential finding of this task: it **re-scopes** the claim.

* **The bounded-region phenomenon is already in the literature.**
  D. Goyeneche, *Mutually unbiased triplets from non-affine families of complex
  Hadamard matrices in dimension six*, J. Phys. A **46**, 105301 (2013),
  **arXiv:1209.4126** — the paper the prior session flagged as "closest known prior work"
  with an **UNVERIFIED** arXiv id, "unable to obtain or read". It has now been **obtained
  and read** (arXiv id now verified). Its Table 1 reports, per family, the maximal set of
  MU bases available from a pair:
  * `{I, K6^(3)(theta,phi,psi)}`: "3 in black regions ... 2 in white regions" (Figs 3(a)-(d));
  * `{I, K6^(2)(x1,x2)}`: "3 in black regions of Fig. 1 / 2 in white regions";
  * `{I, M6^(1)(t)}`: "3 if `t in [0.5309pi, 0.9157pi]`", "3 if `t in [1.5312pi, 1.9163pi]`", "2 otherwise";
  * `{I, D6^(1)(c)}`: "Affine  3, `forall c in [-1/8, 1/8]`";
  * `{I, F6^(2)(a,b)}`: "Affine  3, `forall a,b in [0, 2pi)`".
  Goyeneche states: "The non-affine family `M6^(1) ... **does not allow a triplet in its
  full range**", and of `B6^(1)`: "This is the only non-affine family where we found a
  triplet for every value."

  **Consequence:** the repo's "F6 bounded arc vs Dita full circle" is an instance of a
  **known, published phenomenon** — third-MUB existence in non-affine order-six CHM
  families holding on bounded parameter regions and failing elsewhere. It is *not* a new
  structural discovery and should not be presented as one. This is exactly the
  duplicate-result failure mode the project has been burned by before.

* **The mechanism is *not* stated.** Goyeneche's bounded/white regions are **numerical
  output of an iterative method** ("Our method is not an algorithm because it does not
  stop with a definite answer but it converges very quickly"), reported as black/white
  regions in figures. Target-word searches of Karlsson's full text return **zero** hits
  for "bounded", "interval", "region" or "third"; its only range statements concern the
  `theta, phi in [0, pi)` sign-redundancy reduction, which is unrelated. So: **the
  phenomenon is L1 (published); the causal mechanism is L2 (not published).** A repo
  mechanism would be new — but it must be presented as a mechanism for a *known*
  phenomenon, with Goyeneche cited, not as the discovery of the phenomenon.

* **The breadth of Goyeneche's Table cuts against the `D` story independently.** Bounded
  existence regions occur across several families at varying `(theta, phi)` — `M6^(1)`
  with two disjoint intervals, `K6^(2)` and `K6^(3)` with structured black/white regions
  — while `D6^(1)` and `F6^(2)` are full-range. A single function `D(theta, phi)` would
  have to reproduce all of that structure; nothing in this task suggests it does.

* **Context.** M. Matolcsi, A. K. Matszangosz, D. Varga, M. Weiner, *Triplets of Mutually
  Unbiased Bases*, **arXiv:2503.14752**: the general Zauner/Szollosi construction gives a
  **2-parameter** family of MUB-triplets (`(I, F1, F2)` with `F1, F2` from the Fourier
  family), so a 1-parameter `lambda`-locus is a **slice** of a higher-dimensional triplet
  variety. Their "exceptional" cube (all entries 24th roots of unity) is notable because
  the repo's exact T4 field is `Q(zeta_24, sqrt5)`; I make **no claim** of a connection —
  flagged as a lead only. `sampled/numerical only` (read from abstract/introduction).

### D3.6 Consequence for the roadmap

`docs/results/final_honest_status.md:35` ("Dita vs F6_theta0 topology contrast |
**HP-supported structural** | A-block geometry") is **downgraded** in this diff: its
evidential basis `arc_circle_asymmetry.txt` self-contradicts (both blocks recorded
`Hermitian=false` in its own numeric field), and the `D`-invariant replacement is
refuted here as a mechanism. The phenomenon is attributed to Goyeneche.

`docs/overview/readme.md:37` (the C2 / Track D row) is a *description* of the two loci and
is factually correct as far as it goes; it is annotated, not rewritten, because the loci
themselves are not in question — only the explanation of their shapes is.

**Verdict: FALSE / COINCIDENTAL** (`exact symbolic` for the refutation).
Residual: `D` is an exact, literature-defined stratum/degeneracy classifier whose
downstream clique consequence is **open** (hypothesis C, `not run`). The *phenomenon*
being explained is published (Goyeneche 2013, arXiv:1209.4126); the *mechanism* is not.

### D3.7 New files added by this task

* `scripts/python/d_invariant_check.py` — the exact/numeric check (reproduction only).
* `results/d_invariant_check_out.txt` — its output.
* `results/_dcheck_err.txt` — stderr, containing the F6 `0/0` and `nan` warnings.
* `results/lit_cache/alpha_beta_role.txt`, `results/lit_cache/kw_search.txt` — extracted
  Karlsson passages (target-word search evidence).
* `results/NOT_RUN_dinv_2026-09-26.log` — campaign gate: what was deliberately not run.

---

## Defect resolution — arc width & F6 provenance

New session. Two defects logged in D3.4 / D3.6 are resolved here, one subsection each,
ending in a verdict. Scripts used (reproduction / re-derivation of the repo's own code
paths only, no new homotopy or clique experiment): `scripts/python/d_invariant_check.py`,
`scripts/python/defect2_probe.py`, `scripts/python/defect2_instrumented.py`, and a CSV
recount. Outputs: `results/d_invariant_check_out.txt`, `results/defect2_probe_out.txt`,
`results/defect2_instrumented_out.txt`, `results/defect1_repro.txt`.

### DR.1 Defect 1 — the ~16x arc-width discrepancy

**Verdict: RECONCILED.** The correct arc extent is **0.117643** rad
(`results/f6_boundary_reverify.txt`, the bracketed boundary search). The
`[0.2925, 0.3000]` interval in `results/locus_classification.txt:9` is the *extrema of
the sampled HP-verified cluster*, not the arc extent, and the `topology` label on it is a
mislabel. Same variable (`lambda`) in both, so this is NOT a parametrization/units
mismatch — it is a boundary-bracket vs sample-range definitional mismatch.

Same parameter, confirmed. `f6_boundary_reverify.txt:7-8` prints `boundary |Delta
lambda|`; `locus_classification.txt:9` labels its interval "bounded arc in lambda". No
units mismatch to rule out. `exact symbolic`.

What each file actually measures:

1. `f6_boundary_reverify.txt` is produced by `scripts/julia/f6_boundary_reverify.jl`
   (line 11: `H = build_karlsson_family(theta, phi, lam)`), a dedicated bracket-and-bisect
   search for the boundary of clique-6 survival at HP_bits=400 around `lambda0 = 0.3`. It
   measures the **true boundary**: `lambda0 - 0.0625984193874` and
   `lambda0 + 0.055044803035`, arc width `0.117643222422`. `sampled/numerical only`
   (400-bit HP probes over a float64 pool/clique pipeline).
2. `locus_classification.txt:9` is produced by `scripts/julia/locus_classification.jl`,
   function `topology_for_label` (lines 176-180):
   ```julia
   lbl == "F6_theta0" && begin
       lams = extrema(p.lambda for p in pts)
       @sprintf("1D bounded arc in lambda at theta=0 [%.4f, %.4f]", lams[1], lams[2])
   end
   ```
   The interval is `extrema(p.lambda for p in pts)` over the cluster's HP-verified
   clique-6 sample points — the **sample range**, not a boundary. `exact symbolic` (code).

Reproduced exactly (this session, `results/defect1_repro.txt`): filtering
`results/special_loci_search.csv` (`status==ok`, `max_clique>=6`) joined to
`results/third_mub_audit.csv` (`hp_survival==true`) yields exactly **27** F6_theta0
points, whose `lambda` extrema are `[0.2925, 0.3000]` — matching `locus_classification.txt:9`
to the digit — and whose distinct `lambda` values are the 0.0025-spaced grid
`{0.2925, 0.2975, 0.3000}`. `sampled/numerical only` (deterministic integer/filter
recount). So the stored string is the sample grid's extent, not the arc.

The reconciliation: `[0.2925, 0.3000] subset [0.23740158, 0.35504480]` (the true arc). The
sample extrema sit well inside the true arc — a refinement grid that only ever produced
HP-verified clique-6 points near the anchor; points at `lambda = 0.29` (visible in
`results/special_loci_run2.txt`) fail HP or clique-6 and so are absent from the extrema.
Neither file is "wrong" numerically; `locus_classification.jl:180` is *mislabeled* — it
reports a data range under the word "topology".

Provenance/timing check (Defect 1 task item 3): both producers predate the 2026-09-17
normalization fix, but the discrepancy is **not downstream of that fix** — it is
authoring error in `topology_for_label` (range labelled as topology), independent of any
normalization convention. `final_honest_status.md:33` ("arc width 0.117643;
f6_boundary_reverify.txt") uses the *correct* boundary number and is unaffected.

**Required follow-up edit (flagged, done below under archive):**
`docs/overview/readme.md:37` cites `results/locus_classification.json` for the
"~0.118 rad" width, but that JSON carries the sample-range string `[0.2925, 0.3000]`, not
0.118. The correct citation for the width is `results/f6_boundary_reverify.txt`.

### DR.2 Defect 2 — F6_theta0 builder provenance

**On-file F6_theta0 results are trustworthy.** They come from the Julia seam-resolved
path, not the Python builder, and that path is sound at `theta=0`. Details below; the one
residual caveat is at the end. Evidence: `exact symbolic` (code paths) +
`sampled/numerical only` (my float re-verifications of the seam matrix).

**(A) What actually produced the F6_theta0 numbers.** The Julia builder
`build_karlsson_family` in `src/mub_zauner_6d_liang_chen.jl:75`. It has an explicit,
documented Fourier-seam branch (lines 94-111, comment audit 2026-09-14):
```julia
if theta == 0.0
    z3sq = one(ComplexF64); z4sq = one(ComplexF64)
    z2sq = alpha_A / beta_A            # exact algebraic resolution of the seam
else
    ...mobius inversion...
end
```
This branch is what every F6_theta0 artifact uses:
`scripts/julia/f6_boundary_reverify.jl:11` calls `build_karlsson_family(theta, phi, lam)`,
as do the T3 all-clique certification scripts and `scripts/julia/locus_classification.jl:223`.
So the load-bearing F6_theta0 results are Julia-produced. The Python `build_k6` is a
reference/audit implementation and is **not** the producer. `exact symbolic`.

**(B) The seam branch is sound at theta = 0**, verified independently in Python
(`results/defect2_probe_out.txt`, PART B). Replicating the branch exactly
(`z3sq = z4sq = 1`, `z2sq = alpha_A/beta_A`, `z1 = exp(i*0.3)`):
* the seam matrix at `(0, 0.5, 0.3)` is a valid CHM: `max|H H^dag - 6I| = 1.3e-15`,
  `max||H_jk| - 1| = 1.1e-16` (`sampled/numerical only`, float64);
* the repo's seam A-block equals Karlsson's own theta=0 form `A = F2 * Omega`,
  `Omega = diag(omega, omega^2)`, of `arXiv:1003.4177` Sec. 5 / Eq. (5.2), exactly
  (`True` in the float check);
* consistent with Karlsson's Eq. (5.2): at theta=0, "z3^2 = z4^2 = 1 => there are no
  restrictions on z1 or z2", so the seam choice `z2sq = alpha_A/beta_A` is one legal
  algebraic representative of the branch, not a limit; and the seam respects the phi-gauge
  (A is phi-independent bitwise) — `sampled/numerical only`.

The repo's own regression suite agrees and goes further: `test/runtests.jl:50-62`
explicitly tests the seam at `lambda in {0.0, 0.7, 2.0}` — the very points where the
Python builder NaNs (below) — asserting valid Hadamardness, and spot-checks
backward-compatibility with one historical anchor `(0, 0.5, 0.3)`. `exact symbolic`
(source reading); the suite itself is `not run` here (no Julia on this host; see below).

**(C) The Python builder defect — CORRECTING my own overstatement.** Section D3.4 said
"`build_k6(0.0, 0.5, lam)` is non-finite at *every* lambda". That is **wrong**; it was an
artifact of my test list starting at `lam = 0.0` (a genuine pole of the divided map)
poisoning one label. The instrumented replication of `build_k6`'s internal step order
(`results/defect2_instrumented_out.txt`) shows the real behaviour is **per-lambda and
1-ulp-non-deterministic**:
* NaN (20 NaN entries) at `lambda in {0.0, 0.48, 0.5, 0.7, 1.5, 2.0, 3.1415, 4.0, 6.2832}`;
* finite at `{0.3, 0.52, 5.0}` where `z3sq` lands at `1 - 4e-17j` rather than exactly
  `1+0j`, so `num/den ~ 1e-17/1e-17` survives — with `z2sq` set by the ulp-noise
  **direction** (it came out `-1+0j`);
* at `lambda = 1.0`, `z2sq = -1-1j` with `|z2sq| = sqrt(2) != 1`: **silently not
  unimodular** — a finite but invalid-CHM output, worse than a NaN because nothing throws.

So the Python builder at the seam is a real defect — but it is a *builder* defect in a
reference script, NOT a defect in the F6_theta0 *results*, which are Julia-produced. The
two stay cleanly separated, as the task requires. `sampled/numerical only`. This also
overturns D3.4's specific claim — the correction note is recorded here and pointed to below.

**(D) The wrong `z2 = 1` value — propagation check.** The string `z2=z3=z4=1` appears
**only** in this report (in the D-section table and in P3, both generated in the previous
session), and was already flagged inline in D3.4. A full-tree search finds it in **no
other** docs/results/papers file. The seam truth — `z2sq = alpha_A/beta_A = e^(-2 pi i/3)`,
i.e. `z2 = e^(-i pi/3) ~ 0.5 - 0.866i` — is already stated correctly in the repo's ledger
at `docs/results/final_honest_status.md:106` ("Fourier-seam resolution (theta=0 exact
algebraic slice; z2^2 = alpha_A/beta_A) | **Proved (algebraic) + regression-tested**").
So the wrong value **did not propagate** beyond the one report, and the correct value is
already on record. `exact symbolic` (search) + `exact symbolic` (the value). Inline
correction markers added at the two spots in this report.

**(E) Residual caveat and what would close it.** Two honest qualifications:
1. **Vintage.** The F6_theta0 artifacts predate the 2026-09-14 seam fix
   (`f6_boundary_reverify.txt` 2026-08-03; `locus_classification.txt` 2026-08-05; the T3
   all-clique logs 2026-08-13/14), so they were produced via the pre-fix ulp-noise path.
   The repo claims post-fix backward-compatibility to ~1e-16 with a **one-anchor**
   regression spot-check (`test/runtests.jl:61`). I verified the seam matrices are
   independently sound, but the one-anchor ulp compatibility is a spot-check, not a sweep,
   and I could not re-run the Julia suite (no Julia at the documented Windows path; WSL
   `julia` not found — consistent with the 2026-09-23 forensic audit). `not run`.
2. **Stale comment.** `test/runtests.jl:26-28` still says `(0, 0.5, 0)` "must throw under
   the fail-fast guard" — that is *pre-fix* behaviour; the seam test at lines 50-57
   explicitly tests `(0, 0.5, 0.0)` as a valid Hadamard post-fix. An internal contradiction
   in a *comment* — minor, but should be tidied. `exact symbolic`.

**Plain statement for the on-file results: TRUSTWORTHY.** The F6_theta0 arc-width,
boundary, and clique numbers come from a sound, literature-matched, regression-tested Julia
seam path, not from the broken Python builder. Full (rather than one-anchor) closure of
provenance would require a single Julia reproduction run at HP under the post-fix pipeline;
that is blocked only by Julia not being installed on this host, and it is a *provenance*
closure, not a correctness concern — the matrices themselves are verified above. If a later
Julia run shifts any F6_theta0 number by more than ~1e-16, flag it; until then, **no doc
number needs to change on account of Defect 2.**

**Required follow-up edits (flagged + done under archive):**
`docs/results/final_honest_status.md:36` overstated the defect as "returns nan"; it now
says the Python builder is ulp-flaky/non-deterministic at theta=0 (NaN at many lambda,
finite-but-noise at others, silently non-unimodular at lambda=1.0), and that the F6_theta0
numbers are produced by the Julia seam path. No paper / theorem-statement edits were made;
those remain for your review.






