# Exact incidence systems for Karlsson one-parameter slices

## Scope and exact parameter

Use an exact unnormalised Hadamard \(H(z,t)\), where \(z=e^{i\lambda}\),
\(t=z^{-1}\), and impose \(zt-1=0\). The incidence equations below apply to
either of these two exact specializations.

For the Diţă-circle problem, use the canonical Diţă matrix
\(D(z,t)\), specified by

\[
D(z,t)=\begin{pmatrix}
1&1&1&1&1&1\\
1&-1&i&-i&-i&i\\
1&i&-1&iz&-iz&-i\\
1&-i&-it&-1&i&-it\\
1&-i&-it&i&-1&it\\
1&i&-i&-iz&iz&-1
\end{pmatrix}.
\]

This is the same family as `dita_D(x)` in `src/brierley_weigert_notes.jl`,
with \(z=e^{2\pi i x}\).

For the selected F6 \(\theta=0\) endpoint problem, use the exact Karlsson
specialization. Let \(\omega=e^{2\pi i/3}\),
\[
A=\begin{pmatrix}\omega&\omega^2\\\omega&-\omega^2\end{pmatrix},
\quad B=-F_2-A,\quad
Z_L(u)=\begin{pmatrix}1&1\\u&-u\end{pmatrix},\quad
Z_R(u)=\begin{pmatrix}1&u\\1&-u\end{pmatrix},
\quad F_2=\begin{pmatrix}1&1\\1&-1\end{pmatrix}.
\]
The exact seam convention in `src/mub_zauner_6d_liang_chen.jl` gives
\(z_1=z,\ z_2=\omega,\ z_3=z_4=1\): at \(\theta=0\),
\(\alpha_A/\beta_A=\omega^2\), and the principal square root is \(\omega\).
Thus the target matrix is
\[
H_{F6}(z,t)=
\begin{pmatrix}
F_2&Z_L(z)&Z_L(-\omega)\\
Z_R(1)&\tfrac12Z_R(1)AZ_L(z)&\tfrac12Z_R(1)BZ_L(-\omega)\\
Z_R(1)&\tfrac12Z_R(1)BZ_L(z)&\tfrac12Z_R(1)AZ_L(-\omega)
\end{pmatrix}.
\]
It is independent of \(\phi\) at \(\theta=0\). Its coefficients are in
\(\mathbb Q(\omega)\), with parameter \(z,t\). This is the exact family to
use when studying the reported bounded arc; the Diţă specialization is a
separate incidence problem. If another Karlsson representative is used,
record and verify the exact transformation of the pair of bases rather than
substituting floating-point matrix entries into an exact ideal.

## Variables and physical real locus

Use ordered columns \(v_a\), \(a=1,\ldots,6\), for the third basis. Fix each
column's global phase by writing

\[
v_a=6^{-1/2}(1,x_{a1},\ldots,x_{a5})^T,\qquad
y_{aj}=x_{aj}^{-1}\quad (j=1,\ldots,5).
\]

For the possible fourth-MU vector write

\[
w=6^{-1/2}(1,q_1,\ldots,q_5)^T,\qquad
r_j=q_j^{-1}.
\]

For \(J_3\), use the complex algebraic ambient ring
\[
\mathbb Q(i,\omega)[z,t,x_{aj},y_{aj}],
\]
with the real structure \(i\mapsto-i,\ \omega\leftrightarrow\omega^2\),
\(z\leftrightarrow t\),
\(x_{aj}\leftrightarrow y_{aj}\). For \(J_4\), adjoin
\(q_j,r_j\) and extend the real structure by \(q_j\leftrightarrow r_j\).
The physical locus is its fixed locus: \(t=\bar z\), \(y_{aj}=\bar
x_{aj}\), and \(r_j=\bar q_j\). Equivalently, use real and imaginary parts
and unit-circle equations. The reciprocal equations alone define a complex
torus; they do **not** mean that arbitrary complex points are physical.

In formulas below let \(x_{a0}=y_{a0}=q_0=r_0=1\). Write \(h_{jk}\) for the
entries of the selected \(H(z,t)\), and \(h^\star_{jk}\) for their image
under the coefficient/parameter involution.

## Generators

Let \(J_3\) denote the third-basis incidence ideal and \(J_4\) its extension
by a fourth-MU-vector witness. The first three groups of generators define
\(J_3\); the final group extends it to \(J_4\).

**Parameter and third-basis phase constraints**
\[
zt-1,\qquad x_{aj}y_{aj}-1\ (1\le a\le6,\ 1\le j\le5).
\]

**Unbiasedness of every third-basis vector to \(H(z,t)\)** For each
\(a=1,\ldots,6\) and \(k=0,\ldots,5\), impose
\[
\left(\sum_{j=0}^5h^\star_{jk}x_{aj}\right)
\left(\sum_{j=0}^5h_{jk}y_{aj}\right)-6=0.
\]
The right-hand side is \(6\), since \(H\) is unnormalised and the vectors
have entries scaled by \(1/\sqrt6\). This is the polynomial form of
\(\left|\sum_j\overline{h_{jk}}x_{aj}\right|^2=6\). It is not
\(6(H^\ast v_a)_k(H^T\bar v_a)_k-1=0\); that expression has a factor-of-six
normalisation error with this convention.

**Orthonormality of the third basis** For each \(1\le a<b\le6\), impose both
\[
1+\sum_{j=1}^5 y_{aj}x_{bj}=0,\qquad
1+\sum_{j=1}^5 x_{aj}y_{bj}=0.
\]
On the physical locus these are complex conjugates of one another; retaining
both gives the polynomial complexification of the real and imaginary parts
of the orthogonality condition. Six vectors are required: five pairwise
orthogonal flat vectors alone do not establish a third basis.

**Fourth-vector witness equations** Add
\[
q_jr_j-1=0\quad(j=1,\ldots,5),
\]
and impose
\[
\left(\sum_{j=0}^5h^\star_{jk}q_j\right)
\left(\sum_{j=0}^5h_{jk}r_j\right)-6=0\quad(k=0,\ldots,5),
\]
and, for each \(a=1,\ldots,6\),
\[
\left(1+\sum_{j=1}^5y_{aj}q_j\right)
\left(1+\sum_{j=1}^5x_{aj}r_j\right)-6=0.
\]
Together with flatness, the last equations say that \(w\) is unbiased to
each vector of the third basis. Thus the physical points of \(V(J_4)\) are
exactly witnesses \(w\) extending some third basis at the selected
parameter.

The displayed generating set has 62 variables and 97 equations for \(J_3\),
and 72 variables and 114 equations for \(J_4\), before
removing proven dependencies. The count is deliberately not interpreted as
a dimension or a complete-intersection claim. The six MU equations for a
flat vector against a unitary basis, and the six witness equations against
an orthonormal basis, each have a sum identity; remove one only after proving
the identity in the chosen normalization. Keep the redundant equations for
verification if useful, but do not use their raw Jacobian rank as a
codimension test.

## Symmetry reduction

1. The six independent column phases have already been removed by setting
   coordinate zero to \(1\). This gauge is valid because every vector is
   flat and hence has a nonzero first coordinate.
2. The third-basis column relabeling is a finite \(S_6\) action. It is not
   removed by polynomial gauge equations. Either retain ordered bases or
   enumerate/canonicalize this finite action with an explicit, exact
   lexicographic real-algebraic rule.
3. Further row/column permutations, diagonal phases, conjugation, or
   parameter maps may be quotiented only when an exact transformation is
   verified to preserve the pair \(\{I,H(z,t)\}\) and map the full witness
   constraints to themselves. The symmetry group is the stabilizer of this
   pair, not the full Hadamard equivalence group.
4. Saturate by any denominators introduced by a chosen chart. The present
   phase chart has no zero-coordinate denominator beyond the coordinates
   already set to \(1\).

## Elimination and real certification

First compute the exact symmetry stabilizer and its action on the selected
incidence ideal. For the F6 arc endpoints, the relevant projection is
\(\pi_3:V(J_3)_{\rm phys}\to S^1\); the fourth-vector ideal \(J_4\) answers
the distinct extension question and has its own projection \(\pi_4\). Then
use separate elimination orders: for \(J_3\), eliminate \((x,y)\) before
\((z,t)\); for \(J_4\), eliminate \((q,r)\), then \((x,y)\), before
\((z,t)\). Retain factorization and saturation data. A nonzero elimination
polynomial in \(z\) can reduce a continuum question to finitely many
candidate parameters; it does not by itself prove that any candidate has,
or lacks, a physical witness.

For a no-witness theorem on a parameter arc, certify the **real** projected
incidence set. Valid endpoints include:

* a unit ideal over the exact coefficient field (which excludes complex and
  therefore physical solutions);
* exact real quantifier elimination / a Positivstellensatz certificate; or
* complete real-root isolation of a zero-dimensional reduction, with exact
  checks of the eliminated fibers and interval bounds for the parameter.

If elimination leaves a positive-dimensional variety, complex Gröbner
elimination alone does not decide whether its physical real locus is empty.
The existing fixed-\(H,B_3\) \(W_1\) computation is a different problem: it
does not eliminate over all third bases or over the parameter.

## Critical locus and fold certificate

The word “fold” must refer to a specified projection and branch, not just a
numerically ill-conditioned solve. For the projection
\(\pi:V(J_3)_{\rm phys}\to S^1\) at the F6 arc, first determine the exact dimension and
local codimension of each relevant component after symmetry reduction.
On its smooth locus, a critical point is where the differential of \(\pi\)
restricted to the incidence variety loses rank. For a local complete
intersection with \(c\) independent equations \(F_1,\ldots,F_c\) and fiber
coordinates \(u\), this is detected by
\[
\operatorname{rank}\left(\frac{\partial(F_1,\ldots,F_c)}{\partial u}\right)<c.
\]
With redundant generators, use a minimal local presentation or the
cotangent-space ideal; do not append a determinant of an arbitrary
overdetermined Jacobian.

Separate the incidence singular locus from smooth critical points. Determine
the dimension \(d\) of the source component. For a map from that component
to the one-dimensional parameter circle, a fold critical point has
\(d\pi=0\) and a nondegenerate Hessian of a local parameter coordinate
restricted to the \(d\)-dimensional tangent space (the one-dimensional
branch case reduces to a nonzero second derivative). Rank loss alone proves
neither a fold nor that no other branches exist. The certificate must
include:

1. exact equations for the critical locus and its singular sublocus;
2. exact isolation of all critical parameter values on the target arc;
3. a fold nondegeneracy check at each claimed fold;
4. continuation/real-root certificates on every complementary parameter
   interval, with all components and symmetry orbits accounted for; and
5. an exact fiber check at each critical value.

The stored F6 endpoint artifact does not record a continued third-basis
branch. `f6_boundary_reverify.jl` bisects the Boolean predicate
“the fresh float64 pool has a six-clique” and verifies one selected clique at
400-bit precision. Its brackets,
\([0.0550448030233,0.0550448030466]\) and
\([0.0625984193757,0.062598419399]\), are therefore not certified fold
values or certified endpoints of \(V(J_3)_{\rm phys}\). No exact branch
coordinate or algebraic critical value is currently recorded. The rational
angle \(\lambda=\pi/10\) is a useful exact-anchor candidate:
\(z=e^{i\pi/10}\) is a twentieth root of unity and lies inside the reported
numerical bracket, but existence of a physical third basis there still needs
an exact or interval certificate. The next legitimate fold task is to
certify an exact third basis at such an algebraic parameter, continue and
classify its component in the exact \(H_{F6}(z,t)\) incidence, and then test
the projection criteria above. Do not call the existing clique-search
threshold crossing a fold.

## First exact anchor computation

The exact F6 seam matrix at \(\lambda=\pi/10\) is now checked over
\(\mathbb Q(\zeta_{60})\), using \(z=\zeta_{60}^3\) and
\(\omega=\zeta_{60}^{20}\): all 36 entries have norm one and the exact Gram
matrix is \(6I\). The reproducible exporter is
`scripts/julia/export_f6_pi10_pool_exact.jl`; it writes the exact
single-vector MU-pool ideal to `symbolic_export/f6_pi10_pool_exact.m2`.
That ideal has ten variables and ten generators (five reciprocal equations
and five MU equations; the sixth MU equation follows from the unitary sum
identity).

The existing numerical pool pipeline at this parameter tracked 252 paths,
recovered 48 vectors, and found a six-clique. This is a useful candidate
anchor only: its coefficients/vectors have not yet been exactly reconstructed,
and it does not certify the incidence fiber or either arc endpoint. The
exact Gröbner computation for the single-vector pool ideal is currently
running; until it completes and the output is checked, no dimension or
solution-count claim is made from that ideal. Even a completed pool
calculation would not replace the six-vector \(J_3\) incidence computation.
