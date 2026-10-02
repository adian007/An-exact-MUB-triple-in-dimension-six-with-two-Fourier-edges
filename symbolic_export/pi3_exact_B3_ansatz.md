# Exact π/3 B3 ansatz

This file records an exact algebraic candidate reconstructed from the
100-digit refinement of the first recovered π/3 clique.

## Cyclotomic parameter

Let (q) be the upper-half-plane primitive 36th root of unity satisfying

\[
q^{12}-q^6+1=0,
\qquad q=e^{i\pi/18}.
\]

The unnormalised Dita matrix at λ=π/3 is obtained from the repository
formula with (z=q^6), (i=q^9), and (t=q^{-6}).

## Common phase polynomial

Each free phase (uin\{a,b,c\}) satisfies

\[
\begin{aligned}
P(u)={}&104329u^{24}-193800u^{22}+192924u^{20}-299406u^{18}
 +331452u^{16}-381612u^{14}\\
&+493595u^{12}-381612u^{10}+331452u^8-299406u^6
 +192924u^4-193800u^2+104329=0.
\end{aligned}
\]

The polynomial is irreducible over ℚ in the current symbolic calculation.
Each phase has its own exact linear q-embedding relation; these relations
are recorded in `w1_pi3_exact_abc.m2` and must not be replaced by one common
q-relation. Choose the three roots by their unit-circle arguments (degrees):

- (a): (-2.235 < \arg(a) < -2.234);
- (b): (116.356 < \arg(b) < 116.357);
- (c): (135.254 < \arg(c) < 135.255).

These intervals select the roots matching the refined numerical clique.

## Exact B3 columns

The six unnormalised columns are

\[
\begin{array}{ll}
v_0=(1,a,-qa,q^{-10},q^{-16}a,q^{-8}), &
v_1=(1,-a,qa,q^{-10},q^2a,q^{-8}),\\
v_2=(1,b,q^7b,q^2,q^{-4}b,q^{16}), &
v_3=(1,-b,-q^7b,q^2,-q^{-4}b,q^{16}),\\
v_4=(1,c,q^{-5}c,q^{14},q^8c,q^4), &
v_5=(1,-c,-q^{-5}c,q^{14},-q^8c,q^4).
\end{array}
\]

The corresponding unitary basis is (B_3=[v_0\ \cdots\ v_5]/\sqrt6).

## Exact checks obtained symbolically

Using (q^{12}-q^6+1=0), all 15 pairwise inner products
(v_i^\dagger v_j) reduce identically to zero. The six MU equations for
each of the three pair templates eliminate to the same polynomial (P(u)).
The three templates have the same eliminated polynomial, but distinct linear
q-embedding relations. The full phase-ansatz basis is recorded in
`results/pi3_phase_groebner_basis.txt`, and the corrected W1 export is
`symbolic_export/w1_pi3_exact_abc.m2`.

## Remaining certification work

This is an exact algebraic candidate, not yet the final W1 theorem. The next
steps are to formalize the complex root isolations for (a,b,c), verify the
MU equations over the specified number-field embeddings, and substitute this
exact B3 into the W1 ideal before running the Gröbner elimination.
