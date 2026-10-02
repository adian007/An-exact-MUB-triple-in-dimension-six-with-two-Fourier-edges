# Exact Fourier certificate for the pi/3 algebraic triple — 2026-10-02

## Result

The exact algebraic basis (B_3) recorded in
`symbolic_export/pi3_exact_B3_ansatz.md` is monomially equivalent to a member
of the order-six Fourier family. Consequently, the specified exact triple
\(\{I,H,B_3\}\) cannot be extended to a fourth MUB. The conclusion is local
to this exact triple; it does not rule out all MUB quartets in dimension six.

The exact check is recorded in
`results/phase2_pi3_fourier_exact_direct.json` and can be reproduced with:

```powershell
python scripts/python/prove_pi3_fourier_exact_direct.py
```

The run used Python 3.13.14 and SymPy 1.14.0. It checked all 36 chart entries
and all three phase-modulus identities by polynomial remainder in
\(\mathbb Q[q]/(q^{12}-q^6+1)\); every remainder is zero.

## Exact input and MUB validity

Let \(q\) be the primitive 36th root with the exact complex embedding
\(q=e^{i\pi/18}\), so \(q^{12}-q^6+1=0\). Let \(a,b,c\) be roots of the
three irreducible quadratic factors \(f_a,f_b,f_c\) displayed in
`phase2_w1_coefficient_field_2026-10-01.md`. These factors are the exact
degree-two gcds of the phase polynomial \(P(u)\) and the corresponding
phase-specific relation \(L_i(q,u)\). Their joint quotient is the recorded
degree-96 field.

For the unordered basis, the sign of each root is immaterial: replacing
\(a\), \(b\), or \(c\) by its negative swaps the corresponding column pair
\((v_0,v_1)\), \((v_2,v_3)\), or \((v_4,v_5)\), using \(q^{18}=-1\). Thus
the exact embedding with this fixed \(q\) defines one unordered basis without
choosing a sign for each square root.

The three squared phases are \(d_i(q)=-f_i(0)\). The exact checks in the
certificate establish
\[
d_i(q)d_i(q^{-1})=1\qquad(i=a,b,c).
\]
Since \(q\) is on the unit circle, this implies \(|a|=|b|=|c|=1\). The
corrected exact remainder report records zero for all 15 pairwise inner
products of the six columns and all 36 cleared unbiasedness equations against
the exact \(H\). The latter equations were reduced using \(P\) and \(L_i\),
which vanish for every root of the corresponding \(f_i\). The exact \(H\)
is the \(\lambda=\pi/3\) Diță matrix in the ansatz and in
`results/exact_karlsson_dita_H.txt`.

## Exact Fourier chart

Write \(M=[v_0\ \cdots\ v_5]\), the unnormalised matrix in the exact ansatz,
and use zero-based row and column orders
\[
p=(4,3,1,0,2,5),\qquad s=(2,0,4,3,1,5).
\]
For the permuted matrix \(N_{ij}=M_{p_i,s_j}\), dephase at the upper-left
entry:
\[
C_{ij}=\frac{N_{ij}N_{00}}{N_{i0}N_{0j}}.
\]
Put \(r=b/a\), \(t=b/c\), \(\zeta=q^6\), and \(\omega=q^{12}\). Reduction
using \(q^{12}-q^6+1=0\) gives the exact matrix
\[
C=\begin{pmatrix}
1&1&1&1&1&1\\
1&r&t&-1&-r&-t\\
1&\omega&\omega^2&1&\omega&\omega^2\\
1&\omega r&\omega^2t&-1&-\omega r&-\omega^2t\\
1&\omega^2&\omega&1&\omega^2&\omega\\
1&\omega^2r&\omega t&-1&\zeta r&\zeta^{-1}t
\end{pmatrix}.
\]
This is the unnormalised repository convention for
\(\widetilde F_6^{(2)}(x,y)\) with
\[
x=\frac{b}{q^6a}=q^{30}r,\qquad
y=\frac{b}{q^{12}c}=q^{24}t.
\]
The parameters have unit modulus because \(a,b,c,q\) do. Dividing the chart
by \(\sqrt6\) gives the paper's normalized family \(F(x,y)\); equivalently,
the normalized basis matrix \(B_3=M/\sqrt6\) has the exact form
\(D_1P_1F_6^{(2)}(x,y)P_2D_2\), where \(F_6^{(2)}\) here denotes the
normalized family and the diagonal factors are the row and column phases
removed by the displayed dephasing. No fitting or phase-root approximation
is used in this identity.

## Theorem application and scope

The paper's matrix in Eq. (3) is the transpose placement of the repository
family: \(F_6^{(2)}(x,y)^T\) is exactly its \(F(x,y)\) convention. The
repository's common-unitary/conjugation argument shows that a pair
\((I,B_3)\) extends to a quartet if and only if its transpose-oriented pair
does. Theorem 1.4 of [Jaming, Matolcsi, Móra, Szöllősi, and Weiner,
arXiv:0902.0882v2](https://arxiv.org/html/0902.0882v2) rules out a quartet
containing \((I,F(x,y))\) for every pair of unit-modulus parameters.
Therefore \((I,B_3)\) cannot lie in any MUB quartet. Any quartet extending
\(\{I,H,B_3\}\) would contain that forbidden pair, which is a contradiction.

This establishes only that this exact \(\pi/3\) triple has no fourth basis.
It does not show that every third basis for the fixed \(H\) is represented in
the recovered 72-vector pool, nor that no four MUBs exist anywhere in
\(\mathbb C^6\). The exact ansatz was reconstructed from the stored clique;
the separate identification of that algebraic object with the stored
high-precision vectors remains the numerical 100-digit comparison in
`phase2_pi3_distinct_embedding_audit_2026-10-01.md`. The nonextension proof
itself concerns the exact algebraic triple defined above and does not rely on
that numerical comparison.

## Reproducibility

- Exact chart/modulus checker: `scripts/python/prove_pi3_fourier_exact_direct.py`
- Machine-readable output: `results/phase2_pi3_fourier_exact_direct.json`
- Exact MUB-identity checker: `scripts/python/verify_pi3_exact_remainders.py`
- Machine-readable MUB checks: `results/phase2_pi3_exact_remainder_verification.json`
- Exact B3 ansatz: `symbolic_export/pi3_exact_B3_ansatz.md`
- Exact MUB remainders: `docs/research/reports/phase2_pi3_exact_remainder_verification_2026-10-01.md`
- Common phase field: `docs/research/reports/phase2_w1_coefficient_field_2026-10-01.md`
- Theorem and transpose-convention audit: `docs/results/k3_family_membership_2026-09-29.md`, C0 and C5

SHA-256 hashes for the direct certificate bundle:

| File | SHA-256 |
|---|---|
| `scripts/python/prove_pi3_fourier_exact_direct.py` | `DBD7782452A35710B2E225ACD0A5FF62441EA124C748B53E398A022C397621C2` |
| `results/phase2_pi3_fourier_exact_direct.json` | `ED5A567D59FE3A67E94038FD1A46C9A4ABC83018A295F60933A622127BC0BBC7` |
| `symbolic_export/pi3_exact_B3_ansatz.md` | `8516E91D0E34F8C2669326BF980B14A889CD45965480334A320E99D6815F0BDE` |
| `docs/research/reports/phase2_w1_coefficient_field_2026-10-01.md` | `A552240E36C9FB2A3DA41F44E3CB504D84CE8CAEC7C66AA67E842B943AB23017` |
| `scripts/python/verify_pi3_exact_remainders.py` | `D38A84E4748BAA7C16E9B2485D2B3D4F9E64004D1DEA1DFC9F3409B5EE58D3F9` |
| `results/phase2_pi3_exact_remainder_verification.json` | `316AA01276A25FEC73C86F17E501756D69C8A86E70B793E5CFEFC45E0B4C4F60` |
| `docs/research/reports/phase2_pi3_exact_remainder_verification_2026-10-01.md` | `606D317C4D2102CD6CFA0F228F2FF6BB38C3D264FBEFD31480D90AF6BA182B24` |
