# K3 Fourier-structure diagnostic

## Scope

The diagnostic runner analyzed the three archived pools at
\(\lambda=0.4,\pi/3,2\pi/3\). Each pool contains 72 recovered MU vectors and
the graph search found four six-vector cliques. For each clique it formed
\[
 B_3,\qquad H_D^\dagger B_3,\qquad B_3^\dagger H_D,
\]
checked transition flatness/unitarity, and searched every pivot-row,
pivot-column, and target-column dephasing chart for the three-\(-1\)
criterion.

## Observed numerical pattern

For all 12 supplied cliques:

- `B3` and `H_D^dagger_B3` had a Fourier-family column witness.
- `B3_dagger_H_D` had no Fourier-column witness at the configured tolerance,
  but its transpose had a Fourier-column witness (equivalently, it passed the
  transposed-Fourier row test).
- The largest reported transition flatness and unitarity residuals were below
  \(5\times10^{-15}\); detected dephased \(-1\) residuals were below
  \(7\times10^{-15}\).

The complete charts and residuals are preserved in
[`diagnostic.json`](../../results/campaigns/k3_fourier_structure/diagnostic.json),
including the numerical transition matrices as real/imaginary pairs.

## Exact chart systems

For each of the 12 detected witnesses, the exporter generated an exact
symbolic system containing:

- the exact symbolic transition matrix \(T=B_3^\dagger H_D(z)\), or its
  transpose for the transposed-family chart;
- the generic parameter relation \(zt-1=0\);
- 30 inverse-variable equations for the six dephased flat vectors;
- 36 MU equations against \(H_D(z)\);
- 30 pairwise basis-orthogonality equations;
- 3 polynomial equations imposing the selected three-\(-1\) chart.

The equations and chart metadata are in
[`exact_chart_equations.json`](../../results/campaigns/k3_fourier_structure/exact_chart_equations.json).
Each system has 62 variables and 100 equations. No Groebner elimination or
generic-\(z\) component-containment result has yet been computed.
Substitution of each sample into its three exported chart equations gave a
maximum unscaled polynomial residual of \(3.8\times10^{-14}\).

## Interpretation

This confirms that the simple invariant is a useful way to expose candidate
charts and produce exact equations for subsequent algebraic work. It does not
prove that the sampled numerical clique lies on an exact branch with the same
chart, that the supplied pools are complete, or that every relevant
Karlsson-family component has the transition structure. In particular,
\(B_3^\dagger H_D\) is numerically identified in the transposed Fourier family
for these samples; its reverse transition \(H_D^\dagger B_3\) is numerically
identified in the Fourier family.
