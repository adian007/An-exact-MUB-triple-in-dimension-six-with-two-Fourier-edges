# Addendum: pi/3 Fourier fit update (2026-09-29)

This new addendum supplements, and does not replace, any existing package file.

The earlier `Fourier_transition_verified_samples.csv` entry for \(\lambda=1.047197551197=\pi/3\) described the B3 representative as not passing its best fixed-pattern test and called for another transition to be tested. That fixed-pattern diagnostic is not the exhaustive two-parameter family-equivalence fit performed in the K3 family-membership audit.

In the new exhaustive numerical fit, clique 0's B3 transition fits the repository direct \(F_6^{(2)}(a,b)\) orientation with:

- consistency error: `3.60e-15`
- total maximum residual: `3.58e-15`
- fixed-entry maximum residual: `3.18e-15`
- row permutation `p=(4,3,1,0,2,5)`
- column permutation `q=(2,0,4,3,1,5)`
- fitted `a=0.5211411612207759+0.8534704974874417i`
- fitted `b=-0.7535391848223343-0.6574029943172545i`
- fit tolerance: `1e-8`

The separate transition \(H_D^\dagger B3\) for the same clique also fits the direct repository orientation:

- consistency error: `3.33e-15`
- total maximum residual: `3.69e-15`
- fixed-entry maximum residual: `3.69e-15`
- row permutation `p=(1,0,3,4,5,2)`
- column permutation `q=(1,2,4,0,3,5)`
- fitted `a=-0.959701-0.281022i`, `b=-0.890371+0.455235i`

The earlier pi/3 sample status is therefore superseded as a statement about this selected sample and this exhaustive numerical fitter. These floating-point fits are not exact family-membership proofs, do not imply family- or component-level membership, do not imply pool completeness, and do not establish MUB nonextendibility. Full per-clique and per-transition data are in [`docs/results/k3_family_membership_2026-09-29.md`](../docs/results/k3_family_membership_2026-09-29.md).
