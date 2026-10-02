# Scripts

- `refine_singular_augmented.py`: refines the 30 initial SVD candidates using F=0, J*u=0, ||u||^2=1 and retains the 24 good roots.
- `a4_orbit_refined.py`: applies the two verified A4 monomial generators and computes the induced permutation/orbit structure on the 24 refined roots.
- `fold_diagnostics.py`: computes reduced-chart fold diagnostics using the first five MU equations.

All scripts are exploratory numerical scripts. They should be rerun with higher precision before being used for a formal proof claim.
