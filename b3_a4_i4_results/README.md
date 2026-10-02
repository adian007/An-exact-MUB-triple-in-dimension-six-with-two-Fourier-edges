# B3 A4 quotient + reduced fourth-MUB computation

This package continues the fold/A4 computation. It recovers the third-MUB MU-vector pool on the Dita circle, constructs the complete size-6 orthogonality cliques, quotients those cliques by the verified fixed-z A4 monomial symmetry, and probes the reduced fourth-vector witness.

## Run
Windows/Anaconda:

```bat
conda activate myenv
python -m pip install numpy scipy
python b3_orbit_i4.py
```

The script writes `b3_a4_i4_results/SUMMARY.json` and pool `.npz` files.

## Important convention
The Hadamard matrix is the exact Dita-circle Karlsson representative used in the project. MU equations use the **columns** of H_D. The matrix is unnormalized in the witness equations, so the target is 6.

## Numerical status
The results are numerical. The pool recovery uses random starts and a least-squares/root solve; the clique search uses a numerical orthogonality threshold. The A4 action is checked by direct vector matching. No claim of global completeness or exact fourth-MUB nonexistence is made here.

## Main result
For lambda = 0.4, pi/3, and 2pi/3:

- 72 recovered MU vectors;
- 4 recovered complete B3 cliques;
- all 4 cliques form one A4 orbit;
- one representative per parameter therefore suffices for the reduced fourth-vector numerical probe.

## Relation to the fourth-MUB problem
The reduced witness asks whether a vector w can be simultaneously unbiased to I, H_D(lambda), and the six vectors of the representative B3. The current numerical probes found no zero residual. The next stage is exact/validated incidence algebra on this reduced system.
