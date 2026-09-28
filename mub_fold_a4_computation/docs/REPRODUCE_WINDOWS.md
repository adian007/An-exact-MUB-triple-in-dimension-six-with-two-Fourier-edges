# Windows / Anaconda reproduction

From the package directory:

```bat
conda create -n mub-fold python=3.11 -y
conda activate mub-fold
python -m pip install numpy scipy
python run_fold_a4.py
```

For the full algebraic/certification pipeline, use Julia and HomotopyContinuation.jl in a separate environment. The package deliberately does not assume Julia is installed by the Python script.

Recommended Julia packages for the next stage include `HomotopyContinuation`, `LinearAlgebra`, `Random`, `CSV`, and `DataFrames`. Certification of nonsingular solutions can be performed with `certify`; singular solutions require the dedicated augmented/deflation workflow described above.
