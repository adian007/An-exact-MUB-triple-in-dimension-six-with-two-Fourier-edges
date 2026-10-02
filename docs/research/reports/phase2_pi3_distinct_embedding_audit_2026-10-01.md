# pi/3 distinct-phase exact embedding audit — 2026-10-01

## Result

The corrected distinct-phase export passes a high-precision numerical audit.
The audit uses the six vectors in `results/phase2_pi3_B3_mp100.json` and the
active export `symbolic_export/w1_pi3_exact_abc.m2`; it does not test the
deprecated `a=b=c` shortcut.

| Check | Maximum residual |
|---|---:|
| Exact vector-template reconstruction | `6.76e-101` |
| Pairwise orthogonality | `1.50e-100` |
| Unit-modulus entries / norm | `9.58e-102` |
| Hadamard unbiasedness | `1.42e-101` |
| `P(a), P(b), P(c)` | `4.08e-95` |
| `L0(q,a), L1(q,b), L2(q,c)` | `1.60e-72` |

The larger `L` residual is expected from evaluating degree-22 integer
polynomials with coefficients near `1e28` on finite 100-digit input data. The
corrected `L2` coefficient is `5354443114714978692181958184`; an earlier
export contained a transcription error in that coefficient.

This is numerical evidence for the exact candidate and a reproducible CAS
input, not yet a Gröbner-basis proof. The next proof step remains running the
corrected export in Macaulay2/Singular once a CAS runtime is available.
