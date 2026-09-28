"""defect2_instrumented.py

Instrumented replication of build_k6's INTERNAL step order for theta=0,
so the float-level seam behaviour matches exactly what builds H. The point:
is the Python builder NaN "at every lambda" (last session's claim) or does it
flip per-lambda on 1-ulp rounding? Reproduction only; no new experiment.
Writes results/defect2_instrumented_out.txt.
"""
import io
import sys

import numpy as np

sys.path.insert(0, "scripts/python")
from karlsson_k6_3 import build_A, build_k6, mobius  # noqa: E402

OUT = io.open("results/defect2_instrumented_out.txt", "w", encoding="utf-8")


def P(*a):
    s = " ".join(str(x) for x in a)
    OUT.write(s + "\n")
    print(s, flush=True)


np.seterr(all="ignore")
th, ph = 0.0, 0.5
A = build_A(th, ph)
F2 = np.array([[1, 1], [1, -1]], dtype=complex)
B = -F2 - A
alpha_A, beta_A = A[0, 1] ** 2, A[0, 0] ** 2
alpha_B, beta_B = B[0, 1] ** 2, B[0, 0] ** 2

P("Internal-order replication of build_k6 at (theta=0, phi=0.5):")
P("lam | z3sq            | num (z2)       | den (z2)       | z2sq | H finite | nan")
for lam in (0.0, 0.3, 0.48, 0.5, 0.52, 0.7, 1.0, 1.5, 2.0, 3.1415, 4.0, 5.0, 6.2832):
    z1 = np.exp(1j * lam)           # EXACT internal order: exp(1j*lam), then **2
    z1sq = z1 ** 2
    z3sq = mobius(z1sq, alpha_A, beta_A)
    num = beta_B - z3sq * np.conj(alpha_B)
    den = alpha_B - z3sq * np.conj(beta_B)
    z2sq = num / den
    from hadamard6 import dephase  # noqa: E402
    H = build_k6(th, ph, lam)
    nan_ct = int(np.isnan(H.real).sum() + np.isnan(H.imag).sum())
    P("%.4f | %s | r=%.1e i=%.1e | r=%.1e i=%.1e | %s | %s | %d"
      % (lam, complex(z3sq), num.real, num.imag, den.real, den.imag,
         complex(z2sq), bool(np.all(np.isfinite(H))), nan_ct))
P("")
P("Read: the seam collapse z3^2 == 1 happens only at ULP-exact equality; whether")
P("num/den hit bitwise 0 (then 0/0 = nan) flips per lambda on exp() rounding.")
P("Where it does NOT hit 0/0, z2sq is 1-ulp-noise driven (see lam=0.7 etc).")
P("")

# Same matrices BUT through the explicit seam branch (Julia rule):
P("Seam-branch rule replicated (z3sq=z4sq=1, z2sq=alpha_A/beta_A):")
z2sq_seam = alpha_A / beta_A
P("  z2sq = %s ; |z2sq| = %.15f" % (complex(z2sq_seam), abs(z2sq_seam)))
OUT.close()
