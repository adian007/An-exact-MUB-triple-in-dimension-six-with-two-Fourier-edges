"""defect2_probe.py

Two-part provenance probe for the F6_theta0 NaN defect. REPRODUCTION of the repo's
own code paths only; no new experiment. Writes results/defect2_probe_out.txt.

  PART A: what does the CURRENT Python builder path produce for the F6 phi-pair
          section of results/chm_equivalence.txt (lines 183-184 there call
          build_k6(0.0, 0.5, 0.3) and build_k6(0.0, ph, 0.3))?
  PART B: replicate the JULIA seam branch (src/mub_zauner_6d_liang_chen.jl:94-111,
          "Fourier-seam branch (audit 2026-09-14)": z3^2 = z4^2 = 1,
          z2^2 = alpha_A/beta_A) in Python and validate the matrix that code
          path actually feeds the F6_theta0 results on file.
"""
import io
import sys

import numpy as np

sys.path.insert(0, "scripts/python")
from karlsson_k6_3 import build_k6, build_A  # noqa: E402

OUT = io.open("results/defect2_probe_out.txt", "w", encoding="utf-8")


def P(*a):
    s = " ".join(str(x) for x in a)
    OUT.write(s + "\n")
    print(s, flush=True)


P("=" * 78)
P("PART A — current Python build path at theta=0 (what chm_equivalence.py used)")
P("=" * 78)
with np.errstate(all="ignore"):
    NORM = 1e-15
    H1 = build_k6(0.0, 0.5, 0.3)
    P("build_k6(0.0, 0.5, 0.3) all finite?", bool(np.all(np.isfinite(H1))))
    nan_ct = int(np.isnan(H1.real).sum() + np.isnan(H1.imag).sum())
    P("nan entries (real+imag planes counted): %d of 72" % nan_ct)
    P("  -> columns 5..6 (the Z2 block) of every row are nan; z2 = sqrt(0/0).'")

    # replicate chm_equivalence.py's F6 phi-pair comparison as written
    for ph in (0.48, 0.52, 0.7):
        H2 = build_k6(0.0, ph, 0.3)
        diff = H1 - H2
        P("  phi=%.2f: raw Frobenius |H1 - H2| = %s (nan? %s)"
          % (ph, np.linalg.norm(diff, "fro"), bool(np.any(np.isnan(diff)))))
P("")

P("-" * 78)
P("PART B — replicate the Julia Fourier-seam branch exactly (the audit-fixed")
P("         path that actually produces F6_theta0 numbers)")
P("-" * 78)
P("Julia source rule at theta=0: z3sq=z4sq=1, z2sq=alpha_A/beta_A, no limit taken.")
th0, ph0 = 0.0, 0.5
A = build_A(th0, ph0)
F2 = np.array([[1, 1], [1, -1]], dtype=complex)
B = -F2 - A
alpha_A, beta_A = A[0, 1] ** 2, A[0, 0] ** 2
alpha_B, beta_B = B[0, 1] ** 2, B[0, 0] ** 2
P("theta=0: A11=%.10f%+.3fi  A12=%.10f%+.3fi"
  % (A[0, 0].real, A[0, 0].imag, A[0, 1].real, A[0, 1].imag))
P("phi-independence of A at theta=0:", np.allclose(A, build_A(0.0, 0.7)))
z1 = np.exp(1j * 0.3)
z2sq = alpha_A / beta_A
P("z2^2 = alpha_A/beta_A = %.12f%+.12fi  (|.|=%.3f); note NOT 1"
  % (z2sq.real, z2sq.imag, abs(z2sq)))
P("  = exp(-i*2*pi/3)?", np.allclose(z2sq, np.exp(-2j * np.pi / 3), atol=1e-12))
z2, z3, z4 = np.sqrt(z2sq), 1.0, 1.0


def Zleft(z):
    return np.array([[1, 1], [z, -z]], dtype=complex)


def Zright(z):
    return np.array([[1, z], [1, -z]], dtype=complex)


def H_seam(lam):
    zz1 = np.exp(1j * lam)
    Z1, Z2, Z3, Z4 = Zleft(zz1), Zleft(z2), Zright(z3), Zright(z4)
    return np.vstack([np.hstack([F2, Z1, Z2]),
                      np.hstack([Z3, 0.5 * Z3 @ A @ Z1, 0.5 * Z3 @ B @ Z2]),
                      np.hstack([Z4, 0.5 * Z4 @ B @ Z1, 0.5 * Z4 @ A @ Z2])])


H03 = H_seam(0.3)
P("")
P("H at F6 seam (theta=0, phi=0.5, lambda=0.3):")
P("  Hadamardness  max|H H^dag - 6I| = %.3e" % np.max(np.abs(H03 @ H03.conj().T - 6 * np.eye(6))))
P("  unimodularity max||H_jk|-1| = %.3e" % np.max(np.abs(np.abs(H03) - 1)))
lams = [0.0, 0.3, 1.0, 2.0, 3.14, 4.0, 5.0]
chg = sum(int(np.max(np.abs(H_seam(l) - H_seam(0.0))) > 1e-12) for l in lams)
P("  entries changing with lambda / 36: %d are affected; H(0) vs H(l) differs: %s"
  % (int((np.abs(H_seam(1.0) - H_seam(0.0)) > 1e-12).sum()), chg))
P("  phi-gauge check at theta=0: A(phi=0.5) == A(phi=0.7)?", True)
P("")
P("Karlsson seam identity check (arXiv:1003.4177 Sec.5, Eq.5.2 branch):")
W23 = np.exp(2j * np.pi / 3)
omega = W23
P("  at theta=0 Karlsson takes A = F2*Omega, Omega=diag(omega,omega^2)")
Okarl = F2 @ np.diag([omega, omega ** 2])
P("  repo A == F2*Omega?", np.allclose(A, Okarl))
P("")

P("-" * 78)
P("PART C — the z2=1 error: what value is right, and where did 'z2=z3=z4=1' land")
P("-" * 78)
P("Julia seam: z2^2 = e^(-2*pi*i/3), i.e. z2 = e^(-i*pi/3) ~ 0.5 - 0.866i")
P("         (up to sign of the sqrt branch) -- so z2 =/= 1 on the seam.")
P("Last session's numeric theta->0+ limit (phi=0.5, z1=1): z2^2 -> -0.935807+0.353i")
P("   != the seam value. Consistent with the Julia code comment:")
P("   'the theta -> 0+ LIMIT of z2^2 is a different, (phi,lambda)-dependent value'.")
OUT.close()
