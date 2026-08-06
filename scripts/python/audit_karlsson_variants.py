"""
audit_karlsson_variants.py  (Phase A1 of the research plan)

Numerically compare TWO transcriptions of the A-matrix in Karlsson's
three-parameter family K6(theta, phi, lambda):

  Variant "karlsson2011" (current code, karlsson_k6_3.py):
      A11 = -1/2 + i*sqrt(3)/2 * ( cos(theta) + exp(-i*phi) sin(theta))
      A12 = -1/2 + i*sqrt(3)/2 * (-cos(theta) + exp(+i*phi) sin(theta))
      A   = [[A11, A12], [conj(A12), -conj(A11)]]

  Variant "mw_review" (McNulty & Weigert arXiv:2410.23997, Eqs. 7.3-7.5,
  as printed):
      A11 = +1/2 + i*sqrt(3)/2 * ( cos(theta) + exp(-i*phi) sin(theta))
      A12 = -1/2 + i*sqrt(3)/2 * (-cos(theta) + exp(+i*phi) sin(theta))
      A   = [[A11, A12], [A12, -A11]]

For each variant we check, on a grid of (theta, phi, lambda):
  1. ||A A^dag - 2 I||           (A must be sqrt(2)-scaled unitary, Prop. 4)
  2. ||B B^dag - 2 I||           with B = -F2 - A
  3. is_hadamard(K6)             (unimodular entries + K K^dag = 6 I)
  4. z4^2 Mobius consistency     |M_A(z2^2) - z4^2|

A variant PASSES a grid point if 1-3 hold to 1e-8 and 4 to 1e-6.
"""

import numpy as np
from hadamard6 import is_hadamard, hadamard_defect

F2 = np.array([[1, 1], [1, -1]], dtype=complex)


def build_A_karlsson2011(theta, phi):
    A11 = -0.5 + 1j * (np.sqrt(3) / 2) * (np.cos(theta) + np.exp(-1j * phi) * np.sin(theta))
    A12 = -0.5 + 1j * (np.sqrt(3) / 2) * (-np.cos(theta) + np.exp(1j * phi) * np.sin(theta))
    return np.array([[A11, A12], [np.conj(A12), -np.conj(A11)]], dtype=complex)


def build_A_mw_review(theta, phi):
    A11 = +0.5 + 1j * (np.sqrt(3) / 2) * (np.cos(theta) + np.exp(-1j * phi) * np.sin(theta))
    A12 = -0.5 + 1j * (np.sqrt(3) / 2) * (-np.cos(theta) + np.exp(1j * phi) * np.sin(theta))
    return np.array([[A11, A12], [A12, -A11]], dtype=complex)


def mobius(z, alpha, beta):
    return (alpha * z - beta) / (np.conj(beta) * z - np.conj(alpha))


def build_k6_from_A(A, lam):
    """Assemble K6 given a 2x2 block A and unimodular parameter z1=e^{i lam}.
    Returns (K6, diagnostics dict). Does NOT raise on non-unitary A/B, so we
    can audit broken variants."""
    B = -F2 - A
    diag = {
        "A_unitary_err": float(np.max(np.abs(A @ A.conj().T - 2 * np.eye(2)))),
        "B_unitary_err": float(np.max(np.abs(B @ B.conj().T - 2 * np.eye(2)))),
    }

    alpha_A, beta_A = A[0, 1] ** 2, A[0, 0] ** 2
    alpha_B, beta_B = B[0, 1] ** 2, B[0, 0] ** 2

    z1 = np.exp(1j * lam)
    z1sq = z1 ** 2
    z3sq = mobius(z1sq, alpha_A, beta_A)
    z4sq = mobius(z1sq, alpha_B, beta_B)

    num = beta_B - z3sq * np.conj(alpha_B)
    den = alpha_B - z3sq * np.conj(beta_B)
    z2sq = num / den

    diag["z4_consistency"] = float(np.abs(z4sq - mobius(z2sq, alpha_A, beta_A)))
    diag["z_moduli_err"] = float(max(abs(abs(z) - 1) for z in (z2sq, z3sq, z4sq)))

    z2, z3, z4 = np.sqrt(z2sq), np.sqrt(z3sq), np.sqrt(z4sq)

    def Zleft(z):
        return np.array([[1, 1], [z, -z]], dtype=complex)

    def Zright(z):
        return np.array([[1, z], [1, -z]], dtype=complex)

    Z1, Z2, Z3, Z4 = Zleft(z1), Zleft(z2), Zright(z3), Zright(z4)

    top = np.hstack([F2, Z1, Z2])
    mid = np.hstack([Z3, 0.5 * Z3 @ A @ Z1, 0.5 * Z3 @ B @ Z2])
    bot = np.hstack([Z4, 0.5 * Z4 @ B @ Z1, 0.5 * Z4 @ A @ Z2])
    K6 = np.vstack([top, mid, bot])
    return K6, diag


def audit_variant(name, build_A, grid):
    n_pass = n_fail = 0
    worst = {"A_unitary_err": 0.0, "B_unitary_err": 0.0,
             "hadamard_unitary_err": 0.0, "hadamard_modulus_err": 0.0,
             "z4_consistency": 0.0}
    first_fail = None
    for (theta, phi, lam) in grid:
        A = build_A(theta, phi)
        K6, diag = build_k6_from_A(A, lam)
        u_err, m_err = hadamard_defect(K6)
        worst["A_unitary_err"] = max(worst["A_unitary_err"], diag["A_unitary_err"])
        worst["B_unitary_err"] = max(worst["B_unitary_err"], diag["B_unitary_err"])
        worst["hadamard_unitary_err"] = max(worst["hadamard_unitary_err"], u_err)
        worst["hadamard_modulus_err"] = max(worst["hadamard_modulus_err"], m_err)
        worst["z4_consistency"] = max(worst["z4_consistency"], diag["z4_consistency"])

        ok = (diag["A_unitary_err"] < 1e-8 and diag["B_unitary_err"] < 1e-8
              and is_hadamard(K6) and diag["z4_consistency"] < 1e-6)
        if ok:
            n_pass += 1
        else:
            n_fail += 1
            if first_fail is None:
                first_fail = (theta, phi, lam, diag, u_err, m_err)

    print(f"--- variant: {name} ---")
    print(f"  grid points: {n_pass + n_fail}   PASS: {n_pass}   FAIL: {n_fail}")
    for k, v in worst.items():
        print(f"  worst {k}: {v:.3e}")
    if first_fail is not None:
        th, ph, lm, diag, u_err, m_err = first_fail
        print(f"  first failure at (theta,phi,lambda)=({th:.3f},{ph:.3f},{lm:.3f})")
        print(f"    A_unitary={diag['A_unitary_err']:.3e} B_unitary={diag['B_unitary_err']:.3e} "
              f"K_unitary={u_err:.3e} K_modulus={m_err:.3e} z4={diag['z4_consistency']:.3e}")
    print()
    return n_pass, n_fail


def main():
    rng = np.random.default_rng(20260802)
    # Structured grid + random samples; avoid exact endpoints.
    thetas = np.linspace(0.05, np.pi - 0.05, 7)
    phis = np.linspace(0.05, np.pi - 0.05, 7)
    lams = np.linspace(0.05, 2 * np.pi - 0.05, 5)
    grid = [(t, p, l) for t in thetas for p in phis for l in lams]
    grid += [tuple(rng.uniform(0.01, 3.1, 2)) + (rng.uniform(0, 2 * np.pi),)
             for _ in range(100)]
    # Named special points from the review (Sec 7.1)
    grid += [
        (np.arccos(1 / np.sqrt(3)), np.pi / 4, 0.4),   # Dita family point
        (0.3, 0.5, 0.2), (1.0, 2.0, 0.7), (0.1, 0.1, 3.0),
    ]

    print(f"Auditing {len(grid)} parameter points per variant\n")
    r1 = audit_variant("karlsson2011 (current code)", build_A_karlsson2011, grid)
    r2 = audit_variant("mw_review (Eqs. 7.3-7.5 as printed)", build_A_mw_review, grid)

    print("=== VERDICT ===")
    for name, (np_, nf) in [("karlsson2011", r1), ("mw_review", r2)]:
        status = "VALID CHM family" if nf == 0 else ("INVALID (all fail)" if np_ == 0 else "MIXED")
        print(f"  {name}: {status}")


if __name__ == "__main__":
    main()
