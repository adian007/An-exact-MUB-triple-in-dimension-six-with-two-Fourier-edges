"""
karlsson_k6_3.py

Implements Karlsson's three-parameter family of complex Hadamard
matrices K6(theta, phi, lambda), following Eqs. (7.2)-(7.6) of the
McNulty & Weigert review (which cites Karlsson 2011, arXiv:1003.4177).

Structure (Eq. 7.2):

    K6(theta,phi,lambda) =
        [ F2      Z1        Z2       ]
        [ Z3   (1/2) Z3 A Z1   (1/2) Z3 B Z2 ]
        [ Z4   (1/2) Z4 B Z1   (1/2) Z4 A Z2 ]

where F2 is the 2x2 Fourier (Hadamard) matrix, A is given by (7.3)-(7.5),
B = -F2 - A, and the "left" submatrices Z1, Z2 and "right" submatrices
Z3, Z4 (Eq. 7.6) are

    Zj = [[1, 1], [zj, -zj]]     (j=1,2, "left"/column-type)
    Zk = [[1, zk], [1, -zk]]     (k=3,4, "right"/row-type)

with unimodular z1,...,z4 related by Mobius transforms
M(z) = (alpha z - beta)/(beta_bar z - alpha_bar):

    z3^2 = M_A(z1^2),  z3^2 = M_B(z2^2),
    z4^2 = M_A(z2^2),  z4^2 = M_B(z1^2),

where alpha_A = A12^2, beta_A = A11^2, alpha_B = B12^2, beta_B = B11^2.

Free parameters: theta, phi in [0, pi) determine A (and hence B); z1 is
then a free unimodular parameter z1 = exp(i*lambda); z2, z3, z4 are
DETERMINED from z1 via the Mobius maps above (the review states: "By
choosing z1 = e^{i lambda}, the remaining three z-parameters are
uniquely determined").
"""

import numpy as np

TOL = 1e-8


def build_A(theta, phi):
    """
    Karlsson's Theorem 11 (arXiv:1003.4177), verified against the
    ORIGINAL paper (not just the review, which had a transcription
    slip): A is a 2x2 UNITARY matrix

        A = [[A11, A12], [conj(A12), -conj(A11)]]

    with
        A11 = -1/2 + i*sqrt(3)/2 * (cos(theta) + exp(-i*phi) sin(theta))
        A12 = -1/2 + i*sqrt(3)/2 * (-cos(theta) + exp(i*phi) sin(theta))

    NOTE the leading MINUS 1/2 (not +1/2), and the (2,1),(2,2) entries
    are conj(A12), -conj(A11) (not A12, -A11). Both details matter for
    A to actually be unitary (Prop. 4: A = F2*(-e/2 + i*sqrt(3)/2*Lambda)
    with Lambda Hermitian unitary).
    """
    A11 = -0.5 + 1j * (np.sqrt(3) / 2) * (np.cos(theta) + np.exp(-1j * phi) * np.sin(theta))
    A12 = -0.5 + 1j * (np.sqrt(3) / 2) * (-np.cos(theta) + np.exp(1j * phi) * np.sin(theta))
    A = np.array([[A11, A12], [np.conj(A12), -np.conj(A11)]], dtype=complex)
    return A


def mobius(z, alpha, beta):
    """M(z) = (alpha*z - beta) / (conj(beta)*z - conj(alpha))."""
    return (alpha * z - beta) / (np.conj(beta) * z - np.conj(alpha))


def build_k6(theta, phi, lam, verbose=False):
    """
    Build one member of Karlsson's K6^(3) family.
    Returns the 6x6 matrix (entries of modulus 1, i.e. UNnormalised,
    consistent with our convention elsewhere in this project).
    """
    F2 = np.array([[1, 1], [1, -1]], dtype=complex)
    A = build_A(theta, phi)
    B = -F2 - A

    # Early sanity gate: A and B must themselves be unitary (Prop 4).
    A_unitary_err = np.max(np.abs(A @ A.conj().T - 2 * np.eye(2)))
    B_unitary_err = np.max(np.abs(B @ B.conj().T - 2 * np.eye(2)))
    if verbose:
        print(f"  ||A A^dag - 2I|| = {A_unitary_err:.2e}   ||B B^dag - 2I|| = {B_unitary_err:.2e}")
    if A_unitary_err > 1e-6 or B_unitary_err > 1e-6:
        raise ValueError(f"A or B not unitary (errs {A_unitary_err:.2e}, {B_unitary_err:.2e}) "
                          f"-- construction is wrong, refusing to proceed.")

    alpha_A, beta_A = A[0, 1] ** 2, A[0, 0] ** 2
    alpha_B, beta_B = B[0, 1] ** 2, B[0, 0] ** 2

    z1 = np.exp(1j * lam)
    z1sq = z1 ** 2

    # z3^2 = M_A(z1^2); z4^2 = M_B(z1^2)
    z3sq = mobius(z1sq, alpha_A, beta_A)
    z4sq = mobius(z1sq, alpha_B, beta_B)

    # need z2 such that z3^2 = M_B(z2^2) AND z4^2 = M_A(z2^2) consistently.
    # invert: z2^2 = M_B^{-1}(z3^2) using M_B^{-1} = M with (alpha,beta) -> ?
    # Mobius M(z) = (alpha z - beta)/(betabar z - alphabar) is an
    # involution-type map when |alpha|^2-|beta|^2 = 1 (standard SU(1,1)
    # Mobius); for our purposes we solve directly: from z3^2 = M_A(z1^2)
    # and z3^2 = M_B(z2^2), invert the latter for z2^2.
    # M_B(z2^2) = z3^2  =>  alpha_B z2^2 - beta_B = z3^2 (betabar_B z2^2 - alphabar_B)
    # => z2^2 (alpha_B - z3^2 * conj(beta_B)) = beta_B - z3^2*conj(alpha_B)
    num = beta_B - z3sq * np.conj(alpha_B)
    den = alpha_B - z3sq * np.conj(beta_B)
    z2sq = num / den

    # cross-check via z4^2 = M_A(z2^2) should match the z4sq computed above
    z4sq_check = mobius(z2sq, alpha_A, beta_A)

    z2 = np.sqrt(z2sq)
    z3 = np.sqrt(z3sq)
    z4 = np.sqrt(z4sq)

    if verbose:
        print(f"  |z1|={abs(z1):.6f} |z2|={abs(z2):.6f} |z3|={abs(z3):.6f} |z4|={abs(z4):.6f}")
        print(f"  z4^2 consistency check: {abs(z4sq - z4sq_check):.2e} (should be ~0)")

    def Zleft(z):
        return np.array([[1, 1], [z, -z]], dtype=complex)

    def Zright(z):
        return np.array([[1, z], [1, -z]], dtype=complex)

    Z1 = Zleft(z1)
    Z2 = Zleft(z2)
    Z3 = Zright(z3)
    Z4 = Zright(z4)

    top = np.hstack([F2, Z1, Z2])
    mid = np.hstack([Z3, 0.5 * Z3 @ A @ Z1, 0.5 * Z3 @ B @ Z2])
    bot = np.hstack([Z4, 0.5 * Z4 @ B @ Z1, 0.5 * Z4 @ A @ Z2])

    K6 = np.vstack([top, mid, bot])
    return K6


if __name__ == "__main__":
    from hadamard6 import is_hadamard, dephase

    print("=== Testing Karlsson K6(theta, phi, lambda) construction ===\n")
    test_points = [
        (0.3, 0.5, 0.2),
        (1.0, 2.0, 0.7),
        (0.1, 0.1, 3.0),
        (np.arccos(1 / np.sqrt(3)), np.pi / 4, 0.4),  # should reduce to Dita family (review: Sec 7.1)
    ]
    for (theta, phi, lam) in test_points:
        K6 = build_k6(theta, phi, lam, verbose=True)
        ok = is_hadamard(K6)
        print(f"theta={theta:.3f} phi={phi:.3f} lambda={lam:.3f}  -> is_hadamard: {ok}\n")
