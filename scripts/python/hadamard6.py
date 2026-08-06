"""
hadamard6.py

Shared helpers for 6x6 complex Hadamard matrix (CHM) work.
Conventions match the Julia driver mub_zauner_6d_liang_chen.jl:

- A CHM H is UNnormalised: |H_ij| = 1 and H H^dag = n I (n = matrix order).
- The corresponding unitary / MUB transition matrix is H / sqrt(n).
"""

import numpy as np

TOL = 1e-8


def is_hadamard(H, tol=TOL):
    """True if H has unimodular entries and H H^dag = n I (n = order)."""
    H = np.asarray(H, dtype=complex)
    n = H.shape[0]
    if H.shape != (n, n):
        return False
    unitary_err = np.linalg.norm(H @ H.conj().T - n * np.eye(n))
    modulus_err = np.max(np.abs(np.abs(H) - 1.0))
    return unitary_err < tol and modulus_err < tol


def hadamard_defect(H):
    """Continuous defect: (||H H^dag - n I||_F, max entrywise | |H_ij| - 1 |)."""
    H = np.asarray(H, dtype=complex)
    n = H.shape[0]
    unitary_err = np.linalg.norm(H @ H.conj().T - n * np.eye(n))
    modulus_err = np.max(np.abs(np.abs(H) - 1.0))
    return unitary_err, modulus_err


def dephase(H, tol=TOL):
    """
    Return the dephased form of H: first row and first column made real
    positive (entries = 1 for a CHM) by multiplying columns/rows by phases.
    Same algorithm as the Julia `dephase`.
    """
    H = np.array(H, dtype=complex, copy=True)
    for j in range(1, H.shape[1]):
        if abs(H[0, j]) > tol:
            H[:, j] *= np.conj(H[0, j]) / abs(H[0, j])
    for i in range(1, H.shape[0]):
        if abs(H[i, 0]) > tol:
            H[i, :] *= np.conj(H[i, 0]) / abs(H[i, 0])
    return H
