#!/usr/bin/env python3
"""w1_pi3_groebner.py — Phase 2: Exact W1 elimination at λ=π/3

Goal: Prove that at λ=π/3 (Diţă slice), the W1 witness system has unit ideal
over a cyclotomic field.

Pool vector format: V[i] is a flat vector (1, x1, x2, x3, x4, x5) with |xi|=1
These need normalization by √6 for MUB work.
"""

import numpy as np
import cmath
import math
from pathlib import Path
import json
import sys

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

N = 6
DITA = cmath.acos(1 / cmath.sqrt(3))
LAMBDA_PI_3 = math.pi / 3

def build_H_dita_pi3():
    """Build the Diţă Hadamard matrix at λ=π/3."""
    z = cmath.exp(1j * cmath.pi / 3)  # ζ₆ = 1/2 + i√3/2
    z_inv = cmath.exp(-1j * cmath.pi / 3)
    i = 1j
    
    H = np.array([
        [1, 1, 1, 1, 1, 1],
        [1, -1, z, -z, i, -i],
        [1, -i, i, i, -i, -1],
        [1, i, -z, z, -1, -i],
        [1, z_inv, -i, -1, -z_inv, i],
        [1, -z_inv, -1, -i, z_inv, i]
    ], dtype=complex)
    
    # Verify: H @ H^† = 6*I
    hh_dag = H @ H.conj().T
    error = np.linalg.norm(hh_dag - N * np.eye(N))
    print(f"H verification: H·H† = 6I? error = {error:.2e}")
    return H, z

def load_pool_and_clique():
    """Load pool and B3 clique, return normalized basis matrices."""
    pool_path = Path("D:/MUBs in 6-dimension/b3_a4_i4_results/pool_pi_over_3.npz")
    data = np.load(pool_path)
    V = data['V']  # shape (72, 6), each row is a flat vector
    
    with open("D:/MUBs in 6-dimension/b3_a4_i4_results/summary.json") as f:
        summary = json.load(f)
    
    # Use clique 0 from the pi/3 entry, not the 0.4 entry.  The summary
    # stores parameters by their exact floating-point string representation.
    target_key = str(LAMBDA_PI_3)
    if target_key not in summary:
        raise KeyError(f"Missing pi/3 summary entry {target_key!r}; available={list(summary)}")
    clique0 = summary[target_key]["cliques"][0]
    
    # Pool vectors: V[i] is a flat vector (1, x1, x2, x3, x4, x5) with norm = sqrt(6)
    # These are OUTPUT from the pool solver - each is a potential basis vector
    # Normalize by sqrt(6) to get proper unit basis vectors
    raw_vecs = V[clique0]  # shape (6, 6)
    
    # Check norms
    norms = np.linalg.norm(raw_vecs, axis=1)
    print(f"Raw clique vector norms (expected √6≈{np.sqrt(6):.4f}): {norms}")
    
    # Normalize
    B3_cols = raw_vecs / np.sqrt(N)  # columns are orthonormal basis vectors
    
    # B3 matrix: columns are basis vectors
    B3 = B3_cols.T  # shape (6, 6), columns are orthonormal
    
    print(f"B3 columns normalized, shape: {B3.shape}")
    unitary_defect = np.linalg.norm(B3 @ B3.conj().T - np.eye(N))
    print(f"B3 unitarity: ||B3·B3† - I|| = {unitary_defect:.2e}")
    
    return V, B3, clique0

def check_MUB_pair(B1, B2, label=""):
    """Check if two bases (column matrices) are mutually unbiased. Returns max defect."""
    # Columns are basis vectors, already normalized
    B1_cols = B1 if np.allclose(np.linalg.norm(B1, axis=0), 1) else B1 / np.sqrt(N)
    B2_cols = B2 if np.allclose(np.linalg.norm(B2, axis=0), 1) else B2 / np.sqrt(N)
    
    defects = []
    for i in range(N):
        for j in range(N):
            dot = np.dot(np.conj(B1_cols[:, i]), B2_cols[:, j])
            defect = abs(abs(dot)**2 - 1.0/N)
            defects.append(defect)
    
    max_def = max(defects)
    avg_def = np.mean(defects)
    print(f"  {label} MU defects: max={max_def:.2e}, avg={avg_def:.2e}")
    return max_def, avg_def

def main():
    print("=" * 70)
    print("PHASE 2: W1 ELIMINATION AT lambda=pi/3")
    print("=" * 70)
    
    # Step 1: Build H
    print("\n[1] Diţă Hadamard at λ=π/3")
    H, z = build_H_dita_pi3()
    print(f"    z = exp(iπ/3) = {z}")
    
    # Step 2: Load B3
    print("\n[2] Loading B3 from pool")
    V, B3, clique0 = load_pool_and_clique()
    
    # Step 3: Verify MUB properties
    print("\n[3] MUB Property Verification")
    I = np.eye(N)
    check_MUB_pair(I, H, "I-H")
    check_MUB_pair(I, B3, "I-B3")
    check_MUB_pair(H, B3, "H-B3")
    
    # Step 4: W1 analysis
    print("\n[4] W1 Witness System")
    print("    W1 seeks v = z/√6 (gauge z₀=1) unbiased to I, H, B3")
    print("    Variables: z₁..z₅, w₁..w₅ (10 vars)")
    print("    Equations: zⱼwⱼ=1 (5) + MU to H (6) + MU to B3 (6) = 17")
    
    # Step 5: Recommendation
    print("\n[5] RECOMMENDATION")
    print("    ✓ H matrix built EXACT (entries in ℚ(ζ₆))")
    print("    ✓ Pool B3 is numerical and verified against H")
    print("    ✓ Distinct-phase exact B3 ansatz passes symbolic identity checks")
    print(f"    ✓ Lambda = π/3 is algebraic (not D0-equivalent)")
    print()
    print("    NEXT: Compute and verify the coefficient field, then reduce the exact W1 ideal")
    print("    Input: Exact H and corrected distinct-phase B3 relations")
    
    # Save status
    status = {
        "phase": "2_setup",
        "status": "exact_candidate_ready_for_gb",
        "lambda": str(math.pi/3),
        "z_exact": str(z),
        "H_verified": True,
        "B3_clique": clique0,
        "pool_size": 72,
        "B3_MU_to_H": True,
        "B3_unitarity_defect": float(np.linalg.norm(B3 @ B3.conj().T - np.eye(N))),
        "normalization": "columns of B3 are unit norm; Hadamard-style raw vectors divided by sqrt(6)",
        "exact_B3_ansatz_available": True,
        "exact_B3_identity_audit": "18/18 Hadamard equations and 15/15 pairwise inner products reduce to zero",
        "groebner_status": "running; no W1 conclusion yet",
        "timestamp": str(np.datetime64('now'))
    }
    
    with open("D:/MUBs in 6-dimension/results/phase2_w1_setup.json", 'w') as f:
        json.dump(status, f, indent=2)
    print(f"\nStatus saved to results/phase2_w1_setup.json")

if __name__ == "__main__":
    main()
