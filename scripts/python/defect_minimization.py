"""
defect_minimization.py — Prong A: Global defect minimization over (H, B3, B4)

Defect = sum over all cross-basis pairs | |<v,w>|^2 - 1/6 |
       + sum over intra-basis pairs |<v,w>|^2 (should be 0 for same basis)

H is from Karlsson family K6^(3) (theta, phi, lambda).
B3, B4 are 6×6 unitary matrices parameterized via Cayley transform of
skew-Hermitian matrices.

Uses scipy.optimize with multiple random restarts (at least 100).
"""

import sys
import time
import json
from pathlib import Path

import numpy as np
from scipy.optimize import differential_evolution, minimize, Bounds

# ── Add repo scripts to path ──────────────────────────────────────────────────
SCRIPT_DIR = Path(r"D:/MUBs in 6-dimension/scripts/python")
sys.path.insert(0, str(SCRIPT_DIR))
from karlsson_k6_3 import build_k6
from hadamard6 import is_hadamard

# ── Constants ──────────────────────────────────────────────────────────────────
N = 6                          # dimension
TARGET_CROSS = 1.0 / N         # |<v,w>|^2 should be 1/6 for different bases
TARGET_INTRA = 0.0             # |<v,w>|^2 should be 0 for same basis
NP_RESTARTS = 150              # number of random restarts
SEED = 42

# ── Unitary parameterization ──────────────────────────────────────────────────
# We parameterize a unitary U via the Cayley transform of a skew-Hermitian S:
#   U = (I - S)(I + S)^{-1}
# where S is skew-Hermitian (S^† = -S).
# This maps the set of skew-Hermitian matrices onto the set of unitary matrices
# with no eigenvalue at -1 (which is fine for our purposes; if -1 is hit we
# perturb slightly).
# A 6×6 skew-Hermitian matrix has 6*6 = 36 real parameters (diagonal purely
# imaginary, off-diagonal complex with S_ji = -conj(S_ij)).
# For efficiency we parameterize only the independent entries.

def skew_hermitian_to_unitary(S):
    """Cayley transform: U = (I - S)(I + S)^{-1}."""
    I = np.eye(N, dtype=complex)
    try:
        U = (I - S) @ np.linalg.inv(I + S)
    except np.linalg.LinAlgError:
        # Perturb and retry
        U = (I - S) @ np.linalg.inv(I + S + 1e-10 * np.eye(N, dtype=complex))
    return U


def unitary_defect(U):
    """How far is U from being unitary? ||U U^† - I||_F."""
    return np.linalg.norm(U @ U.conj().T - np.eye(N))


def params_to_B(params):
    """Convert a 36-dimensional real parameter vector to a 6×6 unitary."""
    # params are 36 real numbers: diag (6 imag), upper triangle (15 complex = 30 real)
    S = np.zeros((N, N), dtype=complex)
    # Diagonal: purely imaginary
    for i in range(N):
        S[i, i] = 1j * params[i]
    # Upper triangle: complex
    idx = N
    for i in range(N):
        for j in range(i + 1, N):
            S[i, j] = params[idx] + 1j * params[idx + 1]
            S[j, i] = -np.conj(S[i, j])
            idx += 2
    U = skew_hermitian_to_unitary(S)
    return U


def B_to_params(B):
    """Convert a 6×6 unitary to a 36-dimensional real parameter vector."""
    # Inverse Cayley: S = (I - U)(I + U)^{-1}
    I = np.eye(N, dtype=complex)
    try:
        S = (I - B) @ np.linalg.inv(I + B)
    except np.linalg.LinAlgError:
        S = (I - B) @ np.linalg.inv(I + B + 1e-10 * np.eye(N, dtype=complex))
    # Skew-Hermitian perturbations may accumulate; force skew-Hermitian
    S = 0.5 * (S - S.conj().T)
    params = np.zeros(36)
    for i in range(N):
        params[i] = S[i, i].imag
    idx = N
    for i in range(N):
        for j in range(i + 1, N):
            params[idx] = S[i, j].real
            params[idx + 1] = S[i, j].imag
            idx += 2
    return params


# ── Defect computation ─────────────────────────────────────────────────────────

def compute_defect(H, B3, B4):
    """
    Compute total MU defect.
    H, B3, B4 are 6×6 matrices (columns are basis vectors).
    For MUBs, we normalize: columns of H/√6, B3/√6, B4/√6 are orthonormal.
    Cross-basis: |<col_i, col_j>|^2 should be 1/6.
    Intra-basis: |<col_i, col_j>|^2 should be 0 for i≠j.
    """
    # Normalize columns
    H_cols = H / np.sqrt(N)
    B3_cols = B3 / np.sqrt(N)
    B4_cols = B4 / np.sqrt(N)

    defect = 0.0

    # Intra-basis orthogonality defects (should be 0)
    for B_name, B_cols in [("H", H_cols), ("B3", B3_cols), ("B4", B4_cols)]:
        for i in range(N):
            for j in range(i + 1, N):
                dot = np.dot(np.conj(B_cols[:, i]), B_cols[:, j])
                defect += abs(dot) ** 2

    # Cross-basis unbiasedness defects (should be 1/6)
    bases = [H_cols, B3_cols, B4_cols]
    for a in range(3):
        for b in range(a + 1, 3):
            for i in range(N):
                for j in range(N):
                    dot = np.dot(np.conj(bases[a][:, i]), bases[b][:, j])
                    defect += abs(abs(dot) ** 2 - TARGET_CROSS)

    return defect


def defect_from_params(theta_phi_lambda, b3_params, b4_params):
    """Total defect from combined parameter vector."""
    theta, phi, lam = theta_phi_lambda
    H = build_k6(theta, phi, lam)
    B3 = params_to_B(b3_params)
    B4 = params_to_B(b4_params)

    # Penalize non-unitarity of B3, B4
    defect = compute_defect(H, B3, B4)
    defect += 100.0 * unitary_defect(B3)
    defect += 100.0 * unitary_defect(B4)

    return defect


# ── Objective wrappers for scipy ───────────────────────────────────────────────

def objective_differential_evolution(x):
    """x = [theta, phi, lambda, b3_params..., b4_params...]"""
    theta_phi_lambda = x[:3]
    b3_params = x[3:3 + 36]
    b4_params = x[3 + 36:3 + 36 + 36]
    val = defect_from_params(theta_phi_lambda, b3_params, b4_params)
    return val


def objective_lbfgsb(x):
    return objective_differential_evolution(x)


# ── Main search ────────────────────────────────────────────────────────────────

def run_search():
    print("=" * 70)
    print("PHASE 5 — PRONG A: Defect Minimization for 4th MUB")
    print("=" * 70)
    print(f"Dimension N = {N}")
    print(f"Random restarts: {NP_RESTARTS}")
    print(f"Target: defect = 0 (4th MUB exists)")
    print()

    # Bounds for theta, phi, lambda
    # theta in [0, pi], phi in [0, 2*pi], lambda in [0, 2*pi]
    theta_bounds = (0.0, np.pi)
    phi_bounds = (0.0, 2.0 * np.pi)
    lam_bounds = (0.0, 2.0 * np.pi)

    # Bounds for unitary parameters: reasonable range
    unitary_bounds = (-5.0, 5.0)

    bounds = [theta_bounds, phi_bounds, lam_bounds] + \
             [unitary_bounds] * 36 + \
             [unitary_bounds] * 36

    best_defect = np.inf
    best_params = None
    near_miss_counts = {1e-6: 0, 1e-3: 0, 1e-2: 0, 1e-1: 0, 1.0: 0}

    all_results = []

    np.random.seed(SEED)

    start_time = time.time()

    for restart in range(NP_RESTARTS):
        if restart % 20 == 0:
            elapsed = time.time() - start_time
            print(f"  Restart {restart}/{NP_RESTARTS} (elapsed {elapsed:.1f}s) "
                  f"best defect so far: {best_defect:.6e}")

        # Random initial point
        theta0 = np.random.uniform(*theta_bounds)
        phi0 = np.random.uniform(*phi_bounds)
        lam0 = np.random.uniform(*lam_bounds)
        b3_init = np.random.uniform(*unitary_bounds, size=36)
        b4_init = np.random.uniform(*unitary_bounds, size=36)
        x0 = np.concatenate([[theta0, phi0, lam0], b3_init, b4_init])

        # Local optimization with L-BFGS-B
        try:
            result = minimize(
                objective_lbfgsb,
                x0,
                method='L-BFGS-B',
                bounds=bounds,
                options={'maxiter': 500, 'ftol': 1e-15, 'gtol': 1e-12}
            )
            defect = result.fun
            params = result.x
        except Exception as e:
            # fallback: evaluate at initial point
            defect = objective_differential_evolution(x0)
            params = x0

        all_results.append({
            'restart': restart,
            'defect': float(defect),
            'success': result.success if hasattr(result, 'success') else False,
        })

        if defect < best_defect:
            best_defect = defect
            best_params = params.copy()

        # Near-miss classification
        for threshold in sorted(near_miss_counts.keys()):
            if defect < threshold:
                near_miss_counts[threshold] += 1
                break

        # If we hit a very low defect, report immediately
        if defect < 1e-10:
            print(f"\n*** CANDIDATE FOUND at restart {restart} ***")
            print(f"    defect = {defect:.6e}")
            print(f"    theta={params[0]:.6f}, phi={params[1]:.6f}, lambda={params[2]:.6f}")

    elapsed = time.time() - start_time
    print()
    print(f"Search complete. Total time: {elapsed:.1f}s")
    print()

    # ── Report best result ────────────────────────────────────────────────────
    print("─" * 70)
    print("BEST RESULT")
    print("─" * 70)
    print(f"  Best defect: {best_defect:.6e}")

    if best_params is not None:
        theta, phi, lam = best_params[0], best_params[1], best_params[2]
        b3_params_vec = best_params[3:3 + 36]
        b4_params_vec = best_params[3 + 36:3 + 36 + 36]

        H_best = build_k6(theta, phi, lam)
        B3_best = params_to_B(b3_params_vec)
        B4_best = params_to_B(b4_params_vec)

        print(f"  theta={theta:.10f}")
        print(f"  phi={phi:.10f}")
        print(f"  lambda={lam:.10f}")

        print(f"  H defect (Hadamardness): {np.linalg.norm(H_best @ H_best.conj().T - N * np.eye(N)):.6e}")
        print(f"  B3 unitarity defect: {unitary_defect(B3_best):.6e}")
        print(f"  B4 unitarity defect: {unitary_defect(B4_best):.6e}")
        print(f"  Intra-basis orthogonality defects:")
        print(f"    H: {np.linalg.norm(np.dot(H_best.conj().T, H_best) - N * np.eye(N)):.6e}")
        print(f"    B3: {np.linalg.norm(np.dot(B3_best.conj().T, B3_best) - N * np.eye(N)):.6e}")
        print(f"    B4: {np.linalg.norm(np.dot(B4_best.conj().T, B4_best) - N * np.eye(N)):.6e}")

        # For cross-basis check on best result, compute per-pair
        H_cols = H_best / np.sqrt(N)
        B3_cols = B3_best / np.sqrt(N)
        B4_cols = B4_best / np.sqrt(N)
        bases = [("H", H_cols), ("B3", B3_cols), ("B4", B4_cols)]
        print(f"  Cross-basis |<v,w>|^2 - 1/6 statistics:")
        for a_name, a_cols in bases:
            for b_name, b_cols in bases:
                if a_name >= b_name:
                    continue
                diffs = []
                for i in range(N):
                    for j in range(N):
                        dot = np.dot(np.conj(a_cols[:, i]), b_cols[:, j])
                        diffs.append(abs(abs(dot) ** 2 - TARGET_CROSS))
                diffs = np.array(diffs)
                print(f"    {a_name}-{b_name}: max={diffs.max():.6e}, "
                      f"mean={diffs.mean():.6e}, median={np.median(diffs):.6e}")

    print()
    print("─" * 70)
    print("NEAR-MISS STATISTICS")
    print("─" * 70)
    for threshold in sorted(near_miss_counts.keys()):
        print(f"  Restarts with defect < {threshold:.0e}: {near_miss_counts[threshold]} / {NP_RESTARTS}")

    print()
    print("─" * 70)
    print("INTERPRETATION")
    print("─" * 70)
    if best_defect < 1e-12:
        print("  *** CANDIDATE 4TH MUB FOUND — needs exact verification ***")
    elif best_defect < 1e-6:
        print("  Near-defect found: best defect is < 1e-6 but > 1e-12.")
        print("  This could indicate a 4th MUB very close to existing parameters.")
    elif best_defect < 0.01:
        print("  Moderate defect: best defect < 0.01. No clear 4th MUB.")
    else:
        print(f"  High defect: best defect = {best_defect:.4f}.")
        print("  No 4th MUB found within Karlsson family with this method.")

    # ── Save results ──────────────────────────────────────────────────────────
    report = {
        'best_defect': float(best_defect),
        'best_theta': float(theta) if best_params is not None else None,
        'best_phi': float(phi) if best_params is not None else None,
        'best_lambda': float(lam) if best_params is not None else None,
        'near_miss': near_miss_counts,
        'n_restarts': NP_RESTARTS,
        'n_candidates': int(best_defect < 1e-12),
        'elapsed_seconds': float(elapsed),
    }

    out_path = Path(r"D:/MUBs in 6-dimension/results/phase5_prongA_defect_minimization.json")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(report, f, indent=2)
    print(f"\nResults saved to {out_path}")

    return report


if __name__ == "__main__":
    run_search()
