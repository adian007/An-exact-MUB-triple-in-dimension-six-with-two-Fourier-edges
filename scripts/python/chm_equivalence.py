"""
Phase 1.1b: CHM equivalence search — column permutation + diagonal phases.

Two 6x6 CHMs H1, H2 are equivalent if
    H1 ≈ D_r @ H2[:, perm] @ D_c
for diagonal unitary D_r, D_c and column permutation perm.

Also tests H.T, conj(H), conj(H).T variants (4 total per pair).

Dita D0 vs Karlsson (theta_D, pi/4, lambda): scripts/python/identify_d0_in_karlsson.py
"""

import itertools
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "python"))

from hadamard6 import dephase  # noqa: E402
from karlsson_k6_3 import build_k6  # noqa: E402

TOL_EQUIV = 1e-10
N = 6


def _phase(z):
    if abs(z) < 1e-15:
        return 1.0 + 0j
    return z / abs(z)


def residual_for_perm(H1, H2, perm):
    """
    Given dephased H1, H2 and column permutation perm (tuple of length N),
    find optimal diagonal phases D_r, D_c minimizing ||H1 - D_r @ H2[:,perm] @ D_c||_F.
    Fix D_c[0] = 1; set D_r[i] = 1/conj(H2[i, perm[0]]) so (i,0) matches.
    """
    P = np.array(perm)
    H2p = H2[:, P]

    D_r = np.ones(N, dtype=complex)
    for i in range(N):
        if abs(H2p[i, 0]) > 1e-15:
            D_r[i] = _phase(H1[i, 0] / H2p[i, 0])
        else:
            D_r[i] = 1.0

    D_c = np.ones(N, dtype=complex)
    for j in range(1, N):
        ratios = []
        for i in range(N):
            if abs(H2p[i, j]) > 1e-15 and abs(D_r[i]) > 1e-15:
                ratios.append(H1[i, j] / (D_r[i] * H2p[i, j]))
        if ratios:
            r0 = ratios[0]
            D_c[j] = np.mean([_phase(r / r0) for r in ratios]) * r0
            D_c[j] = _phase(D_c[j]) if abs(D_c[j]) > 1e-15 else 1.0
        else:
            D_c[j] = 1.0

    H_fit = D_r[:, None] * H2p * D_c[None, :]
    res = np.linalg.norm(H1 - H_fit, "fro")
    return res, perm, D_r, D_c


def chm_equivalence_residual(H1, H2, tol=TOL_EQUIV):
    """
    Minimum Frobenius residual over 720 column permutations and 4 matrix variants.
    Returns dict with best residual and transform info.
    """
    best = {"residual": np.inf, "variant": None, "perm": None, "row_perm": None}

    def variants(H):
        return [
            ("H", H),
            ("H.T", H.T),
            ("conj(H)", np.conj(H)),
            ("conj(H).T", np.conj(H).T),
        ]

    for vname, H2v in variants(H2):
        D1 = dephase(H1)
        D2 = dephase(H2v)
        for perm in itertools.permutations(range(N)):
            res, p, _, _ = residual_for_perm(D1, D2, perm)
            if res < best["residual"]:
                best = {
                    "residual": res,
                    "variant": vname,
                    "perm": p,
                    "row_perm": tuple(range(N)),
                    "equivalent": res < tol,
                }
    return best


def chm_equivalence_residual_full(H1, H2, tol=TOL_EQUIV):
    """
    Full monomial CHM equivalence: row perm + column perm + diagonal phases + variants.
    Used to check D0 ≈ D_bc (literature: block-circulant form of D(0)).
    """
    best = {
        "residual": np.inf,
        "variant": None,
        "perm": None,
        "row_perm": None,
        "equivalent": False,
    }

    def variants(H):
        return [
            ("H", H),
            ("H.T", H.T),
            ("conj(H)", np.conj(H)),
            ("conj(H).T", np.conj(H).T),
        ]

    D1 = dephase(H1)
    for vname, H2v in variants(H2):
        for row_perm in itertools.permutations(range(N)):
            H2r = H2v[np.array(row_perm), :]
            D2 = dephase(H2r)
            for col_perm in itertools.permutations(range(N)):
                res, p, _, _ = residual_for_perm(D1, D2, col_perm)
                if res < best["residual"]:
                    best = {
                        "residual": res,
                        "variant": vname,
                        "perm": p,
                        "row_perm": row_perm,
                        "equivalent": res < tol,
                    }
                    if res < tol:
                        return best
    return best


def build_pair_report(label, theta, phi, lam1, lam2):
    H1 = build_k6(theta, phi, lam1)
    H2 = build_k6(theta, phi, lam2)
    r = chm_equivalence_residual(H1, H2)
    return {
        "label": label,
        "params": (theta, phi, lam1, lam2),
        **r,
    }


def main():
    out = ROOT / "results" / "chm_equivalence.txt"
    dita_theta = float(np.arccos(1 / np.sqrt(3)))

    lines = ["=== Phase 1.1b: CHM equivalence (720 perm + diag phases) ===", ""]

    # Dita lambda pairs
    lines.append("--- Dita (theta=arccos(1/sqrt(3)), phi=pi/4): lambda pairs ---")
    lam_dita = [0.4, 0.41, 0.45, 0.5, 0.9, 1.3, 2.0]
    dita_results = []
    for i, la in enumerate(lam_dita):
        if i == 0:
            continue
        r = build_pair_report(
            f"Dita lambda 0.4 vs {la}",
            dita_theta, np.pi / 4, 0.4, la,
        )
        dita_results.append(r)
        lines.append(
            f"  {r['label']}: residual={r['residual']:.6e}  variant={r['variant']}  "
            f"equiv={r['equivalent']}"
        )

    # F6_theta0 phi pairs
    lines.append("")
    lines.append("--- F6_theta0 (theta=0, lambda=0.3): phi pairs ---")
    phi_f6 = [0.48, 0.5, 0.52, 0.7]
    f6_results = []
    for ph in phi_f6:
        if abs(ph - 0.5) < 1e-12:
            continue
        H1 = build_k6(0.0, 0.5, 0.3)
        H2 = build_k6(0.0, ph, 0.3)
        r = chm_equivalence_residual(H1, H2)
        r["label"] = f"F6 phi 0.5 vs {ph}"
        f6_results.append(r)
        lines.append(
            f"  {r['label']}: residual={r['residual']:.6e}  variant={r['variant']}  "
            f"equiv={r.get('equivalent', r['residual'] < TOL_EQUIV)}"
        )

    # Branch decision
    dita_equiv = all(r["residual"] < TOL_EQUIV for r in dita_results)
    f6_equiv = all(r["residual"] < TOL_EQUIV for r in f6_results)
    min_dita = min(r["residual"] for r in dita_results) if dita_results else np.inf
    max_dita = max(r["residual"] for r in dita_results) if dita_results else 0.0
    max_f6 = max(r["residual"] for r in f6_results) if f6_results else 0.0

    lines.extend([
        "",
        "=== Summary ===",
        f"Dita lambda pairs: min_residual={min_dita:.6e}  max_residual={max_dita:.6e}  all_equiv={dita_equiv}",
        f"F6 phi pairs: max_residual={max_f6:.6e}  all_equiv={f6_equiv}",
        "",
    ])
    if dita_equiv and f6_equiv:
        branch = "1.2 (gauge / parametrization artifact)"
    elif not dita_equiv and f6_equiv:
        branch = "1.3 (genuine 1D lambda family at Dita; phi gauge at theta=0)"
    elif dita_equiv and not f6_equiv:
        branch = "1.2 partial (unexpected: lambda equiv but phi not)"
    else:
        branch = "1.3 (both inequivalent — unlikely)"
    lines.append(f"Branch decision: -> {branch}")

    text = "\n".join(lines) + "\n"
    out.write_text(text, encoding="utf-8")
    print(text)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
