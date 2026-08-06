"""
degeneracy_scan.py — emit targeted (theta, phi, lambda) for special Karlsson loci.

Detects:
  1. Mobius z4^2 consistency violation (z4_dev > 1e-6)
  2. A/B block eigenvalue degeneracy or near-singular det
  3. Dephased H matching a pure circulant pattern (Björck-type)

Output: results/degeneracy_candidates.csv (+ .json) for Julia search_special_loci.jl

Usage: python degeneracy_scan.py
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import numpy as np

from hadamard6 import dephase, is_hadamard
from karlsson_k6_3 import build_A, build_k6, mobius

ROOT = Path(__file__).resolve().parent.parent.parent
OUT_CSV = ROOT / "results" / "degeneracy_candidates.csv"
OUT_JSON = ROOT / "results" / "degeneracy_candidates.json"

F2 = np.array([[1, 1], [1, -1]], dtype=complex)


def block_eigenvalue_disc(M: np.ndarray) -> float:
    """Squared discriminant of 2x2 eigenvalues; 0 => repeated eigenvalues."""
    tr = np.trace(M)
    det = np.linalg.det(M)
    return float(abs(tr * tr - 4 * det))


def mobius_diagnostics(theta: float, phi: float, lam: float) -> dict:
    A = build_A(theta, phi)
    B = -F2 - A
    alpha_A, beta_A = A[0, 1] ** 2, A[0, 0] ** 2
    alpha_B, beta_B = B[0, 1] ** 2, B[0, 0] ** 2
    z1sq = np.exp(2j * lam)
    z3sq = mobius(z1sq, alpha_A, beta_A)
    z4sq = mobius(z1sq, alpha_B, beta_B)
    num = beta_B - z3sq * np.conj(alpha_B)
    den = alpha_B - z3sq * np.conj(beta_B)
    z2sq = num / den
    z4sq_check = mobius(z2sq, alpha_A, beta_A)
    z4_dev = float(abs(z4sq - z4sq_check))
    return {
        "A_disc": block_eigenvalue_disc(A),
        "B_disc": block_eigenvalue_disc(B),
        "A_det": float(abs(np.linalg.det(A))),
        "B_det": float(abs(np.linalg.det(B))),
        "z4_dev": z4_dev,
    }


def circulant_match_score(H: np.ndarray, tol: float = 1e-6) -> float:
    """
    Max phase mismatch between dephased H and its circulant template (row k = roll(row0, k)).
    0 = exact circulant up to global column phases removed by dephase.
    """
    D = dephase(H)
    n = D.shape[0]
    ref = D[0, :]
    worst = 0.0
    for k in range(1, n):
        shifted = np.roll(ref, k)
        for j in range(n):
            a, b = D[k, j], shifted[j]
            if abs(a) < tol or abs(b) < tol:
                continue
            phase_diff = abs(np.angle(a * np.conj(b)))
            worst = max(worst, phase_diff)
    return worst


def scan_point(theta: float, phi: float, lam: float) -> dict | None:
    try:
        diag = mobius_diagnostics(theta, phi, lam)
        H = build_k6(theta, phi, lam)
    except (ValueError, FloatingPointError):
        return None
    if not is_hadamard(H):
        return None
    circ_err = circulant_match_score(H)
    flags = []
    if diag["z4_dev"] > 1e-6:
        flags.append("mobius_z4")
    if diag["A_disc"] < 1e-8 or diag["B_disc"] < 1e-8:
        flags.append("block_eig_degen")
    if diag["A_det"] < 1e-6 or diag["B_det"] < 1e-6:
        flags.append("block_det_small")
    if circ_err < 1e-8:
        flags.append("circulant_match")
    if not flags:
        return None
    return {
        "theta": theta,
        "phi": phi,
        "lambda": lam,
        "flags": flags,
        "z4_dev": diag["z4_dev"],
        "A_disc": diag["A_disc"],
        "B_disc": diag["B_disc"],
        "A_det": diag["A_det"],
        "B_det": diag["B_det"],
        "circulant_err": circ_err,
    }


def candidate_grid() -> list[tuple[float, float, float]]:
    pts: list[tuple[float, float, float]] = []
    # Mobius audit-fail neighborhood
    base = (0.9553166181245092, 0.7853981633974483, 0.4)
    for d in np.linspace(-0.05, 0.05, 11):
        pts.append((base[0] + d, base[1] + d * 0.5, base[2] + d * 0.3))
    # theta=0 Fourier slice
    for phi in np.linspace(0, np.pi / 2, 9):
        for lam in np.linspace(0, 2 * np.pi, 12, endpoint=False):
            pts.append((0.0, float(phi), float(lam)))
    # Dita neighborhood
    dita_t = float(np.arccos(1 / np.sqrt(3)))
    for d in np.linspace(-0.03, 0.03, 7):
        pts.append((dita_t + d, np.pi / 4 + d, 0.4 + d))
    # Structured coarse grid (not the 512 sweep — different resolution)
    for t in np.linspace(0.05, np.pi - 0.05, 15):
        for p in np.linspace(0.05, np.pi - 0.05, 15):
            for l in np.linspace(0.05, 2 * np.pi - 0.05, 8):
                pts.append((float(t), float(p), float(l)))
    # Named special points
    pts.extend([
        (0.0, 0.0, 0.0),
        (0.0, np.pi / 4, np.pi / 6),
        (dita_t, np.pi / 4, 0.4),
        base,
        (np.pi / 6, np.pi / 3, np.pi / 4),
    ])
    # Deduplicate
    seen = set()
    out = []
    for t, p, l in pts:
        key = (round(t, 10), round(p, 10), round(l, 10))
        if key not in seen:
            seen.add(key)
            out.append(key)
    return out


def main() -> None:
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    grid = candidate_grid()
    candidates = []
    for theta, phi, lam in grid:
        row = scan_point(theta, phi, lam)
        if row is not None:
            candidates.append(row)

    # Sort by severity (Mobius dev first)
    candidates.sort(key=lambda r: (-r["z4_dev"], r["circulant_err"], r["A_disc"]))

    fieldnames = [
        "theta", "phi", "lambda", "flags", "z4_dev",
        "A_disc", "B_disc", "A_det", "B_det", "circulant_err",
    ]
    def _csv_row(r: dict) -> dict:
        # Full precision: 10-digit CSV coords can miss sharp third-MUB loci (Dita anchor).
        out = {k: r[k] for k in fieldnames if k != "flags"}
        for k in ("theta", "phi", "lambda"):
            out[k] = format(r[k], ".17g")
        out["flags"] = "|".join(r["flags"])
        return out

    with OUT_CSV.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in candidates:
            w.writerow(_csv_row(r))

    with OUT_JSON.open("w") as f:
        json.dump(candidates, f, indent=2)

    print(f"Scanned {len(grid)} points, found {len(candidates)} degeneracy candidates")
    print(f"  CSV:  {OUT_CSV}")
    print(f"  JSON: {OUT_JSON}")
    by_flag: dict[str, int] = {}
    for r in candidates:
        for fl in r["flags"]:
            by_flag[fl] = by_flag.get(fl, 0) + 1
    for fl, n in sorted(by_flag.items()):
        print(f"    {fl}: {n}")


if __name__ == "__main__":
    main()
