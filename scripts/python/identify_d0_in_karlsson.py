"""
E0: which Karlsson (theta_D, pi/4, lambda) is CHM-equivalent to textbook Dita D0?

Textbook D0 is Bengtsson quant-ph/0610161 eq. (11) at x=0 (fourth roots).
CHM distance uses the same dephase + 720 column permutations + 4 variants
as scripts/python/chm_equivalence.py.

Usage: python scripts/python/identify_d0_in_karlsson.py
Output: results/identify_d0_in_karlsson.txt
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "python"))

from chm_equivalence import (  # noqa: E402
    chm_equivalence_residual,
    chm_equivalence_residual_full,
)
from hadamard6 import dephase, is_hadamard  # noqa: E402
from karlsson_k6_3 import build_k6  # noqa: E402

OUT = ROOT / "results" / "identify_d0_in_karlsson.txt"
DITA_THETA = float(np.arccos(1 / np.sqrt(3)))
DITA_PHI = float(np.pi / 4)
LAM_CERT = (0.0, 0.4, np.pi / 3, 2 * np.pi / 3, np.pi, 4 * np.pi / 3, 5 * np.pi / 3)
TOL_EQUIV = 1e-8


def textbook_dita_d0():
    """Unnormalised D(0), Bengtsson eq. (11) with z=1."""
    i = 1j
    return np.array(
        [
            [1, 1, 1, 1, 1, 1],
            [1, -1, i, -i, -i, i],
            [1, i, -1, i, -i, -i],
            [1, -i, i, -1, i, -i],
            [1, -i, -i, i, -1, i],
            [1, i, -i, -i, i, -1],
        ],
        dtype=complex,
    )


def dita_D_bc():
    """Block-circulant equivalent of D(0), Bengtsson eqs. (79)–(80)."""
    w = np.exp(2j * np.pi / 24)
    C3 = np.array(
        [[1, w**6, w**6], [w**6, 1, w**6], [w**6, w**6, 1]], dtype=complex
    )
    C4 = np.array(
        [[w**15, w**3, w**3], [w**3, w**15, w**3], [w**3, w**3, w**15]],
        dtype=complex,
    )
    return np.block([[C3, C4], [C4, -1j * C3.conj().T]])


def lambda_grid():
    algebraic = [
        0.0,
        np.pi / 12,
        np.pi / 8,
        np.pi / 6,
        np.pi / 4,
        np.pi / 3,
        np.pi / 2,
        2 * np.pi / 3,
        3 * np.pi / 4,
        np.pi,
        4 * np.pi / 3,
        5 * np.pi / 3,
        2 * np.pi,
    ]
    dense = list(np.linspace(0.0, 2 * np.pi, 97, endpoint=True))
    cert = list(LAM_CERT)
    vals = sorted(set(round(float(x), 12) for x in algebraic + dense + cert))
    return vals


def fmt_perm(perm):
    if perm is None:
        return "None"
    return "(" + " ".join(str(p) for p in perm) + ")"


def main() -> int:
    D0 = textbook_dita_d0()
    Dbc = dita_D_bc()
    lines = [
        "=== E0: identify textbook Dita D0 inside Karlsson Dita slice ===",
        "theta_D = arccos(1/sqrt(3)), phi = pi/4",
        "D0 = Bengtsson quant-ph/0610161 eq. (11) at x=0 (unnormalised)",
        "CHM residual: dephase + 720 column perms + {H, H.T, conj, conj.T}",
        f"equivalence tolerance = {TOL_EQUIV:g}",
        "",
    ]

    d0_had = is_hadamard(D0)
    dbc_had = is_hadamard(Dbc)
    lines.append(f"D0 is_hadamard={d0_had}  D_bc is_hadamard={dbc_had}")
    r_dbc = chm_equivalence_residual(D0, Dbc)
    lines.append(
        f"D0 vs D_bc (col-only): residual={r_dbc['residual']:.6e}  "
        f"variant={r_dbc['variant']}  perm={fmt_perm(r_dbc['perm'])}  "
        f"equiv={r_dbc['equivalent']}"
    )
    r_dbc_full = chm_equivalence_residual_full(D0, Dbc, tol=TOL_EQUIV)
    lines.append(
        f"D0 vs D_bc (row+col): residual={r_dbc_full['residual']:.6e}  "
        f"variant={r_dbc_full['variant']}  row_perm={fmt_perm(r_dbc_full['row_perm'])}  "
        f"col_perm={fmt_perm(r_dbc_full['perm'])}  equiv={r_dbc_full['equivalent']}"
    )
    lines.append("")

    best = {"residual": np.inf, "lambda": None}
    cert_rows = []
    matches = []
    lines.append("--- residual(Karlsson(theta_D, pi/4, lambda), D0) ---")
    for lam in lambda_grid():
        H = build_k6(DITA_THETA, DITA_PHI, lam)
        r = chm_equivalence_residual(H, D0)
        rec = {
            "lambda": lam,
            "residual": r["residual"],
            "variant": r["variant"],
            "perm": r["perm"],
            "equivalent": r["residual"] < TOL_EQUIV,
        }
        if rec["residual"] < best["residual"]:
            best = rec
        in_cert = any(abs(lam - c) < 1e-12 for c in LAM_CERT)
        mark = "  LAMBDA_CERT" if in_cert else ""
        if rec["equivalent"] or in_cert or rec["residual"] < 1e-3:
            lines.append(
                f"  lambda={lam:.12g}  residual={rec['residual']:.6e}  "
                f"variant={rec['variant']}  perm={fmt_perm(rec['perm'])}  "
                f"equiv={rec['equivalent']}{mark}"
            )
        if rec["equivalent"]:
            matches.append(rec)
        if in_cert:
            cert_rows.append(rec)

    lines.append("")
    lines.append("=== Best match ===")
    lines.append(
        f"lambda={best['lambda']:.12g}  residual={best['residual']:.6e}  "
        f"variant={best['variant']}  perm={fmt_perm(best['perm'])}  "
        f"equiv={best['equivalent']}"
    )
    lines.append(f"n_matches (residual < {TOL_EQUIV:g}): {len(matches)}")
    for m in matches:
        lines.append(
            f"  match lambda={m['lambda']:.12g}  variant={m['variant']}  "
            f"perm={fmt_perm(m['perm'])}"
        )

    lines.append("")
    lines.append("=== T3 Lambda_cert vs D0 ===")
    in_cert_match = False
    for rec in cert_rows:
        lines.append(
            f"  lambda={rec['lambda']:.12g}  residual={rec['residual']:.6e}  "
            f"equiv={rec['equivalent']}"
        )
        if rec["equivalent"]:
            in_cert_match = True

    lines.append("")
    if in_cert_match:
        lines.append(
            "VERDICT: a Lambda_cert point is CHM-equivalent to textbook D0. "
            "T3 at that point is a corollary of Brierley–Weigert PRA 79, 052316 "
            "(arXiv:0901.4051), not a new exact-coefficient theorem."
        )
    else:
        lines.append(
            "VERDICT: no Lambda_cert point matched D0 at the stated tolerance. "
            "T3 is not identified as a BW 2009 corollary via this CHM test."
        )

    deph_note = dephase(D0)
    lines.append("")
    lines.append(
        f"dephased D0 first row (should be 1s): {np.round(deph_note[0], 6).tolist()}"
    )

    text = "\n".join(lines) + "\n"
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(text, encoding="utf-8")
    print(text)
    print(f"Wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
