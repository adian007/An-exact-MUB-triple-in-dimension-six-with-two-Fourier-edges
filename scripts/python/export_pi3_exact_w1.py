#!/usr/bin/env python3
"""Deprecated exploratory export for the invalid a=b=c specialization.

Use export_pi3_exact_abc_w1.py for the corrected distinct-root system.
"""

from __future__ import annotations

from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "symbolic_export" / "w1_pi3_exact_ansatz.m2"


def qpow(q: sp.Symbol, exponent: int) -> sp.Expr:
    return q ** (exponent % 36)


def conjugate_matrix(matrix: list[list[sp.Expr]], q: sp.Symbol, u: sp.Symbol, ui: sp.Symbol):
    return [[sp.sympify(entry).xreplace({q: q**35, u: ui}) for entry in row] for row in matrix]


def m2(expr: sp.Expr) -> str:
    return sp.sstr(sp.expand(expr)).replace("**", "^")


def main() -> None:
    q, u, ui = sp.symbols("q u ui")
    z1, z2, z3, z4, z5 = sp.symbols("z1 z2 z3 z4 z5")
    w1, w2, w3, w4, w5 = sp.symbols("w1 w2 w3 w4 w5")
    x = [sp.Integer(1), z1, z2, z3, z4, z5]
    y = [sp.Integer(1), w1, w2, w3, w4, w5]

    H = [
        [1, 1, 1, 1, 1, 1],
        [1, -1, qpow(q, 6), -qpow(q, 6), qpow(q, 9), -qpow(q, 9)],
        [1, -qpow(q, 9), qpow(q, 9), qpow(q, 9), -qpow(q, 9), -1],
        [1, qpow(q, 9), -qpow(q, 6), qpow(q, 6), -1, -qpow(q, 9)],
        [1, qpow(q, 30), -qpow(q, 9), -1, -qpow(q, 30), qpow(q, 9)],
        [1, -qpow(q, 30), -1, -qpow(q, 9), qpow(q, 30), qpow(q, 9)],
    ]
    B3 = [
        [1, 1, 1, 1, 1, 1],
        [u, -u, u, -u, u, -u],
        [-qpow(q, 1) * u, qpow(q, 1) * u, qpow(q, 7) * u, -qpow(q, 7) * u, qpow(q, 31) * u, -qpow(q, 31) * u],
        [qpow(q, 26), qpow(q, 26), qpow(q, 2), qpow(q, 2), qpow(q, 14), qpow(q, 14)],
        [qpow(q, 20) * u, qpow(q, 2) * u, qpow(q, 32) * u, -qpow(q, 32) * u, qpow(q, 8) * u, -qpow(q, 8) * u],
        [qpow(q, 28), qpow(q, 28), qpow(q, 16), qpow(q, 16), qpow(q, 4), qpow(q, 4)],
    ]
    Hbar = conjugate_matrix(H, q, u, ui)
    Bbar = conjugate_matrix(B3, q, u, ui)

    equations = [z1 * w1 - 1, z2 * w2 - 1, z3 * w3 - 1, z4 * w4 - 1, z5 * w5 - 1]
    for matrix, conjugate in ((H, Hbar), (B3, Bbar)):
        for column in range(6):
            left = sum(y[row] * matrix[row][column] for row in range(6))
            right = sum(x[row] * conjugate[row][column] for row in range(6))
            equations.append(sp.expand(left * right - 6))

    P = (
        104329 * u**24 - 193800 * u**22 + 192924 * u**20 - 299406 * u**18
        + 331452 * u**16 - 381612 * u**14 + 493595 * u**12
        - 381612 * u**10 + 331452 * u**8 - 299406 * u**6
        + 192924 * u**4 - 193800 * u**2 + 104329
    )
    # This is the q-linear member of results/pi3_phase_groebner_basis.txt.
    L = (
        610818468996621564036607682 * q
        + 1514595399406991554715184449 * u**22
        - 1756307929240743070046703434 * u**20
        + 989839703306762314045034317 * u**18
        - 2568263831168312266729485648 * u**16
        + 1531431372534934342773255463 * u**14
        - 2644550445045953692323054841 * u**12
        + 3958308486166522559203865220 * u**10
        - 1177917192554426766591014887 * u**8
        + 2139096543575902440554317075 * u**6
        - 1591397118696353756136381769 * u**4
        + 801861270977657426086402045 * u**2
        - 1811714204362686491473791534
    )

    lines = [
        "-- Exact pi/3 W1 ansatz export; no Groebner result is asserted.",
        "-- B3 uses the exact a=b=c=u specialization of the phase ansatz.",
        "kk = toField(QQ[q,u,ui] / ideal(q^12-q^6+1, u*ui-1, P, L));",
        "R = kk[z1,z2,z3,z4,z5,w1,w2,w3,w4,w5];",
        "witness = ideal(",
        "  " + ",\n  ".join(m2(e) for e in equations),
        ");",
        "print(\"numgens witness = \" | toString numgens witness);",
        "elapsedTime G = gens gb witness;",
        "print(\"Groebner basis = \" | toString G);",
        "print(\"1 % witness = \" | toString (1_R % witness));",
        "print(\"dim = \" | toString dim witness);",
    ]
    # Put field polynomials after the ring declaration so M2 can parse them.
    lines.insert(2, f"P = {m2(P)};")
    lines.insert(3, f"L = {m2(L)};")
    OUTPUT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {OUTPUT} with {len(equations)} witness equations")


if __name__ == "__main__":
    main()
