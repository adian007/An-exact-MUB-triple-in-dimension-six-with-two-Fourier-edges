#!/usr/bin/env python3
"""Check exact Fourier-family membership for the pi/3 phase candidate."""

from __future__ import annotations

import json
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[2]


def qp(q: sp.Symbol, exponent: int) -> sp.Expr:
    return q ** (exponent % 36)


def conjugate(expr: sp.Expr, q: sp.Symbol, u: sp.Symbol) -> sp.Expr:
    return sp.together(expr.subs({q: 1 / q, u: 1 / u}))


def main() -> None:
    q, u = sp.symbols("q u")
    z = qp(q, 6)
    I = qp(q, 9)
    H = sp.Matrix(
        [
            [1, 1, 1, 1, 1, 1],
            [1, -1, z, -z, I, -I],
            [1, -I, I, I, -I, -1],
            [1, I, -z, z, -1, -I],
            [1, qp(q, 30), -I, -1, -qp(q, 30), I],
            [1, -qp(q, 30), -1, -I, qp(q, 30), I],
        ]
    )
    B = sp.Matrix(
        [
            [1, 1, 1, 1, 1, 1],
            [u, -u, u, -u, u, -u],
            [-q * u, q * u, qp(q, 7) * u, -qp(q, 7) * u, qp(q, 31) * u, -qp(q, 31) * u],
            [qp(q, 26), qp(q, 26), qp(q, 2), qp(q, 2), qp(q, 14), qp(q, 14)],
            [qp(q, 20) * u, qp(q, 2) * u, qp(q, 32) * u, -qp(q, 32) * u, qp(q, 8) * u, -qp(q, 8) * u],
            [qp(q, 28), qp(q, 28), qp(q, 16), qp(q, 16), qp(q, 4), qp(q, 4)],
        ]
    )
    # T = H^dagger B; normalization cancels in cross-ratios.
    Hbar = H.applyfunc(lambda entry: conjugate(entry, q, u))
    T = Hbar.T * B
    rows = (1, 2, 5, 4, 3, 0)
    cols = (1, 5, 3, 0, 4, 2)
    chart = sp.Matrix(6, 6, lambda i, j: T[rows[i], cols[j]])
    cross = sp.Matrix(6, 6, lambda i, j: sp.together(chart[i, j] * chart[0, 0] / (chart[i, 0] * chart[0, j])))
    r = qp(q, 6)
    base = sp.Matrix(
        [
            [1, 1, 1, 1, 1, 1],
            [1, r, r**2, r**3, r**4, r**5],
            [1, r**2, r**4, 1, r**2, r**4],
            [1, r**3, 1, r**3, 1, r**3],
            [1, r**4, r**2, 1, r**4, r**2],
            [1, r**5, r**4, r**3, r**2, r],
        ]
    )
    A = sp.together(cross[1, 1] / base[1, 1])
    C = sp.together(cross[1, 2] / base[1, 2])
    predicted = base.copy()
    for row, col in ((1, 1), (1, 4), (3, 1), (3, 4), (5, 1), (5, 4)):
        predicted[row, col] *= A
    for row, col in ((1, 2), (1, 5), (3, 2), (3, 5), (5, 2), (5, 5)):
        predicted[row, col] *= C

    P = 104329 * u**24 - 193800 * u**22 + 192924 * u**20 - 299406 * u**18 + 331452 * u**16 - 381612 * u**14 + 493595 * u**12 - 381612 * u**10 + 331452 * u**8 - 299406 * u**6 + 192924 * u**4 - 193800 * u**2 + 104329
    L = 610818468996621564036607682 * q + 1514595399406991554715184449 * u**22 - 1756307929240743070046703434 * u**20 + 989839703306762314045034317 * u**18 - 2568263831168312266729485648 * u**16 + 1531431372534934342773255463 * u**14 - 2644550445045953692323054841 * u**12 + 3958308486166522559203865220 * u**10 - 1177917192554426766591014887 * u**8 + 2139096543575902440554317075 * u**6 - 1591397118696353756136381769 * u**4 + 801861270977657426086402045 * u**2 - 1811714204362686491473791534
    q_coefficient = sp.Poly(L, q).coeff_monomial(q)
    q_remainder = sp.Poly(L, q).coeff_monomial(1)
    q_as_function_of_u = -q_remainder / q_coefficient
    P_poly = sp.Poly(P, u)
    failures = []
    for i in range(6):
        for j in range(6):
            numerator = sp.together(cross[i, j] - predicted[i, j]).as_numer_denom()[0]
            substituted = sp.together(numerator.subs(q, q_as_function_of_u)).as_numer_denom()[0]
            remainder = sp.rem(sp.Poly(substituted, u), P_poly).as_expr()
            if remainder != 0:
                failures.append((i, j, sp.factor(remainder)))
    result = {
        "status": "EXACT_FOURIER_IDENTITY" if not failures else "FAILED",
        "chart_rows": list(rows),
        "chart_columns": list(cols),
        "remainder_failures": len(failures),
        "field_relations": ["q^12-q^6+1", "P(u)", "L(q,u)"],
    }
    print(json.dumps(result, indent=2))
    output = ROOT / "results/phase2_exact_pi3_fourier_symbolic.json"
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {output}")


if __name__ == "__main__":
    main()
