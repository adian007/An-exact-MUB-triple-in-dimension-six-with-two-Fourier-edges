#!/usr/bin/env python3
"""Exact Fourier-family check using arithmetic in Q[u]/(P)."""

from __future__ import annotations

import json
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[2]


def main() -> None:
    u = sp.Symbol("u")
    P = sp.Poly(
        104329 * u**24 - 193800 * u**22 + 192924 * u**20 - 299406 * u**18
        + 331452 * u**16 - 381612 * u**14 + 493595 * u**12
        - 381612 * u**10 + 331452 * u**8 - 299406 * u**6
        + 192924 * u**4 - 193800 * u**2 + 104329,
        u,
        domain=sp.QQ,
    )
    C = sp.Integer(610818468996621564036607682)
    A = (
        1514595399406991554715184449 * u**22
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

    def red(value: sp.Expr) -> sp.Expr:
        return sp.rem(sp.Poly(sp.cancel(value), u, domain=sp.QQ), P).as_expr()

    def add(left: sp.Expr, right: sp.Expr) -> sp.Expr:
        return red(left + right)

    def neg(value: sp.Expr) -> sp.Expr:
        return red(-value)

    def mul(left: sp.Expr, right: sp.Expr) -> sp.Expr:
        return red(left * right)

    def inv(value: sp.Expr) -> sp.Expr:
        return sp.invert(sp.Poly(red(value), u, domain=sp.QQ), P).as_expr()

    one = sp.Integer(1)
    zero = sp.Integer(0)
    q = red(-A / C)
    qi = inv(q)
    ui = inv(u)

    def qp(exponent: int) -> sp.Expr:
        result = one
        base = q if exponent >= 0 else qi
        for _ in range(abs(exponent) % 36):
            result = mul(result, base)
        return result

    def matmul(left, right):
        return [
            [red(sum((mul(left[i][k], right[k][j]) for k in range(6)), zero)) for j in range(6)]
            for i in range(6)
        ]

    H = [
        [one, one, one, one, one, one],
        [one, neg(one), qp(6), neg(qp(6)), qp(9), neg(qp(9))],
        [one, neg(qp(9)), qp(9), qp(9), neg(qp(9)), neg(one)],
        [one, qp(9), neg(qp(6)), qp(6), neg(one), neg(qp(9))],
        [one, qp(-6), neg(qp(9)), neg(one), neg(qp(-6)), qp(9)],
        [one, neg(qp(-6)), neg(one), neg(qp(9)), qp(-6), qp(9)],
    ]
    Hbar = [
        [one, one, one, one, one, one],
        [one, neg(one), qp(-6), neg(qp(-6)), qp(-9), neg(qp(-9))],
        [one, neg(qp(-9)), qp(-9), qp(-9), neg(qp(-9)), neg(one)],
        [one, qp(-9), neg(qp(-6)), qp(-6), neg(one), neg(qp(-9))],
        [one, qp(6), neg(qp(-9)), neg(one), neg(qp(6)), qp(-9)],
        [one, neg(qp(6)), neg(one), neg(qp(-9)), qp(6), qp(-9)],
    ]
    B = [
        [one, one, one, one, one, one],
        [u, neg(u), u, neg(u), u, neg(u)],
        [neg(mul(qp(1), u)), mul(qp(1), u), mul(qp(7), u), neg(mul(qp(7), u)), mul(qp(31), u), neg(mul(qp(31), u))],
        [qp(26), qp(26), qp(2), qp(2), qp(14), qp(14)],
        [mul(qp(20), u), mul(qp(2), u), mul(qp(32), u), neg(mul(qp(32), u)), mul(qp(8), u), neg(mul(qp(8), u))],
        [qp(28), qp(28), qp(16), qp(16), qp(4), qp(4)],
    ]
    T = matmul(list(map(list, zip(*Hbar))), B)
    rows = (1, 2, 5, 4, 3, 0)
    cols = (1, 5, 3, 0, 4, 2)
    chart = [[T[rows[i]][cols[j]] for j in range(6)] for i in range(6)]
    cross = [[mul(mul(chart[i][j], chart[0][0]), inv(mul(chart[i][0], chart[0][j]))) for j in range(6)] for i in range(6)]
    r = qp(6)
    base = [
        [one, one, one, one, one, one],
        [one, r, mul(r, r), qp(18), qp(24), qp(30)],
        [one, mul(r, r), qp(24), one, mul(r, r), qp(24)],
        [one, qp(18), one, qp(18), one, qp(18)],
        [one, qp(24), mul(r, r), one, qp(24), mul(r, r)],
        [one, qp(30), qp(24), qp(18), mul(r, r), r],
    ]
    param_a = mul(cross[1][1], inv(base[1][1]))
    param_b = mul(cross[1][2], inv(base[1][2]))
    predicted = [row[:] for row in base]
    for row, col in ((1, 1), (1, 4), (3, 1), (3, 4), (5, 1), (5, 4)):
        predicted[row][col] = mul(predicted[row][col], param_a)
    for row, col in ((1, 2), (1, 5), (3, 2), (3, 5), (5, 2), (5, 5)):
        predicted[row][col] = mul(predicted[row][col], param_b)
    failures = [(i, j) for i in range(6) for j in range(6) if red(cross[i][j] - predicted[i][j]) != 0]
    result = {
        "status": "EXACT_FOURIER_IDENTITY" if not failures else "FAILED",
        "failures": failures,
        "field_degree": P.degree(),
        "chart_rows": list(rows),
        "chart_columns": list(cols),
    }
    print(json.dumps(result, indent=2))
    output = ROOT / "results/phase2_exact_pi3_fourier_symbolic.json"
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {output}")


if __name__ == "__main__":
    main()
