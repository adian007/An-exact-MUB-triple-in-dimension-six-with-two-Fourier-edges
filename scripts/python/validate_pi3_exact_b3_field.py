#!/usr/bin/env python3
"""Validate the a=b=c pi/3 B3 candidate in Q[u]/(P)."""

from __future__ import annotations

import json
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[2]


def main() -> None:
    u = sp.Symbol("u")
    P = sp.Poly(104329*u**24 - 193800*u**22 + 192924*u**20 - 299406*u**18 + 331452*u**16 - 381612*u**14 + 493595*u**12 - 381612*u**10 + 331452*u**8 - 299406*u**6 + 192924*u**4 - 193800*u**2 + 104329, u, domain=sp.QQ)
    C = sp.Integer(610818468996621564036607682)
    A = 1514595399406991554715184449*u**22 - 1756307929240743070046703434*u**20 + 989839703306762314045034317*u**18 - 2568263831168312266729485648*u**16 + 1531431372534934342773255463*u**14 - 2644550445045953692323054841*u**12 + 3958308486166522559203865220*u**10 - 1177917192554426766591014887*u**8 + 2139096543575902440554317075*u**6 - 1591397118696353756136381769*u**4 + 801861270977657426086402045*u**2 - 1811714204362686491473791534

    def red(x):
        return sp.rem(sp.Poly(sp.cancel(x), u, domain=sp.QQ), P).as_expr()

    def mul(x, y):
        return red(x * y)

    def inv(x):
        return sp.invert(sp.Poly(red(x), u, domain=sp.QQ), P).as_expr()

    one = sp.Integer(1)
    q = red(-A / C)
    qi = inv(q)
    ui = inv(u)

    def qp(k):
        base = q if k >= 0 else qi
        result = one
        for _ in range(abs(k) % 36):
            result = mul(result, base)
        return result

    H = [
        [one, one, one, one, one, one],
        [one, -one, qp(6), -qp(6), qp(9), -qp(9)],
        [one, -qp(9), qp(9), qp(9), -qp(9), -one],
        [one, qp(9), -qp(6), qp(6), -one, -qp(9)],
        [one, qp(-6), -qp(9), -one, -qp(-6), qp(9)],
        [one, -qp(-6), -one, -qp(9), qp(-6), qp(9)],
    ]
    Hbar = [[red(entry.subs({q if False else u: u})) for entry in row] for row in H]
    # Construct the conjugate H explicitly by q -> q^-1.
    Hbar = [
        [one, one, one, one, one, one],
        [one, -one, qp(-6), -qp(-6), qp(-9), -qp(-9)],
        [one, -qp(-9), qp(-9), qp(-9), -qp(-9), -one],
        [one, qp(-9), -qp(-6), qp(-6), -one, -qp(-9)],
        [one, qp(6), -qp(-9), -one, -qp(6), qp(-9)],
        [one, -qp(6), -one, -qp(-9), qp(6), qp(-9)],
    ]
    B = [
        [one, one, one, one, one, one],
        [u, -u, u, -u, u, -u],
        [-qp(1)*u, qp(1)*u, qp(7)*u, -qp(7)*u, qp(31)*u, -qp(31)*u],
        [qp(26), qp(26), qp(2), qp(2), qp(14), qp(14)],
        [qp(20)*u, qp(2)*u, qp(32)*u, -qp(32)*u, qp(8)*u, -qp(8)*u],
        [qp(28), qp(28), qp(16), qp(16), qp(4), qp(4)],
    ]
    Bbar = [
        [one, one, one, one, one, one],
        [ui, -ui, ui, -ui, ui, -ui],
        [-qp(-1)*ui, qp(-1)*ui, qp(-7)*ui, -qp(-7)*ui, qp(-31)*ui, -qp(-31)*ui],
        [qp(-26), qp(-26), qp(-2), qp(-2), qp(-14), qp(-14)],
        [qp(-20)*ui, qp(-2)*ui, qp(-32)*ui, -qp(-32)*ui, qp(-8)*ui, -qp(-8)*ui],
        [qp(-28), qp(-28), qp(-16), qp(-16), qp(-4), qp(-4)],
    ]
    failures = []
    for col in range(6):
        left = sum(mul(Bbar[row][0], H[row][col]) if False else 0 for row in range(6))
        inner = sum(mul(Bbar[row][0], H[row][col]) for row in range(6))
        # Bbar[row][0] is not the vector column; use Bbar[row][basis].
        inner = sum(mul(Bbar[row][0], H[row][col]) for row in range(6))
        if red(mul(inner, sum(mul(B[row][0], Hbar[row][col]) for row in range(6))) - 6) != 0:
            failures.append(("H", 0, col))
    for basis in range(6):
        for col in range(6):
            inner = sum(mul(Bbar[row][basis], H[row][col]) for row in range(6))
            check = sum(mul(B[row][basis], Hbar[row][col]) for row in range(6))
            if red(mul(inner, check) - 6) != 0:
                failures.append(("H", basis, col))
    for left in range(6):
        for right in range(left + 1, 6):
            inner = sum(mul(Bbar[row][left], B[row][right]) for row in range(6))
            if red(inner) != 0:
                failures.append(("B", left, right))
    result = {"status": "EXACT_B3_VALIDATED" if not failures else "FAILED", "failures": failures, "field_degree": P.degree()}
    print(json.dumps(result, indent=2))
    (ROOT / "results/phase2_exact_pi3_b3_validation.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
