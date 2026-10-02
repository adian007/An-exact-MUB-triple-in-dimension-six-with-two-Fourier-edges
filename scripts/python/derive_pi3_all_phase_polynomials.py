#!/usr/bin/env python3
"""Derive exact phase polynomials for the three observed pi/3 vector pairs."""

from __future__ import annotations

import sympy as sp


def conj(expr: sp.Expr, q: sp.Symbol, u: sp.Symbol) -> sp.Expr:
    return sp.expand(expr.xreplace({q: 1 / q, u: 1 / u}))


def phase_polynomial(q: sp.Symbol, u: sp.Symbol, H: sp.Matrix, vector: sp.Matrix):
    phi = sp.cyclotomic_poly(36, q)
    equations = []
    for column in range(6):
        inner = sum(conj(vector[row], q, u) * H[row, column] for row in range(6))
        equations.append(sp.Poly(sp.together(inner * conj(inner, q, u) - 6).as_numer_denom()[0], q, u).as_expr())
    basis = sp.groebner([phi, *equations], q, u, order="lex")
    return [sp.factor(poly.as_expr()) for poly in basis.polys]


def main() -> None:
    q, u = sp.symbols("q u")
    z = q**6
    I = q**9
    H = sp.Matrix(
        [
            [1, 1, 1, 1, 1, 1],
            [1, -1, z, -z, I, -I],
            [1, -I, I, I, -I, -1],
            [1, I, -z, z, -1, -I],
            [1, q**-6, -I, -1, -q**-6, I],
            [1, -q**-6, -1, -I, q**-6, I],
        ]
    )
    ansatz = {
        "pair_0": [1, u, -q * u, q**-10, q**-16 * u, q**-8],
        "pair_1": [1, u, q**7 * u, q**2, q**-4 * u, q**16],
        "pair_2": [1, u, q**-5 * u, q**14, q**8 * u, q**4],
    }
    for label, entries in ansatz.items():
        print(label)
        for polynomial in phase_polynomial(q, u, H, sp.Matrix(entries)):
            print(polynomial)


if __name__ == "__main__":
    main()
