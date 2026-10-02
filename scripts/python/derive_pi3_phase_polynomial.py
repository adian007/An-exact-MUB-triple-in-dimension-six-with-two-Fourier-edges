#!/usr/bin/env python3
"""Derive exact MU equations for the observed pi/3 phase ansatz."""

from __future__ import annotations

import sympy as sp


def conjugate(expr: sp.Expr, q: sp.Symbol, u: sp.Symbol) -> sp.Expr:
    return sp.expand(expr.xreplace({q: 1 / q, u: 1 / u}))


def main() -> None:
    q, u = sp.symbols("q u")
    phi36 = sp.cyclotomic_poly(36, q)
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
    # Observed first vector phases:
    # [1, u, -q*u, q^-10, q^-16*u, q^-8], q=exp(i*pi/18).
    v = sp.Matrix([1, u, -q * u, q**-10, q**-16 * u, q**-8])
    equations = []
    for column in range(6):
        inner = sum(conjugate(v[row], q, u) * H[row, column] for row in range(6))
        expression = sp.together(sp.expand(inner * conjugate(inner, q, u) - 6))
        numerator = sp.Poly(expression.as_numer_denom()[0], q, u).as_expr()
        equations.append(sp.factor(numerator))

    print("phi36(q)=", phi36)
    for index, equation in enumerate(equations):
        print(f"MU{index}: terms={len(sp.Poly(equation, q, u).terms())}")
        print(sp.factor(equation))
    # Eliminate q for each equation; factorization reveals candidate u-fields.
    resultants = []
    for index, equation in enumerate(equations):
        resultant = sp.resultant(equation, phi36, q)
        resultants.append(sp.Poly(resultant, u, domain=sp.QQ))
        print(f"RESULTANT{index}:")
        print(sp.factor(resultant))

    common = resultants[0]
    for candidate in resultants[1:]:
        common = sp.gcd(common, candidate)
    print("COMMON_GCD:")
    print(sp.factor(common.as_expr()))
    print("JOINT_ELIMINATION:")
    ideal = sp.groebner([phi36, *equations], q, u, order="lex")
    print("FULL_BASIS:")
    for polynomial in ideal.polys:
        print(sp.factor(polynomial.as_expr()))
    from pathlib import Path
    output_path = Path(__file__).resolve().parents[2] / "results" / "pi3_phase_groebner_basis.txt"
    output_path.write_text(
        "Exact pi/3 phase ansatz Groebner basis\n"
        + "\n\n".join(str(sp.factor(polynomial.as_expr())) for polynomial in ideal.polys)
        + "\n",
        encoding="utf-8",
    )
    print("wrote", output_path)
    for polynomial in ideal.polys:
        if not polynomial.as_expr().has(q):
            print(sp.factor(polynomial.as_expr()))


if __name__ == "__main__":
    main()
