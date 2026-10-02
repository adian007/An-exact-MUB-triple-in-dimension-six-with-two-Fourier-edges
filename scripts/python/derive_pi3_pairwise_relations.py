#!/usr/bin/env python3
"""Inspect exact pairwise relations for the pi/3 phase-ansatz B3."""

from __future__ import annotations

import sympy as sp


def bar(expr: sp.Expr, q: sp.Symbol, *variables: sp.Symbol) -> sp.Expr:
    substitutions = {q: 1 / q}
    substitutions.update({variable: 1 / variable for variable in variables})
    return sp.expand(expr.xreplace(substitutions))


def main() -> None:
    q, a, b, c = sp.symbols("q a b c")
    phi = sp.cyclotomic_poly(36, q)
    vectors = [
        sp.Matrix([1, a, -q * a, q**-10, q**-16 * a, q**-8]),
        sp.Matrix([1, -a, q * a, q**-10, q**2 * a, q**-8]),
        sp.Matrix([1, b, q**7 * b, q**2, q**-4 * b, q**16]),
        sp.Matrix([1, -b, -q**7 * b, q**2, -q**-4 * b, q**16]),
        sp.Matrix([1, c, q**-5 * c, q**14, q**8 * c, q**4]),
        sp.Matrix([1, -c, -q**-5 * c, q**14, -q**8 * c, q**4]),
    ]
    for left in range(6):
        for right in range(left + 1, 6):
            expression = sp.together(sum(bar(vectors[left][row], q, a, b, c) * vectors[right][row] for row in range(6)))
            numerator = expression.as_numer_denom()[0]
            reduced = sp.rem(sp.Poly(numerator, q), sp.Poly(phi, q)).as_expr()
            print(f"O{left}{right}: {sp.factor(reduced)}")


if __name__ == "__main__":
    main()
