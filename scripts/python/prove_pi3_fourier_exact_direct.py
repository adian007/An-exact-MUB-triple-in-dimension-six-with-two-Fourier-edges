#!/usr/bin/env python3
"""Exact direct-chart certificate for the pi/3 B3 Fourier membership claim.

The matrix entries are symbolic over Q(q,a,b,c), with q^12-q^6+1=0.
The chart identity is checked after clearing its nonzero monomial
denominators. The same cyclotomic quotient verifies |a|=|b|=|c|=1 from
the exact quadratic coefficient-field relations recorded in the W1 audit.
"""

from __future__ import annotations

import json
import platform
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[2]


q, a, b, c = sp.symbols("q a b c")
PHI36 = sp.Poly(q**12 - q**6 + 1, q, domain=sp.QQ)


def qp(exponent: int) -> sp.Expr:
    return q ** (exponent % 36)


def reduce_q(expression: sp.Expr) -> sp.Expr:
    polynomial = sp.Poly(sp.expand(expression), q, domain=sp.QQ)
    return sp.rem(polynomial, PHI36).as_expr()


def reduce_q_over_phase_field(expression: sp.Expr) -> sp.Expr:
    numerator = sp.together(expression).as_numer_denom()[0]
    polynomial = sp.Poly(
        sp.expand(numerator), q, domain=sp.QQ.frac_field(a, b, c)
    )
    return sp.rem(polynomial, PHI36).as_expr()


def conjugate_q(expression: sp.Expr) -> sp.Expr:
    polynomial = sp.Poly(expression, q, domain=sp.QQ)
    return sp.expand(
        sum(coefficient * q ** ((-monomial[0]) % 36)
            for monomial, coefficient in polynomial.terms())
    )


def main() -> None:
    # The exact, unnormalised columns v_0,...,v_5 in pi3_exact_B3_ansatz.md.
    matrix = sp.Matrix(
        [
            [1, 1, 1, 1, 1, 1],
            [a, -a, b, -b, c, -c],
            [-qp(1) * a, qp(1) * a, qp(7) * b, -qp(7) * b,
             qp(31) * c, -qp(31) * c],
            [qp(26), qp(26), qp(2), qp(2), qp(14), qp(14)],
            [qp(20) * a, qp(2) * a, qp(32) * b, -qp(32) * b,
             qp(8) * c, -qp(8) * c],
            [qp(28), qp(28), qp(16), qp(16), qp(4), qp(4)],
        ]
    )

    # One exact Fourier chart. Permutations are zero-based and act on rows
    # and columns of the unnormalised matrix.
    rows = (4, 3, 1, 0, 2, 5)
    columns = (2, 0, 4, 3, 1, 5)
    chart = matrix.extract(rows, columns)

    # Repository convention F_6^(2)(x,y), where q^6 is exp(i*pi/3).
    x = b / (qp(6) * a)
    y = b / (qp(12) * c)
    fourier = sp.Matrix(
        [
            [1, 1, 1, 1, 1, 1],
            [1, qp(6) * x, qp(12) * y, qp(18), qp(24) * x, qp(30) * y],
            [1, qp(12), qp(24), 1, qp(12), qp(24)],
            [1, qp(18) * x, y, qp(18), x, qp(18) * y],
            [1, qp(24), qp(12), 1, qp(24), qp(12)],
            [1, qp(30) * x, qp(24) * y, qp(18), qp(12) * x, qp(6) * y],
        ]
    )

    # C_ij = chart_ij chart_00 / (chart_i0 chart_0j). Check equality
    # entry-by-entry without division, in Q(q,a,b,c)/(Phi_36(q)).
    chart_failures = []
    for i in range(6):
        for j in range(6):
            numerator = (
                chart[i, j] * chart[0, 0]
                - fourier[i, j] * chart[i, 0] * chart[0, j]
            )
            remainder = reduce_q_over_phase_field(numerator)
            if remainder != 0:
                chart_failures.append([i, j, str(remainder)])

    # Exact coefficient-field data: a^2=d_a(q), b^2=d_b(q), c^2=d_c(q).
    coeff_a = (
        16*q**11 - 35*q**10 + 170*q**9 - 18*q**8 + 120*q**7 + 32*q**6
        - 166*q**5 - 157*q**4 + 8*q**3 - 46*q**2 - 48*q - 66
    )
    coeff_b = (
        -166*q**11 + 192*q**10 + 170*q**9 - 46*q**8 - 72*q**7
        + 32*q**6 + 150*q**5 - 35*q**4 + 8*q**3 + 64*q**2
        + 120*q - 66
    )
    coeff_c = (
        150*q**11 - 157*q**10 + 170*q**9 + 64*q**8 - 48*q**7
        + 32*q**6 + 16*q**5 + 192*q**4 + 8*q**3 - 18*q**2
        - 72*q - 66
    )
    square_roots = {
        "a": -coeff_a / 323,
        "b": -coeff_b / 323,
        "c": -coeff_c / 323,
    }
    modulus_failures = {
        name: str(reduce_q(value * conjugate_q(value) - 1))
        for name, value in square_roots.items()
        if reduce_q(value * conjugate_q(value) - 1) != 0
    }

    result = {
        "status": "PASS" if not chart_failures and not modulus_failures else "FAIL",
        "python_version": platform.python_version(),
        "sympy_version": sp.__version__,
        "inputs": [
            "symbolic_export/pi3_exact_B3_ansatz.md",
            "docs/research/reports/phase2_w1_coefficient_field_2026-10-01.md",
        ],
        "field_relation": "q^12-q^6+1=0",
        "row_permutation": list(rows),
        "column_permutation": list(columns),
        "fourier_parameters": {
            "x": "b/(q^6*a)",
            "y": "b/(q^12*c)",
        },
        "chart_entry_checks": 36,
        "chart_failures": chart_failures,
        "phase_modulus_checks": 3,
        "phase_modulus_failures": modulus_failures,
        "scope": "Exact Fourier membership and phase-modulus checks only; MUB identities are recorded in the exact remainder report.",
    }
    print(json.dumps(result, indent=2))
    if result["status"] != "PASS":
        raise SystemExit(1)
    output = ROOT / "results/phase2_pi3_fourier_exact_direct.json"
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {output}")


if __name__ == "__main__":
    main()
