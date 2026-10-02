#!/usr/bin/env python3
"""Exact symbolic check of D(z) unitarity and Hadamard 2x2 counts.

Only SymPy is required. The computation takes place in Q(i)[z,z^{-1}].
"""

from __future__ import annotations

from itertools import combinations

import sympy as sp


z = sp.symbols("z")
I = sp.I

D = sp.Matrix(
    [
        [1, 1, 1, 1, 1, 1],
        [1, -1, z, -z, I, -I],
        [1, -I, I, I, -I, -1],
        [1, I, -z, z, -1, -I],
        [1, 1 / z, -I, -1, -1 / z, I],
        [1, -1 / z, -1, -I, 1 / z, I],
    ]
)


def conjugate_on_unit_circle(value: sp.Expr) -> sp.Expr:
    """Conjugate in Q(i)[z,z^-1] under z-bar = z^-1."""
    return sp.expand(sp.conjugate(value).xreplace({sp.conjugate(z): 1 / z}))


def is_zero_laurent(value: sp.Expr) -> bool:
    numerator = sp.cancel(value).as_numer_denom()[0]
    return sp.Poly(sp.expand(numerator), z, domain=sp.QQ_I).is_zero


def is_zero_at_one(value: sp.Expr) -> bool:
    return sp.simplify(value.subs(z, 1)) == 0


def main() -> None:
    gram = D * D.applyfunc(conjugate_on_unit_circle).T
    assert gram == 6 * sp.eye(6), "D(z)D(z)^* is not identically 6I"

    conditions = []
    for rows in combinations(range(6), 2):
        r, s = rows
        for columns in combinations(range(6), 2):
            c, d = columns
            orthogonality = sp.expand(
                D[r, c] * conjugate_on_unit_circle(D[s, c])
                + D[r, d] * conjugate_on_unit_circle(D[s, d])
            )
            conditions.append(orthogonality)

    identically_hadamard = sum(is_zero_laurent(condition) for condition in conditions)
    at_one = sum(is_zero_at_one(condition) for condition in conditions)
    assert len(conditions) == 225
    assert identically_hadamard == 51
    assert at_one == 75

    print("PASS: D(z)D(z)^* = 6I in Q(i)[z,z^-1]")
    print(f"PASS: {identically_hadamard}/225 submatrices are identically Hadamard")
    print(f"PASS: {at_one}/225 submatrices are Hadamard at z=1")
    print("For all |z|=1, at least the 51 identically vanishing conditions persist;")
    print("all other conditions can add only finitely many exceptional parameter values.")


if __name__ == "__main__":
    main()
