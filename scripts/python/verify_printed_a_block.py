from __future__ import annotations

import sympy as sp

# Real coordinates for sin/cos(theta) and sin/cos(phi).
s_t, c_t, s_p, c_p = sp.symbols("s_t c_t s_p c_p", real=True)
I = sp.I
sqrt3 = sp.sqrt(3)

A11 = sp.Rational(1, 2) + I * sqrt3 / 2 * (c_t + (c_p - I * s_p) * s_t)
A12 = -sp.Rational(1, 2) + I * sqrt3 / 2 * (-c_t + (c_p + I * s_p) * s_t)
A_printed = sp.Matrix([[A11, A12], [A12, -A11]])
residual = sp.expand(A_printed * A_printed.conjugate().T - 2 * sp.eye(2))

relations = [s_t**2 + c_t**2 - 1, s_p**2 + c_p**2 - 1]
basis = sp.groebner(
    relations,
    s_t, c_t, s_p, c_p,
    domain=sp.QQ.algebraic_field(sqrt3, I),
)
expected = sp.Matrix([
    [sqrt3 * s_t * s_p, -I * s_t * c_p * (sqrt3 + 3 * s_t * s_p)],
    [ I * s_t * c_p * (sqrt3 + 3 * s_t * s_p), sqrt3 * s_t * s_p],
])
difference = (residual - expected).applyfunc(lambda entry: sp.expand(basis.reduce(entry)[1]))
assert all(entry == 0 for entry in difference)
print("PASS: full residual agrees with the displayed Delta/Gamma matrix modulo the two trigonometric circle relations.")
print(expected)
