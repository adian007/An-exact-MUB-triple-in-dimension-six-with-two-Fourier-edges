#!/usr/bin/env python3
"""Export the corrected exact pi/3 W1 system with distinct a,b,c roots."""

from __future__ import annotations

from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "symbolic_export" / "w1_pi3_exact_abc.m2"


def qp(q: sp.Symbol, exponent: int) -> sp.Expr:
    return q ** (exponent % 36)


def bar(entry: sp.Expr, q, a, b, c, ai, bi, ci):
    return sp.sympify(entry).xreplace({q: q**35, a: ai, b: bi, c: ci})


def fmt(expr: sp.Expr) -> str:
    return sp.sstr(sp.expand(expr)).replace("**", "^")


def polynomial(u: sp.Symbol) -> sp.Expr:
    return 104329*u**24 - 193800*u**22 + 192924*u**20 - 299406*u**18 + 331452*u**16 - 381612*u**14 + 493595*u**12 - 381612*u**10 + 331452*u**8 - 299406*u**6 + 192924*u**4 - 193800*u**2 + 104329


def main() -> None:
    q, a, b, c, ai, bi, ci = sp.symbols("q a b c ai bi ci")
    z1, z2, z3, z4, z5 = sp.symbols("z1 z2 z3 z4 z5")
    w1, w2, w3, w4, w5 = sp.symbols("w1 w2 w3 w4 w5")
    x = [1, z1, z2, z3, z4, z5]
    y = [1, w1, w2, w3, w4, w5]
    H = [
        [1, 1, 1, 1, 1, 1],
        [1, -1, qp(q, 6), -qp(q, 6), qp(q, 9), -qp(q, 9)],
        [1, -qp(q, 9), qp(q, 9), qp(q, 9), -qp(q, 9), -1],
        [1, qp(q, 9), -qp(q, 6), qp(q, 6), -1, -qp(q, 9)],
        [1, qp(q, 30), -qp(q, 9), -1, -qp(q, 30), qp(q, 9)],
        [1, -qp(q, 30), -1, -qp(q, 9), qp(q, 30), qp(q, 9)],
    ]
    B = [
        [1, 1, 1, 1, 1, 1],
        [a, -a, b, -b, c, -c],
        [-qp(q, 1)*a, qp(q, 1)*a, qp(q, 7)*b, -qp(q, 7)*b, qp(q, 31)*c, -qp(q, 31)*c],
        [qp(q, 26), qp(q, 26), qp(q, 2), qp(q, 2), qp(q, 14), qp(q, 14)],
        [qp(q, 20)*a, qp(q, 2)*a, qp(q, 32)*b, -qp(q, 32)*b, qp(q, 8)*c, -qp(q, 8)*c],
        [qp(q, 28), qp(q, 28), qp(q, 16), qp(q, 16), qp(q, 4), qp(q, 4)],
    ]
    Hbar = [[bar(entry, q, a, b, c, ai, bi, ci) for entry in row] for row in H]
    Bbar = [[bar(entry, q, a, b, c, ai, bi, ci) for entry in row] for row in B]
    equations = [z1*w1-1, z2*w2-1, z3*w3-1, z4*w4-1, z5*w5-1]
    for matrix, conjugate in ((H, Hbar), (B, Bbar)):
        for column in range(6):
            left = sum(y[row] * matrix[row][column] for row in range(6))
            right = sum(x[row] * conjugate[row][column] for row in range(6))
            equations.append(sp.expand(left * right - 6))

    P_a, P_b, P_c = polynomial(a), polynomial(b), polynomial(c)
    L0 = 610818468996621564036607682*q + 1514595399406991554715184449*a**22 - 1756307929240743070046703434*a**20 + 989839703306762314045034317*a**18 - 2568263831168312266729485648*a**16 + 1531431372534934342773255463*a**14 - 2644550445045953692323054841*a**12 + 3958308486166522559203865220*a**10 - 1177917192554426766591014887*a**8 + 2139096543575902440554317075*a**6 - 1591397118696353756136381769*a**4 + 801861270977657426086402045*a**2 - 1811714204362686491473791534
    L1 = 610818468996621564036607682*q + 3352989164513563728488605330*b**22 - 2770857433192841326203024039*b**20 + 3167065685867132714755475197*b**18 - 6270773258597842323743219399*b**16 + 4579744692264779374914258781*b**14 - 8051612621828886127519440680*b**12 + 7416627722299525522595911037*b**10 - 4176525922160551925590943297*b**8 + 6180409088066584397938320628*b**6 - 2947627835384158355404461405*b**4 + 3079541734716012913628182138*b**2 - 2942616298328776548850611985
    L2 = 610818468996621564036607682*q - 4867584563920555283203789779*c**22 + 4527165362433584396249727473*c**20 - 4156905389173895028800509514*c**18 + 8839037089766154590472705047*c**16 - 6111176064799713717687514244*c**14 + 10696163066874839819842495521*c**12 - 11374936208466048081799776257*c**10 + 5354443114714978692181958184*c**8 - 8319505631642486838492637703*c**6 + 4539024954080512111540843174*c**4 - 3881403005693670339714584183*c**2 + 4754330502691463040324403519

    lines = [
        "-- Exact pi/3 W1 export with distinct algebraic phases a,b,c.",
        "S = QQ[q,a,b,c,ai,bi,ci];",
        f"Pa = {fmt(P_a)};",
        f"Pb = {fmt(P_b)};",
        f"Pc = {fmt(P_c)};",
        f"L0 = {fmt(L0)};",
        f"L1 = {fmt(L1)};",
        f"L2 = {fmt(L2)};",
        "kk = toField(S / ideal(q^12-q^6+1, a*ai-1, b*bi-1, c*ci-1, Pa, Pb, Pc, L0, L1, L2));",
        "R = kk[z1,z2,z3,z4,z5,w1,w2,w3,w4,w5];",
        "witness = ideal(\n  " + ",\n  ".join(fmt(e) for e in equations) + "\n);",
        "print(\"numgens witness = \" | toString numgens witness);",
        "elapsedTime G = gens gb witness;",
        "print(\"Groebner basis = \" | toString G);",
        "print(\"1 % witness = \" | toString (1_R % witness));",
        "print(\"dim = \" | toString dim witness);",
    ]
    OUTPUT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {OUTPUT} with {len(equations)} equations")


if __name__ == "__main__":
    main()
