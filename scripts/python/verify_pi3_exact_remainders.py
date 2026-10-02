"""Exact remainder verification for the three distinct pi/3 phase embeddings.

This is a compact SymPy proof check for the corrected export.  It verifies
that every cleared Hadamard-unbiasedness equation lies in the ideal generated
by the cyclotomic relation, the common degree-24 phase polynomial, and the
phase-specific q-linear relation.  Pairwise B3 orthogonality is checked
separately by the existing symbolic derivation.
"""

from __future__ import annotations

import json
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "phase2_pi3_exact_remainder_verification.json"


def conjugate(expr: sp.Expr, q: sp.Symbol, u: sp.Symbol) -> sp.Expr:
    return sp.expand(sp.sympify(expr).xreplace({q: 1 / q, u: 1 / u}))


def hadamard(q: sp.Symbol) -> sp.Matrix:
    z = q**6
    ii = q**9
    return sp.Matrix([
        [1, 1, 1, 1, 1, 1],
        [1, -1, z, -z, ii, -ii],
        [1, -ii, ii, ii, -ii, -1],
        [1, ii, -z, z, -1, -ii],
        [1, q**-6, -ii, -1, -q**-6, ii],
        [1, -q**-6, -1, -ii, q**-6, ii],
    ])


def phase_ideal(q: sp.Symbol, u: sp.Symbol, relation: sp.Expr) -> sp.GroebnerBasis:
    phi = sp.cyclotomic_poly(36, q)
    p = (104329*u**24 - 193800*u**22 + 192924*u**20 - 299406*u**18
         + 331452*u**16 - 381612*u**14 + 493595*u**12
         - 381612*u**10 + 331452*u**8 - 299406*u**6
         + 192924*u**4 - 193800*u**2 + 104329)
    return sp.groebner([phi, p, relation], q, u, order="grevlex", method="f5b")


def main() -> None:
    q, u = sp.symbols("q u")
    H = hadamard(q)
    ansatz = {
        "pair_0_plus": [1, u, -q*u, q**-10, q**-16*u, q**-8],
        "pair_0_minus": [1, -u, q*u, q**-10, q**2*u, q**-8],
        "pair_1_plus": [1, u, q**7*u, q**2, q**-4*u, q**16],
        "pair_1_minus": [1, -u, -q**7*u, q**2, -q**-4*u, q**16],
        "pair_2_plus": [1, u, q**-5*u, q**14, q**8*u, q**4],
        "pair_2_minus": [1, -u, -q**-5*u, q**14, -q**8*u, q**4],
    }
    relations = {
        "pair_0": (610818468996621564036607682*q
                   + 1514595399406991554715184449*u**22
                   - 1756307929240743070046703434*u**20
                   + 989839703306762314045034317*u**18
                   - 2568263831168312266729485648*u**16
                   + 1531431372534934342773255463*u**14
                   - 2644550445045953692323054841*u**12
                   + 3958308486166522559203865220*u**10
                   - 1177917192554426766591014887*u**8
                   + 2139096543575902440554317075*u**6
                   - 1591397118696353756136381769*u**4
                   + 801861270977657426086402045*u**2
                   - 1811714204362686491473791534),
        "pair_1": (610818468996621564036607682*q
                   + 3352989164513563728488605330*u**22
                   - 2770857433192841326203024039*u**20
                   + 3167065685867132714755475197*u**18
                   - 6270773258597842323743219399*u**16
                   + 4579744692264779374914258781*u**14
                   - 8051612621828886127519440680*u**12
                   + 7416627722299525522595911037*u**10
                   - 4176525922160551925590943297*u**8
                   + 6180409088066584397938320628*u**6
                   - 2947627835384158355404461405*u**4
                   + 3079541734716012913628182138*u**2
                   - 2942616298328776548850611985),
        "pair_2": (610818468996621564036607682*q
                   - 4867584563920555283203789779*u**22
                   + 4527165362433584396249727473*u**20
                   - 4156905389173895028800509514*u**18
                   + 8839037089766154590472705047*u**16
                   - 6111176064799713717687514244*u**14
                   + 10696163066874839819842495521*u**12
                   - 11374936208466048081799776257*u**10
                   + 5354443114714978692181958184*u**8
                   - 8319505631642486838492637703*u**6
                   + 4539024954080512111540843174*u**4
                   - 3881403005693670339714584183*u**2
                   + 4754330502691463040324403519),
    }

    phase_bases = {
        label: phase_ideal(q, u, relation)
        for label, relation in relations.items()
    }
    result = {"status": "PASS", "checks": {}}
    for label, vector in ansatz.items():
        phase_label = label.rsplit("_", 1)[0]
        G = phase_bases[phase_label]
        remainders = []
        for col in range(6):
            inner = sum(conjugate(vector[row], q, u) * H[row, col]
                        for row in range(6))
            expr = sp.together(inner * conjugate(inner, q, u) - 6)
            numerator = sp.Poly(expr.as_numer_denom()[0], q, u).as_expr()
            remainders.append(sp.expand(G.reduce(numerator)[1]) == 0)
        result["checks"][label] = {
            "all_six_hadamard_equations_zero": all(remainders),
            "zero_equations": sum(remainders),
        }
        if not all(remainders):
            result["status"] = "FAIL"

    OUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
