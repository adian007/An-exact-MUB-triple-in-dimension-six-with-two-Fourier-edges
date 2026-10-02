"""Exact Groebner consistency check for the three corrected phase embeddings."""

from __future__ import annotations

import json
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "phase2_pi3_phase_field_ideal.json"


def main() -> None:
    q, a, b, c = sp.symbols("q a b c")
    phi = sp.cyclotomic_poly(36, q)
    P = lambda u: (104329*u**24 - 193800*u**22 + 192924*u**20
                   - 299406*u**18 + 331452*u**16 - 381612*u**14
                   + 493595*u**12 - 381612*u**10 + 331452*u**8
                   - 299406*u**6 + 192924*u**4 - 193800*u**2
                   + 104329)
    L0 = (610818468996621564036607682*q
          + 1514595399406991554715184449*a**22
          - 1756307929240743070046703434*a**20
          + 989839703306762314045034317*a**18
          - 2568263831168312266729485648*a**16
          + 1531431372534934342773255463*a**14
          - 2644550445045953692323054841*a**12
          + 3958308486166522559203865220*a**10
          - 1177917192554426766591014887*a**8
          + 2139096543575902440554317075*a**6
          - 1591397118696353756136381769*a**4
          + 801861270977657426086402045*a**2
          - 1811714204362686491473791534)
    L1 = (610818468996621564036607682*q
          + 3352989164513563728488605330*b**22
          - 2770857433192841326203024039*b**20
          + 3167065685867132714755475197*b**18
          - 6270773258597842323743219399*b**16
          + 4579744692264779374914258781*b**14
          - 8051612621828886127519440680*b**12
          + 7416627722299525522595911037*b**10
          - 4176525922160551925590943297*b**8
          + 6180409088066584397938320628*b**6
          - 2947627835384158355404461405*b**4
          + 3079541734716012913628182138*b**2
          - 2942616298328776548850611985)
    L2 = (610818468996621564036607682*q
          - 4867584563920555283203789779*c**22
          + 4527165362433584396249727473*c**20
          - 4156905389173895028800509514*c**18
          + 8839037089766154590472705047*c**16
          - 6111176064799713717687514244*c**14
          + 10696163066874839819842495521*c**12
          - 11374936208466048081799776257*c**10
          + 5354443114714978692181958184*c**8
          - 8319505631642486838492637703*c**6
          + 4539024954080512111540843174*c**4
          - 3881403005693670339714584183*c**2
          + 4754330502691463040324403519)

    # Eliminate the independent phase variable from each P(u), L_i(q,u).
    # Their resultant is the exact set of possible q-values. The three phase
    # variables are independent, so a common cyclotomic q root is precisely
    # the consistency condition over the algebraic closure.
    resultants = [sp.resultant(P(u), relation, u)
                  for u, relation in ((a, L0), (b, L1), (c, L2))]
    common = sp.Poly(phi, q)
    for resultant in resultants:
        common = common.gcd(sp.Poly(resultant, q))
    is_unit = sp.degree(common, q) == 0
    q_root = sp.exp(sp.I * sp.pi / 18)
    field = sp.QQ.algebraic_field(q_root)
    phase_factors = {}
    for label, variable, relation in zip(
            ("a", "b", "c"), (a, b, c), (L0, L1, L2)):
        phase_poly = sp.Poly(P(variable), variable, domain=field)
        relation_poly = sp.Poly(relation.subs(q, q_root), variable, domain=field)
        factor = phase_poly.gcd(relation_poly).monic()
        phase_factors[label] = {
            "degree_over_Q_q": int(factor.degree()),
            "factor": sp.sstr(factor.as_expr()),
        }
    result = {
        "status": "INCONSISTENT" if is_unit else "CONSISTENT",
        "unit_ideal": bool(is_unit),
        "common_q_factor_degree": int(common.degree()),
        "resultant_degrees": [int(sp.degree(r, q)) for r in resultants],
        "intended_embedding_q": "exp(i*pi/18)",
        "phase_factors_over_Q_q": phase_factors,
        "elimination": "exact univariate resultants and gcd with cyclotomic polynomial",
        "scope": "phase embedding ideal only; excludes the W1 witness variables",
    }
    OUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
