"""High-precision audit of the corrected, distinct a/b/c pi/3 ansatz.

This is a numerical audit of the exact export in
``symbolic_export/w1_pi3_exact_abc.m2``.  It deliberately does not test the
deprecated a=b=c shortcut.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import mpmath as mp


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "results" / "phase2_pi3_B3_mp100.json"
OUT = ROOT / "results" / "phase2_pi3_distinct_embedding_audit.json"

P = [
    104329, 0, -193800, 0, 192924, 0, -299406, 0, 331452, 0,
    -381612, 0, 493595, 0, -381612, 0, 331452, 0, -299406, 0,
    192924, 0, -193800, 0, 104329,
]

# Coefficients are indexed by the even powers u^0,u^2,...,u^22.
L0 = [-1811714204362686491473791534, 0,
      801861270977657426086402045, 0,
      -1591397118696353756136381769, 0,
      2139096543575902440554317075, 0,
      -1177917192554426766591014887, 0,
      3958308486166522559203865220, 0,
      -2644550445045953692323054841, 0,
      1531431372534934342773255463, 0,
      -2568263831168312266729485648, 0,
      989839703306762314045034317, 0,
      -1756307929240743070046703434, 0,
      1514595399406991554715184449]
L1 = [-2942616298328776548850611985, 0,
      3079541734716012913628182138, 0,
      -2947627835384158355404461405, 0,
      6180409088066584397938320628, 0,
      -4176525922160551925590943297, 0,
      7416627722299525522595911037, 0,
      -8051612621828886127519440680, 0,
      4579744692264779374914258781, 0,
      -6270773258597842323743219399, 0,
      3167065685867132714755475197, 0,
      -2770857433192841326203024039, 0,
      3352989164513563728488605330]
L2 = [4754330502691463040324403519, 0,
      -3881403005693670339714584183, 0,
      4539024954080512111540843174, 0,
      -8319505631642486838492637703, 0,
      5354443114714978692181958184, 0,
      -11374936208466048081799776257, 0,
      10696163066874839819842495521, 0,
      -6111176064799713717687514244, 0,
      8839037089766154590472705047, 0,
      -4156905389173895028800509514, 0,
      4527165362433584396249727473, 0,
      -4867584563920555283203789779]
QCOEF_STR = "610818468996621564036607682"


def z(pair: list[str]) -> mp.mpc:
    return mp.mpc(mp.mpf(pair[0]), mp.mpf(pair[1]))


def poly_even(coeffs: list[int], u: mp.mpc) -> mp.mpc:
    return sum(mp.mpf(c) * u ** k for k, c in enumerate(coeffs))


def residuals() -> dict:
    mp.mp.dps = 120
    raw = json.loads(SOURCE.read_text(encoding="utf-8"))
    vv = [[z(x) for x in row] for row in raw["vectors"]]
    q = mp.e ** (mp.j * mp.pi / 18)
    qcoef = mp.mpf(QCOEF_STR)
    a, b, c = vv[0][1], vv[2][1], vv[4][1]

    phase_residuals = {
        "P(a)": abs(poly_even(P, a)),
        "P(b)": abs(poly_even(P, b)),
        "P(c)": abs(poly_even(P, c)),
        "L0(q,a)": abs(poly_even(L0, a) + qcoef * q),
        "L1(q,b)": abs(poly_even(L1, b) + qcoef * q),
        "L2(q,c)": abs(poly_even(L2, c) + qcoef * q),
        "q^12-q^6+1": abs(q**12 - q**6 + 1),
    }

    templates = [
        [1, a, -q*a, q**-10, q**-16*a, q**-8],
        [1, -a, q*a, q**-10, q**2*a, q**-8],
        [1, b, q**7*b, q**2, q**-4*b, q**16],
        [1, -b, -q**7*b, q**2, -q**-4*b, q**16],
        [1, c, q**-5*c, q**14, q**8*c, q**4],
        [1, -c, -q**-5*c, q**14, -q**8*c, q**4],
    ]
    template_error = max(abs(templates[i][j] - vv[i][j])
                         for i in range(6) for j in range(6))

    pair_error = max(abs(sum(mp.conj(vv[i][k]) * vv[j][k] for k in range(6)))
                     for i in range(6) for j in range(i + 1, 6))
    norm_error = max(abs(sum(abs(x)**2 for x in row) - 6) for row in vv)
    zeta = q**6
    ii = q**9
    H = [
        [1, 1, 1, 1, 1, 1],
        [1, -1, zeta, -zeta, ii, -ii],
        [1, -ii, ii, ii, -ii, -1],
        [1, ii, -zeta, zeta, -1, -ii],
        [1, q**-6, -ii, -1, -q**-6, ii],
        [1, -q**-6, -1, -ii, q**-6, ii],
    ]
    h_mu_error = mp.mpf(0)
    for row in vv:
        for col in range(6):
            inner = sum(mp.conj(row[r]) * H[r][col] for r in range(6))
            h_mu_error = max(h_mu_error, abs(abs(inner) - mp.sqrt(6)))

    all_residuals = {k: float(v) for k, v in phase_residuals.items()}
    all_residuals.update({
        "template_max_abs": float(template_error),
        "pairwise_max_abs": float(pair_error),
        "unit_norm_max_abs": float(norm_error),
        "H_MU_max_abs": float(h_mu_error),
    })
    return {
        # The degree-22 relations amplify the ~100-digit input rounding to
        # about 1e-72; the structural vector identities remain ~1e-100.
        "status": "PASS" if max(all_residuals.values()) < 1e-60 else "FAIL",
        "source": str(SOURCE.relative_to(ROOT)),
        "export": "symbolic_export/w1_pi3_exact_abc.m2",
        "precision_dps": mp.mp.dps,
        "shortcut_tested": False,
        "max_residual": max(all_residuals.values()),
        "residuals": all_residuals,
    }


if __name__ == "__main__":
    result = residuals()
    OUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
