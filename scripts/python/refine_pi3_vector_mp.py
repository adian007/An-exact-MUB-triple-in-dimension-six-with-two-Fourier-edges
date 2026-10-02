#!/usr/bin/env python3
"""High-precision refinement of one stored pi/3 MU vector.

This refines a numerical seed against ten independent equations. It is not
an exact algebraic reconstruction; the output is intended as input for a
later relation-recovery step.
"""

from __future__ import annotations

import json
from pathlib import Path

import mpmath as mp


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "results/campaigns/k3_fourier_structure/diagnostic_2026-09-28.json"


def main() -> None:
    mp.mp.dps = 100
    data = json.loads(SOURCE.read_text(encoding="utf-8"))
    run = next(item for item in data["runs"] if abs(item["lambda"] - mp.pi / 3) < 1e-12)
    raw = run["cliques"][0]["transitions"]["B3"]["matrix"]
    zeta = mp.mpc("0.5", mp.sqrt(3) / 2)
    I = mp.j
    H = [
        [1, 1, 1, 1, 1, 1],
        [1, -1, zeta, -zeta, I, -I],
        [1, -I, I, I, -I, -1],
        [1, I, -zeta, zeta, -1, -I],
        [1, 1 / zeta, -I, -1, -1 / zeta, I],
        [1, -1 / zeta, -1, -I, 1 / zeta, I],
    ]
    def equations(*args):
        vec = [mp.mpc(1)] + [mp.mpc(args[2 * j], args[2 * j + 1]) for j in range(5)]
        out = [mp.re(vec[j] * mp.conj(vec[j])) - 1 for j in range(1, 6)]
        for col_index in range(5):
            inner = sum(mp.conj(vec[row]) * H[row][col_index] for row in range(6))
            out.append(mp.re(inner * mp.conj(inner)) - 6)
        return tuple(out)

    refined = []
    for col in range(6):
        seed = [complex(raw[row][col][0], raw[row][col][1]) for row in range(6)]
        seed = [value / seed[0] for value in seed]
        seed_real = []
        for value in seed[1:]:
            seed_real.extend([mp.mpf(value.real), mp.mpf(value.imag)])
        root = mp.findroot(equations, tuple(seed_real), tol=mp.mpf("1e-80"), maxsteps=100)
        refined.append([mp.mpc(1)] + [mp.mpc(root[2 * j], root[2 * j + 1]) for j in range(5)])
        print(f"vector {col}: max_equation_residual {mp.nstr(max(abs(value) for value in equations(*root)), 8)}")

    pairwise = []
    for left in range(6):
        for right in range(left + 1, 6):
            pairwise.append(abs(sum(mp.conj(refined[left][row]) * refined[right][row] for row in range(6))))
    h_defects = []
    for vector in refined:
        for col in range(6):
            inner = sum(mp.conj(vector[row]) * H[row][col] for row in range(6))
            h_defects.append(abs(abs(inner) ** 2 - 6))

    output = {
        "parameter": "pi/3",
        "precision_bits_approx": 332,
        "source": str(SOURCE.relative_to(ROOT)).replace("\\", "/"),
        "method": "mpmath.findroot; five unit-circle and five MU equations per vector",
        "max_pairwise_orthogonality_residual": mp.nstr(max(pairwise), 30),
        "max_H_MU_residual_all_six_columns": mp.nstr(max(h_defects), 30),
        "vectors": [
            [[mp.nstr(value.real, 110), mp.nstr(value.imag, 110)] for value in vector]
            for vector in refined
        ],
    }
    output_path = ROOT / "results/phase2_pi3_B3_mp100.json"
    output_path.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print("max_pairwise_orthogonality_residual", mp.nstr(max(pairwise), 8))
    print("max_H_MU_residual_all_six_columns", mp.nstr(max(h_defects), 8))
    print("wrote", output_path)


if __name__ == "__main__":
    main()
