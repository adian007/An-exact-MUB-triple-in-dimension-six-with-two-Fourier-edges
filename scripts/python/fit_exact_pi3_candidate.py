#!/usr/bin/env python3
"""Numerical family-fit diagnostic for the exact phase-ansatz candidate."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from k3_family_membership import fit_fourier_both_orientations, fit_as_dict


ROOT = Path(__file__).resolve().parents[2]


def main() -> None:
    data = json.loads((ROOT / "results/phase2_pi3_B3_mp100.json").read_text(encoding="utf-8"))
    vectors = np.asarray(
        [[complex(float(pair[0]), float(pair[1])) for pair in vector] for vector in data["vectors"]],
        dtype=complex,
    ).T
    lam = np.pi / 3
    z = np.exp(1j * lam)
    H = np.array(
        [
            [1, 1, 1, 1, 1, 1],
            [1, -1, z, -z, 1j, -1j],
            [1, -1j, 1j, 1j, -1j, -1],
            [1, 1j, -z, z, -1, -1j],
            [1, np.conj(z), -1j, -1, -np.conj(z), 1j],
            [1, -np.conj(z), -1, -1j, np.conj(z), 1j],
        ],
        dtype=complex,
    )
    transition = (H.conj().T / np.sqrt(6)) @ (vectors / np.sqrt(6))
    direct, transpose = fit_fourier_both_orientations(transition)
    result = {"direct": fit_as_dict(direct), "transpose": fit_as_dict(transpose)}
    print(json.dumps(result, indent=2, default=lambda value: [value.real, value.imag]))
    output = ROOT / "results/phase2_exact_pi3_fourier_fit.json"
    output.write_text(json.dumps(result, indent=2, default=lambda value: [value.real, value.imag]) + "\n", encoding="utf-8")
    print(f"wrote {output}")


if __name__ == "__main__":
    main()
