#!/usr/bin/env python3
"""Fit the exact a=b=c=u specialization, separately from the recovered clique."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from k3_family_membership import fit_fourier_both_orientations, fit_as_dict


ROOT = Path(__file__).resolve().parents[2]


def main() -> None:
    data = json.loads((ROOT / "results/phase2_pi3_B3_mp100.json").read_text(encoding="utf-8"))
    u = complex(*[float(value) for value in data["vectors"][0][1]])
    q = np.exp(1j * np.pi / 18)
    B = np.array(
        [
            [1, 1, 1, 1, 1, 1],
            [u, -u, u, -u, u, -u],
            [-q * u, q * u, q**7 * u, -q**7 * u, q**31 * u, -q**31 * u],
            [q**26, q**26, q**2, q**2, q**14, q**14],
            [q**20 * u, q**2 * u, q**32 * u, -q**32 * u, q**8 * u, -q**8 * u],
            [q**28, q**28, q**16, q**16, q**4, q**4],
        ],
        dtype=complex,
    )
    z = np.exp(1j * np.pi / 3)
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
    transition = (H.conj().T @ B) / 6
    direct, transpose = fit_fourier_both_orientations(transition)
    result = {"direct": fit_as_dict(direct), "transpose": fit_as_dict(transpose)}
    print(json.dumps(result, indent=2, default=lambda value: [value.real, value.imag]))
    (ROOT / "results/phase2_abc_fourier_fit.json").write_text(json.dumps(result, indent=2, default=lambda value: [value.real, value.imag]) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
