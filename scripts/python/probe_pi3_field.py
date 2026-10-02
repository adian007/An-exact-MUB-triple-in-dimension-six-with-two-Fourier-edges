#!/usr/bin/env python3
"""Probe small real number fields for the stored pi/3 B3 representative.

This is a recognition diagnostic only. It never upgrades floating data to an
exact algebraic representation.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import sympy as sp


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "results/campaigns/k3_fourier_structure/diagnostic_2026-09-28.json"


def main() -> None:
    data = json.loads(SOURCE.read_text(encoding="utf-8"))
    run = next(item for item in data["runs"] if abs(item["lambda"] - np.pi / 3) < 1e-12)
    matrix = np.asarray(run["cliques"][0]["transitions"]["B3"]["matrix"], dtype=float)
    constants = {
        "Q(sqrt2,sqrt3)": [sp.sqrt(2), sp.sqrt(3), sp.sqrt(6)],
    }
    print("pi/3 B3 field-recognition probe; diagnostic only")
    for name, basis in constants.items():
        errors = []
        nontrivial = 0
        for component in matrix.reshape(-1, 2).flat:
            candidate = sp.nsimplify(float(component), basis, full=True, tolerance=1e-12)
            error = abs(float(candidate) - float(component))
            errors.append(error)
            if error > 1e-12:
                nontrivial += 1
        print(f"{name}: max_error={max(errors):.3e}, mean_error={np.mean(errors):.3e}, unresolved={nontrivial}/72")
    print("Conclusion: small-field recognition is not an exact certificate.")


if __name__ == "__main__":
    main()
