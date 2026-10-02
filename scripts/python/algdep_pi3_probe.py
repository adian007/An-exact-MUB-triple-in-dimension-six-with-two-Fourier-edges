#!/usr/bin/env python3
"""Bounded integer-relation probe for the 100-digit pi/3 B3 refinement."""

from __future__ import annotations

import json
from pathlib import Path

import mpmath as mp


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "results/phase2_pi3_B3_mp100.json"


def relation(value: mp.mpf, degree: int):
    powers = [mp.mpf(1)]
    for _ in range(degree):
        powers.append(powers[-1] * value)
    return mp.pslq(mp.matrix(powers), tol=mp.mpf("1e-85"), maxcoeff=10**10, maxsteps=10000)


def main() -> None:
    mp.mp.dps = 100
    data = json.loads(SOURCE.read_text(encoding="utf-8"))
    found = []
    unresolved = 0
    for vector_index, vector in enumerate(data["vectors"]):
        for component_index, pair in enumerate(vector):
            for part, text in (("re", pair[0]), ("im", pair[1])):
                value = mp.mpf(text)
                if abs(value) < mp.mpf("1e-95"):
                    found.append((vector_index, component_index, part, {"degree": 1, "coefficients": [0, 1]}))
                    continue
                answer = None
                for degree in range(1, 13):
                    candidate = relation(value, degree)
                    if candidate is not None:
                        answer = {"degree": degree, "coefficients": candidate}
                        break
                if answer is None:
                    unresolved += 1
                else:
                    found.append((vector_index, component_index, part, answer))
    print(f"recognized={len(found)} unresolved={unresolved} of 72 real components")
    for item in found:
        print(item)


if __name__ == "__main__":
    main()
