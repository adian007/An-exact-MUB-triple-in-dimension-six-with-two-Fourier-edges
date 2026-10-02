#!/usr/bin/env python3
"""Audit whether augmented fold singular candidates occur in the Julia pool."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
SINGULAR = ROOT / "provenance/fold_a4_extracted/mub_fold_a4_computation/results/refined_singular_roots.json"
SWEEP = ROOT / "results/campaigns/i3_singular_locus/fold_sweep_postnorm_20260929_reconciled.json"
ORBITS = ROOT / "provenance/fold_a4_extracted/mub_fold_a4_computation/results/a4_refined_orbits.json"


def circular_distance(left: np.ndarray, right: np.ndarray) -> float:
    return float(np.max(np.abs(np.angle(np.exp(1j * (left - right))))))


def main() -> None:
    singular = json.loads(SINGULAR.read_text(encoding="utf-8"))
    good = [item for item in singular if item["F_norm"] < 1e-7 and item["J_u_norm"] < 1e-7 and item["norm_err"] < 1e-8]
    sweep = json.loads(SWEEP.read_text(encoding="utf-8"))
    critical = next(item for item in sweep["records"] if item["label"] == "fold_critical")
    pool = np.asarray(critical["pool_vectors_re_im"], dtype=float)
    pool_complex = pool[:, :, 0] + 1j * pool[:, :, 1]
    pool_phases = np.angle(pool_complex[:, 1:] / pool_complex[:, :1])
    matches = []
    for candidate in good:
        phase = np.asarray(candidate["a"], dtype=float)
        distances = np.asarray([circular_distance(phase, row) for row in pool_phases])
        matches.append(float(np.min(distances)))
    orbit_data = json.loads(ORBITS.read_text(encoding="utf-8"))
    orbit_sizes = [len(orbit["indices"]) for orbit in orbit_data["orbits"]]
    result = {
        "singular_records": len(singular),
        "filtered_singular_candidates": len(good),
        "critical_julia_pool": int(len(pool_complex)),
        "candidate_to_pool_min_phase_distance_max": max(matches),
        "candidate_to_pool_matches_at_1e-7": sum(distance < 1e-7 for distance in matches),
        "candidate_to_pool_matches_at_1e-5": sum(distance < 1e-5 for distance in matches),
        "reported_a4_orbit_sizes": orbit_sizes,
        "reported_a4_max_match_error": max(max(item["matches"][i][1] for i in range(len(item["matches"]))) for item in orbit_data["orbits"]),
    }
    print(json.dumps(result, indent=2))
    output = ROOT / "results/phase2_fold_critical_match_audit.json"
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {output}")


if __name__ == "__main__":
    main()
