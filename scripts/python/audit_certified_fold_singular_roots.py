"""Match critical-parameter physical roots to interval-certified roots."""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "results/campaigns/i3_singular_locus"
SOURCE = BASE / "augmented_singular_roots.json"
CERT = BASE / "root_certification.json"
OUT = ROOT / "results/phase2_fold_certified_singular_roots.json"
REFINED = ROOT / "provenance/fold_a4_extracted/mub_fold_a4_computation/results/refined_singular_roots.json"
ORBITS = ROOT / "provenance/fold_a4_extracted/mub_fold_a4_computation/results/a4_refined_orbits.json"
LAMBDA_STAR = 0.1114802243779665542913031975274172717685818097497045


def complex_list(re_im) -> np.ndarray:
    return np.asarray([complex(*pair) for pair in re_im], dtype=np.complex128)


def main() -> None:
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    certified = json.loads(CERT.read_text(encoding="utf-8"))
    roots = certified["roots"]
    cert_vectors = [np.asarray([complex(*pair) for pair in item["approximation"]])
                   for item in roots]
    physical = source["physical_roots"]
    refined_raw = json.loads(REFINED.read_text(encoding="utf-8"))
    refined = [item for item in refined_raw if item["F_norm"] < 1e-7 and
               item["J_u_norm"] < 1e-7 and item["norm_err"] < 1e-8]
    orbit_data = json.loads(ORBITS.read_text(encoding="utf-8"))

    near = []
    all_match_dists = []
    used = set()
    used_all = set()
    near_vectors = []
    for candidate in physical:
        z = complex(*candidate["z_re_im"])
        angle = math.atan2(z.imag, z.real)
        delta = abs(math.atan2(math.sin(angle-LAMBDA_STAR),
                               math.cos(angle-LAMBDA_STAR)))
        x = complex_list(candidate["x_re_im"])
        y = complex_list(candidate["y_re_im"])
        u = complex_list(candidate["u_re_im"])
        packed = np.concatenate([x, y, [z, 1/z], u])
        distances = [float(np.max(np.abs(packed - vec))) for vec in cert_vectors]
        index = int(np.argmin(distances))
        dist = distances[index]
        all_match_dists.append(dist)
        used_all.add(index)
        if delta < 1e-8:
            near_vectors.append(np.angle(x))
            near.append({
                "lambda": angle,
                "delta_lambda_from_target": delta,
                "physical_system_residual": candidate["system_residual"],
                "physical_defect": candidate["physical_defect"],
                "certified_root_index": roots[index]["root_index"],
                "certified": bool(roots[index]["certified"]),
                "max_coordinate_match_error": dist,
                "unique_certificate_match": index not in used,
            })
            used.add(index)

    refined_phases = [np.asarray(item["a"], dtype=float) for item in refined]
    refined_matches = []
    refined_to_near = []
    matched_refined = set()
    for phase in refined_phases:
        if not near_vectors:
            break
        distances = [float(np.max(np.abs(np.angle(np.exp(1j * (phase - target))))))
                     for target in near_vectors]
        index = int(np.argmin(distances))
        refined_matches.append(distances[index])
        refined_to_near.append(index)
        matched_refined.add(index)

    orbit_hits = []
    for orbit in orbit_data["orbits"]:
        hits = [refined_to_near[index] for index in orbit["indices"]]
        orbit_hits.append({
            "recorded_orbit_size": orbit["orbit_size"],
            "certified_root_hits": len(set(hits)),
            "max_match_error_in_source_orbit": orbit["max_match_error"],
        })

    result = {
        "status": "PASS" if len(near) == 24 and all(x["certified"] and
                  x["unique_certificate_match"] and
                  x["max_coordinate_match_error"] < 1e-7 for x in near) else "FAIL",
        "source_physical_root_count": len(physical),
        "certified_augmented_root_count": certified["certified_root_count"],
        "uncertified_augmented_root_count": certified["uncertified_root_count"],
        "near_lambda_star_count_within_1e-8": len(near),
        "lambda_star": LAMBDA_STAR,
        "maximum_delta_lambda_for_near_roots": max(
            (item["delta_lambda_from_target"] for item in near), default=None),
        "independent_refined_singular_candidates": len(refined),
        "refined_candidate_to_certified_root_max_phase_error": max(refined_matches, default=None),
        "refined_candidate_matches_within_1e-8": sum(x < 1e-8 for x in refined_matches),
        "distinct_certified_roots_hit_by_refined_candidates": len(matched_refined),
        "recorded_a4_orbits_mapped_to_certified_roots": orbit_hits,
        "all_physical_roots_uniquely_matched": len(used_all) == len(physical),
        "all_physical_roots_max_match_error": max(all_match_dists, default=None),
        "near_roots": near,
        "scope": "isolated roots of the augmented system near lambda_star; not a proof of exact equality with the rounded target parameter or of global pool completeness",
    }
    OUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in result.items() if k != "near_roots"}, indent=2))


if __name__ == "__main__":
    main()
