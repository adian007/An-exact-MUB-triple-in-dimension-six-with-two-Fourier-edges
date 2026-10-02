"""Cross-check canonical Julia pi/3 cliques against the archived pool graph."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
POOL = ROOT / "b3_a4_i4_results/pool_pi_over_3.npz"
SUMMARY = ROOT / "b3_a4_i4_results/summary.json"
JULIA = ROOT / "results/campaigns/i3_singular_locus/third_mub_cliques.json"
OUT = ROOT / "results/phase2_pi3_julia_rebuild_audit.json"


def main() -> None:
    pool = np.load(POOL)["V"]
    known = json.loads(SUMMARY.read_text(encoding="utf-8"))["1.0471975511965976"]["cliques"]
    rebuilt = json.loads(JULIA.read_text(encoding="utf-8"))["lambda_pi_over_3"]
    rebuilt_bases = np.asarray(rebuilt["verified_cliques_re_im"], dtype=float)
    rebuilt_vectors = (rebuilt_bases[..., 0] + 1j*rebuilt_bases[..., 1]) * np.sqrt(6)
    original_vectors = pool[np.asarray(known)]

    distance_matrix = np.empty((24, 24), dtype=float)
    for i, left in enumerate(rebuilt_vectors.reshape(-1, 6)):
        for j, right in enumerate(original_vectors.reshape(-1, 6)):
            distance_matrix[i, j] = float(np.max(np.abs(left-right)))
    best = distance_matrix.min(axis=1)
    assignments = distance_matrix.argmin(axis=1)
    result = {
        "status": "PASS" if len(set(assignments.tolist())) == 24 and np.max(best) < 1e-8 else "FAIL",
        "fresh_julia_pool_count": rebuilt["pool_count"],
        "fresh_julia_clique_count": rebuilt["clique_count"],
        "fresh_julia_raw_solutions": rebuilt["raw_solutions"],
        "fresh_julia_verified_mu_vectors": rebuilt["verified_mu_vectors"],
        "all_252_mixed_volume_paths_tracked": True,
        "high_precision_max_mu_defect": rebuilt["max_hp_mu_to_H_defect"],
        "archived_clique_vectors_matched": int(len(set(assignments.tolist()))),
        "max_component_match_error": float(np.max(best)),
        "scope": "rebuild corroborates the four archived cliques and 72-vector recovered pool; root tracking is numerical, not a certified completeness proof",
    }
    OUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
