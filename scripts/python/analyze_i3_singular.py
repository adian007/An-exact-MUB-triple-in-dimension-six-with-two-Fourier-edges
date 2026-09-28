"""Analyze A4 orbits of isolated I3 roots and recovered third-MUB cliques."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from dita_i3 import (
    _same_projective_vector,
    a4_clique_orbits,
    a4_vector_orbits,
)

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CAMPAIGN = ROOT / "results" / "campaigns" / "i3_singular_locus"
LAMBDA_STAR = 0.1114802243779665542913031975274172717685818097497045


def complex_pair(pair: list[float]) -> complex:
    return complex(pair[0], pair[1])


def physical_singular_records(data: dict[str, object]) -> list[dict[str, object]]:
    records = data.get("all_finite_roots", data.get("physical_roots", []))
    physical = []
    for record in records:
        zvalue = complex_pair(record["z_re_im"])
        xvalues = [complex_pair(pair) for pair in record["x_re_im"]]
        yvalues = [complex_pair(pair) for pair in record["y_re_im"]]
        uvalues = [complex_pair(pair) for pair in record["u_re_im"]]
        defect = max(
            abs(abs(zvalue) - 1),
            max(abs(abs(value) - 1) for value in xvalues),
            max(abs(yvalues[index] - xvalues[index].conjugate()) for index in range(5)),
            abs(sum(value * value for value in uvalues) - 1),
            float(record["system_residual"]),
        )
        if "t_re_im" in record:
            defect = max(
                defect,
                abs(complex_pair(record["t_re_im"]) - zvalue.conjugate()),
            )
        if defect >= 1e-6:
            continue
        vector = np.asarray([1.0 + 0j, *xvalues])
        duplicate = next(
            (
                previous
                for previous in physical
                if abs(complex_pair(previous["z_re_im"]) - zvalue) < 1e-6
                and _same_projective_vector(
                    vector,
                    np.asarray(
                        [
                            1.0 + 0j,
                            *(complex_pair(pair) for pair in previous["x_re_im"]),
                        ]
                    ),
                    1e-7,
                )
            ),
            None,
        )
        if duplicate is None:
            physical.append(record)
        elif record.get("certified") is True:
            duplicate["certified"] = True
    return physical


def analyze_singular_roots(
    path: Path,
    lambda_tolerance: float,
    parameter_cluster_tolerance: float,
) -> dict[str, object]:
    data = json.loads(path.read_text(encoding="utf-8"))
    certification_path = path.parent / "root_certification.json"
    if certification_path.is_file():
        certification = json.loads(certification_path.read_text(encoding="utf-8"))
        certified_by_index = {
            int(root["root_index"]): bool(root["certified"])
            for root in certification["roots"]
        }
        for index, record in enumerate(data.get("all_finite_roots", []), start=1):
            record["certified"] = certified_by_index.get(index, False)
    all_physical = physical_singular_records(data)
    data["physical_candidates_raw"] = data.get("physical_candidates")
    data["physical_candidates"] = len(all_physical)
    data["physical_roots"] = all_physical
    path.write_text(json.dumps(data, separators=(",", ":")) + "\n", encoding="utf-8")
    target_z = np.exp(1j * LAMBDA_STAR)
    records = all_physical
    selected = [
        record
        for record in records
        if abs(complex_pair(record["z_re_im"]) - target_z) < lambda_tolerance
    ]
    clusters: list[list[dict[str, object]]] = []
    for record in sorted(
        selected, key=lambda row: complex_pair(row["z_re_im"]).real
    ):
        zvalue = complex_pair(record["z_re_im"])
        cluster = next(
            (
                group
                for group in clusters
                if abs(complex_pair(group[0]["z_re_im"]) - zvalue)
                < parameter_cluster_tolerance
            ),
            None,
        )
        if cluster is None:
            clusters.append([record])
        else:
            cluster.append(record)

    cluster_results = []
    for cluster in clusters:
        parameter = sum(
            complex_pair(record["z_re_im"]) for record in cluster
        ) / len(cluster)
        vectors = [
            [1.0 + 0j, *(complex_pair(pair) for pair in record["x_re_im"])]
            for record in cluster
        ]
        orbits = a4_vector_orbits(vectors, parameter)
        sizes = [orbit["orbit_size"] for orbit in orbits]
        cluster_results.append(
            {
                "parameter_z_re_im": [parameter.real, parameter.imag],
                "lambda": float(np.angle(parameter)),
                "root_count": len(vectors),
                "certified_root_count": sum(
                    record.get("certified") is True for record in cluster
                ),
                "certification_available": all(
                    "certified" in record for record in cluster
                ),
                "all_roots_certified": all(
                    record.get("certified") is True for record in cluster
                ),
                "max_system_residual": max(
                    float(record["system_residual"]) for record in cluster
                ),
                "orbit_sizes": sizes,
                "stabilizer_orders": [
                    orbit["stabilizer_order"] for orbit in orbits
                ],
                "orbit_decomposition": orbits,
                "matches_two_12_orbit_hypothesis": sizes == [12, 12],
            }
        )
    try:
        source = str(path.resolve().relative_to(ROOT.resolve()))
    except ValueError:
        source = str(path.resolve())
    results: dict[str, object] = {
        "source": source,
        "target_lambda": LAMBDA_STAR,
        "lambda_tolerance_in_z": lambda_tolerance,
        "parameter_cluster_tolerance_in_z": parameter_cluster_tolerance,
        "physical_candidate_count": len(selected),
        "critical_parameter_clusters": cluster_results,
        "augmented_finite_roots": data.get("finite_augmented_solutions"),
        "augmented_mixed_volume_estimate": data.get("mixed_volume"),
        "claim_scope": "numerical roots in the supplied H_D(z) gauge",
    }
    counts_path = path.parent / "transition_counts.json"
    if counts_path.is_file() and cluster_results:
        counts = json.loads(counts_path.read_text(encoding="utf-8"))
        regular = int(counts["critical"]["distinct_physical_vectors"])
        singular = sum(item["root_count"] for item in cluster_results)
        total = regular + singular
        results["transition_reconciliation"] = {
            "below_regular_count": int(counts["below"]["distinct_physical_vectors"]),
            "critical_regular_count": regular,
            "critical_singular_root_count": singular,
            "critical_total_after_augmentation": total,
            "above_regular_count": int(counts["above"]["distinct_physical_vectors"]),
            "expected_counts": [120, 96, 72],
            "observed_counts": [
                int(counts["below"]["distinct_physical_vectors"]),
                total,
                int(counts["above"]["distinct_physical_vectors"]),
            ],
            "matches_expected": (
                int(counts["below"]["distinct_physical_vectors"]) == 120
                and total == 96
                and int(counts["above"]["distinct_physical_vectors"]) == 72
            ),
            "interpretation": (
                "The critical pool solve counts regular roots; the augmented "
                "physical singular roots are added for the critical total."
            ),
        }
    else:
        results["transition_reconciliation"] = None
    return results


def analyze_cliques(path: Path) -> dict[str, object]:
    data = json.loads(path.read_text(encoding="utf-8"))
    output: dict[str, object] = {}
    for label, entry in data.items():
        if label == "_provenance":
            continue
        lambda_value = float(entry["lambda"])
        cliques = [
            [[complex_pair(pair) for pair in vector] for vector in clique]
            for clique in entry["verified_cliques_re_im"]
        ]
        orbits = a4_clique_orbits(cliques, np.exp(1j * lambda_value)) if cliques else []
        output[label] = {
            "lambda": lambda_value,
            "pool_count": entry["pool_count"],
            "verified_clique_count": len(cliques),
            "orbit_sizes": [orbit["orbit_size"] for orbit in orbits],
            "stabilizer_orders": [
                12 // int(orbit["orbit_size"]) for orbit in orbits
            ],
            "orbit_decomposition": orbits,
            "claim_scope": "A4 invariance of verified cliques in the supplied recovered pool",
        }
    return output


def export_roots_for_certification(campaign_dir: Path) -> Path:
    source = campaign_dir / "augmented_singular_roots.json"
    data = json.loads(source.read_text(encoding="utf-8"))
    records = data.get("all_finite_roots", [])
    if not records:
        raise ValueError(f"No finite augmented roots found in {source}.")
    output = campaign_dir / "roots_to_certify.csv"
    header = ["root_index"]
    for prefix, count in (("x", 5), ("y", 5)):
        for index in range(1, count + 1):
            header.extend((f"{prefix}{index}r", f"{prefix}{index}i"))
    header.extend(("zr", "zi", "tr", "ti"))
    for index in range(1, 11):
        header.extend((f"u{index}r", f"u{index}i"))

    with output.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(header)
        for root_index, record in enumerate(records, start=1):
            row: list[float | int] = [root_index]
            for prefix in ("x_re_im", "y_re_im"):
                for pair in record[prefix]:
                    row.extend(pair)
            zvalue = complex_pair(record["z_re_im"])
            tvalue = (
                complex_pair(record["t_re_im"])
                if "t_re_im" in record
                else 1 / zvalue
            )
            row.extend((zvalue.real, zvalue.imag, tvalue.real, tvalue.imag))
            for pair in record["u_re_im"]:
                row.extend(pair)
            writer.writerow(row)
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "mode", choices=("singular", "cliques", "all", "prepare-certification")
    )
    parser.add_argument("--campaign-dir", type=Path, default=DEFAULT_CAMPAIGN)
    parser.add_argument("--lambda-tolerance", type=float, default=1e-5)
    parser.add_argument("--parameter-cluster-tolerance", type=float, default=1e-7)
    args = parser.parse_args()

    results: dict[str, object] = {}
    if args.mode == "prepare-certification":
        output_path = export_roots_for_certification(args.campaign_dir)
        print(f"Wrote {output_path}")
        return
    if args.mode in ("singular", "all"):
        results["singular"] = analyze_singular_roots(
            args.campaign_dir / "augmented_singular_roots.json",
            args.lambda_tolerance,
            args.parameter_cluster_tolerance,
        )
    if args.mode in ("cliques", "all"):
        results["cliques"] = analyze_cliques(
            args.campaign_dir / "third_mub_cliques.json"
        )
    output_path = args.campaign_dir / "a4_orbit_analysis.json"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
    if "singular" in results:
        singular = results["singular"]
        print(f"Physical singular roots near lambda*: {singular['physical_candidate_count']}")
        for cluster in singular["critical_parameter_clusters"]:
            print(
                f"  lambda={cluster['lambda']:.17g}: "
                f"orbits={cluster['orbit_sizes']}, "
                f"stabilizers={cluster['stabilizer_orders']}, "
                f"certified={cluster['certified_root_count']}/"
                f"{cluster['root_count']}"
            )
        reconciliation = singular["transition_reconciliation"]
        if reconciliation is not None:
            print(
                "Transition counts (below, critical+singular, above): "
                f"{reconciliation['observed_counts']} "
                f"match={reconciliation['matches_expected']}"
            )
    if "cliques" in results:
        for label, entry in results["cliques"].items():
            print(
                f"{label}: {entry['verified_clique_count']} verified cliques, "
                f"A4 orbit sizes={entry['orbit_sizes']}, "
                f"stabilizers={entry['stabilizer_orders']}"
            )
    print(f"Wrote {output_path}")


if __name__ == "__main__":
    main()
