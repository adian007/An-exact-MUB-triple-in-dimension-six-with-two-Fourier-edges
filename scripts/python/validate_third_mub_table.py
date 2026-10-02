#!/usr/bin/env python3
"""Steps 1-3 table + clustering from special_loci_search CSV (no Julia deps)."""
import csv
import math
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results"
CSV_PATH = RESULTS / "special_loci_search.csv"
DITA = math.acos(1 / math.sqrt(3))


def load_rows(path):
    with open(path, encoding="utf-8") as f:
        r = csv.DictReader(f)
        return list(r)


def third_candidates(rows, exclude_deep_ref=True):
    out = []
    for row in rows:
        if row.get("status") != "ok":
            continue
        try:
            mc = int(row.get("max_clique") or 0)
        except ValueError:
            continue
        if mc < 6:
            continue
        if exclude_deep_ref and "ref_ref" in row.get("name", ""):
            continue
        out.append(row)
    return out


def locus_info(theta, phi, lam):
    d_theta0 = abs(float(theta))
    d_dita = math.hypot(float(theta) - DITA, float(phi) - math.pi / 4)
    dists = {"Fourier_theta0": d_theta0, "Dita": d_dita, "Mobius_z4_dev": float("nan")}
    nearest = min(dists, key=dists.get)
    return d_theta0, nearest, dists[nearest]


def exact_match(theta, phi, lam, atol=1e-6):
    rels = []
    labels = [(0, "0"), (math.pi / 6, "pi/6"), (math.pi / 4, "pi/4"), (DITA, "arccos(1/sqrt(3))")]
    for val, lab in labels:
        if abs(float(theta) - val) < atol:
            rels.append(f"theta={lab}")
        if abs(float(phi) - val) < atol:
            rels.append(f"phi={lab}")
        if abs(float(lam) - val) < atol:
            rels.append(f"lambda={lab}")
    return "; ".join(rels)


def cluster_rank(candidates):
    if len(candidates) < 2:
        return 0, []
    xs = [[float(c["theta"]), float(c["phi"]), float(c["lambda"])] for c in candidates]
    n = len(xs)
    # center
    mean = [sum(x[i] for x in xs) / n for i in range(3)]
    xc = [[x[i] - mean[i] for i in range(3)] for x in xs]
    # SVD via covariance
    import numpy as np
    M = np.array(xc)
    _, s, _ = np.linalg.svd(M, full_matrices=False)
    rank = int(np.sum(s > 1e-6 * s[0]))
    return rank, s.tolist()


def main():
    rows = load_rows(CSV_PATH)
    cand = third_candidates(rows)
    print(f"CSV rows: {len(rows)}")
    print(f"Third-MUB candidates (exclude ref_ref): {len(cand)}\n")
    print(f"{'idx':>3}  {'name':<28} {'theta':>10} {'phi':>10} {'lambda':>10}  {'near':<16} {'d':>10}  exact")
    audit_path = RESULTS / "third_mub_candidates_audit_py.csv"
    with open(audit_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["index", "name", "theta", "phi", "lambda", "d_theta0", "nearest_locus",
                    "nearest_dist", "exact_confirmed", "found_fourth"])
        for i, c in enumerate(cand, 1):
            th, ph, lm = c["theta"], c["phi"], c["lambda"]
            d0, near, nd = locus_info(th, ph, lm)
            ex = exact_match(th, ph, lm)
            print(f"{i:3d}  {c['name']:<28} {float(th):10.6f} {float(ph):10.6f} {float(lm):10.6f}  "
                  f"{near:<16} {nd:10.2e}  {ex or '—'}")
            w.writerow([i, c["name"], th, ph, lm, d0, near, nd, ex, c.get("found_fourth")])
    rank, sv = cluster_rank(cand)
    print(f"\nCluster SVD rank: {rank}  singular values: {[round(x, 6) for x in sv]}")
    summary = RESULTS / "validation_summary_py.txt"
    summary.write_text(
        f"n_csv_rows={len(rows)}\nn_candidates={len(cand)}\ncluster_rank={rank}\n"
        f"cluster_sv={','.join(str(x) for x in sv)}\nfound_fourth_any="
        f"{any(c.get('found_fourth')=='true' for c in cand)}\n",
        encoding="utf-8",
    )
    print(f"Wrote {audit_path}")
    print(f"Wrote {summary}")


if __name__ == "__main__":
    main()
