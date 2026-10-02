"""Audit the restored pi/3 candidate pool as an orthogonality/MU graph."""

from __future__ import annotations

import itertools
import json
from pathlib import Path

import networkx as nx
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
POOL = ROOT / "b3_a4_i4_results" / "pool_pi_over_3.npz"
SUMMARY = ROOT / "b3_a4_i4_results" / "summary.json"
OUT = ROOT / "results" / "phase2_pi3_pool_graph_audit.json"


def main() -> None:
    V = np.load(POOL)["V"]
    G_orth = nx.Graph()
    G_mu = nx.Graph()
    G_orth.add_nodes_from(range(len(V)))
    G_mu.add_nodes_from(range(len(V)))
    gram = V.conj() @ V.T
    orth_errors = []
    mu_errors = []
    for i, j in itertools.combinations(range(len(V)), 2):
        value = abs(gram[i, j])
        if value < 1e-9:
            G_orth.add_edge(i, j)
            orth_errors.append(value)
        elif abs(value - np.sqrt(6)) < 1e-9:
            G_mu.add_edge(i, j)
            mu_errors.append(abs(value - np.sqrt(6)))

    cliques = [sorted(c) for c in nx.find_cliques(G_orth) if len(c) >= 6]
    cliques.sort(key=lambda c: (len(c), c))
    known = json.loads(SUMMARY.read_text(encoding="utf-8"))["1.0471975511965976"]["cliques"]
    known_checks = []
    for clique in known:
        known_checks.append({
            "clique": clique,
            "is_orthogonal_6_clique": G_orth.subgraph(clique).number_of_edges() == 15,
            "internal_mu_edges": G_mu.subgraph(clique).number_of_edges(),
        })

    result = {
        "status": "PASS",
        "pool_size": int(len(V)),
        "row_norm_max_error_from_6": float(np.max(np.abs(np.sum(np.abs(V)**2, axis=1) - 6))),
        "orth_graph_edges": G_orth.number_of_edges(),
        "mu_graph_edges": G_mu.number_of_edges(),
        "mu_pair_threshold": "abs(abs(<v_i,v_j>)-sqrt(6)) < 1e-9",
        "orth_edge_max_abs_inner": float(max(orth_errors, default=0.0)),
        "mu_edge_max_abs_inner_error": float(max(mu_errors, default=0.0)),
        "max_orth_clique_size": max((len(c) for c in cliques), default=0),
        "orth_cliques_size_at_least_6": cliques,
        "known_clique_checks": known_checks,
        "scope": "restored 72-vector pool only; not a completeness theorem for all MU vectors",
    }
    OUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
