"""Audit B3 cliques and Fourier-family transition charts.

This is the diagnostic stage of the k3_fourier_structure campaign. It consumes
the archived pool format (V has shape n x 6 and row norm sqrt(6)) and writes
machine-readable evidence. It does not certify pool completeness or prove
Fourier-family containment.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np


LAMBDAS = {
    "0.400000000000": 0.4,
    "1.047197551197": math.pi / 3,
    "2.094395102393": 2 * math.pi / 3,
}


def dita_hadamard(lam: float) -> np.ndarray:
    z = np.exp(1j * lam)
    return np.array(
        [
            [1, 1, 1, 1, 1, 1],
            [1, -1, z, -z, 1j, -1j],
            [1, -1j, 1j, 1j, -1j, -1],
            [1, 1j, -z, z, -1, -1j],
            [1, z.conjugate(), -1j, -1, -z.conjugate(), 1j],
            [1, -z.conjugate(), -1, -1j, z.conjugate(), 1j],
        ],
        dtype=complex,
    ) / np.sqrt(6)


def six_cliques(vectors: np.ndarray, tolerance: float = 1e-7) -> list[tuple[int, ...]]:
    gram = np.abs(vectors @ vectors.conj().T)
    adjacency = [
        set(np.flatnonzero((gram[i] < tolerance) & (np.arange(len(vectors)) != i)))
        for i in range(len(vectors))
    ]
    found: set[tuple[int, ...]] = set()

    def visit(current: list[int], candidates: set[int]) -> None:
        if len(current) == 6:
            found.add(tuple(sorted(current)))
            return
        while candidates:
            index = candidates.pop()
            visit(current + [index], candidates & adjacency[index])

    visit([], set(range(len(vectors))))
    return sorted(found)


def three_minus_one_test(
    matrix: np.ndarray, tolerance: float = 1e-8
) -> dict[str, object]:
    """Test the three-minus-one criterion over all dephasing charts.

    Row and column permutations followed by dephasing are represented by
    choosing a pivot row, pivot column, and target column. The dephased entry
    at row i and target column j is -1 exactly when
    M[i,j] M[pivot_row,pivot_column] +
    M[i,pivot_column] M[pivot_row,j] = 0.
    """
    if matrix.shape != (6, 6):
        raise ValueError(f"Expected a 6x6 transition matrix, got {matrix.shape}")

    flatness_error = float(np.max(np.abs(np.sqrt(6) * np.abs(matrix) - 1)))
    unitarity_error = float(np.max(np.abs(matrix.conj().T @ matrix - np.eye(6))))
    witnesses: list[dict[str, object]] = []
    best_count = 0

    for pivot_row in range(6):
        other_rows = [row for row in range(6) if row != pivot_row]
        for pivot_column in range(6):
            for target_column in range(6):
                if target_column == pivot_column:
                    continue
                pivot = matrix[pivot_row, pivot_column]
                target_pivot = matrix[pivot_row, target_column]
                residuals = [
                    (
                        matrix[row, target_column] * pivot
                        + matrix[row, pivot_column] * target_pivot
                    )
                    / (
                        matrix[row, pivot_column] * target_pivot
                    )
                    for row in other_rows
                ]
                hits = [
                    index
                    for index, residual in enumerate(residuals)
                    if abs(residual) <= tolerance
                ]
                best_count = max(best_count, len(hits))
                if len(hits) >= 3:
                    witnesses.append(
                        {
                            "pivot_row": pivot_row,
                            "pivot_column": pivot_column,
                            "target_column": target_column,
                            "minus_one_rows": [other_rows[index] for index in hits],
                            "max_dephased_minus_one_error": float(
                                max(abs(residuals[index]) for index in hits)
                            ),
                        }
                    )
                    break
            if witnesses:
                break

    is_hadamard = flatness_error <= tolerance and unitarity_error <= tolerance
    return {
        "criterion": "dephased column contains at least three entries equal to -1",
        "flatness_error": flatness_error,
        "unitarity_error": unitarity_error,
        "criterion_applicable": is_hadamard,
        "witness_found": bool(witnesses),
        "witness": witnesses[0] if witnesses else None,
        "max_minus_one_count": best_count,
        "tolerance": tolerance,
        "status": "numerical_candidate" if is_hadamard and witnesses else (
            "not_detected_numerically" if is_hadamard else "input_not_verified_hadamard"
        ),
    }


def analyze_pool(path: Path, lam: float, tolerance: float) -> dict[str, object]:
    vectors = np.load(path)["V"]
    cliques = six_cliques(vectors)
    h = dita_hadamard(lam)
    representatives = []
    for clique in cliques:
        basis = vectors[list(clique)].T / np.sqrt(6)
        transitions = {
            "B3": basis,
            "H_D_dagger_B3": h.conj().T @ basis,
            "B3_dagger_H_D": basis.conj().T @ h,
        }
        representatives.append(
            {
                "clique": [int(index) for index in clique],
                "transitions": {
                    name: {
                        "matrix": [
                            [[float(value.real), float(value.imag)] for value in row]
                            for row in matrix
                        ],
                        "Fourier_family_column_test": three_minus_one_test(
                            matrix, tolerance
                        ),
                        "transposed_Fourier_family_row_test": three_minus_one_test(
                            matrix.T, tolerance
                        ),
                    }
                    for name, matrix in transitions.items()
                },
            }
        )
    return {
        "pool": str(path),
        "lambda": lam,
        "pool_count": int(len(vectors)),
        "clique_count": len(cliques),
        "cliques": representatives,
        "status": "diagnostic_numerical_stage",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--pool-dir",
        type=Path,
        default=Path("mub_solution_attack_package_v2/01_results"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/campaigns/k3_fourier_structure/diagnostic.json"),
    )
    parser.add_argument("--pool", choices=sorted(LAMBDAS))
    parser.add_argument("--tolerance", type=float, default=1e-8)
    args = parser.parse_args()

    selected = [args.pool] if args.pool else sorted(LAMBDAS)
    results = []
    for stem in selected:
        path = args.pool_dir / f"pool_{stem}.npz"
        if not path.is_file():
            raise FileNotFoundError(f"Missing pool input: {path}")
        results.append(analyze_pool(path, LAMBDAS[stem], args.tolerance))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps({"campaign_id": "k3_fourier_structure", "runs": results}, indent=2))
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
