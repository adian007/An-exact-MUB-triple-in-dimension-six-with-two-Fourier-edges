"""Negative controls for the legacy one-column/three-row diagnostic.

Run from the repository root:
    python test/test_three_minus_one_negative_controls.py

The report distinguishes the implemented single-column/three-row predicate
from Theorem 1's three-distinct-columns predicate. These are numerical
diagnostics, not classifications or certificates.
"""

from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "python"))

from analyze_k3_fourier_structure import three_minus_one_test
from karlsson_k6_3 import build_k6


TOLERANCES = (1e-10, 1e-8, 1e-6)
PERTURBATIONS = (1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2)


def tao_s6_unitary() -> np.ndarray:
    omega = np.exp(2j * np.pi / 3)
    omega2 = omega**2
    matrix = np.array(
        [
            [1, 1, 1, 1, 1, 1],
            [1, 1, omega, omega, omega2, omega2],
            [1, omega, 1, omega2, omega2, omega],
            [1, omega, omega2, 1, omega, omega2],
            [1, omega2, omega2, omega, 1, omega],
            [1, omega2, omega, omega2, omega, 1],
        ],
        dtype=complex,
    )
    return matrix / np.sqrt(6)


def minimum_single_column_residual(matrix: np.ndarray) -> float:
    """Minimum third-smallest -1 residual over every dephasing chart."""
    best = float("inf")
    for pivot_row in range(6):
        for pivot_column in range(6):
            for target_column in range(6):
                if target_column == pivot_column:
                    continue
                errors = []
                for row in range(6):
                    if row == pivot_row:
                        continue
                    denominator = (
                        matrix[row, pivot_column] * matrix[pivot_row, target_column]
                    )
                    if abs(denominator) == 0:
                        continue
                    error = (
                        matrix[row, target_column] * matrix[pivot_row, pivot_column]
                        + denominator
                    ) / denominator
                    errors.append(abs(error))
                if len(errors) >= 3:
                    best = min(best, float(sorted(errors)[2]))
    return best


def theorem_three_distinct_column_residual(matrix: np.ndarray) -> float:
    """Minimum threshold residual for three columns each containing a -1."""
    best = float("inf")
    for pivot_row in range(6):
        other_rows = [row for row in range(6) if row != pivot_row]
        for pivot_column in range(6):
            per_column = []
            for target_column in range(6):
                if target_column == pivot_column:
                    continue
                errors = []
                for row in other_rows:
                    denominator = (
                        matrix[row, pivot_column] * matrix[pivot_row, target_column]
                    )
                    if abs(denominator) == 0:
                        continue
                    error = (
                        matrix[row, target_column] * matrix[pivot_row, pivot_column]
                        + denominator
                    ) / denominator
                    errors.append(abs(error))
                per_column.append(min(errors, default=float("inf")))
            if len(per_column) >= 3:
                best = min(best, float(sorted(per_column)[2]))
    return best


def hadamard_errors(matrix: np.ndarray) -> tuple[float, float]:
    return (
        float(np.max(np.abs(np.sqrt(6) * np.abs(matrix) - 1))),
        float(np.max(np.abs(matrix.conj().T @ matrix - np.eye(6)))),
    )


def control_matrix_summary(matrix: np.ndarray, tolerance: float) -> dict[str, object]:
    implementation = three_minus_one_test(matrix, tolerance)
    theorem_residual = theorem_three_distinct_column_residual(matrix)
    return {
        "implementation_witness_found": implementation["witness_found"],
        "implementation_min_dephased_residual": minimum_single_column_residual(matrix),
        "theorem_three_distinct_columns_candidate": theorem_residual <= tolerance,
        "theorem_condition_minimum_residual": theorem_residual,
        "flatness_error": implementation["flatness_error"],
        "unitarity_error": implementation["unitarity_error"],
        "criterion_applicable": implementation["criterion_applicable"],
    }


class ThreeMinusOneNegativeControls(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.controls: list[tuple[str, np.ndarray]] = []
        for index, (theta, phi, lam) in enumerate(
            ((0.31, 0.73, 0.41), (0.82, 1.19, 2.07))
        ):
            matrix = build_k6(theta, phi, lam) / np.sqrt(6)
            cls.controls.append((f"karlsson_off_dita_{index + 1}", matrix))
        cls.controls.append(("tao_S6", tao_s6_unitary()))

        diagnostic_path = (
            ROOT / "results" / "campaigns" / "k3_fourier_structure" / "diagnostic.json"
        )
        data = json.loads(diagnostic_path.read_text(encoding="utf-8"))
        cls.cliques = []
        for run in data["runs"]:
            pool = np.load(ROOT / run["pool"])["V"]
            for clique_index, clique in enumerate(run["cliques"]):
                vectors = pool[clique["clique"]]
                basis = vectors.T / np.sqrt(6)
                cls.cliques.append(
                    (f"lambda_{run['lambda']:.12f}_clique_{clique_index}", basis)
                )

    def test_negative_controls_and_perturbation_sensitivity(self) -> None:
        report: dict[str, object] = {
            "test_scope": "diagnostic only; no theorem-level classification",
            "tolerances": list(TOLERANCES),
            "perturbations": list(PERTURBATIONS),
            "negative_controls": {},
            "perturbations_by_clique": {},
            "accepted_rejected_residual_gap": {},
        }

        for name, matrix in self.controls:
            flatness, unitarity = hadamard_errors(matrix)
            self.assertLess(flatness, 1e-10, name)
            self.assertLess(unitarity, 1e-10, name)
            report["negative_controls"][name] = {
                str(tol): control_matrix_summary(matrix, tol) for tol in TOLERANCES
            }

        report["unperturbed_b3_cliques"] = {
            name: {
                str(tol): control_matrix_summary(matrix, tol) for tol in TOLERANCES
            }
            for name, matrix in self.cliques
        }

        for name, matrix in self.cliques:
            per_clique: dict[str, object] = {}
            for delta in PERTURBATIONS:
                perturbed = matrix.copy()
                perturbed[1, 1] *= np.exp(1j * delta)
                per_clique[str(delta)] = {
                    str(tol): control_matrix_summary(perturbed, tol)
                    for tol in TOLERANCES
                }
            report["perturbations_by_clique"][name] = per_clique

        for tol in TOLERANCES:
            accepted_residuals = []
            rejected_residuals = []
            eligible_cases = [*self.controls, *self.cliques]
            for _, matrix in eligible_cases:
                summary = control_matrix_summary(matrix, tol)
                if not summary["criterion_applicable"]:
                    continue
                target = (
                    accepted_residuals
                    if summary["implementation_witness_found"]
                    else rejected_residuals
                )
                target.append(summary["implementation_min_dephased_residual"])
            for _, matrix in self.cliques:
                for delta in PERTURBATIONS:
                    perturbed = matrix.copy()
                    perturbed[1, 1] *= np.exp(1j * delta)
                    summary = control_matrix_summary(perturbed, tol)
                    if not summary["criterion_applicable"]:
                        continue
                    target = (
                        accepted_residuals
                        if summary["implementation_witness_found"]
                        else rejected_residuals
                    )
                    target.append(summary["implementation_min_dephased_residual"])
            report["accepted_rejected_residual_gap"][str(tol)] = {
                "eligible_accepted_case_count": len(accepted_residuals),
                "eligible_rejected_case_count": len(rejected_residuals),
                "largest_accepted_minimum_residual": max(accepted_residuals, default=None),
                "smallest_rejected_minimum_residual": min(rejected_residuals, default=None),
                "gap": (
                    min(rejected_residuals) - max(accepted_residuals)
                    if accepted_residuals and rejected_residuals
                    else None
                ),
            }

        print(json.dumps(report, indent=2))


if __name__ == "__main__":
    unittest.main(verbosity=2)
