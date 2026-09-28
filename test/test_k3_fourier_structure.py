"""Tests for the one-column diagnostic and Matszangosz-Szollosi Theorem 1.

Theorem 1 (Designs, Codes and Cryptography 92 (2024), Theorem 1):
for a normalized order-six CHM, three distinct columns each containing
at least one -1 characterize membership in the transposed Fourier family
or the 2-circulant family. The row analogue is inferred by transposing.
The 2-circulant X_6 matrix and parameter constraint are equations (2)-(3).
"""

import json
import sys
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "python"))

from analyze_k3_fourier_structure import (
    analyze_pool,
    three_distinct_columns_theorem_test,
    three_minus_one_test,
)
from karlsson_k6_3 import build_k6


TOLERANCES = (1e-10, 1e-8, 1e-6)
THEOREM_URL = "https://link.springer.com/article/10.1007/s10623-024-01503-w"


def fourier6() -> np.ndarray:
    index = np.arange(6)
    return np.exp(2j * np.pi * np.outer(index, index) / 6) / np.sqrt(6)


def tao_s6() -> np.ndarray:
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


def dita_d6(parameter: float) -> np.ndarray:
    """Repository Dita D(x), Bengtsson quant-ph/0610161, Eq. (11)."""
    z = np.exp(2j * np.pi * parameter)
    zb = np.conjugate(z)
    i = 1j
    return np.array(
        [
            [1, 1, 1, 1, 1, 1],
            [1, -1, i, -i, -i, i],
            [1, i, -1, i * z, -i * z, -i],
            [1, -i, i * zb, -1, i, -i * zb],
            [1, -i, -i * zb, i, -1, i * zb],
            [1, i, -i, -i * z, i * z, -1],
        ],
        dtype=complex,
    ) / np.sqrt(6)


def two_circulant_constraint(beta, gamma, epsilon, phi) -> complex:
    return (
        beta * gamma * epsilon**2
        + beta * gamma * phi
        + beta**2 * epsilon * phi
        + gamma * epsilon * phi
        + beta * gamma**2 * epsilon * phi
        + beta * gamma * epsilon * phi**2
    )


def two_circulant_x6(beta, gamma, epsilon, phi) -> np.ndarray:
    """Equation (2) of the linked paper; use only with Eq. (3) satisfied."""
    return np.array(
        [
            [1, 1, 1, 1, 1, 1],
            [1, -1, -1 / (gamma * epsilon), -1 / (beta * phi),
             1 / (gamma * epsilon), 1 / (beta * phi)],
            [1, -epsilon / beta, -1, -epsilon / (gamma * phi),
             epsilon / (gamma * phi), epsilon / beta],
            [1, -phi / gamma, -phi / (beta * epsilon),
             -1, phi / gamma, phi / (beta * epsilon)],
            [1, epsilon / beta, phi / (beta * epsilon), 1 / (beta * phi),
             1 / (beta * gamma), gamma / beta**2],
            [1, phi / gamma, 1 / (gamma * epsilon), epsilon / (gamma * phi),
             beta / gamma**2, 1 / (beta * gamma)],
        ],
        dtype=complex,
    ) / np.sqrt(6)


def random_equivalent(matrix: np.ndarray, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    rows = rng.permutation(6)
    columns = rng.permutation(6)
    row_phases = np.exp(1j * rng.uniform(-np.pi, np.pi, 6))
    column_phases = np.exp(1j * rng.uniform(-np.pi, np.pi, 6))
    return (
        row_phases[:, None]
        * matrix[rows, :][:, columns]
        * column_phases[None, :]
    )


def chm_errors(matrix: np.ndarray) -> tuple[float, float]:
    return (
        float(np.max(np.abs(np.sqrt(6) * np.abs(matrix) - 1))),
        float(np.max(np.abs(matrix.conj().T @ matrix - np.eye(6)))),
    )


def theorem_residual(matrix: np.ndarray, tolerance: float) -> float | None:
    result = three_distinct_columns_theorem_test(matrix, tolerance)
    return result["minimum_worst_of_three_residual"]


def make_two_circulant_control(
    phase: complex | None = None,
) -> tuple[np.ndarray, dict[str, complex]]:
    phase = np.exp(2j * np.pi / 3) if phase is None else phase
    parameters = {
        "beta": 1 + 0j,
        "gamma": phase,
        "epsilon": 1 + 0j,
        "phi": phase,
    }
    residual = two_circulant_constraint(**parameters)
    return two_circulant_x6(**parameters), {**parameters, "constraint_residual": residual}


def control_summary(matrix: np.ndarray, tolerance: float) -> dict[str, object]:
    direct = three_distinct_columns_theorem_test(matrix, tolerance)
    transposed = three_distinct_columns_theorem_test(matrix.T, tolerance)
    return {
        "criterion_applicable": direct["criterion_applicable"],
        "witness_found": direct["witness_found"],
        "minimum_worst_of_three_residual": direct["minimum_worst_of_three_residual"],
        "pivot_choice": direct["pivot_choice"],
        "columns_achieving_minimum": direct["columns_achieving_minimum"],
        "transpose_application": {
            "criterion_applicable": transposed["criterion_applicable"],
            "witness_found": transposed["witness_found"],
            "minimum_worst_of_three_residual": transposed[
                "minimum_worst_of_three_residual"
            ],
        },
    }


def compact_summary(matrix: np.ndarray, tolerance: float) -> dict[str, object]:
    direct = three_distinct_columns_theorem_test(matrix, tolerance)
    transposed = three_distinct_columns_theorem_test(matrix.T, tolerance)
    return {
        "criterion_applicable": direct["criterion_applicable"],
        "witness_found": direct["witness_found"],
        "residual": direct["minimum_worst_of_three_residual"],
        "transpose_criterion_applicable": transposed["criterion_applicable"],
        "transpose_witness_found": transposed["witness_found"],
        "transpose_residual": transposed["minimum_worst_of_three_residual"],
    }


class K3FourierStructureTests(unittest.TestCase):
    def test_single_column_helper_is_named_as_diagnostic_not_family_test(self):
        matrix = fourier6()
        result = three_minus_one_test(matrix)
        self.assertTrue(result["criterion_applicable"])
        self.assertTrue(result["witness_found"])
        self.assertIn("not Theorem 1", result["criterion"])
        self.assertIn("criterion_applicable", list(result))
        self.assertEqual(list(result).index("witness_found"), list(result).index("criterion_applicable") + 1)
        theorem_result = three_distinct_columns_theorem_test(matrix)
        self.assertEqual(len(theorem_result["pivot_chart_results"]), 36)
        for chart in theorem_result["pivot_chart_results"]:
            self.assertEqual(
                abs(
                    list(chart).index("witness_found")
                    - list(chart).index("criterion_applicable")
                ),
                1,
            )

    def test_theorem_predicate_rejects_unverified_hadamard_input_without_verdict(self):
        rng = np.random.default_rng(1729)
        matrix = np.exp(1j * rng.uniform(-np.pi, np.pi, (6, 6)))
        result = three_distinct_columns_theorem_test(matrix)
        self.assertFalse(result["criterion_applicable"])
        self.assertFalse(result["witness_found"])
        self.assertEqual(result["status"], "not_verified_hadamard")
        self.assertIsNone(result["verdict"])

    def test_theorem_positive_controls_survive_equivalence(self):
        omega = np.exp(2j * np.pi / 3)
        two_circulant = [
            ("X6_omega", *make_two_circulant_control(omega)),
            ("X6_omega_squared", *make_two_circulant_control(omega**2)),
        ]
        self.assertTrue(
            all(abs(parameters["constraint_residual"]) < 1e-12
                for _, _, parameters in two_circulant)
        )
        controls = [
            ("F6", fourier6()),
            ("F6_transpose", fourier6().T),
            *[(name, matrix) for name, matrix, _ in two_circulant],
        ]
        two_circulant_parameters = {
            name: {
                key: [float(value.real), float(value.imag)]
                for key, value in parameters.items()
                if key != "constraint_residual"
            }
            | {"constraint_residual": float(abs(parameters["constraint_residual"]))}
            for name, _, parameters in two_circulant
        }
        report = {}
        for name, matrix in controls:
            variants = [("base", matrix)] + [
                (f"random_equivalence_{seed}", random_equivalent(matrix, seed))
                for seed in (104, 205, 306)
            ]
            report[name] = {}
            for variant_name, variant in variants:
                flatness, unitarity = chm_errors(variant)
                self.assertLess(flatness, 1e-10, (name, variant_name))
                self.assertLess(unitarity, 1e-10, (name, variant_name))
                for tolerance in TOLERANCES:
                    result = three_distinct_columns_theorem_test(variant, tolerance)
                    self.assertTrue(
                        result["criterion_applicable"],
                        (name, variant_name, tolerance),
                    )
                    self.assertTrue(result["witness_found"], (name, variant_name, tolerance))
                    transpose_result = three_distinct_columns_theorem_test(
                        variant.T, tolerance
                    )
                    self.assertTrue(
                        transpose_result["criterion_applicable"],
                        (name, variant_name, tolerance, "transpose"),
                    )
                    self.assertTrue(
                        transpose_result["witness_found"],
                        (name, variant_name, tolerance, "transpose"),
                    )
                report[name][variant_name] = compact_summary(variant, 1e-8)
        print(
            "POSITIVE_CONTROL_SUMMARY",
            json.dumps(
                {
                    "2_circulant_parameters_and_constraint_residuals": two_circulant_parameters,
                    "cases": report,
                },
                indent=2,
            ),
        )
        print("POSITIVE_CONTROL_REPORT", json.dumps(report, indent=2))

    def test_negative_controls_and_genuine_chm_family_sweep(self):
        x6, x6_parameters = make_two_circulant_control()
        controls = [
            ("Tao_S6", tao_s6()),
            ("D6_Dita_x_0.137", dita_d6(0.137)),
            ("Karlsson_0.31_0.73_0.41", build_k6(0.31, 0.73, 0.41) / np.sqrt(6)),
            ("Karlsson_0.82_1.19_2.07", build_k6(0.82, 1.19, 2.07) / np.sqrt(6)),
        ]
        for name, matrix in controls:
            flatness, unitarity = chm_errors(matrix)
            self.assertLess(flatness, 1e-10, name)
            self.assertLess(unitarity, 1e-10, name)

        negative_results = {}
        for name, matrix in controls:
            negative_results[name] = {
                str(tolerance): compact_summary(matrix, tolerance)
                for tolerance in TOLERANCES
            }
        self.assertTrue(
            all(
                not three_distinct_columns_theorem_test(controls[0][1], tolerance)[
                    "witness_found"
                ]
                for tolerance in TOLERANCES
            ),
            "Tao_S6",
        )
        self.assertTrue(
            all(
                three_distinct_columns_theorem_test(controls[1][1], tolerance)[
                    "witness_found"
                ]
                for tolerance in TOLERANCES
            ),
            "D6 is a theorem-positive control, not a negative control",
        )
        for name, matrix in controls[2:]:
            self.assertTrue(
                all(
                    not three_distinct_columns_theorem_test(matrix, tolerance)["witness_found"]
                    for tolerance in TOLERANCES
                ),
                name,
            )

        d6_sweep = []
        for parameter in np.linspace(0, 1, 21):
            matrix = dita_d6(float(parameter))
            flatness, unitarity = chm_errors(matrix)
            self.assertLess(flatness, 1e-10)
            self.assertLess(unitarity, 1e-10)
            self.assertTrue(
                three_distinct_columns_theorem_test(matrix, 1e-8)["witness_found"]
            )
            d6_sweep.append(
                {
                    "family": "Diță D6(x)",
                    "parameter_x": float(parameter),
                    "result": compact_summary(matrix, 1e-8),
                }
            )

        theta, phi = 0.82, 1.19
        broad_lambdas = np.linspace(0, 2 * np.pi, 65)
        local_lambdas = np.arange(1.95, 2.1001, 0.005)
        lambda_sweep = []
        for lam in np.unique(np.concatenate((broad_lambdas, local_lambdas, [2.07]))):
            matrix = build_k6(theta, phi, float(lam)) / np.sqrt(6)
            flatness, unitarity = chm_errors(matrix)
            self.assertLess(flatness, 1e-10)
            self.assertLess(unitarity, 1e-10)
            lambda_sweep.append(
                {
                    "lambda": float(lam),
                    "minimum_worst_of_three_residual": theorem_residual(matrix, 1e-8),
                }
            )

        minimum = min(
            lambda_sweep,
            key=lambda item: item["minimum_worst_of_three_residual"],
        )
        lambda_2_matrix = build_k6(theta, phi, 2.0) / np.sqrt(6)
        lambda_2_07_matrix = build_k6(theta, phi, 2.07) / np.sqrt(6)
        tolerance_sensitivity = {
            "lambda_2.0": {
                str(tolerance): compact_summary(lambda_2_matrix, tolerance)
                for tolerance in (1e-8, 1e-6, 1e-5)
            },
            "lambda_2.07": {
                str(tolerance): compact_summary(lambda_2_07_matrix, tolerance)
                for tolerance in (1e-8, 1e-2, 2e-2)
            },
        }
        d6_residuals = [
            point["result"]["residual"]
            for point in d6_sweep
        ]
        print(
            "CONTROL_SUMMARY",
            json.dumps(
                {
                    "source": THEOREM_URL,
                    "tolerances": TOLERANCES,
                    "2-circulant_constraint_residual": float(
                        abs(x6_parameters["constraint_residual"])
                    ),
                    "2-circulant_parameter_values": {
                        key: [float(value.real), float(value.imag)]
                        for key, value in x6_parameters.items()
                        if key != "constraint_residual"
                    },
                    "negative_controls": {
                        name: {
                            "by_tolerance": negative_results[name],
                            "flatness_and_unitarity": chm_errors(matrix),
                        }
                        for name, matrix in controls
                    },
                    "D6_Dita_family_sweep": {
                        "family": "repository dita_D(x), Bengtsson quant-ph/0610161, Eq. (11)",
                        "parameter_values": [point["parameter_x"] for point in d6_sweep],
                        "all_theorem_witnesses": all(
                            point["result"]["witness_found"] for point in d6_sweep
                        ),
                        "residual_range": [min(d6_residuals), max(d6_residuals)],
                    },
                    "Karlsson_fixed_theta_phi_lambda_sweep": {
                        "theta": theta,
                        "phi": phi,
                        "tolerance": 1e-8,
                        "minimum": minimum,
                        "tolerance_sensitivity": tolerance_sensitivity,
                        "at_lambda_2_07": next(
                            point
                            for point in lambda_sweep
                            if abs(point["lambda"] - 2.07) < 1e-12
                        ),
                        "local_table": [
                            point
                            for point in lambda_sweep
                            if 1.95 <= point["lambda"] <= 2.1000001
                        ],
                    },
                    "B6_status": "not_constructed; no B6/Bjoerck CHM constructor found in repository source",
                },
                indent=2,
            ),
        )

    def test_transition_theorem_audit_covers_all_supplied_cliques_and_orientations(self):
        source_path = ROOT / "results" / "campaigns" / "k3_fourier_structure" / "diagnostic.json"
        legacy = json.loads(source_path.read_text(encoding="utf-8"))
        report = []
        grouped: dict[tuple[float, str, str], list[float]] = {}
        for run in legacy["runs"]:
            analyzed = analyze_pool(ROOT / run["pool"], run["lambda"], 1e-8)
            self.assertEqual(analyzed["clique_count"], 4)
            for clique_index, clique in enumerate(analyzed["cliques"]):
                for transition_name, transition in clique["transitions"].items():
                    for application_name in (
                        "theorem1_column_application_to_matrix",
                        "theorem1_row_application_inferred_via_transpose",
                    ):
                        result = transition[application_name]
                        self.assertTrue(result["criterion_applicable"])
                        report.append(
                            {
                                "lambda": run["lambda"],
                                "clique": clique_index,
                                "transition": transition_name,
                                "application": application_name,
                                "witness_found": result["witness_found"],
                                "minimum_worst_of_three_residual": result[
                                    "minimum_worst_of_three_residual"
                                ],
                            }
                        )
                        grouped.setdefault(
                            (run["lambda"], transition_name, application_name), []
                        ).append(
                            result["minimum_worst_of_three_residual"]
                        )
        self.assertEqual(len(report), 3 * 4 * 3 * 2)
        grouped_report = [
            {
                "lambda": key[0],
                "transition": key[1],
                "application": key[2],
                "clique_count": len(residuals),
                "minimum_residual": min(residuals),
                "maximum_residual": max(residuals),
                "spread": max(residuals) - min(residuals),
            }
            for key, residuals in sorted(grouped.items())
        ]
        print("TRANSITION_THEOREM_AUDIT_SUMMARY", json.dumps(grouped_report, indent=2))


if __name__ == "__main__":
    unittest.main(verbosity=2)
