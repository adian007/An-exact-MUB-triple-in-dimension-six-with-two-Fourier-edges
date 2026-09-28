"""Regression controls for exhaustive order-six family-equivalence fits."""

import sys
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "python"))

from k3_family_membership import (
    FIT_TOLERANCE,
    build_stage_b_controls,
    fit_fourier_both_orientations,
    fit_fourier_family,
    fit_two_circulant,
    fourier_family,
    random_equivalent,
    two_circulant_constraint,
    two_circulant_family,
)


class FamilyMembershipFitTests(unittest.TestCase):
    def test_fourier_family_parameter_placement_matches_cited_convention(self):
        a = np.exp(0.371j)
        b = np.exp(-1.117j)
        matrix = fourier_family(a, b)
        fit = fit_fourier_family(matrix)
        self.assertEqual(fit.status, "FIT")
        self.assertLessEqual(fit.consistency_error, FIT_TOLERANCE)
        self.assertLessEqual(fit.total_max_residual, FIT_TOLERANCE)
        self.assertEqual(fit.charts_checked, 720 * 720)
        np.testing.assert_allclose(abs(fit.a), 1, atol=FIT_TOLERANCE)
        np.testing.assert_allclose(abs(fit.b), 1, atol=FIT_TOLERANCE)

    def test_fourier_fit_keeps_direct_and_transposed_families_separate(self):
        base = fourier_family(np.exp(0.28j), np.exp(-0.91j))
        equivalent = random_equivalent(base, 4201)
        direct, transpose_of_equivalent = fit_fourier_both_orientations(equivalent)
        self.assertEqual(direct.status, "FIT")
        self.assertEqual(transpose_of_equivalent.status, "NONE")
        self.assertEqual(fit_fourier_family(base).status, "FIT")
        self.assertEqual(fit_fourier_family(base.T).status, "NONE")
        transpose_direct, transpose_via_transpose = fit_fourier_both_orientations(
            base.T
        )
        self.assertEqual(transpose_direct.status, "NONE")
        self.assertEqual(transpose_via_transpose.status, "FIT")

    def test_two_circulant_fit_checks_constraint_and_permutation_equivalence(self):
        omega = np.exp(2j * np.pi / 3)
        parameters = (1 + 0j, omega, 1 + 0j, omega)
        self.assertLess(abs(two_circulant_constraint(*parameters)), FIT_TOLERANCE)
        matrix = random_equivalent(two_circulant_family(*parameters), 4217)
        fit = fit_two_circulant(matrix)
        self.assertEqual(fit.status, "FIT")
        self.assertLessEqual(fit.matrix_max_residual, FIT_TOLERANCE)
        self.assertLessEqual(fit.constraint_3_residual, FIT_TOLERANCE)
        self.assertEqual(fit.charts_checked, 720 * 720)
        self.assertGreater(fit.solver_starts, 0)

    def test_positive_controls_are_members_of_their_expected_fit_family(self):
        for case in build_stage_b_controls():
            if case["expected"] == "none":
                continue
            matrix = case["matrix"]
            if case["expected"] == "F":
                self.assertEqual(fit_fourier_family(matrix).status, "FIT", case["name"])
            elif case["expected"] == "F^T":
                self.assertEqual(fit_fourier_family(matrix.T).status, "FIT", case["name"])
            else:
                self.assertEqual(fit_two_circulant(matrix).status, "FIT", case["name"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
