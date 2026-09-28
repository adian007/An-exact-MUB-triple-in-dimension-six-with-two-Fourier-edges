import sys
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "python"))

from analyze_k3_fourier_structure import three_minus_one_test


class K3FourierStructureTests(unittest.TestCase):
    def test_fourier_matrix_is_detected(self):
        indices = np.arange(6)
        fourier = np.exp(2j * np.pi * np.outer(indices, indices) / 6) / np.sqrt(6)
        result = three_minus_one_test(fourier)
        self.assertTrue(result["criterion_applicable"])
        self.assertTrue(result["witness_found"])
        self.assertLess(result["witness"]["max_dephased_minus_one_error"], 1e-12)

    def test_equivalent_fourier_matrix_is_detected_after_dephasing(self):
        indices = np.arange(6)
        fourier = np.exp(2j * np.pi * np.outer(indices, indices) / 6) / np.sqrt(6)
        rows = (2, 4, 0, 5, 1, 3)
        columns = (3, 1, 5, 0, 4, 2)
        row_phases = np.exp(1j * np.arange(6) * 0.37)
        column_phases = np.exp(1j * np.arange(6) * -0.29)
        equivalent = fourier[list(rows), :][:, list(columns)]
        equivalent = row_phases[:, None] * equivalent * column_phases[None, :]
        result = three_minus_one_test(equivalent)
        self.assertTrue(result["criterion_applicable"])
        self.assertTrue(result["witness_found"])

    def test_non_hadamard_input_is_not_classified(self):
        rng = np.random.default_rng(1729)
        matrix = np.exp(1j * rng.uniform(-np.pi, np.pi, (6, 6)))
        result = three_minus_one_test(matrix)
        self.assertFalse(result["criterion_applicable"])
        self.assertEqual(result["status"], "input_not_verified_hadamard")


if __name__ == "__main__":
    unittest.main()
