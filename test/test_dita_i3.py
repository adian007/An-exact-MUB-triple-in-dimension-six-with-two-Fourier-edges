import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "python"))

import dita_i3
import analyze_i3_singular


class DiitaI3Tests(unittest.TestCase):
    def test_exact_hadamard_and_single_vector_system(self):
        hbar = dita_i3.H_D.xreplace(
            {dita_i3.z: dita_i3.t, dita_i3.t: dita_i3.z, dita_i3.I: -dita_i3.I}
        )
        defect = (dita_i3.H_D * hbar.T - 6 * sp.eye(6)).applyfunc(
            dita_i3.reduce_parameter_relation
        )
        self.assertEqual(defect, sp.zeros(6))
        self.assertEqual(len(dita_i3.single_vector_equations()), 6)
        self.assertEqual(len(dita_i3.single_vector_relations()), 6)

    def test_single_vector_symmetries(self):
        equations = dita_i3.single_vector_equations()
        for equation in equations:
            conjugate = dita_i3.conjugate_expression(equation)
            self.assertIn(0, [sp.simplify(conjugate - candidate) for candidate in equations])

        shift_targets = (1, 2, 4, 3, 5, 6)
        for equation, target in zip(equations, shift_targets):
            self.assertEqual(
                sp.simplify(dita_i3.shift_pi_expression(equation) - equations[target - 1]),
                0,
            )

    def test_fixed_z_symmetry_certificates_and_a4_order(self):
        for U, Q in dita_i3.fixed_z_symmetry_generators():
            self.assertEqual(
                (U * dita_i3.H_D * Q - dita_i3.H_D).applyfunc(
                    dita_i3.reduce_parameter_relation
                ),
                sp.zeros(6),
            )
        self.assertEqual(len(dita_i3.row_permutation_group()), 12)

    def test_branch_ideal_counts_and_orthogonality_equations(self):
        self.assertEqual(dita_i3.branch_ideal_counts(), {
            "variables": 62,
            "generators": 97,
            "orthogonality_equations": 30,
        })
        self.assertEqual(len(dita_i3.branch_orthogonality_equations()), 30)

    def test_numeric_vector_orbit_has_a4_order_at_generic_parameter(self):
        phases = np.array([0.1, 0.7, 1.8, 2.6, 4.0, 5.2])
        vector = np.exp(1j * phases)
        orbit = dita_i3.a4_vector_orbit(vector, np.exp(0.371j))
        self.assertEqual(len(orbit), 12)
        decomposition = dita_i3.a4_vector_orbits(
            orbit, np.exp(0.371j)
        )
        self.assertEqual([entry["orbit_size"] for entry in decomposition], [12])
        self.assertEqual(decomposition[0]["stabilizer_order"], 1)

        with self.assertRaisesRegex(RuntimeError, "outside the supplied set"):
            dita_i3.a4_vector_orbits([vector], np.exp(0.371j))
    def test_clique_orbit_rejects_missing_images_and_groups_complete_orbit(self):
        parameter = np.exp(0.371j)
        generators = dita_i3._numeric_generator_matrices(parameter)
        initial = [
            np.exp(1j * (np.arange(6) * (0.19 + index * 0.07) + index * 0.31))
            for index in range(6)
        ]
        cliques = [initial]
        queue = [initial]
        while queue:
            current = queue.pop()
            for generator in generators:
                image = [dita_i3._phase_canonical(generator @ vector) for vector in current]
                if not any(
                    dita_i3._same_projective_set(image, known, 1e-8)
                    for known in cliques
                ):
                    cliques.append(image)
                    queue.append(image)

        orbits = dita_i3.a4_clique_orbits(cliques, parameter)
        self.assertEqual(len(orbits), 1)
        self.assertEqual(orbits[0]["orbit_size"], len(cliques))
        with self.assertRaisesRegex(RuntimeError, "outside the supplied clique set"):
            dita_i3.a4_clique_orbits([initial], parameter)

    def test_generated_equation_module_is_importable(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir)
            dita_i3._write_artifacts(output_dir)
            spec = importlib.util.spec_from_file_location(
                "generated_i3_equations", output_dir / "I3_equations.py"
            )
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            self.assertEqual(len(module.F), 6)
            self.assertEqual(len(module.relations), 6)
            self.assertTrue((output_dir / "run_manifest.json").is_file())

    def test_singular_analysis_does_not_assume_two_orbits(self):
        parameter = np.exp(1j * analyze_i3_singular.LAMBDA_STAR)
        initial = np.exp(1j * np.array([0.1, 0.7, 1.8, 2.6, 4.0, 5.2]))
        orbit = dita_i3.a4_vector_orbit(initial, parameter)
        roots = [
            {
                "z_re_im": [parameter.real, parameter.imag],
                "t_re_im": [parameter.real, -parameter.imag],
                "x_re_im": [[value.real, value.imag] for value in vector[1:]],
                "y_re_im": [[value.real, -value.imag] for value in vector[1:]],
                "u_re_im": [[1.0, 0.0], *[[0.0, 0.0] for _ in range(9)]],
                "system_residual": 0.0,
            }
            for vector in orbit
        ]
        with tempfile.TemporaryDirectory() as temp_dir:
            source = Path(temp_dir) / "roots.json"
            source.write_text(json.dumps({"physical_roots": roots}), encoding="utf-8")
            result = analyze_i3_singular.analyze_singular_roots(
                source, lambda_tolerance=1e-5, parameter_cluster_tolerance=1e-7
            )
        self.assertEqual(result["physical_candidate_count"], 12)
        self.assertEqual(
            result["critical_parameter_clusters"][0]["orbit_sizes"], [12]
        )
        self.assertFalse(
            result["critical_parameter_clusters"][0][
                "matches_two_12_orbit_hypothesis"
            ]
        )


if __name__ == "__main__":
    unittest.main()
