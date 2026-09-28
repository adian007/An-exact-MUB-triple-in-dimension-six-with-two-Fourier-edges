"""Exact symbolic model and symmetry actions for the Diţă-circle I3 problem."""

from __future__ import annotations

import argparse
from collections import deque
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess
import sys
from typing import Sequence

import numpy as np
import sympy as sp

I = sp.I
z, t = sp.symbols("z t")
x = sp.symbols("x1:6")
y = sp.symbols("y1:6")

H_D = sp.Matrix(
    [
        [1, 1, 1, 1, 1, 1],
        [1, -1, z, -z, I, -I],
        [1, -I, I, I, -I, -1],
        [1, I, -z, z, -1, -I],
        [1, t, -I, -1, -t, I],
        [1, -t, -1, -I, t, I],
    ]
)


def single_vector_equations() -> list[sp.Expr]:
    """Return the six exact unbiasedness equations for v=(1,x1,...,x5)."""
    hbar = H_D.xreplace({z: t, t: z, I: -I})
    v = sp.Matrix([1, *x])
    vbar = sp.Matrix([1, *y])
    equations = []
    for column in range(6):
        inner = (vbar.T * H_D[:, column])[0]
        conjugate_inner = (v.T * hbar[:, column])[0]
        equations.append(sp.expand(inner * conjugate_inner - 6))
    return equations


def single_vector_relations() -> list[sp.Expr]:
    """Return parameter and inverse-variable relations defining the torus."""
    return [z * t - 1, *(xj * yj - 1 for xj, yj in zip(x, y))]


def reduce_parameter_relation(expression: sp.Expr) -> sp.Expr:
    """Reduce a Laurent expression modulo the shared relation z*t=1."""
    return sp.factor(sp.cancel(expression.subs(t, 1 / z)))


def branch_orthogonality_equations() -> list[str]:
    """Describe the 30 equations joining six gauged flat vectors."""
    equations = []
    for a in range(1, 7):
        for b in range(a + 1, 7):
            ab = " + ".join(
                ["1", *(f"y{a}{j}*x{b}{j}" for j in range(1, 6))]
            )
            ba = " + ".join(
                ["1", *(f"x{a}{j}*y{b}{j}" for j in range(1, 6))]
            )
            equations.extend((f"({ab}) = 0", f"({ba}) = 0"))
    return equations


def branch_ideal_counts() -> dict[str, int]:
    """Return the unreduced six-vector branch ideal's generator/variable counts."""
    return {
        "variables": 2 + 6 * 10,
        "generators": 1 + 6 * 5 + 6 * 6 + 30,
        "orthogonality_equations": 30,
    }


def conjugate_expression(expression: sp.Expr) -> sp.Expr:
    """Apply z<->t, x<->y and complex conjugation exactly."""
    substitution = {z: t, t: z, I: -I}
    substitution.update({xj: yj for xj, yj in zip(x, y)})
    substitution.update({yj: xj for xj, yj in zip(x, y)})
    return sp.expand(expression.xreplace(substitution))


def shift_pi_expression(expression: sp.Expr) -> sp.Expr:
    """Apply the exact z->-z, x4<->x5 parameter-shift symmetry."""
    substitution = {
        z: -z,
        t: -t,
        x[3]: x[4],
        x[4]: x[3],
        y[3]: y[4],
        y[4]: y[3],
    }
    return sp.expand(expression.xreplace(substitution))


def fixed_z_symmetry_generators() -> tuple[tuple[sp.Matrix, sp.Matrix], ...]:
    """Return the two exact monomial identities U*H_D*Q = H_D."""
    U2 = sp.zeros(6)
    for row, source in enumerate((0, 3, 2, 1, 5, 4)):
        U2[row, source] = 1
    Q2 = sp.Matrix(
        [
            [1, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 1, 0],
            [0, 0, 0, 1, 0, 0],
            [0, 0, 1, 0, 0, 0],
            [0, 1, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 1],
        ]
    )
    U3 = sp.Matrix(
        [
            [0, -1, 0, 0, 0, 0],
            [0, 0, 0, 0, z, 0],
            [0, 0, 0, -I, 0, 0],
            [0, 0, 0, 0, 0, -z],
            [1, 0, 0, 0, 0, 0],
            [0, 0, I, 0, 0, 0],
        ]
    )
    Q3 = sp.Matrix(
        [
            [0, 0, 0, -1, 0, 0],
            [1, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, -1 / z, 0],
            [0, 1 / z, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, I],
            [0, 0, -I, 0, 0, 0],
        ]
    )
    return ((U2, Q2), (U3, Q3))


def row_permutation_group() -> set[tuple[int, ...]]:
    """Generate the 12-element row-permutation group from the two certificates."""
    identity = tuple(range(6))
    generators = ((0, 3, 2, 1, 5, 4), (1, 4, 3, 5, 0, 2))
    group = {identity}
    queue = deque([identity])
    while queue:
        current = queue.popleft()
        for generator in generators:
            product = tuple(current[generator[i]] for i in range(6))
            if product not in group:
                group.add(product)
                queue.append(product)
    return group


def _numeric_generator_matrices(parameter: complex) -> tuple[np.ndarray, ...]:
    matrices = []
    for U, _ in fixed_z_symmetry_generators():
        evaluated = U.subs({z: parameter, t: 1 / parameter})
        matrices.append(np.asarray(evaluated.evalf(), dtype=np.complex128))
    return tuple(matrices)


def _phase_canonical(vector: np.ndarray) -> np.ndarray:
    vector = np.asarray(vector, dtype=np.complex128).reshape(6)
    pivot = vector[0]
    if abs(pivot) < 1e-14:
        raise ValueError("The vector's first coordinate is zero; gauge is undefined.")
    return vector / pivot


def _same_projective_vector(left: np.ndarray, right: np.ndarray, tol: float) -> bool:
    left_norm = np.linalg.norm(left)
    right_norm = np.linalg.norm(right)
    if left_norm == 0 or right_norm == 0:
        raise ValueError("Zero vectors do not define projective states.")
    left_unit = left / left_norm
    right_unit = right / right_norm
    overlap = np.vdot(left_unit, right_unit)
    if abs(overlap) == 0:
        return False
    phase = overlap / abs(overlap)
    return np.linalg.norm(left_unit - right_unit * np.conj(phase)) <= tol


def a4_vector_orbit(
    vector: Sequence[complex], parameter: complex, *, tol: float = 1e-8
) -> list[np.ndarray]:
    """Return the discovered A4 orbit of a flat vector for the exact H_D gauge.

    This acts on vectors for the packaged H_D(z) representative, not directly
    on vectors returned by the project's Karlsson-family pool solver.
    """
    if abs(parameter) == 0:
        raise ValueError("The Diţă parameter z must be nonzero.")
    initial = _phase_canonical(np.asarray(vector, dtype=np.complex128))
    generators = _numeric_generator_matrices(parameter)
    orbit = [initial]
    queue = deque([initial])
    while queue:
        current = queue.popleft()
        for generator in generators:
            image = _phase_canonical(generator @ current)
            if not any(_same_projective_vector(image, known, tol) for known in orbit):
                orbit.append(image)
                queue.append(image)
                if len(orbit) > 12:
                    raise RuntimeError("Symmetry orbit exceeded the verified A4 order 12.")
    return orbit


def _same_projective_set(
    left: Sequence[np.ndarray], right: Sequence[np.ndarray], tol: float
) -> bool:
    if len(left) != len(right):
        return False
    unmatched = list(right)
    for vector in left:
        match = next(
            (
                index
                for index, candidate in enumerate(unmatched)
                if _same_projective_vector(vector, candidate, tol)
            ),
            None,
        )
        if match is None:
            return False
        unmatched.pop(match)
    return not unmatched


def a4_clique_orbits(
    cliques: Sequence[Sequence[Sequence[complex]]],
    parameter: complex,
    *,
    tol: float = 1e-8,
) -> list[dict[str, object]]:
    """Partition supplied H_D-gauge cliques into A4 orbits.

    An error is raised if a generator maps a supplied clique outside the input
    collection; this avoids reporting an incomplete orbit as complete.
    """
    generators = _numeric_generator_matrices(parameter)
    canonical = [
        [_phase_canonical(np.asarray(vector, dtype=np.complex128)) for vector in clique]
        for clique in cliques
    ]
    if any(len(clique) != 6 for clique in canonical):
        raise ValueError("Every third-basis clique must contain exactly six vectors.")

    def locate(image: Sequence[np.ndarray]) -> int | None:
        matches = [
            index
            for index, candidate in enumerate(canonical)
            if _same_projective_set(image, candidate, tol)
        ]
        if len(matches) > 1:
            raise RuntimeError(f"Ambiguous clique match for transformed basis: {matches}")
        return matches[0] if matches else None

    unassigned = set(range(len(canonical)))
    orbits: list[dict[str, object]] = []
    while unassigned:
        representative = min(unassigned)
        members = {representative}
        queue = deque([representative])
        while queue:
            current = queue.popleft()
            for generator in generators:
                image = [_phase_canonical(generator @ v) for v in canonical[current]]
                target = locate(image)
                if target is None:
                    raise RuntimeError(
                        f"A4 generator maps clique {current} outside the supplied "
                        "clique set; the input may be incomplete or use a different "
                        "Hadamard representative."
                    )
                if target not in members:
                    members.add(target)
                    queue.append(target)
        unassigned.difference_update(members)
        orbits.append(
            {
                "representative": representative,
                "members": sorted(members),
                "orbit_size": len(members),
            }
        )
    return orbits


def a4_vector_orbits(
    vectors: Sequence[Sequence[complex]],
    parameter: complex,
    *,
    tol: float = 1e-8,
) -> list[dict[str, object]]:
    """Partition a supplied H_D(z)-gauge vector set into discovered A4 orbits.

    Missing generator images raise an error rather than being silently
    interpreted as a smaller complete orbit. Stabilizers use the verified
    order-12 row-permutation group and are reported as candidate sizes.
    """
    generators = _numeric_generator_matrices(parameter)
    canonical = [
        _phase_canonical(np.asarray(vector, dtype=np.complex128)) for vector in vectors
    ]

    def locate(image: np.ndarray) -> int | None:
        matches = [
            index
            for index, candidate in enumerate(canonical)
            if _same_projective_vector(image, candidate, tol)
        ]
        if len(matches) > 1:
            raise RuntimeError(f"Ambiguous vector match for transformed root: {matches}")
        return matches[0] if matches else None

    unassigned = set(range(len(canonical)))
    orbits: list[dict[str, object]] = []
    group_order = len(row_permutation_group())
    while unassigned:
        representative = min(unassigned)
        members = {representative}
        queue = deque([representative])
        while queue:
            current = queue.popleft()
            for generator in generators:
                image = _phase_canonical(generator @ canonical[current])
                target = locate(image)
                if target is None:
                    raise RuntimeError(
                        f"A4 generator maps vector {current} outside the supplied "
                        "set; the roots may be incomplete or not form an invariant set."
                    )
                if target not in members:
                    members.add(target)
                    queue.append(target)
        unassigned.difference_update(members)
        orbit_size = len(members)
        if group_order % orbit_size:
            raise RuntimeError(
                f"Observed orbit size {orbit_size} does not divide group order {group_order}."
            )
        orbits.append(
            {
                "representative": representative,
                "members": sorted(members),
                "orbit_size": orbit_size,
                "stabilizer_order": group_order // orbit_size,
            }
        )
    return orbits


def _write_artifacts(output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    repo_root = Path(__file__).resolve().parents[2]

    def display_path(path: Path) -> str:
        resolved = path.resolve()
        try:
            return str(resolved.relative_to(repo_root))
        except ValueError:
            return str(resolved)

    equations = single_vector_equations()
    variables = [z, t, *x, *y]
    equation_text = output_dir / "I3_single_vector.txt"
    equation_text.write_text(
        "Exact single-vector MU ideal I3(z) for the Diita circle\n"
        "Variables: z,t,x1..x5,y1..y5; relations z*t-1 and xj*yj-1.\n\n"
        + "\n\n".join(
            f"F{index} = {sp.sstr(equation)}"
            for index, equation in enumerate(equations, start=1)
        )
        + "\n",
        encoding="utf-8",
    )

    python_text = [
        '"""Generated exact I3 equations; regenerate with scripts/python/dita_i3.py."""',
        "import sympy as sp",
        "I = sp.I",
        'z, t = sp.symbols("z t")',
        'x = sp.symbols("x1:6")',
        'y = sp.symbols("y1:6")',
        "x1, x2, x3, x4, x5 = x",
        "y1, y2, y3, y4, y5 = y",
        "relations = [z*t - 1, *(xj*yj - 1 for xj, yj in zip(x, y))]",
        "F = [",
        *[f"    {sp.sstr(equation)}," for equation in equations],
        "]",
        "",
    ]
    (output_dir / "I3_equations.py").write_text(
        "\n".join(python_text), encoding="utf-8"
    )

    stats = [
        "6 exact MU equations generated.",
        f"equations={len(equations)}",
        f"variables={len(variables)}",
        f"ideal_generators={len(equations) + len(single_vector_relations())}",
    ]
    for index, equation in enumerate(equations, start=1):
        polynomial = sp.Poly(equation, *variables)
        stats.append(
            f"F{index}: terms={len(polynomial.terms())}, "
            f"degree={polynomial.total_degree()}, "
            f"z={equation.has(z)}, t={equation.has(t)}"
        )
    (output_dir / "I3_stats.txt").write_text("\n".join(stats) + "\n", encoding="utf-8")

    symmetry_lines = [
        "Exact conjugation candidate: z<->t, x_j<->y_j, i->-i.",
    ]
    for index, equation in enumerate(equations, start=1):
        image = conjugate_expression(equation)
        matches = [
            target
            for target, candidate in enumerate(equations, start=1)
            if sp.simplify(image - candidate) == 0
        ]
        symmetry_lines.append(f"F{index} -> {matches}")
    symmetry_lines.extend(
        ["", "Exact parameter symmetry z->-z, x4<->x5, y4<->y5:"]
    )
    for index, equation in enumerate(equations, start=1):
        image = shift_pi_expression(equation)
        matches = [
            target
            for target, candidate in enumerate(equations, start=1)
            if sp.simplify(image - candidate) == 0
        ]
        symmetry_lines.append(f"F{index} -> {matches}")
    for index, (U, Q) in enumerate(fixed_z_symmetry_generators(), start=2):
        residual = (U * H_D * Q - H_D).applyfunc(reduce_parameter_relation)
        if residual != sp.zeros(6):
            raise ArithmeticError(f"Exact fixed-z generator G{index} failed.")
        symmetry_lines.append(f"G{index}: U*H_D*Q = H_D (exact modulo z*t-1)")
    symmetry_lines.append(f"row_permutation_group_order={len(row_permutation_group())}")
    (output_dir / "symmetry_test.txt").write_text(
        "\n".join(symmetry_lines) + "\n", encoding="utf-8"
    )
    (output_dir / "i3_branch_definition.md").write_text(
        "# Symbolic I3 branch ideal\n\n"
        "For six gauged vectors `v_a=(1,x[a,1],...,x[a,5])/sqrt(6)`, "
        "use one shared `z*t-1` relation, six copies of the six MU equations, "
        "and inverse relations `x[a,j]*y[a,j]-1`. For each `a<b`, add both "
        "conjugate orthogonality equations "
        "`1 + sum_j y[a,j]*x[b,j] = 0` and "
        "`1 + sum_j x[a,j]*y[b,j] = 0`. Thus the unreduced system has "
        f"{branch_ideal_counts()['variables']} variables and "
        f"{branch_ideal_counts()['generators']} generators, including "
        "30 pairwise orthogonality equations. The count shares the parameter "
        "relation instead of duplicating it six times.\n",
        encoding="utf-8",
    )
    manifest_inputs = [Path(__file__).resolve()]
    manifest_outputs = [
        output_dir / name
        for name in (
            "i3_branch_definition.md",
            "I3_equations.py",
            "I3_single_vector.txt",
            "I3_stats.txt",
            "symmetry_test.txt",
        )
    ]
    manifest = {
        "schema_version": "1.0",
        "runner": display_path(Path(__file__)).replace("\\", "/"),
        "inputs": [
            {
                "path": display_path(path).replace("\\", "/"),
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            }
            for path in manifest_inputs
        ],
        "outputs": [
            {
                "path": display_path(path).replace("\\", "/"),
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            }
            for path in manifest_outputs
        ],
        "python_version": sys.version.split()[0],
        "sympy_version": sp.__version__,
        "git_commit": _git_value("rev-parse", "HEAD"),
        "git_dirty": _git_dirty(),
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "random_seed": None,
        "completion_status": "completed",
        "remarks": "Deterministic exact symbolic generation; timestamp is run metadata.",
    }
    (output_dir / "run_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )


def _git_value(*args: str) -> str:
    try:
        return subprocess.run(
            ["git", "-C", str(Path(__file__).resolve().parents[2]), *args],
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return "unavailable"


def _git_dirty() -> bool | None:
    try:
        result = subprocess.run(
            [
                "git",
                "-C",
                str(Path(__file__).resolve().parents[2]),
                "status",
                "--porcelain",
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        return bool(result.stdout.strip())
    except (OSError, subprocess.CalledProcessError):
        return None


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path(__file__).resolve().parents[2] / "symbolic_export" / "I3",
        help="directory for generated equations and exact symmetry checks",
    )
    args = parser.parse_args()
    _write_artifacts(args.out_dir)
    print(f"Generated exact Diita I3 artifacts in {args.out_dir.resolve()}")


if __name__ == "__main__":
    main()
