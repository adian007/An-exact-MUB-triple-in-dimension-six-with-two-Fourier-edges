"""Numerical equivalence fits for the order-six Fourier and 2-circulant families.

The Fourier convention follows Bengtsson et al., arXiv:quant-ph/0610161,
Eq. (4), as transcribed in ``src/brierley_weigert_notes.jl::fourier_F``.
The 2-circulant convention is Matszangosz-Szollosi, Eqs. (2)-(3).

These routines search all row/column permutation charts after dephasing.
They are floating-point fits, not exact classification algorithms.
"""

from __future__ import annotations

import itertools
import json
from dataclasses import asdict, dataclass
from typing import Any

import numpy as np
from scipy.optimize import least_squares


FIT_TOLERANCE = 1e-8
PERMUTATIONS = np.asarray(list(itertools.permutations(range(6))), dtype=np.intp)
ROOTS_OF_UNITY_3 = np.exp(2j * np.pi * np.arange(3) / 3)
Q6 = np.exp(2j * np.pi / 6)

FOURIER_A_POSITIONS = tuple(
    (row, column)
    for row in (1, 3, 5)
    for column in (1, 4)
)
FOURIER_B_POSITIONS = tuple(
    (row, column)
    for row in (1, 3, 5)
    for column in (2, 5)
)
FOURIER_MODIFIED = np.zeros((6, 6), dtype=bool)
for _row, _column in (*FOURIER_A_POSITIONS, *FOURIER_B_POSITIONS):
    FOURIER_MODIFIED[_row, _column] = True
FOURIER_FIXED = ~FOURIER_MODIFIED
FOURIER_A_ROWS = np.asarray([row for row, _ in FOURIER_A_POSITIONS])
FOURIER_A_COLS = np.asarray([column for _, column in FOURIER_A_POSITIONS])
FOURIER_B_ROWS = np.asarray([row for row, _ in FOURIER_B_POSITIONS])
FOURIER_B_COLS = np.asarray([column for _, column in FOURIER_B_POSITIONS])


@dataclass(frozen=True)
class FourierFit:
    family: str
    status: str
    fit_tolerance: float
    row_permutation: tuple[int, ...] | None
    column_permutation: tuple[int, ...] | None
    a: complex | None
    b: complex | None
    consistency_error: float | None
    fixed_entry_max_residual: float | None
    total_max_residual: float | None
    total_frobenius_residual: float | None
    charts_checked: int
    admissible_fixed_charts: int


@dataclass(frozen=True)
class TwoCirculantFit:
    family: str
    status: str
    fit_tolerance: float
    row_permutation: tuple[int, ...] | None
    column_permutation: tuple[int, ...] | None
    beta: complex | None
    gamma: complex | None
    epsilon: complex | None
    phi: complex | None
    constraint_3_residual: float | None
    matrix_max_residual: float | None
    matrix_frobenius_residual: float | None
    objective_max_residual: float | None
    solver: str
    solver_starts: int
    solver_nfev: int
    charts_checked: int


def _phase(values: np.ndarray, axis: int = -1) -> np.ndarray:
    values = np.asarray(values, dtype=complex)
    magnitudes = np.abs(values)
    safe = np.where(magnitudes > 0, values / np.maximum(magnitudes, 1e-300), 1)
    mean = np.mean(safe, axis=axis)
    mean_magnitude = np.abs(mean)
    fallback = np.take(safe, 0, axis=axis)
    return np.where(
        mean_magnitude > 1e-14,
        mean / np.maximum(mean_magnitude, 1e-300),
        fallback,
    )


def _unit_phase(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=complex)
    magnitude = np.abs(values)
    return np.where(
        magnitude > 0,
        values / np.maximum(magnitude, 1e-300),
        np.ones_like(values, dtype=complex),
    )


def fourier_family(a: complex = 1, b: complex = 1) -> np.ndarray:
    """Return the normalized F6^(2)(a,b) in the cited Eq. (4) convention.

    With q=exp(2*pi*i/6), a multiplies (0-based) rows 1,3,5 x columns 1,4;
    b multiplies rows 1,3,5 x columns 2,5. Both parameters are unimodular.
    """
    q = Q6
    z = complex(a)
    w = complex(b)
    return np.asarray(
        [
            [1, 1, 1, 1, 1, 1],
            [1, q * z, q**2 * w, q**3, q**4 * z, q**5 * w],
            [1, q**2, q**4, 1, q**2, q**4],
            [q**0, q**3 * z, w, q**3, z, q**3 * w],
            [1, q**4, q**2, 1, q**4, q**2],
            [1, q**5 * z, q**4 * w, q**3, q**2 * z, q * w],
        ],
        dtype=complex,
    )


def two_circulant_constraint(
    beta: complex, gamma: complex, epsilon: complex, phi: complex
) -> complex:
    return (
        beta * gamma * epsilon**2
        + beta * gamma * phi
        + beta**2 * epsilon * phi
        + gamma * epsilon * phi
        + beta * gamma**2 * epsilon * phi
        + beta * gamma * epsilon * phi**2
    )


def two_circulant_family(
    beta: complex, gamma: complex, epsilon: complex, phi: complex
) -> np.ndarray:
    """Return the normalized X6 matrix in Matszangosz-Szollosi Eq. (2)."""
    return np.asarray(
        [
            [1, 1, 1, 1, 1, 1],
            [
                1,
                -1,
                -1 / (gamma * epsilon),
                -1 / (beta * phi),
                1 / (gamma * epsilon),
                1 / (beta * phi),
            ],
            [
                1,
                -epsilon / beta,
                -1,
                -epsilon / (gamma * phi),
                epsilon / (gamma * phi),
                epsilon / beta,
            ],
            [
                1,
                -phi / gamma,
                -phi / (beta * epsilon),
                -1,
                phi / gamma,
                phi / (beta * epsilon),
            ],
            [
                1,
                epsilon / beta,
                phi / (beta * epsilon),
                1 / (beta * phi),
                1 / (beta * gamma),
                gamma / beta**2,
            ],
            [
                1,
                phi / gamma,
                1 / (gamma * epsilon),
                epsilon / (gamma * phi),
                beta / gamma**2,
                1 / (beta * gamma),
            ],
        ],
        dtype=complex,
    )


def dephase_chart(matrix: np.ndarray, rows: tuple[int, ...], columns: tuple[int, ...]) -> np.ndarray:
    """Permute and dephase at the top-left entry via cross-ratios."""
    chart = matrix[np.asarray(rows)][:, np.asarray(columns)]
    pivot = chart[0, 0]
    return chart * pivot / (chart[:, :1] * chart[:1, :])


def _validate_matrix(matrix: np.ndarray) -> np.ndarray:
    array = np.asarray(matrix, dtype=complex)
    if array.shape != (6, 6):
        raise ValueError(f"Expected a 6x6 matrix, received {array.shape}")
    if not np.all(np.isfinite(array)):
        raise ValueError("Matrix contains non-finite entries")
    if np.any(np.abs(array) < 1e-14):
        raise ValueError("Cannot dephase a matrix with zero entries")
    return array


def _fit_fourier_family_with_positions(
    matrix: np.ndarray,
    a_rows: np.ndarray,
    a_cols: np.ndarray,
    b_rows: np.ndarray,
    b_cols: np.ndarray,
    tolerance: float = FIT_TOLERANCE,
) -> tuple[dict[str, Any], int]:
    target = _validate_matrix(matrix)
    base = fourier_family(1, 1)
    best: dict[str, Any] | None = None
    admissible = 0

    for row_perm in PERMUTATIONS:
        row_matrix = target[row_perm]
        charts = np.transpose(row_matrix[:, PERMUTATIONS], (1, 0, 2))
        pivot = charts[:, 0, 0]
        dephased = charts * pivot[:, None, None] / (
            charts[:, :, :1] * charts[:, :1, :]
        )

        a_ratios = dephased[:, a_rows, a_cols] / base[a_rows, a_cols]
        b_ratios = dephased[:, b_rows, b_cols] / base[b_rows, b_cols]
        a_values = _phase(a_ratios, axis=1)
        b_values = _phase(b_ratios, axis=1)
        a_spread = np.max(np.abs(a_ratios - a_values[:, None]), axis=1)
        b_spread = np.max(np.abs(b_ratios - b_values[:, None]), axis=1)
        consistency = np.maximum.reduce(
            (
                a_spread,
                b_spread,
                np.abs(np.abs(a_values) - 1),
                np.abs(np.abs(b_values) - 1),
            )
        )

        predicted = np.broadcast_to(base, dephased.shape).copy()
        predicted[:, a_rows, a_cols] *= a_values[:, None]
        predicted[:, b_rows, b_cols] *= b_values[:, None]
        difference = dephased - predicted
        total_max = np.max(np.abs(difference), axis=(1, 2))
        fixed_max = np.max(np.abs(difference[:, FOURIER_FIXED]), axis=1)
        score = np.maximum.reduce((total_max, fixed_max, consistency))
        admissible += int(np.count_nonzero(fixed_max <= tolerance))

        local = int(np.argmin(score))
        if best is None or score[local] < best["score"]:
            best = {
                "score": float(score[local]),
                "row_permutation": tuple(int(value) for value in row_perm),
                "column_permutation": tuple(int(value) for value in PERMUTATIONS[local]),
                "a": complex(a_values[local]),
                "b": complex(b_values[local]),
                "consistency_error": float(consistency[local]),
                "fixed_entry_max_residual": float(fixed_max[local]),
                "total_max_residual": float(total_max[local]),
                "total_frobenius_residual": float(
                    np.linalg.norm(difference[local])
                ),
            }

    if best is None:
        raise RuntimeError("No Fourier fit chart was evaluated")
    return best, admissible


def fit_fourier_family(
    matrix: np.ndarray,
    tolerance: float = FIT_TOLERANCE,
) -> FourierFit:
    """Exhaust all 720x720 permutation charts for F6^(2)(a,b).

    This tests only the direct parameter placement defined by
    ``FOURIER_A_POSITIONS`` and ``FOURIER_B_POSITIONS``. Test membership in
    the transposed family by calling this function with ``matrix.T``.
    """
    best, admissible = _fit_fourier_family_with_positions(
        matrix,
        FOURIER_A_ROWS,
        FOURIER_A_COLS,
        FOURIER_B_ROWS,
        FOURIER_B_COLS,
        tolerance,
    )
    matched = (
        best["consistency_error"] <= tolerance
        and best["fixed_entry_max_residual"] <= tolerance
        and best["total_max_residual"] <= tolerance
    )
    return FourierFit(
        family="F6^(2)(a,b)",
        status="FIT" if matched else "NONE",
        fit_tolerance=tolerance,
        row_permutation=best["row_permutation"],
        column_permutation=best["column_permutation"],
        a=best["a"],
        b=best["b"],
        consistency_error=best["consistency_error"],
        fixed_entry_max_residual=best["fixed_entry_max_residual"],
        total_max_residual=best["total_max_residual"],
        total_frobenius_residual=best["total_frobenius_residual"],
        charts_checked=720 * 720,
        admissible_fixed_charts=admissible,
    )


def _x6_batched(
    beta: np.ndarray, gamma: np.ndarray, epsilon: np.ndarray, phi: np.ndarray
) -> np.ndarray:
    count = len(beta)
    result = np.ones((count, 6, 6), dtype=complex)
    result[:, 1, 1] = -1
    result[:, 1, 2] = -1 / (gamma * epsilon)
    result[:, 1, 3] = -1 / (beta * phi)
    result[:, 1, 4] = 1 / (gamma * epsilon)
    result[:, 1, 5] = 1 / (beta * phi)
    result[:, 2, 1] = -epsilon / beta
    result[:, 2, 2] = -1
    result[:, 2, 3] = -epsilon / (gamma * phi)
    result[:, 2, 4] = epsilon / (gamma * phi)
    result[:, 2, 5] = epsilon / beta
    result[:, 3, 1] = -phi / gamma
    result[:, 3, 2] = -phi / (beta * epsilon)
    result[:, 3, 3] = -1
    result[:, 3, 4] = phi / gamma
    result[:, 3, 5] = phi / (beta * epsilon)
    result[:, 4, 1] = epsilon / beta
    result[:, 4, 2] = phi / (beta * epsilon)
    result[:, 4, 3] = 1 / (beta * phi)
    result[:, 4, 4] = 1 / (beta * gamma)
    result[:, 4, 5] = gamma / beta**2
    result[:, 5, 1] = phi / gamma
    result[:, 5, 2] = 1 / (gamma * epsilon)
    result[:, 5, 3] = epsilon / (gamma * phi)
    result[:, 5, 4] = beta / gamma**2
    result[:, 5, 5] = 1 / (beta * gamma)
    return result


def _constraint_batched(
    beta: np.ndarray, gamma: np.ndarray, epsilon: np.ndarray, phi: np.ndarray
) -> np.ndarray:
    return (
        beta * gamma * epsilon**2
        + beta * gamma * phi
        + beta**2 * epsilon * phi
        + gamma * epsilon * phi
        + beta * gamma**2 * epsilon * phi
        + beta * gamma * epsilon * phi**2
    )


def _x6_solver_residual(angles: np.ndarray, target: np.ndarray) -> np.ndarray:
    beta, gamma, epsilon, phi = np.exp(1j * angles)
    matrix_residual = two_circulant_family(beta, gamma, epsilon, phi) - target
    constraint_residual = two_circulant_constraint(beta, gamma, epsilon, phi)
    combined = np.concatenate(
        (matrix_residual.ravel(), np.asarray([constraint_residual]))
    )
    return np.concatenate((combined.real, combined.imag))


def fit_two_circulant(
    matrix: np.ndarray,
    tolerance: float = FIT_TOLERANCE,
    solver_starts_per_chart: int = 3,
) -> TwoCirculantFit:
    """Exhaust permutation charts, infer phase starts, and polish by least squares.

    For each chart, beta is obtained from the three cube-root branches of
    (beta*gamma)/(gamma/beta^2), while gamma, epsilon and phi follow from
    entries in Eq. (2). The best chart is then refined from all three branches
    with scipy.optimize.least_squares against the full matrix and Eq. (3).
    """
    target = _validate_matrix(matrix)
    best: dict[str, Any] | None = None

    for row_perm in PERMUTATIONS:
        row_matrix = target[row_perm]
        charts = np.transpose(row_matrix[:, PERMUTATIONS], (1, 0, 2))
        pivot = charts[:, 0, 0]
        dephased = charts * pivot[:, None, None] / (
            charts[:, :, :1] * charts[:, :1, :]
        )

        gamma_epsilon = _unit_phase(-1 / dephased[:, 1, 2])
        beta_phi = _unit_phase(-1 / dephased[:, 1, 3])
        beta_gamma = _unit_phase(1 / dephased[:, 4, 4])
        gamma_over_beta2 = _unit_phase(dephased[:, 4, 5])
        beta_cubed = _unit_phase(beta_gamma / gamma_over_beta2)
        beta0 = np.exp(1j * np.angle(beta_cubed) / 3)
        betas = beta0[:, None] * ROOTS_OF_UNITY_3[None, :]
        gammas = beta_gamma[:, None] / betas
        epsilons = gamma_epsilon[:, None] / gammas
        phis = beta_phi[:, None] / betas
        flat_beta = betas.reshape(-1)
        flat_gamma = gammas.reshape(-1)
        flat_epsilon = epsilons.reshape(-1)
        flat_phi = phis.reshape(-1)
        repeated = np.repeat(dephased, 3, axis=0)

        models = _x6_batched(flat_beta, flat_gamma, flat_epsilon, flat_phi)
        matrix_max = np.max(np.abs(models - repeated), axis=(1, 2))
        constraints = np.abs(
            _constraint_batched(flat_beta, flat_gamma, flat_epsilon, flat_phi)
        )
        scores = np.maximum(matrix_max, constraints)
        local = int(np.argmin(scores))
        if best is None or scores[local] < best["score"]:
            chart_index = local // 3
            root_index = local % 3
            best = {
                "score": float(scores[local]),
                "row_permutation": tuple(int(value) for value in row_perm),
                "column_permutation": tuple(
                    int(value) for value in PERMUTATIONS[chart_index]
                ),
                "dephased_matrix": dephased[chart_index],
                "initial_phases": np.asarray(
                    [
                        flat_beta[local],
                        flat_gamma[local],
                        flat_epsilon[local],
                        flat_phi[local],
                    ],
                    dtype=complex,
                ),
                "root_index": root_index,
            }

    assert best is not None
    starts = []
    initial = best["initial_phases"]
    for root_offset in range(solver_starts_per_chart):
        phases = initial.copy()
        phases[0] *= ROOTS_OF_UNITY_3[root_offset % 3]
        phases[1] = (1 / best["dephased_matrix"][4, 4]) / phases[0]
        phases[2] = (
            -1 / best["dephased_matrix"][1, 2]
        ) / phases[1]
        phases[3] = (
            -1 / best["dephased_matrix"][1, 3]
        ) / phases[0]
        starts.append(np.angle(phases))

    solved = []
    total_nfev = 0
    for start in starts:
        result = least_squares(
            _x6_solver_residual,
            start,
            args=(best["dephased_matrix"],),
            method="trf",
            ftol=1e-12,
            xtol=1e-12,
            gtol=1e-12,
            max_nfev=250,
        )
        total_nfev += int(result.nfev)
        parameters = np.exp(1j * result.x)
        beta, gamma, epsilon, phi = parameters
        fitted = two_circulant_family(beta, gamma, epsilon, phi)
        matrix_max = float(np.max(np.abs(fitted - best["dephased_matrix"])))
        matrix_frob = float(np.linalg.norm(fitted - best["dephased_matrix"]))
        constraint = float(
            abs(two_circulant_constraint(beta, gamma, epsilon, phi))
        )
        solved.append(
            {
                "score": max(matrix_max, constraint),
                "beta": complex(beta),
                "gamma": complex(gamma),
                "epsilon": complex(epsilon),
                "phi": complex(phi),
                "constraint": constraint,
                "matrix_max": matrix_max,
                "matrix_frobenius": matrix_frob,
            }
        )
    polished = min(solved, key=lambda item: item["score"])
    matched = (
        polished["matrix_max"] <= tolerance
        and polished["constraint"] <= tolerance
    )
    return TwoCirculantFit(
        family="2-circulant X6",
        status="FIT" if matched else "NONE",
        fit_tolerance=tolerance,
        row_permutation=best["row_permutation"],
        column_permutation=best["column_permutation"],
        beta=polished["beta"],
        gamma=polished["gamma"],
        epsilon=polished["epsilon"],
        phi=polished["phi"],
        constraint_3_residual=polished["constraint"],
        matrix_max_residual=polished["matrix_max"],
        matrix_frobenius_residual=polished["matrix_frobenius"],
        objective_max_residual=polished["score"],
        solver="scipy.optimize.least_squares(method='trf', ftol=xtol=gtol=1e-12, max_nfev=250)",
        solver_starts=len(starts),
        solver_nfev=total_nfev,
        charts_checked=720 * 720,
    )


def fit_fourier_both_orientations(
    matrix: np.ndarray, tolerance: float = FIT_TOLERANCE
) -> tuple[FourierFit, FourierFit]:
    return (
        fit_fourier_family(matrix, tolerance),
        fit_fourier_family(np.asarray(matrix).T, tolerance),
    )


def _jsonable(value: Any) -> Any:
    if isinstance(value, complex):
        return [float(value.real), float(value.imag)]
    if isinstance(value, tuple):
        return list(value)
    if isinstance(value, dict):
        return {key: _jsonable(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_jsonable(item) for item in value]
    return value


def fit_as_dict(result: FourierFit | TwoCirculantFit) -> dict[str, Any]:
    return _jsonable(asdict(result))


def source_fourier_placement() -> dict[str, Any]:
    return {
        "source": "Bengtsson et al., arXiv:quant-ph/0610161, Eq. (4), repository transcription src/brierley_weigert_notes.jl::fourier_F",
        "parameter_a_zero_based": [list(position) for position in FOURIER_A_POSITIONS],
        "parameter_b_zero_based": [list(position) for position in FOURIER_B_POSITIONS],
        "q": [float(Q6.real), float(Q6.imag)],
        "source_convention": "q=exp(2*pi*i/6); a=exp(2*pi*i*x1); b=exp(2*pi*i*x2)",
    }


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


def build_stage_b_controls() -> list[dict[str, Any]]:
    """Create reproducible B4 control cases without writing any result files."""
    rng = np.random.default_rng(29092026)
    cases: list[dict[str, Any]] = []
    fourier_parameters = [
        ("F6_base", 1 + 0j, 1 + 0j),
        ("F6_random_ab_1", np.exp(0.37j), np.exp(-1.11j)),
        ("F6_random_ab_2", np.exp(-2.2j), np.exp(2.41j)),
        ("F6_random_ab_3", np.exp(1.72j), np.exp(0.66j)),
    ]
    for name, a, b in fourier_parameters:
        cases.append({"name": name, "matrix": fourier_family(a, b), "expected": "F"})
    cases.append(
        {
            "name": "F6_random_ab_1_transpose",
            "matrix": fourier_family(*fourier_parameters[1][1:]).T,
            "expected": "F^T",
        }
    )

    omega = np.exp(2j * np.pi / 3)
    x_controls = [
        ("X6_omega", (1 + 0j, omega, 1 + 0j, omega)),
        ("X6_omega_squared", (1 + 0j, omega**2, 1 + 0j, omega**2)),
    ]
    for index, (name, parameters) in enumerate(x_controls):
        matrix = two_circulant_family(*parameters)
        cases.append({"name": name, "matrix": matrix, "expected": "X6"})
        cases.append(
            {
                "name": f"{name}_equivalent_{index}",
                "matrix": random_equivalent(matrix, 9100 + index),
                "expected": "X6",
            }
        )
    random_parameters = {
        "beta": np.exp(1j * rng.uniform(-np.pi, np.pi)),
        "gamma": np.exp(1j * rng.uniform(-np.pi, np.pi)),
        "epsilon": 1 + 0j,
    }
    # Solve Eq. (3) for phi on the unit circle using dense scanning and
    # one-dimensional least squares, so the control is a valid X6 member.
    scan = np.exp(2j * np.pi * np.arange(4096) / 4096)
    constraints = np.asarray(
        [
            abs(
                two_circulant_constraint(
                    random_parameters["beta"],
                    random_parameters["gamma"],
                    random_parameters["epsilon"],
                    phi,
                )
            )
            for phi in scan
        ]
    )
    phi = scan[int(np.argmin(constraints))]
    if constraints.min() > 1e-8:
        random_parameters["gamma"] = omega
        random_parameters["beta"] = 1 + 0j
        phi = omega
    random_x = (
        random_parameters["beta"],
        random_parameters["gamma"],
        random_parameters["epsilon"],
        phi,
    )
    if abs(two_circulant_constraint(*random_x)) <= 1e-8:
        random_x_matrix = two_circulant_family(*random_x)
        cases.append({"name": "X6_random_phases", "matrix": random_x_matrix, "expected": "X6"})
        cases.append(
            {
                "name": "X6_random_phases_rephased",
                "matrix": random_equivalent(random_x_matrix, 9102),
                "expected": "X6",
            }
        )

    for parameter in (0.137, 0.25, 0.5):
        z = np.exp(2j * np.pi * parameter)
        zb = np.conjugate(z)
        i = 1j
        dita = np.asarray(
            [
                [1, 1, 1, 1, 1, 1],
                [1, -1, i, -i, -i, i],
                [1, i, -1, i * z, -i * z, -i],
                [1, -i, i * zb, -1, i, -i * zb],
                [1, -i, -i * zb, i, -1, i * zb],
                [1, i, -i, -i * z, i * z, -1],
            ],
            dtype=complex,
        )
        cases.append({"name": f"D6_x_{parameter}", "matrix": dita, "expected": "X6"})

    omega = np.exp(2j * np.pi / 3)
    tao = np.asarray(
        [
            [1, 1, 1, 1, 1, 1],
            [1, 1, omega, omega, omega**2, omega**2],
            [1, omega, 1, omega**2, omega**2, omega],
            [1, omega, omega**2, 1, omega, omega**2],
            [1, omega**2, omega**2, omega, 1, omega],
            [1, omega**2, omega, omega**2, omega, 1],
        ],
        dtype=complex,
    )
    cases.append({"name": "Tao_S6", "matrix": tao, "expected": "none"})
    cases.append(
        {
            "name": "Karlsson_off_Dita_1",
            "matrix": _build_karlsson(0.31, 0.73, 0.41),
            "expected": "none",
        }
    )
    cases.append(
        {
            "name": "Karlsson_off_Dita_2",
            "matrix": _build_karlsson(0.82, 1.19, 2.07),
            "expected": "none",
        }
    )
    return cases


def _build_karlsson(theta: float, phi: float, lam: float) -> np.ndarray:
    import sys
    from pathlib import Path

    script_dir = str(Path(__file__).resolve().parent)
    if script_dir not in sys.path:
        sys.path.insert(0, script_dir)
    from karlsson_k6_3 import build_k6

    return build_k6(theta, phi, lam)


def run_stage_b_controls(tolerance: float = FIT_TOLERANCE) -> dict[str, Any]:
    report: list[dict[str, Any]] = []
    for case in build_stage_b_controls():
        matrix = case["matrix"]
        fourier = fit_fourier_family(matrix, tolerance)
        transpose_fourier = fit_fourier_family(matrix.T, tolerance)
        two_circulant = fit_two_circulant(matrix, tolerance)
        flatness = float(np.max(np.abs(np.abs(matrix) - 1)))
        unitarity = float(np.max(np.abs(matrix.conj().T @ matrix - 6 * np.eye(6))))
        report.append(
            {
                "name": case["name"],
                "expected": case["expected"],
                "flatness_error": flatness,
                "unnormalized_unitarity_error": unitarity,
                "F_fit": fit_as_dict(fourier),
                "F_transpose_fit": fit_as_dict(transpose_fourier),
                "X6_fit": fit_as_dict(two_circulant),
                "positive_control_pass": (
                    case["expected"] == "none"
                    or (
                        (fourier.status == "FIT" if case["expected"] == "F" else False)
                        or (
                            transpose_fourier.status == "FIT"
                            if case["expected"] == "F^T"
                            else False
                        )
                        or (two_circulant.status == "FIT" if case["expected"] == "X6" else False)
                    )
                ),
            }
        )
    return {
        "tier": "NUMERICAL",
        "fit_tolerance": tolerance,
        "fourier_definition": source_fourier_placement(),
        "fourier_matrix_F6_2_a_b": fourier_family(1, 1).tolist(),
        "two_circulant_source": "Matszangosz-Szollosi, Theorem 1 paper, Eqs. (2)-(3)",
        "cases": report,
    }


def main() -> None:
    print(json.dumps(run_stage_b_controls(), indent=2))


if __name__ == "__main__":
    main()
