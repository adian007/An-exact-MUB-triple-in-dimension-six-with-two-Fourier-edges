"""Export exact B3 incidence and three-minus-one chart equations.

The output is an exact polynomial-system input for each numerical witness
chart. It is not an elimination result: symbolic branch containment still has
to be proved.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import sympy as sp

import dita_i3


def build_b3_system(chart: dict[str, object], transpose_view: bool) -> dict[str, object]:
    x = sp.symbols("x1_1:7 x2_1:7 x3_1:7 x4_1:7 x5_1:7")
    y = sp.symbols("y1_1:7 y2_1:7 y3_1:7 y4_1:7 y5_1:7")
    b = sp.ones(6, 6)
    bbar = sp.ones(6, 6)
    for row in range(1, 6):
        for col in range(6):
            index = (row - 1) * 6 + col
            b[row, col] = x[index]
            bbar[row, col] = y[index]

    hbar = dita_i3.H_D.xreplace(
        {dita_i3.z: dita_i3.t, dita_i3.t: dita_i3.z, dita_i3.I: -dita_i3.I}
    )
    equations: list[sp.Expr] = [dita_i3.z * dita_i3.t - 1]
    equations.extend(
        b[row, col] * bbar[row, col] - 1
        for row in range(1, 6)
        for col in range(6)
    )

    for col in range(6):
        vector = b[:, col]
        vector_bar = bbar[:, col]
        for hcol in range(6):
            inner = (vector_bar.T * dita_i3.H_D[:, hcol])[0]
            conjugate_inner = (vector.T * hbar[:, hcol])[0]
            equations.append(sp.expand(inner * conjugate_inner - 6))

    for left in range(6):
        for right in range(left + 1, 6):
            equations.append(
                sp.expand(sum(bbar[row, left] * b[row, right] for row in range(6)))
            )
            equations.append(
                sp.expand(sum(bbar[row, right] * b[row, left] for row in range(6)))
            )

    transition = bbar.T * dita_i3.H_D
    transpose = ".T" if transpose_view else ""

    def entry(row: int, col: int) -> sp.Expr:
        return transition[col, row] if transpose_view else transition[row, col]

    transition_view = [
        [str(sp.expand(entry(row, col))) for col in range(6)]
        for row in range(6)
    ]
    pivot_row = int(chart["pivot_row"])
    pivot_column = int(chart["pivot_column"])
    target_column = int(chart["target_column"])
    incidence_equations = []
    for row in chart["minus_one_rows"]:
        row = int(row)
        incidence_equations.append(
            sp.expand(
                entry(row, target_column) * entry(pivot_row, pivot_column)
                + entry(row, pivot_column) * entry(pivot_row, target_column)
            )
        )

    equations.extend(incidence_equations)
    variables = [str(dita_i3.z), str(dita_i3.t), *(str(v) for v in x), *(str(v) for v in y)]
    return {
        "variables": variables,
        "transition_definition": f"T = Bbar^T * H_D; chart view T{transpose}",
        "transition_matrix": transition_view,
        "equation_counts": {
            "parameter_and_torus": 31,
            "mutual_unbiasedness": 36,
            "basis_orthogonality": 30,
            "three_minus_one_chart": 3,
            "total": len(equations),
        },
        "equations": [str(sp.expand(equation)) for equation in equations],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("results/campaigns/k3_fourier_structure/diagnostic.json"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/campaigns/k3_fourier_structure/exact_chart_equations.json"),
    )
    args = parser.parse_args()

    diagnostic = json.loads(args.input.read_text())
    charts = []
    for run in diagnostic["runs"]:
        for clique_index, clique in enumerate(run["cliques"]):
            target = clique["transitions"]["B3_dagger_H_D"]
            for family_test, is_transpose in (
                ("Fourier_family_column_test", False),
                ("transposed_Fourier_family_row_test", True),
            ):
                result = target[family_test]
                if not result["witness_found"]:
                    continue
                witness = result["witness"]
                chart_id = (
                    f"lambda_{run['lambda']:.12f}_clique_{clique_index}_"
                    f"{'FT' if is_transpose else 'F'}"
                )
                charts.append(
                    {
                        "chart_id": chart_id,
                        "lambda_sample": run["lambda"],
                        "clique_indices": clique["clique"],
                        "family_test": family_test,
                        "witness": witness,
                        "system": build_b3_system(witness, is_transpose),
                    }
                )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(
            {
                "campaign_id": "k3_fourier_structure",
                "status": "exact_chart_equations_not_eliminated",
                "source_diagnostic": str(args.input),
                "chart_count": len(charts),
                "charts": charts,
            },
            indent=2,
        )
    )
    print(f"Wrote {len(charts)} exact witness-chart systems to {args.output}")


if __name__ == "__main__":
    main()
