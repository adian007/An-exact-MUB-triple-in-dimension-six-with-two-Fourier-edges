# K3 Fourier-structure diagnostic

## Retraction

The earlier interpretation in this file is retracted. The legacy runner
searched for three \(-1\)-like entries in a single dephased column and labelled
those detections as Fourier-family or transposed-Fourier evidence. That is not
the theorem predicate and does not establish either family membership.

The correct result is Theorem 1 of Matszangosz–Szöllősi:
<https://link.springer.com/article/10.1007/s10623-024-01503-w>. For a
normalized order-six complex Hadamard matrix it tests three distinct columns,
each containing at least one \(-1\), and concludes transposed Fourier **or**
2-circulant family membership. The row form is inferred by applying the
theorem to the transpose.

The old files
[`diagnostic.json`](../../results/campaigns/k3_fourier_structure/diagnostic.json)
and
[`exact_chart_equations.json`](../../results/campaigns/k3_fourier_structure/exact_chart_equations.json)
are preserved byte-for-byte as legacy diagnostic data. Their labels
`Fourier_family_column_test`, `transposed_Fourier_family_row_test`,
`family_test`, and chart-ID suffixes `F`/`FT` are incorrect for the data they
encode. See the adjacent
[`RETRACTION_NOTICE.md`](../../results/campaigns/k3_fourier_structure/RETRACTION_NOTICE.md).

Corrected theorem-predicate outputs are regenerated only to dated result
paths. Numerical detections remain sample evidence; they do not prove exact
branch membership, pool completeness, or component containment.
