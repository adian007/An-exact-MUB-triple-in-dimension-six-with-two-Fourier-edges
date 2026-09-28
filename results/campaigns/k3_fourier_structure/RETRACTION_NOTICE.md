# Retraction notice: legacy Fourier-family labels

Date: 2026-09-28

The labels `Fourier_family_column_test`, `transposed_Fourier_family_row_test`,
`family_test`, and the chart-ID suffixes `F`/`FT` in the legacy
`diagnostic.json` and `exact_chart_equations.json` are incorrect for the data
they contain.

Those files encode only the one-column three-row diagnostic: three
\(-1\)-like entries in a single dephased column (or row after transposition).
This must not be read as Fourier-family, transposed-Fourier-family, or
2-circulant-family membership. The `criterion_applicable` flag must always be
read alongside `witness_found`; a detection flag alone is not a verdict when
the input did not pass the Hadamard check.

The correct published criterion is Theorem 1 of Matszangosz–Szöllősi,
“A characterization of complex Hadamard matrices appearing in families of
MUB triplets,” Designs, Codes and Cryptography 92 (2024), Theorem 1:
<https://link.springer.com/article/10.1007/s10623-024-01503-w>. For a
normalized order-six complex Hadamard matrix, three distinct columns must
each contain at least one \(-1\). The theorem concludes that the matrix is in
the transposed Fourier family or the 2-circulant family. The analogous row
statement follows by applying the theorem to the transpose.

The legacy JSON files remain unchanged as provenance. Corrected regeneration
is written only to new dated paths; the corrected outputs distinguish
numerical candidate detections from exact or certified conclusions.
