# Research provenance and hashing

## Why provenance matters

This problem combines symbolic formulas, numerical continuation, reconstructed bases, and literature transcriptions. A result can change because of a source edit, a parameter rounding, a solver version, or a tolerance. Hashes make those causes distinguishable.

## Hash manifest should include

- source files that construct Karlsson, Dita, Fourier, and Tao matrices;
- polynomial-system builders;
- certification scripts;
- exact parameter manifest;
- result JSON/CSV/TXT outputs;
- literature transcription notes;
- environment report.

## Recommended record

For each artifact store: relative path, SHA-256, timestamp UTC, git commit, Julia version, package manifest hash, seed, precision, and tolerance. Hash generated data separately from source code.

## Interpretation

A matching hash proves byte identity, not mathematical correctness. A changed hash is a trigger for rerunning the focused validation, not evidence that the old result was wrong.

The existing `research/hashes/phase0_source_hashes.txt` and benchmark provenance should be extended with the new Markdown reports and with the exact literature version URLs used in this dossier.
