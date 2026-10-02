# Fold critical-root match audit

Date: 2026-09-30

## Result

The archived augmented-root computation contains 30 numerical records, of
which 24 pass the recorded residual filters. The regenerated Julia critical
pool contains 72 normalized vectors.

Direct comparison in the common dephased phase convention found:

- 0/24 singular candidates within (10^{-7}) radians of a Julia-pool vector;
- 0/24 within (10^{-5}) radians;
- maximum nearest-vector phase distance: 1.3935 radians.

The candidate roots independently satisfy the augmented MU equations at about
(10^{-12}) in the archived data, while the Julia pool independently satisfies
the fixed-H MU equations at about (10^{-15}). Thus the failed match is not
caused by an obvious residual failure.

The archived A4 data still reports two 12-element numerical orbits with
maximum matching error (8.51\times10^{-10}). This audit does not certify
those roots or prove that they are missing physical solutions; it shows that
the 24-to-72 discrepancy is a genuine unresolved set-convention/root-set
issue, not an already demonstrated duplicate list.

Artifact: `results/phase2_fold_critical_match_audit.json`.
