# File Guide

## 01_symbolic_I3

- readme.md: original package-level notes.
- I3_single_vector.txt: generated exact single-vector equations.
- I3_equations.py: symbolic construction of the single-vector and branch ideals.
- i3_branch_definition.md: mathematical definition of the six-vector branch ideal.
- I3_stats.txt: generator/variable/term statistics.
- build_i3.py: builder used to generate the symbolic artifacts.
- symmetry_test.txt: exact symmetry test output.
- fixed_symmetry.txt: numerical fixed-z symmetry search output.
- fixed_symmetry_exact.txt: exact symbolic verification of the discovered symmetries.
- search_fixed_sym.py: exhaustive numerical search over row permutations for fixed-z monomial symmetries.
- search_zminus_sym.py: z -> -z symmetry test/search.
- extract_fixed_sym.py: helper for extracting/verifying the fixed-z symmetry data.

## 02_symmetry

Reserved for the next A4 orbit-decomposition scripts and symmetry data. No fabricated results are placed here.

## 03_dita_exact

Reserved for the exact Diţă family representation, parameter map, row/column permutations and rephasing certificate. These formulas are documented in docs/01_research_context.md and should be copied into the main project's exact-family module.

## 04_fold_analysis

Reserved for the numerical fold/singularity computation at lambda*. Existing numerical results are documented but not reconstructed here as source code unless the original script is available.

## 05_fourth_MUB

Reserved for the fourth-vector witness/incidence ideal and certified obstruction scripts. Do not treat this directory as containing a global fourth-MUB proof.

## 06_results

Use this directory for raw CSV/JSON/TXT output from reproducible runs. Keep generated results separate from source code.

## project_integration

Contains integration instructions and a proposed folder layout for merging this package into the main MUB repository.
