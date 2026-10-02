# INTERFACES.md — MUB 4th-Basis Search Contract

**Version**: 2026-09-29  
**Purpose**: Coordinate parallel subagents for the dimension-6 MUB problem.

## Overview

The work partitions by phase. **All phases can run in parallel** since they read from existing data and write to different output files.

## Directory Ownership (Write)

| Phase | Directory | Output Files | Read-only Inputs |
|-------|-----------|--------------|------------------|
| 0 | `scripts/python/` | `ref_ref_resolve.py`, `phase0_manifest.json` | All CSV/JSON/npz files |
| 1A | `scripts/python/` | `pool_exact_pi3.py` | `b3_a4_i4_results/*.npz` |
| 1B | `scripts/python/` | `pool_graph_pi3.py` | `b3_a4_i4_results/pool_pi_over_3.npz` |
| 2 | `scripts/python/` | `w1_pi3_groebner.py`, `phase2_result.md` | `b3_a4_i4_results/pool_pi_over_3.npz`, `b3_a4_i4_results/summary.json` |
| 3 | `scripts/python/` | `fourier_bridge.py` | `docs/results/k3_family_membership_2026-09-29.md` |
| 4 | `scripts/python/` | `full_circle_strat.py` | `provenance/fold_a4_final.zip` |
| 5 | `scripts/python/` | `defect_search.py`, `qcqp_check.py`, `sos_check.py` | Already exists (defect_minimization.py) |

**NO OVERLAP**: Each agent writes to unique files. No agent modifies another's directory.

## Interfaces

### Pool Data Format (`b3_a4_i4_results/pool_pi_over_3.npz`)
- Key `V`: `np.ndarray` shape `(72, 6)`, dtype `complex128`
- Rows are vectors, columns are components
- Vectors are already normalized: `|v|=1`

### Summary JSON Format (`b3_a4_i4_results/summary.json`)
```json
{
  "pool_size": 72,
  "cliques": [[2, 24, 35, 47, 66, 67], ...],
  "generator_clique_permutations": [[2, 3, 0, 1], ...]
}
```

### Output Spec: Phase 2 (W1 GB)
File `phase2_result.md` must contain:
1. H matrix at λ=π/3 (exact, in Q(ζ₂₄))
2. B3 reconstruction (exact or high-precision)
3. GB computation: `{1}` or not, degree, time
4. Cofactor certificate: `1 = Σ c_i f_i` with explicit c_i
5. Independent verification result

### Output Spec: Phase 1B (Pool Graph)
File `phase1b_pool_graph_report.md` must contain:
1. Pool statistics: 72 vectors, each 6-component complex
2. Graph: edges by ORTHOGONAL (0), UNBIASED (1/6), NEITHER
3. Verification of 4 known cliques (6-cliques, mutually unbiased)
4. Automorphism group: size, generators, monomial unitary check
5. SAT/ILP result: UNSAT proves no quartet

## Invariants

1. **Pool completeness**: A "no quartet" result only proves no 4th MUB for the RESTORED pool, not all possible B3.
2. **Reference point**: All paths relative to `D:/MUBs in 6-dimension/`
3. **Python**: Use `C:/Users/adian/AppData/Local/Programs/Python/Python313/python.exe`

## Definition of Done

For each phase, run the agent's script and verify:
1. Output file exists and is non-empty
2. Status field indicates success (no uncaught exceptions)
3. Cross-check at least one numerical result against existing data