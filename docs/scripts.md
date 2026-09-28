# Scripts directory: runnable entry points

The `scripts/` directory contains the main workflows and verification utilities.
Each script is a standalone entry point that:

- Reads from `research/` (parameterizations and protocols)
- Uses code from `src/` (reusable algorithms)
- Writes results to `results/` or `symbolic_export/`
- Records a run manifest for reproducibility

## Structure

```
scripts/
  julia/              # Julia runnable workflows (HomotopyContinuation-based)
  python/             # Python utility, audit, and exact symbolic scripts
  wolfram/            # Wolfram language verification (standalone)
  docker/             # Docker container runners (e.g., Macaulay2)
```

## Script categories (by role)

### Search scripts
- `search_special_loci.jl` — enumerate third-MUB candidates via polynomial degeneracy
- `enumerate_d0_cliques.jl` — exhaustive third-basis and fourth-vector search
- `i3_singular_locus.jl` — transition-pool counts and augmented I₃ Jacobian-null solve

### Validation & certification
- `validate_third_mub_candidates.jl` — per-basis HP verification
- `validate_perH_pipeline.jl` — end-to-end validation pipeline
- `certify_nwit1_all_cliques_*.jl` — W₁ witness certification (high-precision)

### Analysis & inspection
- `locus_phi_sweep_full_circle.jl` — locus geometric survey
- `lambda_periodicity_dita.jl` — parameter-space periodicity
- `locus_classification.jl` — orbit classification (Track D)
- `gauge_equivalence_analysis.jl` — gauge structure verification

### Export & symbolic computation
- `symbolic_elimination.jl` — generate M2 and Groebner systems
- `export_f6_pi10_pool_exact.jl` — exact pool exports
- `python/dita_i3.py` — generate and test the exact Diţă-circle I₃ ideal and symmetries

### Benchmarks & regression
- `run_benchmarks.jl` — regression suite (F6, D0, Tao)

### Audits & diagnostics
- `audit_clique_pipeline.jl` — planted-clique pipeline validation
- `csv_reconciliation.jl` — row-count audit
- `finish_project_audit.py` — final project state summary
- `verify_*_construction.jl` — specific finding verification

## How to run a script

All Julia scripts should:

1. Include `_paths.jl` to set `ROOT`, `RESULTS_DIR`, `SYMBOLIC_EXPORT_DIR`.
2. Include `src/MubSearch.jl` to load all modules.
3. Accept optional command-line arguments (e.g., `--quick`, `--resume`).
4. Call `write_result()` to save outputs with provenance.
5. Generate a `run_manifest.json` with the command, inputs, and outputs.

Example invocation (from project root):

```bash
julia --project=. scripts/julia/run_benchmarks.jl --quick
```

### Checking script status

Before running a long script, check its header comments:

```julia
# ============================================================================
# script_name.jl — Part X: [description]
#
# Objective: [what it does]
# Expected runtime: [estimate]
# Outputs: results/[files]
# Usage: julia --project=. scripts/julia/script_name.jl [options]
# ============================================================================
```

## Fast vs. long-running scripts

- **Fast** (< 1 min): `audit_*.jl`, `formalize_*.jl`, `csv_reconciliation.jl`
- **Medium** (5–30 min): `validate_*.jl`, locus sweeps, exports
- **Long** (30 min – hours): pool searches, W₁ certification, special loci enumeration

Use `MUB_LONG_TESTS=0` environment variable to skip expensive computations in
test scripts. See `docs/overview/reproduce.md` for the full regeneration table.

## Organizing new scripts

For new features or campaigns:

- Keep utility scripts grouped by function (e.g., all validation scripts together).
- Use the same `_paths.jl` and `MubSearch.jl` inclusion pattern.
- Document expected runtime and memory usage in the header.
- Accept command-line arguments for tuning (seed, tolerance, max iterations).
- Write results with `write_result()` and a `run_manifest.json`.

If a new feature has multiple entry points (e.g., `search_v1.jl`, `search_v2.jl`,
`validate_search_output.jl`), consider creating a `search/` subdirectory under
`scripts/julia/` to keep them organized.

The exact Diţă-circle I₃ implementation and its interface boundary with the
Karlsson pool solver are documented in
[`docs/research/dita_i3_workflow.md`](research/dita_i3_workflow.md).
Run the singular-locus campaign with
`julia --project=. scripts/julia/i3_singular_locus.jl [counts|singular|cliques|certify-saved]`,
then run [`analyze_i3_singular.py`](../scripts/python/analyze_i3_singular.py)
to report A₄ orbit sizes and stabilizers for generated roots or cliques.
Use its `prepare-certification` mode before `certify-saved` to certify saved
augmented roots without rerunning the augmented solve. Recorded results and
their completeness caveat are in
[`dita_i3_workflow.md`](research/dita_i3_workflow.md).

## Common debugging patterns

- Set `JULIA_DEBUG=MubSearch` to enable debug prints in modules.
- Save intermediate results to a named file in `results/` for inspection.
- Use small parameter sets for initial test runs (e.g., `--test10`).
- Check `results/*.log` and `results/*.err` for diagnostics.
