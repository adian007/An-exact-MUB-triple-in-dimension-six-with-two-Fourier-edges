# MUB-in-dimension-6

Research code and computational evidence for mutually unbiased bases in
dimension six, with project scope documented in
[docs/scientific_status.md](docs/scientific_status.md).

## Start here

- [Repository organization](docs/organization.md) — directory roles and the
  research-to-results workflow.
- [Scripts guide](docs/scripts.md) — overview and categories of entry points.
- [Documentation index](docs/readme.md) — status, reproduction, claims, and
  research documentation.
- [Reproduction guide](docs/overview/reproduce.md) — environment setup,
  tests, and commands for regenerating results.

## Main directories

| Directory | Role |
|---|---|
| `src/` | Reusable Julia implementations. |
| `scripts/` | Runnable Julia, Python, and Wolfram workflows. |
| `test/` | Regression tests. |
| `research/` | Structured research inputs and protocols. |
| `docs/research/` | Human-readable research notes and reports. |
| `results/` | Generated outputs and run logs. |
| `symbolic_export/` | Exact systems for symbolic algebra tools. |

## Supporting bundles

- [`mub_complete_computational_package/`](mub_complete_computational_package/readme.md) — self-contained computational research package.
- [`mub_fold_a4_computation/`](mub_fold_a4_computation/README.md) — fold and A4 computation bundle.
- [`mub_solution_attack_package_v2/`](mub_solution_attack_package_v2/README_CURRENT.md) — current research handoff and provenance bundle.
- [`b3_a4_i4_results/`](b3_a4_i4_results/README.md) — B3/A4 and reduced fourth-vector result bundle.
- `archives/`, `provenance/`, and `runs/` — dated snapshots, provenance records, and preserved run inputs/outputs.

See [the directory map](docs/directory_map.md) for the full root layout. These bundle paths are referenced by reports and campaign inputs, so keep them stable unless performing a coordinated path migration.

The project uses Julia dependencies declared in `Project.toml` and locked in
`Manifest.toml`. See the reproduction guide before running expensive searches.
# MUB-in-dimension-6
