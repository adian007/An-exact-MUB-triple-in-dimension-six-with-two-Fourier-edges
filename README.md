# MUB-in-dimension-6

Research code and computational evidence for mutually unbiased bases in
dimension six, with project scope documented in
[docs/SCIENTIFIC_STATUS.md](docs/SCIENTIFIC_STATUS.md).

## Start here

- [Repository organization](docs/ORGANIZATION.md) — directory roles and the
  research-to-results workflow.
- [Scripts guide](docs/SCRIPTS.md) — overview and categories of entry points.
- [Documentation index](docs/README.md) — status, reproduction, claims, and
  research documentation.
- [Reproduction guide](docs/overview/REPRODUCE.md) — environment setup,
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

The project uses Julia dependencies declared in `Project.toml` and locked in
`Manifest.toml`. See the reproduction guide before running expensive searches.
