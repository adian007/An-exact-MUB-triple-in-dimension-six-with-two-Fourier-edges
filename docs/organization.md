# Repository organization

The repository contains both reusable code and a long-running research record.
Keep those roles separate: code defines and executes methods; research inputs
define what is being studied; results preserve what a run produced; docs explain
what the evidence supports.

## Directory responsibilities

| Directory | Contents and intended use |
|---|---|
| `src/` | Reusable Julia modules and mathematical implementations. |
| `scripts/` | Runnable analysis, validation, certification, and export entry points, grouped by language. |
| `test/` | Automated regression and structural tests for `src/`. |
| `research/` | Structured or machine-readable research inputs: catalogues, parameterizations, claim data, and protocols. |
| `docs/research/` | Human-readable research notes, literature assessments, and reports. |
| `docs/overview/` | Project status, claims, known results, and reproduction guidance. |
| `docs/results/` | Curated, human-readable summaries of computational findings. |
| `docs/paper/` | Paper drafts, proofs, and publication source. |
| `results/` | Generated numerical outputs, machine-readable results, and execution logs. |
| `symbolic_export/` | Generated exact systems consumed by symbolic algebra tools such as Macaulay2. |
| `runs/` | Snapshots of a particular run when preserving its working inputs and outputs together is useful. |
| `pdf/` | Literature PDFs and their associated source or extraction utilities. |

`research/` and `docs/research/` are deliberately different: put data or
protocols that scripts consume in the former; put explanations written for
researchers in the latter. Do not keep duplicate copies of the same report in
both locations.

## How work moves through the repository

```
Code and protocols (reusable)
├─ src/             (algorithms)
├─ test/            (regression tests)
└─ research/        (structured inputs)
        ↓
Scripts (entry points)
├─ scripts/julia/
├─ scripts/python/
└─ scripts/wolfram/
        ↓
Outputs (raw results and logs)
├─ results/         (JSON, CSV, TXT)
└─ symbolic_export/ (M2, etc.)
        ↓
Interpretation and curation
├─ docs/research/   (analysis)
├─ docs/results/    (summaries)
├─ docs/overview/   (status & claims)
└─ docs/paper/      (proofs & drafts)
```

- **Code** (src/) should be independent of any particular run; tests verify its
  correctness.
- **Protocols and inputs** (research/) guide what a run will do, and are
  machine-readable (CSV, JSON, TOML).
- **Scripts** (scripts/) are entry points that accept command-line arguments,
  call code, consume inputs, and write results.
- **Results** (results/, symbolic_export/) are generated and contain both
  structured outputs and diagnostic logs.
- **Documentation** (docs/) explains what the results mean, not how to run them
  (that goes in the reproduction guide).

## Organizing new computational campaigns

The repository has accumulated legacy flat files in `scripts/julia/` and
`results/`. Preserve those paths unless a coordinated migration updates every
reference. For new or substantially revised campaigns:

### Layout

```
research/
  campaigns/
    <campaign-id>/
      protocol.md         (description of the search or validation)
      parameters.json     (machine-readable inputs)
      config.toml         (runtime settings: solver, tolerances)

results/
  campaigns/
    <campaign-id>/
      <run-id>/
        run_manifest.json (record of how this run was produced)
        metrics.csv       (numerical outputs)
        summary.md        (human interpretation, if needed)
        log.txt           (solver diagnostic and stderr)
```

- Put the protocol, parameter definitions, and machine-readable inputs under
  `research/campaigns/<campaign-id>/`.
- Put outputs under `results/campaigns/<campaign-id>/<run-id>/`, where
  `<run-id>` is timestamped (UTC) or otherwise unique.
- Keep reusable algorithms in `src/`; keep campaign orchestration and
  command-line parsing in `scripts/`.
- Group newly added scripts by role (`search/`, `validate/`, `certify/`, or
  `export/`) only when a group has multiple entry points. Existing script paths
  remain supported.
- Promote only interpreted, claim-relevant conclusions into `docs/results/`
  or `docs/overview/`. Raw logs and interim tables are not themselves a
  scientific conclusion.

### Run manifests

Each significant run should have a `run_manifest.json` next to its outputs.
Record the following:

```json
{
  "schema_version": "1.0",
  "runner": "scripts/julia/my_workflow.jl",
  "arguments": ["--param", "value"],
  "inputs": [
    "scripts/julia/my_workflow.jl",
    "Project.toml",
    "Manifest.toml",
    "research/campaigns/campaign-x/parameters.json"
  ],
  "inputs_hashes": {
    "scripts/julia/my_workflow.jl": "sha256:abcd...",
    "research/campaigns/campaign-x/parameters.json": "sha256:1234..."
  },
  "outputs": [
    "results/campaigns/campaign-x/run-123/metrics.csv"
  ],
  "git_commit": "abc1234...",
  "git_dirty": false,
  "julia_version": "1.12.6",
  "timestamp_utc": "2026-09-28 20:30:00",
  "random_seed": 42,
  "completion_status": "completed",
  "remarks": "optional notes"
}
```

A manifest supports reproducibility; it does not by itself validate the
mathematics or replace a human-readable interpretation. See the benchmark
workflow for an example in `results/benchmarks/run_manifest.json` and
`src/Provenance.jl` for the standard provenance record function.

## Guidelines

1. **Do not move existing results.** Legacy code often hardcodes paths
   like `results/verify_dita_construction.txt`. Keep those paths stable;
   add new campaigns under the `campaigns/` subdirectory instead.

2. **Distinguish content types.**
   - `research/`: machine-readable (CSV, JSON, TOML)
   - `docs/research/`: human-written (Markdown)
   - `results/`: generated outputs (any format)
   - `docs/results/` and `docs/overview/`: curated interpretations for publication

3. **Keep data and narrative separate.** If a script produces both a data table
   and an interpretation, store the table in `results/` and the interpretation
   (if significant) in `docs/results/` with a cross-reference.

4. **Use the same random seed and parameters for reproducible reruns.** Store
   them in the run manifest and preserve them in `research/campaigns/` if
   others may need to repeat the run.

5. **Don't accumulate duplicate evidence.** If a finding is reported in
   `docs/research/`, don't copy it to `docs/overview/` unless it supports
   a specific claim or status update. Link instead.
