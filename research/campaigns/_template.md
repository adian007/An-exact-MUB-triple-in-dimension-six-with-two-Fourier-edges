# Campaign template: structured research inputs

This directory holds the definition of a computational campaign: what it searches
for, how it does it, and what parameters it uses. Outputs live separately under
`results/campaigns/<campaign-id>/`.

## Files in this directory

- **protocol.md** — Human-readable description of the research question, method,
  expected outcomes, and success criteria.
- **parameters.json** — Machine-readable search space, parameter ranges, or
  point definitions that scripts will parse and use.
- **config.toml** — Solver settings, numerical tolerances, computational limits
  (timeout, memory), and other tuning parameters.

## Example: protocol.md

```markdown
# Campaign: [name]

## Objective
Describe what we are computing and why.

## Method
How the script will approach the problem (e.g., pool enumeration, gradient
search, exhaustive certification).

## Parameters
- Dimension: 6
- Parameter space: [description]
- Expected pool size: [estimate]

## Tolerances
- Hadamard orthogonality error: 1e-10
- Solver confidence: 1e-8

## Success criteria
When the campaign will be considered complete (e.g., "all 588 candidates
certified" or "locus verified to 1000 sample points").

## Remarks
Any special considerations (e.g., memory, runtime estimates, dependencies on
other campaigns).
```

## Example: parameters.json

```json
{
  "campaign_id": "f6_arc_theta0",
  "timestamp_created": "2026-09-28T20:30:00Z",
  "search_space": {
    "theta": 0.0,
    "phi": [0.0, 0.5, 1.0],
    "lambda": {
      "type": "linspace",
      "start": -1.0,
      "end": 1.0,
      "count": 100
    }
  },
  "solver": {
    "max_iterations": 1000,
    "hp_bits": 256
  },
  "expected_outcomes": {
    "third_mub_bases": "~10",
    "fourth_mub_vectors": "varies"
  }
}
```

## Example: config.toml

```toml
[solver]
tolerance = 1e-8
hp_bits = 256
max_iterations = 1000

[limits]
timeout_seconds = 3600
memory_mb = 8192

[dedup]
projective_tolerance = 1e-8
clustering_bits = 128

[output]
write_interval = 10
log_level = "info"
```

## Running a campaign

A campaign is typically invoked from a script entry point, which:

1. Loads this protocol and parameters.
2. Instantiates the search or validation.
3. Writes results to `results/campaigns/<campaign-id>/<run-id>/`.
4. Generates a `run_manifest.json` in the output directory.

Example from a hypothetical script:

```julia
include(joinpath(@__DIR__, "_paths.jl"))

protocol_dir = joinpath(ROOT, "research", "campaigns", "campaign_name")
params = JSON.parsefile(joinpath(protocol_dir, "parameters.json"))
config = TOML.parsefile(joinpath(protocol_dir, "config.toml"))

# ... run the campaign ...

results_dir = joinpath(RESULTS_DIR, "campaigns", params["campaign_id"], run_id)
mkpath(results_dir)
write_result(joinpath(results_dir, "metrics.json"), results; provenance = prov)
```

## Preserving campaigns for reproducibility

- **Never edit** parameters.json or config.toml after a run completes; create a
  new version (`parameters_v2.json`) if you need to change them.
- **Record the version** in the `run_manifest.json` so someone can recreate the
  exact run.
- **Link campaigns to code versions** via the `git_commit` field in the manifest.
  If a run used a specific branch or tag, note it.
