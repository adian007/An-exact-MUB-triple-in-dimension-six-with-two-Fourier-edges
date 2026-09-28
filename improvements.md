# Improvements to repository organization and reproducibility

## Summary

This change clarifies the repository's directory structure, roles, and workflows without moving or deleting existing research artifacts. Key improvements:

### 1. **Clear directory semantics** (`docs/organization.md`)
   - Documented the role of each top-level directory.
   - Distinguished between `research/` (machine-readable inputs) and `docs/research/` (human-written notes).
   - Clarified that `results/` holds raw generated output, while `docs/results/` and `docs/overview/` hold curated interpretations.

### 2. **Run manifests for reproducibility**
   - Added `run_manifest.json` generation to `scripts/julia/run_benchmarks.jl`.
   - The manifest records:
     - Runner path and command-line arguments
     - Input files and their roles
     - Output paths and expected artifacts
     - Git commit, Julia version, timestamp, seed, and completion status
   - Other scripts can follow this pattern for future campaigns.

### 3. **Structured campaigns framework**
   - Created `research/campaigns/` directory with a template (`_template.md`).
   - Created corresponding `results/campaigns/` for outputs.
   - New campaigns can organize their inputs, parameters, and protocols clearly.
   - Manifests link inputs to outputs, supporting reproducibility.

### 4. **Updated documentation**
   - Expanded `readme.md` with brief directory overview.
   - Added `docs/scripts.md` to document script categories and roles.
   - Updated `docs/readme.md` to cross-reference the organization guide.
   - Updated reproduction guide to mention the new `run_manifest.json`.

## Non-changes (intentionally preserved)

- **No file moves:** Legacy scripts and results remain in their existing paths.
- **No deletions:** All existing outputs and evidence are preserved.
- **Backward compatibility:** The changes are additive and do not break existing workflows.

## How to use the new structure

### For new computational campaigns:

1. Create a directory under `research/campaigns/<campaign-id>/`.
2. Add `protocol.md`, `parameters.json`, and `config.toml` (see template).
3. Write or adapt a script entry point in `scripts/julia/`.
4. Have the script call `write_result()` with `run_manifest.json` as one of its outputs.
5. Results live under `results/campaigns/<campaign-id>/<run-id>/`.

### For understanding existing results:

- Consult `docs/overview/reproduce.md` for the one-command regeneration table.
- Check the provenance fields in JSON outputs (`_provenance` block).
- For legacy flat files in `results/`, refer to the script header comments.

## Files changed

- `readme.md` — expanded with directory overview and links.
- `docs/readme.md` — added reference to `organization.md`.
- `docs/organization.md` — **new** comprehensive guide.
- `docs/scripts.md` — **new** script categories and usage guide.
- `scripts/julia/run_benchmarks.jl` — added `run_manifest.json` generation.
- `docs/overview/reproduce.md` — updated to mention manifests.
- `research/campaigns/` — **new** directory with `_template.md`.
- `results/campaigns/` — **new** empty directory for future outputs.

## Next steps (optional, for future sessions)

- As new campaigns or searches are added, use the `research/campaigns/` structure.
- Gradually migrate long-running scripts to use `run_manifest.json` for consistency.
- Consider documenting the existing campaigns (e.g., special loci search, Dita
  periodicity) in `research/campaigns/<campaign-id>/protocol.md` for clarity.
