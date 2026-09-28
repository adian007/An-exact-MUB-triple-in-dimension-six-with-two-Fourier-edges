# Directory organization at a glance

## Current state (visual reference)

```
MUBs in 6-dimension/
│
├─ 📋 Project files
│  ├─ README.md (START HERE)
│  ├─ Project.toml (Julia environment)
│  ├─ Manifest.toml (locked dependencies)
│  └─ IMPROVEMENTS.md (summary of changes)
│
├─ 📚 docs/ (human-readable documentation)
│  ├─ README.md (index)
│  ├─ ORGANIZATION.md (THIS GUIDES STRUCTURE)
│  ├─ SCRIPTS.md (script categories)
│  ├─ SCIENTIFIC_STATUS.md (project scope)
│  ├─ EXACT_SYSTEM.md (algebraic setup)
│  │
│  ├─ overview/
│  │  ├─ REPRODUCE.md (commands to regenerate)
│  │  ├─ claim_audit.md
│  │  ├─ known_results_table.md
│  │  └─ ... (status, findings, audit)
│  │
│  ├─ research/
│  │  ├─ definitions/ (mathematical definitions)
│  │  ├─ literature/ (literature review)
│  │  ├─ methods/ (method summaries)
│  │  └─ ... (human-written explanations)
│  │
│  ├─ paper/
│  │  ├─ main_theorems.tex
│  │  ├─ proofs/ (detailed proofs)
│  │  └─ preamble_common.tex
│  │
│  └─ results/ (summaries of findings)
│     └─ benchmarks.md
│
├─ 🔬 src/ (reusable Julia code)
│  ├─ MubSearch.jl (module entry point)
│  ├─ Karlsson.jl (construction)
│  ├─ Pool.jl (pool audit)
│  ├─ Certification.jl (W₁ witness)
│  ├─ Provenance.jl (reproducibility)
│  └─ ... (7 modules total)
│
├─ 🚀 scripts/ (runnable workflows)
│  ├─ julia/
│  │  ├─ run_benchmarks.jl (regression)
│  │  ├─ search_special_loci.jl (search)
│  │  ├─ validate_third_mub_candidates.jl (certification)
│  │  ├─ symbolic_elimination.jl (export)
│  │  ├─ certify_nwit1_all_cliques_*.jl (W₁)
│  │  └─ ... (70+ scripts)
│  │
│  ├─ python/
│  │  ├─ finish_project_audit.py
│  │  ├─ degeneracy_scan.py
│  │  └─ ... (20 utilities)
│  │
│  ├─ wolfram/ (4 verification scripts)
│  └─ docker/
│     └─ run_m2.ps1 (Macaulay2 runner)
│
├─ ✅ test/
│  ├─ runtests.jl (regression suite)
│  └─ norm_detection.jl
│
├─ 🔍 research/ (machine-readable inputs)
│  ├─ claims/ (claim tracking)
│  ├─ definitions/ (parameter definitions)
│  ├─ exact_systems/ (algebraic systems)
│  ├─ literature/ (literature data)
│  ├─ methods/ (method comparison)
│  ├─ parameter_spaces/ (parameterizations)
│  │
│  ├─ campaigns/ ⭐ (NEW: structured inputs for new work)
│  │  ├─ _TEMPLATE.md (how to define a campaign)
│  │  └─ <campaign-id>/
│  │     ├─ protocol.md
│  │     ├─ parameters.json
│  │     └─ config.toml
│  │
│  └─ ... (hashes, logs, reports)
│
├─ 📊 results/ (generated numerical outputs)
│  ├─ benchmarks/
│  │  ├─ benchmarks.json (provenance embedded)
│  │  ├─ run_manifest.json ⭐ (reproducibility record)
│  │  └─ benchmarks.md
│  │
│  ├─ campaigns/ ⭐ (NEW: organized by campaign)
│  │  └─ <campaign-id>/
│  │     └─ <run-id>/
│  │        ├─ run_manifest.json
│  │        ├─ metrics.csv
│  │        └─ ... (outputs)
│  │
│  └─ ... (100+ other logs, CSVs, TXTs)
│
├─ 🔤 symbolic_export/ (exact systems for M2)
│  ├─ pool_F6.m2
│  ├─ pool_Dita.m2
│  └─ ... (17 Macaulay2 system files)
│
├─ 📖 pdf/ (literature)
│  └─ ... (40+ research papers and scripts)
│
├─ 📚 runs/ (saved run snapshots)
│  └─ 2026-09-17T170232Z_phase1_definitions/
│
└─ .julia-depot/ (Julia packages, not tracked)
```

## What goes where: decision tree

```
Does it describe what we're computing or how?
├─ YES, human-readable (English/math)
│  ├─ Is it background, tutorial, or method explanation?
│  │  → docs/research/ (e.g., definitions, literature summary)
│  ├─ Is it a project status update or claim?
│  │  → docs/overview/ or docs/results/
│  └─ Is it publication content (proofs, theorems)?
│     → docs/paper/
│
├─ YES, machine-readable (JSON, CSV, TOML)
│  ├─ Does the script consume it as input?
│  │  → research/campaigns/<id>/ (or research/<category>/)
│  └─ Did a script generate it?
│     → results/ (or results/campaigns/<id>/<run>/)
│
└─ YES, is it reusable code or a test?
   ├─ Does it export functions for use by scripts?
   │  → src/
   └─ Does it test src/ or provide regression checks?
      → test/
```

## Key distinction: research/ vs. docs/research/

| Aspect | `research/` | `docs/research/` |
|---|---|---|
| **Format** | JSON, CSV, TOML, machine scripts | Markdown, human-written |
| **Purpose** | Inputs consumed by scripts; definitions | Explanations, analysis, literature review |
| **Consumption** | Script workflows read these | Humans read these for understanding |
| **Example** | `research/campaigns/f6_arc/parameters.json` | `docs/research/definitions/mub_dimension6_research.md` |
| **Audience** | Scripts and automation | Researchers and reviewers |

✋ **Don't duplicate.** If you write an explanation in `docs/research/`, don't copy it to
`research/` as a machine file. Link instead, or keep data and explanation in separate locations.

## How results move toward publication

```
1. Script executes                     (scripts/julia/)
   ↓
2. Writes raw output with provenance  (results/)
   ↓
3. Researcher analyzes & interprets   (docs/research/)
   ↓
4. Important findings summarized      (docs/results/ + docs/overview/)
   ↓
5. Incorporated into paper            (docs/paper/)
```

Each layer adds clarity and interpretation; raw data is never deleted.

## Running a new campaign

```
Step 1: Define inputs
  research/campaigns/<id>/protocol.md
  research/campaigns/<id>/parameters.json
  research/campaigns/<id>/config.toml

Step 2: Write runner script
  scripts/julia/run_<id>.jl
  (reads from research/campaigns/<id>/)
  (writes to results/campaigns/<id>/<timestamp>/)

Step 3: Generate manifest
  results/campaigns/<id>/<timestamp>/run_manifest.json
  (records what ran, on what inputs, producing what outputs)

Step 4: Analyze & interpret (optional)
  docs/research/... or docs/results/...
  (human-readable summary for publication)
```

See [docs/ORGANIZATION.md](docs/ORGANIZATION.md) for full details.
