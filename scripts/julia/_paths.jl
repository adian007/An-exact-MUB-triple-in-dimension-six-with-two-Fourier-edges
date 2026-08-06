# Shared project paths for scripts in scripts/julia/
const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const RESULTS_DIR = joinpath(ROOT, "results")
const SYMBOLIC_EXPORT_DIR = joinpath(ROOT, "symbolic_export")
include(joinpath(ROOT, "src", "mub_zauner_6d_liang_chen.jl"))
