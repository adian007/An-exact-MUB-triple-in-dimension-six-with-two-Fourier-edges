# ============================================================================
# Provenance.jl — Part XV: standardized metadata for every serious result.
#
# Every exported artifact carries: git commit, worktree dirtiness, Julia
# version, platform, thread count, timestamp, random seed, parameter
# values, solver settings, tolerances, path counts, certification status.
# A result without provenance is not citable.
#
# JSON output uses a tiny built-in serializer (no new dependencies);
# NamedTuples and Dicts are supported.
# ============================================================================

function _git_info()
    out = Dict{String,Any}("git_commit" => "unavailable", "git_dirty" => missing)
    # Prefer the canonical worktree (env override for WSL copies), else cwd.
    dirs = String[]
    haskey(ENV, "MUB_GIT_DIR") && push!(dirs, ENV["MUB_GIT_DIR"])
    push!(dirs, pwd())
    for d in dirs
        try
            commit = strip(String(read(`git -C $d rev-parse HEAD`)))
            dirty = !isempty(strip(String(read(`git -C $d status --porcelain`))))
            out["git_commit"] = commit
            out["git_dirty"] = dirty
            out["git_worktree"] = d
            break
        catch
            continue
        end
    end
    if out["git_commit"] == "unavailable"
        out["git_commit"] = "unavailable (not a git worktree; running from a copy)"
    end
    return out
end

"""
    provenance_record(; seed=nothing, parameters=nothing, solver=Dict(),
                      tolerances=Dict(), extra=Dict{String,Any}())

Build the standard provenance dictionary. Attach it to every result file.
"""
function provenance_record(; seed = nothing, parameters = nothing,
                           solver = Dict{String,Any}(), tolerances = Dict{String,Any}(),
                           extra = Dict{String,Any}())
    info = Dict{String,Any}()
    merge!(info, _git_info())
    info["julia_version"] = string(VERSION)
    info["platform"] = string(Sys.MACHINE)
    info["os"] = Sys.islinux() ? "Linux" : Sys.iswindows() ? "Windows" :
                 Sys.isapple() ? "macOS" : "other"
    info["hostname"] = try
        strip(String(read(`hostname`)))
    catch
        "unknown"
    end
    info["julia_threads"] = Threads.nthreads()
    info["timestamp_utc"] = replace(string(Dates.now(Dates.UTC)), "T" => " ", "Z" => "")
    info["random_seed"] = seed
    info["parameters"] = parameters
    info["solver"] = solver
    info["tolerances"] = tolerances
    info["homotopycontinuation_version"] = string(pkgversion(HomotopyContinuation))
    merge!(info, extra)
    return info
end

# --- minimal JSON writer ----------------------------------------------------

_json_escape(s::String) = (s = replace(s, "\\" => "\\\\"); replace(s, "\"" => "\\\""))

_json(x) = x   # identity default for the recursion below

_json(x::Union{String,Symbol,AbstractString}) = "\"$(_json_escape(string(x)))\""
_json(x::Bool) = x ? "true" : "false"
_json(x::Nothing) = "null"
_json(x::Missing) = "null"
_json(x::Integer) = string(x)
_json(x::AbstractFloat) = isfinite(x) ? string(x) : "\"$x\""
_json(x::Complex) = "\"$x\""
_json(x::Real) = "\"$x\""   # BigFloat etc.
_json(x::Union{Tuple,AbstractVector}) = "[" * join(map(_json, collect(x)), ",") * "]"
_json(p::Pair) = _json(Dict{String,Any}(string(first(p)) => last(p)))
_json(d::NamedTuple) = _json(Dict{String,Any}(string(k) => v for (k, v) in pairs(d)))
_json(d::AbstractDict) =
    "{" * join(["\"$(_json_escape(string(k)))\":$(_json(v))" for (k, v) in d], ",") * "}"
# Generic struct (e.g. PoolAudit): serialize its fields.
_json(x) = _json(Dict{String,Any}(string(f) => getfield(x, f)
                                  for f in fieldnames(typeof(x))))

"""Serialize any nested combination of dicts / named tuples / vectors /
scalars / structs to JSON (defensive; used for result artifacts)."""
json(x) = _json(x)

"""
    write_result(path::String, data; provenance=nothing)

Write `data` as JSON with a `_provenance` block, creating parent dirs.
Every serious computational result in this project must go through this.
"""
function write_result(path::AbstractString, data; provenance = nothing)
    mkpath(dirname(abspath(path)))
    payload = Dict{String,Any}(string(k) => v for (k, v) in
                               (data isa AbstractDict ? data :
                                Dict("result" => data)))
    if provenance !== nothing
        payload["_provenance"] = provenance
    end
    open(path, "w") do io
        write(io, json(payload))
    end
    return path
end
