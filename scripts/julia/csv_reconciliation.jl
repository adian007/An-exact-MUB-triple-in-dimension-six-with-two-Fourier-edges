# Item 3: Reconcile special_loci_search.csv row count vs intended queue.
include(joinpath(@__DIR__, "_paths.jl"))
using CSV, DataFrames, Printf

const OUT = joinpath(RESULTS_DIR, "csv_reconciliation.txt")
const DEGEN_CSV = joinpath(RESULTS_DIR, "degeneracy_candidates.csv")
const SEARCH_CSV = joinpath(RESULTS_DIR, "special_loci_search.csv")

function parse_degen_cap()
    cap = 80
    args = copy(ARGS)
    if "--degen-cap" in args
        i = findfirst(==("--degen-cap"), args)
        i !== nothing && length(args) >= i + 1 && (cap = parse(Int, args[i + 1]))
    end
    haskey(ENV, "DEGEN_CAP") && (cap = parse(Int, ENV["DEGEN_CAP"]))
    return cap
end

function count_degen_raw()
    !isfile(DEGEN_CSV) && return 0, 0
    lines = readlines(DEGEN_CSV)
    n_total = max(0, length(lines) - 1)
    seen = Set{Tuple{Float64,Float64,Float64}}()
    for line in lines[2:end]
        isempty(strip(line)) && continue
        parts = split(line, ',')
        length(parts) < 3 && continue
        key = (round(parse(Float64, parts[1]), digits = 15),
               round(parse(Float64, parts[2]), digits = 15),
               round(parse(Float64, parts[3]), digits = 15))
        push!(seen, key)
    end
    return n_total, length(seen)
end

function expected_primary_queue(; degen_cap = 80)
    dita_theta = acos(1 / sqrt(3))
    n_anchors = 20  # F6_cyclic, F6_theta0, Dita, Mobius_audit_fail, 5 theta_*, 12 lambda_kpi6
    n_degen_raw, n_degen_dedup = count_degen_raw()
    n_degen_loaded = min(n_degen_dedup, degen_cap)
    return (
        n_anchors = n_anchors,
        n_degen_raw = n_degen_raw,
        n_degen_dedup = n_degen_dedup,
        n_degen_cap = degen_cap,
        n_degen_loaded = n_degen_loaded,
        n_primary = n_anchors + n_degen_loaded,
    )
end

const CSV_HEADER = [
    "name", "theta", "phi", "lambda", "status", "z4_dev",
    "n_pool", "n_raw", "n_conjugate", "n_verified",
    "max_clique", "found_third", "found_fourth",
    "found_fourth_1e6", "found_fourth_1e8", "found_fourth_1e10",
    "n_tracked", "n_certified", "n_distinct", "n_distinct_real",
    "pool_complete_flag", "n_dedup_certified",
]

"""Manual UTF-8 CSV parse (robust to malformed rows during concurrent --resume)."""
function load_search_rows(path)
    lines = filter(!isempty, strip.(readlines(path)))
    isempty(lines) && return NamedTuple[], 0
    header = split(replace(lines[1], "\ufeff" => ""), ',')
    ncols = length(CSV_HEADER)
    rows = NamedTuple[]
    bad = 0
    for (i, line) in enumerate(lines[2:end])
        parts = split(line, ',')
        if length(parts) != ncols
            bad += 1
            continue
        end
        g(k) = parts[findfirst(==(k), header)]
        try
            push!(rows, (
                name = g("name"),
                theta = parse(Float64, g("theta")),
                phi = parse(Float64, g("phi")),
                lambda = parse(Float64, g("lambda")),
                status = g("status"),
                max_clique = isempty(g("max_clique")) ? 0 : parse(Int, g("max_clique")),
                pool_complete_flag = g("pool_complete_flag") == "true",
                n_certified = isempty(g("n_certified")) ? missing : parse(Int, g("n_certified")),
                n_tracked = isempty(g("n_tracked")) ? missing : parse(Int, g("n_tracked")),
            ))
        catch
            bad += 1
        end
    end
    return rows, bad
end

function _safe_startswith(s, prefix)
    s isa AbstractString || return false
    return startswith(s, prefix)
end

"""Simulate one-level refinement count from primary queue (same rules as search_special_loci.jl)."""
function simulate_refinements(primary_pts; quick = false)
    quick && return 0, 0
    n_ref_total = 0
    n_ref_unique = 0
    seen = Set{Tuple{Float64,Float64,Float64}}()
    for (name, th, ph, lam) in primary_pts
        key = (round(th, digits = 10), round(ph, digits = 10), round(lam, digits = 10))
        push!(seen, key)
    end
    # We need actual evaluation to know which spawn refinements; use CSV cache instead
    return n_ref_total, n_ref_unique
end

function analyze_csv()
    degen_cap = parse_degen_cap()
    lines = String[]
    push!(lines, "=== Item 3: CSV reconciliation (588 vs actual) ===")
    push!(lines, "Date: $(Dates.now())")
    push!(lines, @sprintf("  degen_cap=%d (see docs/SEARCH_SETS.md for S588 vs S*)", degen_cap))
    push!(lines, "")

    exp = expected_primary_queue(; degen_cap = degen_cap)
    push!(lines, "--- Intended primary queue ---")
    push!(lines, @sprintf("  anchors: %d", exp.n_anchors))
    push!(lines, @sprintf("  degeneracy_candidates.csv rows (raw): %d", exp.n_degen_raw))
    push!(lines, @sprintf("  degeneracy unique (theta,phi,lambda): %d", exp.n_degen_dedup))
    push!(lines, @sprintf("  degeneracy loaded (cap %d, sorted -z4_dev): %d", exp.n_degen_cap, exp.n_degen_loaded))
    push!(lines, @sprintf("  primary total: %d", exp.n_primary))
    push!(lines, "")

    !isfile(SEARCH_CSV) && (push!(lines, "ERROR: $(SEARCH_CSV) missing"); return lines)
    rows, n_bad = load_search_rows(SEARCH_CSV)
    n_csv = length(rows)
    push!(lines, "--- Actual CSV ---")
    push!(lines, @sprintf("  data rows (valid): %d", n_csv))
    n_bad > 0 && push!(lines, @sprintf("  malformed rows skipped: %d (concurrent --resume artifact)", n_bad))
    push!(lines, @sprintf("  expected (verify_csv): 588"))
    push!(lines, @sprintf("  gap vs 588: %d (negative = more rows than target)", 588 - n_csv))
    push!(lines, @sprintf("  historical baseline (phase2): 471 rows, gap=117"))
    push!(lines, "")

    # Row categories
    n_anchor = count(r -> !_safe_startswith(r.name, "degen_") && !_safe_startswith(r.name, "ref_"), rows)
    n_degen = count(r -> _safe_startswith(r.name, "degen_"), rows)
    n_ref = count(r -> _safe_startswith(r.name, "ref_"), rows)
    n_ref1 = count(r -> _safe_startswith(r.name, "ref_") && !occursin("ref_ref", r.name), rows)
    n_ref2 = count(r -> occursin("ref_ref", r.name), rows)
    push!(lines, "--- Row categories ---")
    push!(lines, @sprintf("  anchor/non-ref names: %d", n_anchor))
    push!(lines, @sprintf("  degen_* rows: %d", n_degen))
    push!(lines, @sprintf("  ref_* (level 1): %d", n_ref1))
    push!(lines, @sprintf("  ref_ref_* (level 2+): %d", n_ref2))
    push!(lines, @sprintf("  total ref rows: %d", n_ref))
    push!(lines, "")

    # Dedup keys
    keys = Set{Tuple{Float64,Float64,Float64}}()
    for r in rows
        push!(keys, (round(r.theta, digits = 10), round(r.phi, digits = 10), round(r.lambda, digits = 10)))
    end
    push!(lines, @sprintf("  unique coordinate keys: %d (duplicates: %d)", length(keys), n_csv - length(keys)))
    push!(lines, "")

    # Status breakdown
    push!(lines, "--- Status breakdown ---")
    for st in sort(unique(r.status for r in rows))
        push!(lines, @sprintf("  %s: %d", st, count(r -> r.status == st, rows)))
    end
    push!(lines, "")

    # Refinement spawners: max_clique>=6 or pool incomplete
    spawn = filter(r -> r.status == "ok" && (
        r.max_clique >= 6 ||
        r.pool_complete_flag == false ||
        (r.n_certified !== missing && r.n_tracked !== missing && r.n_certified != r.n_tracked)
    ), rows)
    push!(lines, @sprintf("  rows eligible for refinement spawn: %d", length(spawn)))
    push!(lines, "")

    # Estimate expected refinements from spawn-eligible primaries
    refine_per = 5^3  # n=5, scale=0.005 for clique>=6
    refine_per_lo = 3^3  # n=3, scale=0.01 otherwise
    est_ref_max = 0
    est_ref_typ = 0
    for r in spawn
        _safe_startswith(r.name, "ref_ref") && continue
        if r.max_clique >= 6
            est_ref_max += refine_per
            est_ref_typ += refine_per
        else
            est_ref_max += refine_per_lo
            est_ref_typ += refine_per_lo
        end
    end
    push!(lines, "--- Refinement estimate (theoretical max, no dedup) ---")
    push!(lines, @sprintf("  from %d spawn-eligible rows: up to %d ref points (125 each if clique>=6)", length(spawn), est_ref_max))
    push!(lines, @sprintf("  actual ref rows in CSV: %d", n_ref))
    push!(lines, "")

    # Gap explanation
    gap = 588 - n_csv
    push!(lines, "=== GAP EXPLANATION (588 - $n_csv = $gap) ===")
    reasons = String[]
    if exp.n_degen_dedup > exp.n_degen_cap
        push!(reasons, "degen cap: only $(exp.n_degen_cap)/$(exp.n_degen_dedup) unique degen candidates loaded (−$(exp.n_degen_dedup - exp.n_degen_cap) never queued)")
    end
    missing_degen = exp.n_degen_loaded - n_degen
    missing_degen > 0 && push!(reasons, @sprintf("degen rows missing from CSV: %d loaded but only %d present", exp.n_degen_loaded, n_degen))
    missing_ref = est_ref_typ - n_ref
    if missing_ref > 0
        push!(reasons, @sprintf("refinement shortfall: ~%d ref points not written (crash/incomplete run/dedup skips)", missing_ref))
    end
    n_primary_csv = n_anchor + n_degen
    if n_primary_csv < exp.n_primary
        push!(reasons, @sprintf("primary shortfall: %d/%d primary points in CSV", n_primary_csv, exp.n_primary))
    end
    n_bad > 0 && push!(reasons, @sprintf("malformed rows during concurrent write: %d skipped", n_bad))
    push!(reasons, "second-level ref_ref blocked by can_spawn_refinement (only primary + ref_* level-1 spawn)")
    push!(reasons, "coordinate dedup: refinements overlapping anchors/degen skipped")
    for r in reasons
        push!(lines, "  - $r")
    end
    push!(lines, "")

    # Log crash evidence
    log_path = joinpath(RESULTS_DIR, "special_loci_log.txt")
    if isfile(log_path)
        bytes = read(log_path)
        utf16 = length(bytes) >= 2 && bytes[1] == 0xFF && bytes[2] == 0xFE
        log_txt = utf16 ? String(transcode(String, bytes)) : String(bytes)
        crash_markers = [
            occursin("RemoteException", log_txt),
            occursin("killed", log_txt),
            occursin("Error", log_txt),
            occursin("No such file", log_txt),
        ]
        push!(lines, "--- Log evidence ---")
        push!(lines, @sprintf("  special_loci_log.txt size: %d bytes (utf16=%s)", length(bytes), utf16))
        push!(lines, @sprintf("  contains errors/exceptions: %s", any(crash_markers)))
    end

    open(OUT, "w") do io
        for ln in lines
            println(io, ln)
        end
    end
    for ln in lines
        println(ln)
    end
    println("\nWrote $OUT")
end

using Dates
analyze_csv()
