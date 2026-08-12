# Targeted search near special Karlsson loci + certified root counts on ALL points.
# Does NOT repeat the coarse 512-point uniform sweep.
#
# Usage:
#   julia --project=. search_special_loci.jl [--quick] [--resume] [--no-refine]
#   julia --project=. search_special_loci.jl --degen-cap 200 [--no-refine]
#   julia --project=. search_special_loci.jl --verify-only
#   julia --project=. search_special_loci.jl --degen-cap 500 [--no-refine] [--resume]

include(joinpath(@__DIR__, "_paths.jl"))

using CSV
using DataFrames
using Printf

const MU_TOLS = (1e-6, 1e-8, 1e-10)
const DEGEN_CSV = joinpath(RESULTS_DIR, "degeneracy_candidates.csv")
const OUT_CSV = get(ENV, "SPECIAL_LOCI_CSV", joinpath(RESULTS_DIR, "special_loci_search.csv"))

const CSV_HEADER = [
    "name", "theta", "phi", "lambda", "status", "z4_dev",
    "n_pool", "n_raw", "n_conjugate", "n_verified",
    "max_clique", "found_third", "found_fourth",
    "found_fourth_1e6", "found_fourth_1e8", "found_fourth_1e10",
    "n_tracked", "n_certified", "n_distinct", "n_distinct_real",
    "pool_complete_flag", "n_dedup_certified",
]

"""Mobius z4 consistency deviation at (theta, phi, lambda)."""
function karlsson_z4_dev(theta, phi, lam)
    A = build_A(theta, phi)
    B = -[1 1; 1 -1] - A
    alpha_A, beta_A = A[1, 2]^2, A[1, 1]^2
    alpha_B, beta_B = B[1, 2]^2, B[1, 1]^2
    z1sq = exp(2im * lam)
    z3sq = mobius(z1sq, alpha_A, beta_A)
    z4sq = mobius(z1sq, alpha_B, beta_B)
    num = beta_B - z3sq * conj(alpha_B)
    den = alpha_B - z3sq * conj(beta_B)
    z2sq = num / den
    z4sq_check = mobius(z2sq, alpha_A, beta_A)
    return abs(z4sq - z4sq_check)
end

"""Points where Karlsson meets known special CHMs or parameter degeneracies."""
function special_anchor_points()
    dita_theta = acos(1 / sqrt(3))
    anchors = NamedTuple[]
    push!(anchors, (:name => "F6_cyclic", :theta => 0.0, :phi => 0.0, :lambda => 0.0))
    push!(anchors, (:name => "F6_theta0", :theta => 0.0, :phi => 0.5, :lambda => 0.3))
    push!(anchors, (:name => "Dita", :theta => dita_theta, :phi => pi / 4, :lambda => 0.4))
    push!(anchors, (:name => "Mobius_audit_fail",
                     :theta => 0.9553166181245092, :phi => 0.7853981633974483, :lambda => 0.4))
    for (name, th) in [("theta_0", 0.0), ("theta_pi2", pi / 2), ("theta_pi3", pi / 3),
                       ("theta_pi4", pi / 4), ("theta_pi6", pi / 6)]
        push!(anchors, (:name => name, :theta => th, :phi => pi / 4, :lambda => pi / 6))
    end
    for k in 0:11
        push!(anchors, (:name => "lambda_kpi6_k=$k", :theta => 0.3, :phi => 0.5, :lambda => k * pi / 6))
    end
    return anchors
end

const COORD_DEDUP_DIGITS = 15  # full Float64 coords; 10-digit keys miss sharp loci (Dita)

"""Load degeneracy candidates from Python scan (full Float64 coords, severity sort)."""
function load_degeneracy_candidates(; degen_cap = 80)
    !isfile(DEGEN_CSV) && return NamedTuple[]
    rows = NamedTuple[]
    open(DEGEN_CSV) do io
        header = split(strip(readline(io)), ',')
        for line in eachline(io)
            isempty(strip(line)) && continue
            parts = split(line, ',')
            length(parts) < 3 && continue
            row = Dict(zip(header, parts))
            z4 = haskey(row, "z4_dev") && !isempty(row["z4_dev"]) ? parse(Float64, row["z4_dev"]) : 0.0
            push!(rows, (
                name = "degen_" * replace(get(row, "flags", "cand"), "|" => "_"),
                theta = parse(Float64, row["theta"]),
                phi = parse(Float64, row["phi"]),
                lambda = parse(Float64, row["lambda"]),
                z4_dev = z4,
            ))
        end
    end
    sort!(rows, by = r -> -r.z4_dev)
    seen = Set{Tuple{Float64,Float64,Float64}}()
    out = NamedTuple[]
    for r in rows
        key = (round(r.theta, digits = COORD_DEDUP_DIGITS),
               round(r.phi, digits = COORD_DEDUP_DIGITS),
               round(r.lambda, digits = COORD_DEDUP_DIGITS))
        key in seen && continue
        push!(seen, key)
        push!(out, (;
            name = r.name,
            theta = r.theta,
            phi = r.phi,
            lambda = r.lambda,
        ))
    end
    return length(out) > degen_cap ? out[1:degen_cap] : out
end

function parse_degen_cap()
    cap = 80
    args = copy(ARGS)
    if "--degen-cap" in args
        i = findfirst(==("--degen-cap"), args)
        i !== nothing && length(args) >= i + 1 && (cap = parse(Int, args[i + 1]))
    end
    haskey(ENV, "DEGEN_CAP") && (cap = parse(Int, ENV["DEGEN_CAP"]))
    return clamp(cap, 1, 1845)
end

"""Small neighborhoods around an anchor."""
function refine_around(theta, phi, lam; n = 5, scale = 0.02)
    out = Tuple{Float64,Float64,Float64}[]
    for dth in range(-scale, scale; length = n)
        for dph in range(-scale, scale; length = n)
            for dlm in range(-scale, scale; length = n)
                push!(out, (theta + dth, phi + dph, lam + dlm))
            end
        end
    end
    return unique(out)
end

function evaluate_point(name, theta, phi, lam; certify = true, verbose = false)
    z4_dev = karlsson_z4_dev(theta, phi, lam)
    H = try
        build_karlsson_family(theta, phi, lam)
    catch e
        return (name = name, theta = theta, phi = phi, lambda = lam,
                status = :build_fail, error = string(e), z4_dev = z4_dev)
    end
    if !is_hadamard(H)
        return (name = name, theta = theta, phi = phi, lambda = lam,
                status = :not_hadamard, z4_dev = z4_dev)
    end

    pc = nothing
    cert_info = nothing
    if certify
        try
            pc = pool_completeness_report(H; tol = 1e-8, verbose = verbose)
            cert_info = (
                n_certified = pc.n_certified,
                n_distinct = pc.n_distinct_certified,
                n_tracked = pc.n_tracked,
                n_distinct_real = missing,
                pool_complete_flag = pc.pool_complete_flag,
            )
        catch e
            cert_info = (error = string(e),)
            pc = (pool_complete_flag = false, n_tracked = missing, n_certified = missing,
                  n_distinct_certified = missing, n_dedup_pool = missing)
        end
    end

    pool, stats = generate_candidate_pool_fresh(H; verbose = verbose)
    pool = deduplicate_pool(pool)
    ext = check_four_mub_extension(pool, H)
    by_tol = Dict{Float64, NamedTuple}()
    for mu_tol in MU_TOLS
        by_tol[mu_tol] = check_four_mub_extension(pool, H; mu_tol = mu_tol)
    end

    pool_complete = pc === nothing ? missing : get(pc, :pool_complete_flag, missing)

    return (
        name = name,
        theta = theta,
        phi = phi,
        lambda = lam,
        status = :ok,
        z4_dev = z4_dev,
        n_pool = length(pool),
        n_raw = stats.n_raw,
        n_conjugate = stats.n_conjugate,
        n_verified = stats.n_verified,
        max_clique = ext.max_clique,
        found_third = ext.found_third,
        found_fourth = ext.found_fourth,
        found_fourth_1e6 = by_tol[1e-6].found_fourth,
        found_fourth_1e8 = by_tol[1e-8].found_fourth,
        found_fourth_1e10 = by_tol[1e-10].found_fourth,
        n_tracked = cert_info === nothing ? missing : get(cert_info, :n_tracked, missing),
        n_certified = cert_info === nothing ? missing : get(cert_info, :n_certified, missing),
        n_distinct = cert_info === nothing ? missing : get(cert_info, :n_distinct, missing),
        n_distinct_real = cert_info === nothing ? missing : get(cert_info, :n_distinct_real, missing),
        pool_complete_flag = pool_complete,
        n_dedup_certified = pc === nothing ? missing : get(pc, :n_dedup_pool, missing),
    )
end

function needs_refinement(r; prev_n_pool = nothing)
    r.status != :ok && return false
    r.found_fourth && return false
    r.max_clique >= 6 && return true
    prev_n_pool !== nothing && r.n_pool != prev_n_pool && return true
    if r.pool_complete_flag === false
        return true
    end
    if r.n_certified !== missing && r.n_tracked !== missing && r.n_certified != r.n_tracked
        return true
    end
    return false
end

function _fmt_val(x)
    x === missing && return ""
    x isa Bool && return x ? "true" : "false"
    x isa Symbol && return string(x)
    x isa Number && (isnan(x) || !isfinite(x)) && return "NaN"
    return string(x)
end

function _row_to_dict(r)
    d = Dict{String,String}()
    for k in CSV_HEADER
        sym = Symbol(k)
        val = hasproperty(r, sym) ? getproperty(r, sym) : missing
        d[k] = _fmt_val(val)
    end
    return d
end

function init_csv_writer(; resume = false)
    mkpath(dirname(OUT_CSV))
    if resume && isfile(OUT_CSV)
        return open(OUT_CSV, "a")
    end
    io = open(OUT_CSV, "w")
    println(io, join(CSV_HEADER, ","))
    flush(io)
    return io
end

"""Load completed points from an existing CSV for --resume."""
function load_csv_state(path = OUT_CSV)
    seen = Set{Tuple{Float64,Float64,Float64}}()
    cache = Dict{Tuple{Float64,Float64,Float64}, NamedTuple}()
    !isfile(path) && return seen, cache, 0, 0, 0, 0
    df = CSV.read(path, DataFrame)
    n_third = n_fourth = n_cert = 0
    for r in eachrow(df)
        key = (round(r.theta, digits = 10), round(r.phi, digits = 10), round(r.lambda, digits = 10))
        push!(seen, key)
        status = Symbol(r.status)
        max_clique = coalesce(r.max_clique, 0)
        found_third = coalesce(r.found_third, false)
        found_fourth = coalesce(r.found_fourth, false)
        pool_complete = r.pool_complete_flag === true
        cache[key] = (
            status = status,
            max_clique = max_clique,
            found_fourth = found_fourth,
            n_pool = r.n_pool,
            pool_complete_flag = r.pool_complete_flag,
            n_certified = r.n_certified,
            n_tracked = r.n_tracked,
        )
        status == :ok && found_third && (n_third += 1)
        status == :ok && found_fourth && (n_fourth += 1)
        status == :ok && pool_complete && (n_cert += 1)
    end
    return seen, cache, nrow(df), n_third, n_fourth, n_cert
end

"""Whether this point may spawn adaptive refinements (primary, or first-level ref third-MUB hits)."""
function can_spawn_refinement(name, from_primary; max_clique = 0)
    from_primary && return true
    startswith(name, "ref_") && !occursin("ref_ref", name) && max_clique >= 6 && return true
    return false
end
function _maybe_enqueue_refinements!(refine_queue, seen, name, theta, phi, lam, r)
    needs_refinement(r) || return
    scale = r.max_clique >= 6 ? 0.005 : 0.01
    n_ref = r.max_clique >= 6 ? 5 : 3
    for (th, ph, lm) in refine_around(theta, phi, lam; n = n_ref, scale = scale)
        rk = (round(th, digits = 10), round(ph, digits = 10), round(lm, digits = 10))
        rk in seen && continue
        push!(refine_queue, ("ref_$(name)", th, ph, lm, 1, scale))
    end
end

function append_csv_row!(io, r)
    d = _row_to_dict(r)
    println(io, join((d[k] for k in CSV_HEADER), ","))
    flush(io)
end

"""Verify CSV: UTF-8, expected row count, required columns present."""
function verify_csv(path = OUT_CSV; expected_rows = 588)
    !isfile(path) && return (ok = false, reason = "file missing")
    bytes = read(path)
    if length(bytes) >= 2 && bytes[1] == 0xFF && bytes[2] == 0xFE
        return (ok = false, reason = "UTF-16 BOM detected (corrupted pipeline)")
    end
    if length(bytes) >= 2 && bytes[1] == 0xFE && bytes[2] == 0xFF
        return (ok = false, reason = "UTF-16 BE BOM detected")
    end
    lines = filter(!isempty, strip.(readlines(path)))
    length(lines) < 1 && return (ok = false, reason = "empty file")
    header = split(replace(lines[1], "\ufeff" => ""), ',')
    missing_cols = setdiff(CSV_HEADER, header)
    !isempty(missing_cols) && return (ok = false, reason = "missing columns: $(missing_cols)")
    n = length(lines) - 1
    n != expected_rows && return (ok = false, reason = "expected $expected_rows data rows, got $n", n_rows = n)
    bad = count(i -> length(split(lines[i], ',')) != length(CSV_HEADER), 2:length(lines))
    bad > 0 && return (ok = false, reason = "$bad rows with wrong column count", n_rows = n)
    return (ok = true, n_rows = n, n_cols = length(header))
end

function run_search(; quick = false, resume = false, degen_cap = 80, no_refine = false)
    seen, cache, n_points, n_third, n_fourth, n_certified_set =
        resume ? load_csv_state() : (Set{Tuple{Float64,Float64,Float64}}(), Dict(), 0, 0, 0, 0)
    io = init_csv_writer(; resume = resume)
    resume && println("Resume: $(length(seen)) points already in CSV")
    no_refine && println("No-refine mode: skipping adaptive refinement spawn")

    queue = Vector{NamedTuple}()
    for anchor in special_anchor_points()
        push!(queue, anchor)
    end
    degen_loaded = load_degeneracy_candidates(; degen_cap = degen_cap)
    println("Degeneracy candidates loaded: $(length(degen_loaded)) (cap=$degen_cap, sorted by -z4_dev)")
    for cand in degen_loaded
        push!(queue, cand)
    end

    refine_queue = Tuple{String,Float64,Float64,Float64,Int,Float64}[]
    if resume && !quick && !no_refine
        for pt in queue
            key = (round(pt.theta, digits = 10), round(pt.phi, digits = 10), round(pt.lambda, digits = 10))
            r = get(cache, key, nothing)
            r === nothing && continue
            _maybe_enqueue_refinements!(refine_queue, seen, pt.name, pt.theta, pt.phi, pt.lambda, r)
        end
        df_resume = CSV.read(OUT_CSV, DataFrame)
        for row in eachrow(df_resume)
            startswith(row.name, "ref_") || continue
            occursin("ref_ref", row.name) && continue
            coalesce(row.max_clique, 0) >= 6 || continue
            key = (round(row.theta, digits = 10), round(row.phi, digits = 10), round(row.lambda, digits = 10))
            r = get(cache, key, nothing)
            r === nothing && continue
            _maybe_enqueue_refinements!(refine_queue, seen, row.name, row.theta, row.phi, row.lambda, r)
        end
        println("  rebuilt refine_queue: $(length(refine_queue)) pending points")
    end

    while !isempty(queue) || !isempty(refine_queue)
        from_primary = false
        if !isempty(queue)
            pt = popfirst!(queue)
            name, theta, phi, lam = pt.name, pt.theta, pt.phi, pt.lambda
            from_primary = true
        else
            name, theta, phi, lam, level, scale = popfirst!(refine_queue)
        end
        key = (round(theta, digits = 10), round(phi, digits = 10), round(lam, digits = 10))
        if key in seen
            if !quick && !no_refine && can_spawn_refinement(name, from_primary; max_clique = get(get(cache, key, (; max_clique=0)), :max_clique, 0))
                r = get(cache, key, nothing)
                r !== nothing && _maybe_enqueue_refinements!(refine_queue, seen, name, theta, phi, lam, r)
            end
            continue
        end
        push!(seen, key)

        println("\n--- Point: $name ($(theta), $(phi), $(lam)) ---")
        r = evaluate_point(name, theta, phi, lam; certify = true)
        n_points += 1
        append_csv_row!(io, r)

        if r.status == :ok
            cache[key] = (
                status = r.status,
                max_clique = r.max_clique,
                found_fourth = r.found_fourth,
                n_pool = r.n_pool,
                pool_complete_flag = r.pool_complete_flag,
                n_certified = r.n_certified,
                n_tracked = r.n_tracked,
            )
            r.found_third && (n_third += 1)
            r.found_fourth && (n_fourth += 1)
            r.pool_complete_flag === true && (n_certified_set += 1)

            if r.found_fourth
                println("  *** COUNTEREXAMPLE: found_fourth=true — escalating to investigate_candidate ***")
                close(io)
                include(joinpath(@__DIR__, "investigate_candidate.jl"))
                investigate(r.theta, r.phi, r.lambda)
                return n_points, n_third, n_fourth, n_certified_set
            end

            if !quick && !no_refine && can_spawn_refinement(name, from_primary; max_clique = r.max_clique)
                _maybe_enqueue_refinements!(refine_queue, seen, name, theta, phi, lam, r)
            end
        end
    end
    close(io)
    return n_points, n_third, n_fourth, n_certified_set
end

function main()
    if "--verify-only" in ARGS
        v = verify_csv()
        if v.ok
            println("CSV OK: $(v.n_rows) data rows, $(v.n_cols) columns, UTF-8")
        else
            println("CSV FAIL: $(v.reason)")
        end
        return
    end

    quick = "--quick" in ARGS
    resume = "--resume" in ARGS
    no_refine = "--no-refine" in ARGS
    degen_cap = parse_degen_cap()
    quick && println("Quick mode: anchors + degeneracy only, no adaptive refinement spawn")
    resume && println("Resume mode: append to existing CSV, skip completed points")
    println("=== Targeted Karlsson loci search (certified on ALL points) ===\n")
    println("  degen_cap=$degen_cap  no_refine=$no_refine")

    if !isfile(DEGEN_CSV)
        println("Note: run `python degeneracy_scan.py` first for degeneracy candidates.")
    end

    n_points, n_third, n_fourth, n_cert = run_search(
        quick = quick, resume = resume, degen_cap = degen_cap, no_refine = no_refine)

    println("\n=== Summary ===")
    println("  points searched: $n_points")
    println("  pool_complete (certified set S members): $n_cert")
    println("  third MUB found:  $n_third")
    println("  fourth MUB found: $n_fourth")
    println("  output: $OUT_CSV")

    expected = quick ? n_points : 588
    v = verify_csv(; expected_rows = expected)
    if v.ok
        println("  CSV verify: PASS ($(v.n_rows) rows, UTF-8)")
    else
        println("  CSV verify: FAIL — $(v.reason)")
    end

    println("\nINTERPRETATION: third MUB at special loci (F6, Dita) is EXPECTED.")
    println("  found_fourth=true would be a counterexample to Zauner within Karlsson.")
    println("  Absence on certified set S is NOT a family-wide proof.")
end

main()
