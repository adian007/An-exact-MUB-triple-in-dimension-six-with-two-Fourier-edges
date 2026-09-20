# Path A3: interval-certified fourth-MUB witness exclusion at Dita anchor.
#
# Solves the reduced witness system (35 eq, 20 var, n_wit=2) via HomotopyContinuation
# and applies certify() — same engine as certify_pool_at_H.
#
# Usage:
#   julia --project=. scripts/julia/certify_fourth_mub_witness.jl
#   julia --project=. scripts/julia/certify_fourth_mub_witness.jl --lambda 0.4 --verbose
#
# Output: results/certify_fourth_mub_witness_dita.txt

include(joinpath(@__DIR__, "_paths.jl"))
include(joinpath(ROOT, "src", "dita_third_mub_construction.jl"))

using Printf
using Dates

const OUT = joinpath(RESULTS_DIR, "certify_fourth_mub_witness_dita.txt")

function load_cached_b3(lambda_)
    cache = joinpath(RESULTS_DIR, "reconstruct_b3_B3_lambda$(lambda_).csv")
    isfile(cache) || return nothing
    rows = readlines(cache)
    length(rows) < 8 && return nothing
    B3 = Matrix{ComplexF64}(undef, 6, 6)
    for i in 1:6
        cols = split(rows[2 + i], ';')
        for j in 1:6
            re, im_ = parse.(Float64, split(cols[j], ','))
            B3[i, j] = ComplexF64(re, im_)
        end
    end
    clique_idx = parse.(Int, split(split(rows[2], '=')[2], ';'))
    return (B3=B3, clique_indices=clique_idx)
end

function parse_args(args)
    λ = 0.4
    verbose = false
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--lambda" && i < length(args)
            λ = parse(Float64, args[i+1]); i += 2
        elseif a == "--verbose"
            verbose = true; i += 1
        elseif a in ("-h", "--help")
            println("Usage: certify_fourth_mub_witness.jl [--lambda 0.4] [--verbose]")
            exit(0)
        else
            error("Unknown arg: $a")
        end
    end
    return (lambda=λ, verbose=verbose)
end

function main(args=ARGS)
    opts = parse_args(args)
    λ = opts.lambda
    mkpath(RESULTS_DIR)

    println("="^72)
    println("Fourth-MUB witness certification @ Dita λ=$λ")
    println("="^72)

    H = build_karlsson_family(DITA_THETA, DITA_PHI, λ)
    t0 = time()

    cached = load_cached_b3(λ)
    if cached !== nothing
        println("Using cached B3 from reconstruct_b3_B3_lambda$(λ).csv")
        info = certify_fourth_mub_witness_at_H(H; n_wit=2,
            clique_indices=cached.clique_indices, B3=cached.B3,
            verbose=opts.verbose)
    else
        info = certify_fourth_mub_witness_at_H(H; n_wit=2, verbose=opts.verbose)
    end
    elapsed = round(time() - t0; digits=2)

    # Path statistics from solve result
    if info.result === nothing
        n_paths = 0
        n_finite = 0
        verdict = info.skipped_solve ? "SKIPPED_MV_$(info.mixed_volume)" : "NO_SOLVE"
    else
        res = info.result
        n_paths = npaths(res)
        n_finite = count(s -> s.return_code == :success, paths(res))
        verdict = if info.n_certified == 0 && info.n_tracked == 0
            "NO_FINITE_WITNESS_ROOTS"
        elseif info.n_certified == 0 && info.n_tracked > 0
            "TRACKED_BUT_NONE_CERTIFIED"
        elseif info.n_certified > 0
            "WITNESS_ROOTS_FOUND"
        else
            "INCONCLUSIVE"
        end
    end

    open(OUT, "w") do io
        println(io, "Fourth-MUB witness certification (Path A3)")
        println(io, "timestamp = ", Dates.now())
        println(io, "lambda = ", λ)
        println(io, "theta_D = ", DITA_THETA)
        println(io, "phi = ", DITA_PHI)
        println(io, "clique_indices = ", info.clique_indices)
        println(io, "n_wit = ", info.n_wit)
        println(io, "n_eqs = ", info.n_eqs)
        println(io, "n_vars = ", info.n_vars)
        println(io, "mixed_volume = ", info.mixed_volume)
        if hasproperty(info, :skipped_solve) && info.skipped_solve
            println(io, "skipped_solve = true")
            println(io, "skip_reason = ", info.skip_reason)
        end
        println(io, "n_paths = ", n_paths)
        println(io, "n_finite_success = ", n_finite)
        println(io, "n_tracked_solutions = ", info.n_tracked)
        println(io, "n_certified = ", info.n_certified)
        println(io, "n_distinct_certified = ", info.n_distinct_certified)
        println(io, "n_real_certified = ", info.n_real_certified)
        println(io, "found_fourth_combinatorial = ", info.found_fourth_combinatorial)
        println(io, "elapsed_s = ", elapsed)
        println(io, "verdict = ", verdict)
        println(io)
        println(io, "INTERPRETATION:")
        if startswith(string(verdict), "SKIPPED_MV") || startswith(string(verdict), "INFEASIBLE_MV")
            println(io, "  Mixed volume ", info.mixed_volume, " exceeds tractable homotopy cap.")
            println(io, "  certify() was NEVER CALLED — the solve was skipped before tracking.")
            println(io, "  This is NOT an interval-arithmetic non-existence certificate.")
            println(io, "  HC certify proves existence of roots near found solutions;")
            println(io, "  it does NOT certify ideal emptiness. Interval exclusion remains open.")
            println(io, "  Status fields: skipped_solve=true, n_certified=0, n_paths=0.")
        elseif verdict == "NO_FINITE_WITNESS_ROOTS"
            println(io, "  HomotopyContinuation found 0 finite roots of the reduced witness")
            println(io, "  system (35 eq, 20 var). certify() therefore certifies nothing.")
            println(io, "  Combined with found_fourth=false from clique search, this is")
            println(io, "  stronger than combinatorial-only evidence BUT is NOT a proof of")
            println(io, "  ideal emptiness (paths may miss roots; no exclusion certificate).")
        elseif verdict == "WITNESS_ROOTS_FOUND"
            println(io, "  WARNING: certified witness roots exist — contradicts combinatorial")
            println(io, "  found_fourth=false; investigate clique vs polynomial encoding.")
        else
            println(io, "  See verdict above; HC certify does not prove non-existence globally.")
        end
    end

    println("\nResults:")
    @printf("  mixed_volume = %d\n", info.mixed_volume)
    @printf("  n_paths = %d  finite_success = %d\n", n_paths, n_finite)
    @printf("  n_tracked = %d  n_certified = %d  n_distinct = %d\n",
            info.n_tracked, info.n_certified, info.n_distinct_certified)
    println("  found_fourth (combinatorial) = ", info.found_fourth_combinatorial)
    println("  verdict = ", verdict)
    println("  elapsed = $(elapsed)s")
    println("  wrote: ", OUT)
    println("="^72)
    return info
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
