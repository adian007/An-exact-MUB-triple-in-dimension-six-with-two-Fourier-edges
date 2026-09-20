# Path A3 staged reformulation: n_wit=1 fourth-MUB witness.
#
# Single candidate fourth-basis vector: unitarity + MU to I (via |z_i|=1) +
# MU to H + MU to B3. Strict subset of the n_wit=2 system (no pairwise ortho).
#
# Logical status: a fourth MUB would supply six solutions of this system;
# emptiness ⇒ no fourth MUB at this (H,B3). Not a "weaker" probe.
#
# Usage:
#   julia --project=. scripts/julia/certify_fourth_mub_witness_nwit1.jl
#   julia --project=. scripts/julia/certify_fourth_mub_witness_nwit1.jl --lambda 0.4 --verbose
#
# Output: results/certify_fourth_mub_witness_nwit1_dita.txt

include(joinpath(@__DIR__, "_paths.jl"))
include(joinpath(ROOT, "src", "dita_third_mub_construction.jl"))

using Printf
using Dates

const OUT = joinpath(RESULTS_DIR, "certify_fourth_mub_witness_nwit1_dita.txt")
const N_WIT = 1

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
            println("Usage: certify_fourth_mub_witness_nwit1.jl [--lambda 0.4] [--verbose]")
            exit(0)
        else
            error("Unknown arg: $a")
        end
    end
    return (lambda=λ, verbose=verbose)
end

function interpret_certify(info)
    if hasproperty(info, :skipped_solve) && info.skipped_solve
        return "SKIPPED_MV_$(info.mixed_volume)"
    end
    n_pass = hasproperty(info, :n_certified_pass_full_residual) ?
        info.n_certified_pass_full_residual : 0
    if n_pass > 0
        return "CERTIFIED_ROOTS_FOUND_$(info.n_certified)_fullpass_$(n_pass)"
    elseif info.squared_for_certify && info.n_certified > 0 && n_pass == 0
        # All isolated roots of a generic square subsystem fail full residual
        # ⇒ overdetermined witness has no isolated common roots (Bertini pattern).
        return "EMPTY_FULL_RESIDUAL_AFTER_CERTIFY_SQUARE_$(info.n_certified)"
    elseif info.n_tracked == 0
        return "ZERO_TRACKED_ZERO_CERTIFIED"
    elseif info.n_certified == 0 && info.n_tracked > 0
        return "TRACKED_BUT_NONE_CERTIFIED"
    else
        return "INCONCLUSIVE"
    end
end

function main(args=ARGS)
    opts = parse_args(args)
    λ = opts.lambda
    mkpath(RESULTS_DIR)

    println("="^72)
    println("Fourth-MUB n_wit=1 necessary-condition witness @ Dita λ=$λ")
    println("Cap = $WITNESS_MV_SOLVE_CAP")
    println("="^72)

    H = build_karlsson_family(DITA_THETA, DITA_PHI, λ)
    t0 = time()

    cached = load_cached_b3(λ)
    if cached === nothing
        error("Need cached B3 at results/reconstruct_b3_B3_lambda$(λ).csv")
    end
    println("Using cached B3; clique_indices = ", cached.clique_indices)

    probe = probe_fourth_mub_witness_mv(H, cached.B3; n_wit=N_WIT)
    @printf("Profile: n_eqs=%d n_vars=%d mixed_volume=%d\n",
            probe.n_eqs, probe.n_vars, probe.mixed_volume)

    info = certify_fourth_mub_witness_at_H(H; n_wit=N_WIT,
        clique_indices=cached.clique_indices, B3=cached.B3,
        verbose=opts.verbose)
    elapsed = round(time() - t0; digits=2)
    verdict = interpret_certify(info)

    n_paths = info.mixed_volume  # BKK start-path count
    n_finite = info.n_tracked

    open(OUT, "w") do io
        println(io, "Fourth-MUB n_wit=1 necessary-condition witness (Path A3 staged)")
        println(io, "timestamp = ", Dates.now())
        println(io, "lambda = ", λ)
        println(io, "theta_D = ", DITA_THETA)
        println(io, "phi = ", DITA_PHI)
        println(io, "clique_indices = ", info.clique_indices)
        println(io, "n_wit = ", N_WIT)
        println(io, "n_eqs = ", info.n_eqs)
        println(io, "n_vars = ", info.n_vars)
        println(io, "mixed_volume = ", info.mixed_volume)
        println(io, "solve_cap = ", WITNESS_MV_SOLVE_CAP)
        if info.skipped_solve
            println(io, "skipped_solve = true")
            println(io, "skip_reason = ", info.skip_reason)
            println(io, "certify_called = false")
        else
            println(io, "skipped_solve = false")
            println(io, "certify_called = true")
        end
        println(io, "squared_for_certify = ", info.squared_for_certify)
        println(io, "n_paths_BKK = ", n_paths)
        println(io, "n_tracked_solutions = ", info.n_tracked)
        println(io, "n_certified = ", info.n_certified)
        println(io, "n_distinct_certified = ", info.n_distinct_certified)
        println(io, "n_real_certified = ", info.n_real_certified)
        println(io, "n_pass_full_residual = ", info.n_certified_pass_full_residual)
        println(io, "found_fourth_combinatorial = ", info.found_fourth_combinatorial)
        println(io, "elapsed_s = ", elapsed)
        println(io, "verdict = ", verdict)
        println(io)
        println(io, "SCOPE:")
        println(io, "  System asks for ONE unit vector v = z/√6 with |z_i|=1 (MU to I),")
        println(io, "  MU to all columns of H, and MU to all columns of B3.")
        println(io, "  A fourth MUB would supply six such vectors ⇒ emptiness of this")
        println(io, "  system proves no fourth MUB at this (H,B3).")
        println(io)
        println(io, "INTERPRETATION:")
        if startswith(string(verdict), "SKIPPED_MV")
            println(io, "  Mixed volume still exceeds cap; certify() was NOT called.")
        elseif startswith(string(verdict), "CERTIFIED_ROOTS_FOUND")
            println(io, "  certify() FOUND certified roots of the (squared) system that also")
            println(io, "  pass full 17-equation residual checks. Necessary condition SATISFIABLE")
            println(io, "  at this (H,B3). Does NOT prove a fourth MUB exists.")
        elseif startswith(string(verdict), "EMPTY_FULL_RESIDUAL_AFTER_CERTIFY_SQUARE")
            println(io, "  Method: generic square (10 random linear combinations of 17 eqs),")
            println(io, "  polyhedral homotopy (mv=$(info.mixed_volume)), certify() on all")
            println(io, "  $(info.n_certified) nonsingular roots, then residual-check full system.")
            println(io, "  Result: $(info.n_certified)/$(info.n_certified) certified square-roots;")
            println(io, "  0 pass the full 17-equation residual (tol 1e-8).")
            println(io, "  Claim: no fourth MUB at λ=$λ for this B3 (W1 empty ⇒ no fourth).")
            println(io, "  Stronger than minimal T3: also no single MU-to-{I,H,B3} vector.")
            println(io, "  NOT a symbolic proof for all λ on the circle.")
        elseif startswith(string(verdict), "CERTIFIED_SQUARED_ONLY")
            println(io, "  certify() certified roots of the random-squared system, but they")
            println(io, "  FAIL the full overdetermined residual — discard as false positives")
            println(io, "  of the squared subsystem.")
        elseif verdict == "ZERO_TRACKED_ZERO_CERTIFIED"
            println(io, "  Polyhedral homotopy (mv=$(info.mixed_volume) start paths) returned")
            println(io, "  0 finite tracked solutions; certify() certified 0.")
            println(io, "  Inconclusive without the Bertini-square residual pipeline.")
        elseif verdict == "TRACKED_BUT_NONE_CERTIFIED"
            println(io, "  Tracked $(info.n_tracked) solutions but certify() accepted 0.")
            println(io, "  Inconclusive for emptiness (approximate roots may be singular/junk).")
        else
            println(io, "  See verdict; do not oversell.")
        end
    end

    println("\nResults:")
    @printf("  mixed_volume = %d  (cap %d)\n", info.mixed_volume, WITNESS_MV_SOLVE_CAP)
    @printf("  skipped_solve = %s  certify_called = %s  squared = %s\n",
            info.skipped_solve, !info.skipped_solve, info.squared_for_certify)
    @printf("  tracked = %d  certified = %d  full_residual_pass = %d\n",
            info.n_tracked, info.n_certified, info.n_certified_pass_full_residual)
    println("  verdict = ", verdict)
    println("  elapsed = $(elapsed)s")
    println("  wrote: ", OUT)
    println("="^72)
    return info
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
