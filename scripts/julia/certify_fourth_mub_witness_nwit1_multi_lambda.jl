# Path A3: n_wit=1 witness certification at several Dita λ on the circle.
#
# At each λ: build H(θ_D, φ=π/4, λ), obtain a third-MUB B3 from the fresh MU pool,
# form the n_wit=1 system, Bertini-square, solve+certify, residual-check.
#
# Logical status of emptiness (see FINDINGS §12b):
#   n_wit=1 empty ⇒ no fourth MUB at that (H,B3)
#   (a fourth MUB would supply six solutions of the n_wit=1 system).
#
# Usage:
#   julia --project=. scripts/julia/certify_fourth_mub_witness_nwit1_multi_lambda.jl
#
# Output: results/certify_fourth_mub_witness_nwit1_multi_lambda.txt

include(joinpath(@__DIR__, "_paths.jl"))
include(joinpath(ROOT, "src", "dita_third_mub_construction.jl"))

using Printf
using Dates

const OUT = joinpath(RESULTS_DIR, "certify_fourth_mub_witness_nwit1_multi_lambda.txt")
const N_WIT = 1

# Spread across [0, 2π); 0.4 already certified in the single-λ run.
const LAMBDAS = [
    0.4,
    π / 3,       # ≈ 1.047
    2π / 3,      # ≈ 2.094
    π,           # ≈ 3.142
    4π / 3,      # ≈ 4.189
    5π / 3,      # ≈ 5.236
]

function load_cached_b3(lambda_)
    # Exact string match only for the known cache naming (λ=0.4).
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

function interpret_certify(info)
    if hasproperty(info, :skipped_solve) && info.skipped_solve
        return "SKIPPED_MV_$(info.mixed_volume)"
    end
    n_pass = hasproperty(info, :n_certified_pass_full_residual) ?
        info.n_certified_pass_full_residual : 0
    if n_pass > 0
        return "CERTIFIED_ROOTS_FOUND_$(info.n_certified)_fullpass_$(n_pass)"
    elseif info.squared_for_certify && info.n_certified > 0 && n_pass == 0
        return "EMPTY_FULL_RESIDUAL_AFTER_CERTIFY_SQUARE_$(info.n_certified)"
    elseif info.n_tracked == 0
        return "ZERO_TRACKED_ZERO_CERTIFIED"
    elseif info.n_certified == 0 && info.n_tracked > 0
        return "TRACKED_BUT_NONE_CERTIFIED"
    else
        return "INCONCLUSIVE"
    end
end

function run_one(λ; verbose=false)
    H = build_karlsson_family(DITA_THETA, DITA_PHI, λ)
    cached = load_cached_b3(λ)
    t0 = time()
    if cached !== nothing
        info = certify_fourth_mub_witness_at_H(H; n_wit=N_WIT,
            clique_indices=cached.clique_indices, B3=cached.B3, verbose=verbose)
        b3_source = "cache"
    else
        info = certify_fourth_mub_witness_at_H(H; n_wit=N_WIT, verbose=verbose)
        b3_source = "fresh_pool"
    end
    elapsed = round(time() - t0; digits=2)
    verdict = interpret_certify(info)
    return (
        lambda=λ,
        b3_source=b3_source,
        clique_indices=info.clique_indices,
        n_eqs=info.n_eqs,
        n_vars=info.n_vars,
        mixed_volume=info.mixed_volume,
        skipped_solve=info.skipped_solve,
        certify_called=!info.skipped_solve,
        squared=info.squared_for_certify,
        n_tracked=info.n_tracked,
        n_certified=info.n_certified,
        n_pass_full=info.n_certified_pass_full_residual,
        found_fourth_combinatorial=info.found_fourth_combinatorial,
        verdict=verdict,
        elapsed_s=elapsed,
    )
end

function main()
    mkpath(RESULTS_DIR)
    println("="^72)
    println("n_wit=1 multi-λ certification @ Dita circle")
    println("λ list = ", LAMBDAS)
    println("Cap = ", WITNESS_MV_SOLVE_CAP)
    println("="^72)

    rows = []
    for λ in LAMBDAS
        @printf("\n--- λ = %.10f ---\n", λ)
        r = run_one(λ)
        push!(rows, r)
        @printf("  mv=%d tracked=%d certified=%d fullpass=%d verdict=%s (%.1fs, B3=%s)\n",
                r.mixed_volume, r.n_tracked, r.n_certified, r.n_pass_full,
                r.verdict, r.elapsed_s, r.b3_source)
    end

    n_empty = count(r -> startswith(string(r.verdict), "EMPTY_FULL_RESIDUAL"), rows)
    n_sat = count(r -> startswith(string(r.verdict), "CERTIFIED_ROOTS_FOUND"), rows)
    n_skip = count(r -> r.skipped_solve, rows)
    n_other = length(rows) - n_empty - n_sat - n_skip

    open(OUT, "w") do io
        println(io, "Fourth-MUB n_wit=1 multi-λ certification (Path A3)")
        println(io, "timestamp = ", Dates.now())
        println(io, "theta_D = ", DITA_THETA)
        println(io, "phi = ", DITA_PHI)
        println(io, "n_wit = ", N_WIT)
        println(io, "solve_cap = ", WITNESS_MV_SOLVE_CAP)
        println(io, "n_lambda = ", length(rows))
        println(io, "n_empty = ", n_empty)
        println(io, "n_satisfiable = ", n_sat)
        println(io, "n_skipped = ", n_skip)
        println(io, "n_other = ", n_other)
        println(io)
        println(io, "LOGIC:")
        println(io, "  A fourth MUB at (H,B3) is six orthonormal vectors each MU to")
        println(io, "  {I, H, B3}. Each such vector is a solution of the n_wit=1 system.")
        println(io, "  Therefore: n_wit=1 emptiness ⇒ no fourth MUB at that point.")
        println(io, "  (Converse false in general: candidate vectors may exist without")
        println(io, "  forming a 6-ONB. Emptiness is strictly stronger than minimal T3.)")
        println(io)
        println(io, "PER-λ RESULTS:")
        for r in rows
            @printf(io, "lambda=%.10f  mv=%d  tracked=%d  certified=%d  fullpass=%d  ",
                    r.lambda, r.mixed_volume, r.n_tracked, r.n_certified, r.n_pass_full)
            println(io, "verdict=", r.verdict)
            println(io, "  B3_source=", r.b3_source, "  clique=", r.clique_indices)
            println(io, "  certify_called=", r.certify_called,
                    "  squared=", r.squared,
                    "  found_fourth_comb=", r.found_fourth_combinatorial,
                    "  elapsed_s=", r.elapsed_s)
        end
        println(io)
        println(io, "SUMMARY:")
        if n_empty == length(rows) && n_skip == 0 && n_sat == 0
            println(io, "  All ", length(rows), " λ values: certified empty (full residual 0).")
            println(io, "  T3 upgrade: certified non-existence at ", length(rows),
                    " independent Dita-λ points + existing 628/628 numerical density.")
            println(io, "  Still NOT a symbolic proof for all λ; λ-invariance remains optional.")
        elseif n_sat > 0
            println(io, "  WARNING: n_wit=1 SATISFIABLE at some λ — necessary condition holds;")
            println(io, "  does NOT prove a fourth MUB, but blocks emptiness-based T3 there.")
        else
            println(io, "  Mixed outcomes; see per-λ rows. Do not oversell.")
        end
    end

    println("\n" * "="^72)
    @printf("Done: empty=%d sat=%d skip=%d other=%d / %d\n",
            n_empty, n_sat, n_skip, n_other, length(rows))
    println("wrote: ", OUT)
    println("="^72)
    return rows
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
