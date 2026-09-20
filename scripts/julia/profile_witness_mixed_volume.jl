# Path A3 salvage check: mixed volume of the reduced (n_wit=2) witness
# vs any smaller subsystem — do NOT raise the solve cap; only profile.
#
# Usage: julia --project=. scripts/julia/profile_witness_mixed_volume.jl
#
# Writes: results/witness_mixed_volume_profile.txt

include(joinpath(@__DIR__, "_paths.jl"))
include(joinpath(ROOT, "src", "dita_third_mub_construction.jl"))
using Printf, Dates

const OUT = joinpath(RESULTS_DIR, "witness_mixed_volume_profile.txt")
const CACHE = joinpath(RESULTS_DIR, "reconstruct_b3_B3_lambda0.4.csv")

function load_cached_b3()
    rows = readlines(CACHE)
    B3 = Matrix{ComplexF64}(undef, 6, 6)
    for i in 1:6
        cols = split(rows[2 + i], ';')
        for j in 1:6
            re, im_ = parse.(Float64, split(cols[j], ','))
            B3[i, j] = ComplexF64(re, im_)
        end
    end
    idx = parse.(Int, split(split(rows[2], '=')[2], ';'))
    return B3, idx
end

function main()
    λ = 0.4
    H = build_karlsson_family(DITA_THETA, DITA_PHI, λ)
    B3, idx = load_cached_b3()

    open(OUT, "w") do io
        println(io, "=== A3 witness mixed-volume profile ===")
        println(io, "Date: ", Dates.now())
        println(io, "lambda = ", λ)
        println(io, "clique_indices = ", idx)
        println(io, "")
        println(io, "NOTE: build_numeric_fourth_mub_witness_system only implements n_wit=2")
        println(io, "  (20 vars, 35 eqs) — same size as symbolic_export/fourth_mub_reduced_Dita.m2.")
        println(io, "  There is NO 5-vector (50-var) formulation in the current pipeline.")
        println(io, "")

        sys2 = build_numeric_fourth_mub_witness_system(H, B3; n_wit=2)
        mv2 = mixed_volume(sys2)
        @printf(io, "n_wit=2 reduced witness: n_eqs=%d n_vars=%d mixed_volume=%d\n",
                length(sys2), nvariables(sys2), mv2)
        @printf(io, "  WITNESS_MV_SOLVE_CAP = %d\n", WITNESS_MV_SOLVE_CAP)
        @printf(io, "  tractable for HC solve+certify under current cap? %s\n",
                mv2 <= WITNESS_MV_SOLVE_CAP ? "YES" : "NO")
        println(io, "")
        println(io, "CONCLUSION:")
        if mv2 > WITNESS_MV_SOLVE_CAP
            println(io, "  The ALREADY-REDUCED 2-vector witness has mv=$mv2.")
            println(io, "  A3 is not blocked by using a larger 5-vector system — we are already")
            println(io, "  on the reduced formulation matching fourth_mub_reduced_Dita.m2.")
            println(io, "  Salvage options are reformulation (different equations / gauges),")
            println(io, "  not further n_wit reduction in the current API.")
        else
            println(io, "  Reduced witness is within cap — full solve+certify should be attempted.")
        end
    end
    println("Wrote $OUT")
end

main()
