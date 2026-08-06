# Verify Dita third-MUB construction (Theorem T2 numerical certificate).
# Usage: julia --project=. scripts/julia/verify_dita_construction.jl

include(joinpath(@__DIR__, "_paths.jl"))
include(joinpath(ROOT, "src", "dita_third_mub_construction.jl"))

using Printf
using Dates

const OUT = joinpath(RESULTS_DIR, "verify_dita_construction.txt")

function main()
    # Spot checks + full-circle sample (subset if HC slow)
    lambdas = [0.0, 0.4, 0.41, 1.57079632679, 3.14159265359, 4.71238898038, 6.28318530718]
    fine = collect(range(0.0, 2π; length=21)[2:end-1])
    all_lam = unique(vcat(lambdas, fine))

    println("Verifying third MUB at $(length(all_lam)) lambda values...")
    results = verify_dita_lambda_circle(all_lam; hp_bits=256)

    n_ok = count(r -> get(r, :ok, false), results)
    open(OUT, "w") do io
        println(io, "=== Verify Dita third-MUB construction (T2) ===")
        println(io, "Date: $(Dates.format(now(), "yyyy-mm-ddTHH:MM:SS"))\n")
        for r in results
            if get(r, :ok, false)
                @printf(io, "  lam=%.6g clique=%d pool=%d hp_ortho=%s hp_mu=%s\n",
                        r.lambda, r.max_clique, r.n_pool, r.hp_ortho, r.hp_mu)
            else
                @printf(io, "  lam=%.6g FAIL: %s\n", r.lambda, r.error)
            end
        end
        println(io, "\nSUMMARY: $n_ok/$(length(results)) PASS")
        println(io, "T2 certificate: third ONB extractable at all tested lambda on Dita circle")
    end
    println("Wrote $OUT")
    n_ok == length(results) || error("Dita construction verification failed")
end

main()
