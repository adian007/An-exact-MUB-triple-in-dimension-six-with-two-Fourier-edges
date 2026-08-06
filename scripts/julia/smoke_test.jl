# Smoke test: per-H pipeline for Karlsson family.
# Usage:  julia --project=. smoke_test.jl

include(joinpath(@__DIR__, "_paths.jl"))

function main()
    println("=== Per-H smoke test: does {I, H} extend to four MUBs? ===\n")
    tracker = init_parametric_pool_tracker(verbose=true)

    for (theta, phi, lam) in ((0.3, 0.5, 0.2), (0.0, 0.5, 0.3), (acos(1/sqrt(3)), pi/4, 0.4))
        r = run_search(:karlsson, theta, phi, lam; tracker=tracker, cross_check=true)
        println("  θ=$(round(theta;digits=3)) φ=$(round(phi;digits=3)) λ=$(round(lam;digits=3)):" *
                "  n_pool=$(r.n_pool)  max_clique=$(r.max_clique)" *
                "  third=$(r.found_third)  fourth=$(r.found_fourth)" *
                "  drop=$(r.meta.drop_suspected)")
    end
    println("\nRIGOR: heuristic (homotopy + fresh fallback, floating tolerances).")
end

main()
