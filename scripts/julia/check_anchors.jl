# Quick third/fourth MUB check at anchor points.
# Usage: julia --project=. check_anchors.jl

include(joinpath(@__DIR__, "_paths.jl"))

function report(name, H)
    pool, stats = generate_candidate_pool_fresh(H; verbose = false)
    pool = deduplicate_pool(pool)
    ext = check_four_mub_extension(pool, H)
    println("$name: n_pool=$(length(pool)) raw=$(stats.n_raw) max_clique=$(ext.max_clique) " *
            "third=$(ext.found_third) fourth=$(ext.found_fourth)")
end

report("F6", F6)
report("generic K6(0.3,0.5,0.2)", build_karlsson_family(0.3, 0.5, 0.2))
report("Dita", build_karlsson_family(acos(1 / sqrt(3)), pi / 4, 0.4))
report("theta=0", build_karlsson_family(0.0, 0.5, 0.3))
