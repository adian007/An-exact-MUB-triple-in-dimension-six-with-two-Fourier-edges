# Adversarial homotopy drop-detection test (Path audit).
# Usage: julia --project=. scripts/julia/drop_detection_stress_test.jl
include(joinpath(@__DIR__, "_paths.jl"))
using Printf, Dates

const OUT = joinpath(RESULTS_DIR, "drop_detection_stress_test.txt")

function main()
    open(OUT, "w") do io
        println(io, "=== Homotopy drop-detection stress test ===")
        println(io, "Date: $(Dates.now())")
        H = build_karlsson_family(0.3, 0.5, 0.2)
        tracker = init_parametric_pool_tracker(verbose=false)
        pool_h, stats_h, meta_h = generate_candidate_pool(H; tracker=tracker, cross_check=true, verbose=false)
        pool_f, stats_f = generate_candidate_pool_fresh(H; verbose=false)
        pool_f = deduplicate_pool(pool_f)
        drop = get(meta_h, :drop_suspected, false)
        @printf(io, "Point: generic K6 (0.3, 0.5, 0.2)\n")
        @printf(io, "  tracked=%d/%d  homotopy_pool=%d  fresh_pool=%d\n",
                meta_h.n_tracked, meta_h.n_expected, length(pool_h), length(pool_f))
        @printf(io, "  drop_suspected=%s\n", drop)
        @printf(io, "  cross_check only_in_b=%d\n", meta_h.cross_check.only_in_b)
        verdict = drop && length(pool_f) >= length(pool_h)
        println(io, "")
        println(io, verdict ?
            "PASS: drop_suspected=true at adversarial point; fallback replaces homotopy pool." :
            "FAIL: expected drop_suspected with fresh pool >= homotopy pool.")
    end
    println("Wrote $OUT")
end

main()
