# Task 2: Re-tag 512-point sweep with corrected pool_complete_flag.
# OLD: n_tracked >= 0.9 * mixed_volume (strict mixed-volume threshold from METHODS v1)
# NEW: n_certified == n_tracked && n_certified > 0 && n_verified > 0
#
# Usage: julia --project=. retag_512_completeness.jl [--sample N]

include(joinpath(@__DIR__, "_paths.jl"))
using Printf

const SWEEP_CSV = joinpath(RESULTS_DIR, "karlsson_perH_sweep.csv")
const OUT_CSV = joinpath(RESULTS_DIR, "karlsson_perH_sweep_retagged.csv")
const SUMMARY_TXT = joinpath(RESULTS_DIR, "retag_512_summary.txt")

function old_pool_complete(n_tracked, mv)
    n_tracked >= 0.9 * mv
end

function parse_sweep_csv(path)
    rows = NamedTuple[]
    open(path) do io
        readline(io)  # header
        for line in eachline(io)
            isempty(strip(line)) && continue
            p = split(line, ',')
            length(p) < 11 && continue
            push!(rows, (
                theta = parse(Float64, p[1]),
                phi = parse(Float64, p[2]),
                lambda = parse(Float64, p[3]),
                n_pool = parse(Float64, p[4]),
                max_clique = parse(Float64, p[5]),
                found_third = parse(Float64, p[6]) > 0.5,
                found_fourth = parse(Float64, p[7]) > 0.5,
                homotopy_used = parse(Float64, p[9]) > 0.5,
                fresh_fallback = parse(Float64, p[10]) > 0.5,
                drop_suspected = parse(Float64, p[11]) > 0.5,
            ))
        end
    end
    return rows
end

function write_retag_csv(path, rows)
    hdr = join([
        "theta", "phi", "lambda", "mixed_volume", "n_tracked", "n_certified",
        "n_distinct_certified", "n_verified", "n_dedup_pool",
        "old_pool_complete", "new_pool_complete", "cert_match", "mv_gap",
        "orig_n_pool", "orig_max_clique", "orig_found_third", "orig_found_fourth",
        "homotopy_used", "fresh_fallback", "drop_suspected",
    ], ",")
    open(path, "w") do io
        println(io, hdr)
        for r in rows
            println(io, join([
                r.theta, r.phi, r.lambda, r.mixed_volume, r.n_tracked, r.n_certified,
                r.n_distinct_certified, r.n_verified, r.n_dedup_pool,
                r.old_pool_complete, r.new_pool_complete, r.cert_match, r.mv_gap,
                r.orig_n_pool, r.orig_max_clique, r.orig_found_third, r.orig_found_fourth,
                r.homotopy_used, r.fresh_fallback, r.drop_suspected,
            ], ","))
        end
    end
end

function retag_point(theta, phi, lam; verbose = false)
    H = build_karlsson_family(theta, phi, lam)
    pc = pool_completeness_report(H; verbose = verbose)
    old_flag = old_pool_complete(pc.n_tracked, pc.mixed_volume)
    new_flag = pc.pool_complete_flag
    return (
        theta = theta, phi = phi, lambda = lam,
        mixed_volume = pc.mixed_volume,
        n_tracked = pc.n_tracked,
        n_certified = pc.n_certified,
        n_distinct_certified = pc.n_distinct_certified,
        n_verified = pc.n_verified,
        n_dedup_pool = pc.n_dedup_pool,
        old_pool_complete = old_flag,
        new_pool_complete = new_flag,
        cert_match = pc.cert_match,
        mv_gap = pc.mv_gap,
    )
end

function main()
    sample_n = nothing
    for (i, a) in enumerate(ARGS)
        a == "--sample" && (sample_n = parse(Int, ARGS[i+1]))
    end

    sweep_rows = parse_sweep_csv(SWEEP_CSV)
    n_total = length(sweep_rows)
    if sample_n !== nothing
        sweep_rows = sweep_rows[1:min(sample_n, n_total)]
        println("Sampling first $(length(sweep_rows)) of $n_total points")
    end

    rows = NamedTuple[]
    t0 = time()
    for (i, r) in enumerate(sweep_rows)
        rt = retag_point(r.theta, r.phi, r.lambda)
        push!(rows, merge(rt, (
            orig_n_pool = r.n_pool,
            orig_max_clique = r.max_clique,
            orig_found_third = r.found_third,
            orig_found_fourth = r.found_fourth,
            homotopy_used = r.homotopy_used,
            fresh_fallback = r.fresh_fallback,
            drop_suspected = r.drop_suspected,
        )))
        i % 16 == 0 && @printf("  ... %d/%d  elapsed %.0fs\n", i, length(sweep_rows), time() - t0)
    end

    write_retag_csv(OUT_CSV, rows)

    old_complete = count(r -> r.old_pool_complete, rows)
    new_complete = count(r -> r.new_pool_complete, rows)
    flip_c2i = count(r -> r.old_pool_complete && !r.new_pool_complete, rows)
    flip_i2c = count(r -> !r.old_pool_complete && r.new_pool_complete, rows)
    max_clique_ge3 = count(r -> r.orig_max_clique >= 3, rows)

    c2i_rows = filter(r -> r.old_pool_complete && !r.new_pool_complete, rows)
    c2i_max_clique = [r.orig_max_clique for r in c2i_rows]

    open(SUMMARY_TXT, "w") do io
        println(io, "=== Task 2: 512-point pool completeness re-tag ===")
        @printf(io, "Points processed: %d\n\n", length(rows))
        @printf(io, "OLD complete (tracked >= 0.9*mv): %d\n", old_complete)
        @printf(io, "NEW complete (n_certified==n_tracked): %d\n", new_complete)
        @printf(io, "Flip complete->incomplete: %d\n", flip_c2i)
        @printf(io, "Flip incomplete->complete: %d\n", flip_i2c)
        @printf(io, "Points with max_clique>=3: %d\n\n", max_clique_ge3)

        if !isempty(c2i_rows)
            println(io, "Complete->incomplete points (re-check max_clique trust):")
            for r in c2i_rows[1:min(20, end)]
                @printf(io, "  (%.4f, %.4f, %.4f) tracked=%d cert=%d mv=%d max_clique=%.0f drop=%s\n",
                        r.theta, r.phi, r.lambda, r.n_tracked, r.n_certified, r.mixed_volume,
                        r.orig_max_clique, r.drop_suspected)
            end
            println(io, "\n  FLAG: these need fresh/monodromy re-solve for clique trust.")
        end
    end

    @printf("\n=== Summary ===\n")
    @printf("OLD complete: %d/%d\n", old_complete, length(rows))
    @printf("NEW complete: %d/%d\n", new_complete, length(rows))
    @printf("Flip c->i: %d  i->c: %d\n", flip_c2i, flip_i2c)
    @printf("max_clique>=3: %d\n", max_clique_ge3)
    println("Wrote $OUT_CSV and $SUMMARY_TXT")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
