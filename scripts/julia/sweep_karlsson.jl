# Per-H Karlsson sweep: for each (theta, phi, lambda), solve the pool of vectors
# MU to {I, H(theta,phi,lambda)} and test whether {I, H} extends to four MUBs.
#
# Pool generation: parameter homotopy from F6 anchor with automatic fresh-solve
# fallback when tracking drops solutions (cross_check=true).
#
# Usage:  julia --project=. sweep_karlsson.jl

include(joinpath(@__DIR__, "_paths.jl"))

using Random

const MU_TOLS = (1e-6, 1e-8, 1e-10)

function run_point(tracker, theta, phi, lam; ortho_tol=1e-8)
    H = build_family_matrix(:karlsson, theta, phi, lam)
    had_defect = norm(H * H' - 6 * I(6))
    pool, _, meta = generate_candidate_pool(H; tracker=tracker, cross_check=true, verbose=false)
    results = Dict{Float64, NamedTuple}()
    for mu_tol in MU_TOLS
        ext = check_four_mub_extension(pool, H; ortho_tol=ortho_tol, mu_tol=mu_tol)
        results[mu_tol] = ext
    end
    primary = results[1e-8]
    return (H=H, pool=pool, meta=meta, had_defect=had_defect, primary=primary, by_tol=results)
end

function main()
    println("=== Per-H Karlsson sweep: does {I, H} extend to four MUBs? ===")
    tracker = init_parametric_pool_tracker(verbose=true)

    thetas = range(0.02, pi - 0.02; length=8)
    phis   = range(0.02, pi - 0.02; length=8)
    lams   = range(0.0, 2pi; length=9)[1:8]
    grid = [(t, p, l) for t in thetas for p in phis for l in lams]
    println("Grid: $(length(grid)) points\n")

    mkpath(RESULTS_DIR)
    out = joinpath(RESULTS_DIR, "karlsson_perH_sweep.csv")
    io = open(out, "w")
    println(io, "theta,phi,lambda,n_pool,max_clique,found_third,found_fourth," *
                "mu_defect_min,homotopy_used,fresh_fallback,drop_suspected," *
                "found_fourth_tol1e-6,found_fourth_tol1e-8,found_fourth_tol1e-10")

    n_fourth = 0
    n_third = 0
    n_fresh_fallback = 0
    t0 = time()
    for (i, (theta, phi, lam)) in enumerate(grid)
        r = run_point(tracker, theta, phi, lam)
        p = r.primary
        meta = r.meta
        homotopy_used = meta.method == :homotopy
        fresh_fallback = meta.drop_suspected && length(r.pool) > 0
        fresh_fallback |= meta.drop_suspected && meta.cross_check.only_in_b > 0
        n_fresh_fallback += fresh_fallback ? 1 : 0
        n_third += p.found_third ? 1 : 0
        n_fourth += p.found_fourth ? 1 : 0
        f6 = r.by_tol[1e-6].found_fourth
        f8 = r.by_tol[1e-8].found_fourth
        f10 = r.by_tol[1e-10].found_fourth
        println(io, join([theta, phi, lam, p.n_pool, p.max_clique, p.found_third, p.found_fourth,
                          p.mu_defect_min, homotopy_used, fresh_fallback, meta.drop_suspected,
                          f6, f8, f10], ","))
        if i % 16 == 0
            println("  ... $(i)/$(length(grid))  elapsed $(round(time()-t0;digits=0))s")
        end
    end
    close(io)

    println("\n--- Summary ---")
    println("  points:           $(length(grid))")
    println("  third MUB found:  $(n_third)")
    println("  fourth MUB found: $(n_fourth)")
    println("  fresh fallback:   $(n_fresh_fallback)/$(length(grid))")
    println("  output:           $(out)")
    println("\nRIGOR: heuristic. Per-H pools via homotopy+fresh fallback; not certified complete.")
end

main()
