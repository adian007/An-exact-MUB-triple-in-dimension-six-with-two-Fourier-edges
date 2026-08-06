include(joinpath(@__DIR__, "_paths.jl"))
th = acos(1 / sqrt(3))
ph = pi / 4
lam = 0.4
for (label, t) in [("exact Dita", th), ("CSV rounded", 0.9553166181)]
    H = build_karlsson_family(t, ph, lam)
    pool, _ = generate_candidate_pool_fresh(H; verbose = false)
    pool = deduplicate_pool(pool)
    ext = check_four_mub_extension(pool, H; ortho_tol = 1e-12, mu_tol = 1e-12)
    println("$label: max_clique=$(ext.max_clique) found_third=$(ext.found_third) found_fourth=$(ext.found_fourth) n_pool=$(length(pool))")
end
