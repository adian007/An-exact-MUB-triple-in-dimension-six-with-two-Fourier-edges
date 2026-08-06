include(joinpath(@__DIR__, "_paths.jl"))
using Printf

anchors = [
    ("F6_theta0", 0.0, 0.5, 0.3),
    ("Dita", acos(1/sqrt(3)), pi/4, 0.4),
]

for (name, th, ph, lm) in anchors
    println("\n=== $name ===")
    H = build_karlsson_family(th, ph, lm)
    pool, = generate_candidate_pool_fresh(H; verbose=false)
    pool = deduplicate_pool(pool)
    ext = check_four_mub_extension(pool, H)
    @printf("n_pool=%d max_clique=%d third=%s fourth=%s\n",
            length(pool), ext.max_clique, ext.found_third, ext.found_fourth)
    for rep in fourth_mub_per_basis_report(pool, H)
        @printf("  basis %d: n_filtered=%d max_clique=%d found_fourth=%s\n",
                rep.basis_idx, rep.n_filtered, rep.max_clique_filtered, rep.found_fourth)
    end
    g = _orthogonality_graph(pool; ortho_tol=1e-12)
    cliques = [c for c in maximal_cliques(g) if length(c) >= 6]
    best_c = argmax(c -> length(c), cliques)
    hp = verify_clique_hp([pool[i] for i in best_c])
    @printf("  HP clique: ortho_max=%.3e mu_max=%.3e ok=%s\n",
            hp.ortho_max, hp.mu_norm_max, hp.ortho_ok && hp.mu_ok)
end
