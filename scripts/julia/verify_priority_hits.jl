# Re-verify priority clique>=6 hits from third_mub_locus_sampling at higher precision.
include(joinpath(@__DIR__, "_paths.jl"))
using Printf

const DITA_THETA = acos(1 / sqrt(3))

function full_probe(label, theta, phi, lam; mu_tol = 1e-12, ortho_tol = 1e-12, hp_bits = 400)
    println("\n=== $label ===")
    @printf("  (θ,φ,λ) = (%.17g, %.17g, %.17g)\n", theta, phi, lam)
    H = build_karlsson_family(theta, phi, lam)
    @printf("  Hadamard: %s\n", is_hadamard(H; tol = 1e-8))
    pc = pool_completeness_report(H; verbose = false)
    @printf("  pool: tracked=%d certified=%d complete=%s dedup=%d\n",
            pc.n_tracked, pc.n_certified, pc.pool_complete_flag, pc.n_dedup_pool)
    pool, _ = generate_candidate_pool_fresh(H; verbose = false)
    pool = deduplicate_pool(pool; tol = 1e-10)
    ext = check_four_mub_extension(pool, H; mu_tol = mu_tol, ortho_tol = ortho_tol)
    @printf("  max_clique=%d found_third=%s found_fourth=%s mu_min=%.2e\n",
            ext.max_clique, ext.found_third, ext.found_fourth, ext.mu_defect_min)
    if ext.found_third
        g = _orthogonality_graph(pool; ortho_tol = ortho_tol)
        cliques = [c for c in maximal_cliques(g) if length(c) >= 6]
        if !isempty(cliques)
            best_c = argmax(c -> length(c), cliques)
            hp = verify_clique_hp([pool[i] for i in best_c]; bits = hp_bits)
            @printf("  HP(%d-bit): ortho_max=%.3e mu_max=%.3e\n", hp_bits, hp.ortho_max, hp.mu_norm_max)
        end
    end
    return ext.max_clique
end

# Priority hits at |delta|=0.01
full_probe("Dita +lambda +0.01", DITA_THETA, pi / 4, 0.4 + 0.01)
full_probe("Dita +lambda -0.01", DITA_THETA, pi / 4, 0.4 - 0.01)
full_probe("F6_theta0 +phi +0.01", 0.0, 0.5 + 0.01, 0.3)
full_probe("F6_theta0 +phi -0.01", 0.0, 0.5 - 0.01, 0.3)
full_probe("F6_theta0 +lambda +0.01", 0.0, 0.5, 0.3 + 0.01)
full_probe("F6_theta0 +lambda -0.01", 0.0, 0.5, 0.3 - 0.01)
full_probe("F6_theta0 phi+lambda (+0.01,+0.01)", 0.0, 0.5 + 0.01, 0.3 + 0.01)

# Controls: isotropic-scale points that should be clique=2
full_probe("Dita +theta +0.01", DITA_THETA + 0.01, pi / 4, 0.4)
full_probe("F6_theta0 +theta +0.01", 0.01, 0.5, 0.3)

# Finer sweep on lambda axis near Dita where hit was reported
for d in [1e-5, 1e-4, 1e-3, 5e-3, 1e-2, 2e-2]
    full_probe("Dita +lambda +$d", DITA_THETA, pi / 4, 0.4 + d)
end
