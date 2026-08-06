# investigate_candidate core logic without CSV dependency (Task 1b)
include(joinpath(@__DIR__, "_paths.jl"))
using Printf

const THETA = 0.9553166181
const PHI = 0.7853981634
const LAM = 0.4

function pslq_params(theta, phi, lam; atol = 1e-6)
    relations = String[]
    for denom in 1:12, num in -24:24
        v = num * pi / denom
        abs(theta - v) < atol && push!(relations, "theta ≈ $(num)*pi/$denom")
        abs(phi - v) < atol && push!(relations, "phi ≈ $(num)*pi/$denom")
        abs(lam - v) < atol && push!(relations, "lambda ≈ $(num)*pi/$denom")
    end
    dita = acos(1 / sqrt(3))
    abs(theta - dita) < atol && push!(relations, "theta = arccos(1/sqrt(3)) [Dita]")
    abs(theta) < atol && push!(relations, "theta = 0 [Fourier slice]")
    return unique(relations)
end

println("=== investigate_candidate (CSV-free) ===")
@printf("  (θ,φ,λ) = (%.16f, %.16f, %.16f)\n", THETA, PHI, LAM)

H = build_karlsson_family(THETA, PHI, LAM)
pc = pool_completeness_report(H; verbose = true)

println("\n[2] Re-certify max_precision=128")
system = build_numeric_pool_system(H)
res = solve(system; show_progress = false)
cert = certify(system, res; show_progress = false, threading = false, max_precision = 128)
@printf("  certified=%d tracked=%d\n", ncertified(cert), length(solutions(res)))

pool, stats = generate_candidate_pool_fresh(H; tol = 1e-12, verbose = true)
pool = deduplicate_pool(pool; tol = 1e-10)
ext = check_four_mub_extension(pool, H; mu_tol = 1e-12, ortho_tol = 1e-12)
@printf("  max_clique=%d found_third=%s found_fourth=%s mu_min=%.2e\n",
        ext.max_clique, ext.found_third, ext.found_fourth, ext.mu_defect_min)

println("\n[HP] BigFloat bits=400 verify")
if ext.found_third
    g = _orthogonality_graph(pool; ortho_tol = 1e-12)
    cliques = [c for c in maximal_cliques(g) if length(c) >= 6]
    best_c = argmax(c -> length(c), cliques)
    hp = verify_clique_hp([pool[i] for i in best_c]; bits = 400)
    @printf("  ortho_max=%.3e mu_max=%.3e\n", hp.ortho_max, hp.mu_norm_max)
else
    println("  (no 6-clique to HP-verify)")
end

println("\n[PSLQ]")
for r in pslq_params(THETA, PHI, LAM)
    println("  $r")
end

# Also test exact Dita for comparison
println("\n--- Exact Dita reference ---")
th = acos(1/sqrt(3)); ph = pi/4; lam = 0.4
H2 = build_karlsson_family(th, ph, lam)
pool2, _ = generate_candidate_pool_fresh(H2; verbose = false)
pool2 = deduplicate_pool(pool2)
ext2 = check_four_mub_extension(pool2, H2; mu_tol = 1e-12, ortho_tol = 1e-12)
@printf("  exact: max_clique=%d found_third=%s found_fourth=%s n_pool=%d\n",
        ext2.max_clique, ext2.found_third, ext2.found_fourth, length(pool2))
