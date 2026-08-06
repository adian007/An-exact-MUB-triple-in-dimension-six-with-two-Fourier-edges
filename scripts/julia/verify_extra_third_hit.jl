# One-off: verify third-MUB hit from special_loci_search.csv row 78 (post subagent run).
include(joinpath(@__DIR__, "_paths.jl"))
using Printf

const TH = 1.5707963268
const PH = 0.05
const LAM = 6.2331853072

H = build_karlsson_family(TH, PH, LAM)
pc = pool_completeness_report(H; verbose = false)
pool, _ = generate_candidate_pool_fresh(H; tol = 1e-12, verbose = false)
pool = deduplicate_pool(pool)
ext = check_four_mub_extension(pool, H; mu_tol = 1e-12, ortho_tol = 1e-12)

@printf("point: theta=%.10f phi=%.10f lambda=%.10f\n", TH, PH, LAM)
@printf("n_tracked=%d n_certified=%d pool_complete=%s\n",
        pc.n_tracked, pc.n_certified, pc.pool_complete_flag)
@printf("max_clique=%d found_third=%s found_fourth=%s n_pool=%d\n",
        ext.max_clique, ext.found_third, ext.found_fourth, ext.n_pool)
