using Dates
println("start=$(now()) Julia=$(VERSION)"); flush(stdout)
using LinearAlgebra, Random
const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(ROOT,"src","MubSearch.jl"))
using .MubSearch
println("module_loaded=$(now())"); flush(stdout)
H = build_karlsson_family(DITA_THETA, DITA_PHI, pi/3)
println("H_defect=$(maximum(abs,H*H'-6I))"); flush(stdout)
Random.seed!(20260917)
out = solve_pool_audited(H)
println("audit=$(out.audit); worst_MU_defect=$(out.worst_defect)"); flush(stdout)
for tol in (1e-8,1e-12,1e-14)
    c = enumerate_all_third_mub_bases(out.pool,H; ortho_tol=tol,mu_tol=tol)
    println("tol=$tol third_bases=$(c.n_verified) raw=$(c.n_raw)")
end
println("complete=$(now()); numerical recovered-pool test only")
