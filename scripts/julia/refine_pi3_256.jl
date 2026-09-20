# ============================================================================
# refine_pi3_256.jl — bounded precision experiment at the Dita point λ = π/3.
#
# Stage 1: Float64 polyhedral HC solve (the audited pipeline, unchanged).
# Stage 2: genuine 256-bit Newton refinement of every recovered pool vector
#          against H(θ_D, π/4, π/3) computed in Complex{BigFloat} (256-bit).
#
# NOT 256-bit homotopy path tracking: discovery is Float64; only refinement
# is high precision. No W1 witness, no exact Groebner, no completeness claim.
#
# Question answered: does the recovered-pool six-clique structure at λ = π/3
# survive 256-bit refinement and orthogonality thresholds far below 1e-12?
#
# Output: results/pi3_precision_256.txt
# ============================================================================
using Dates, Printf, LinearAlgebra, Random
const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(ROOT, "src", "MubSearch.jl"))
using .MubSearch
using Graphs

println("start=$(now())"); flush(stdout)

# ---------- 256-bit Karlsson matrix at (θ_D, π/4, π/3) ----------
function dita_pi3_matrix()
    ii = Complex{BigFloat}(0, 1)
    one_c = one(Complex{BigFloat})
    θ = acos(BigFloat(1) / sqrt(BigFloat(3)))
    φ = big(pi) / 4
    λ = big(pi) / 3
    s, c = sin(θ), cos(θ)
    A11 = -one_c / 2 + ii * (sqrt(BigFloat(3)) / 2) * (c + exp(-ii * φ) * s)
    A12 = -one_c / 2 + ii * (sqrt(BigFloat(3)) / 2) * (-c + exp(ii * φ) * s)
    A = [A11 A12; conj(A12) -conj(A11)]
    F2 = Complex{BigFloat}.(BigFloat[1 1; 1 -1])
    Bm = -F2 - A
    αA = A[1, 2]^2;  βA = A[1, 1]^2
    αB = Bm[1, 2]^2; βB = Bm[1, 1]^2
    mob(z, al, be) = (al * z - be) / (conj(be) * z - conj(al))
    z1sq = exp(ii * 2 * λ)
    z3sq = mob(z1sq, αA, βA)
    z4sq = mob(z1sq, αB, βB)
    z2sq = (βB - z3sq * conj(αB)) / (αB - z3sq * conj(βB))
    z1 = exp(ii * λ)
    z2, z3, z4 = sqrt(z2sq), sqrt(z3sq), sqrt(z4sq)
    L(z) = [one_c one_c; z -z]
    R(z) = [one_c z; one_c -z]
    half = one_c / 2
    top = hcat(F2, L(z1), L(z2))
    mid = hcat(R(z3), half .* (R(z3) * A * L(z1)), half .* (R(z3) * Bm * L(z2)))
    bot = hcat(R(z4), half .* (R(z4) * Bm * L(z1)), half .* (R(z4) * A * L(z2)))
    return vcat(top, mid, bot)
end

# ---------- 256-bit Newton refinement on the phase gauge ----------
function eqs_jac(a, H)
    z = vcat(one(Complex{BigFloat}), cis.(a))
    f = zeros(BigFloat, 5)
    J = zeros(BigFloat, 5, 5)
    for k in 1:5
        s = sum(conj(H[j, k]) * z[j] for j in 1:6)
        f[k] = abs2(s) - 6
        for j in 1:5
            J[k, j] = 2 * real(conj(s) * (Complex{BigFloat}(0, 1) * conj(H[j+1, k]) * z[j+1]))
        end
    end
    return f, J
end

function refine(v, H; target = big"1e-70", maxit = 60)
    a = angle.(Complex{BigFloat}.(v[2:6] ./ v[1]))
    for it in 1:maxit
        f, J = eqs_jac(a, H)
        maximum(abs, f) < target &&
            return vcat(one(Complex{BigFloat}), cis.(a)) ./ sqrt(BigFloat(6)), it - 1
        a -= J \ f
    end
    error("Newton did not reach $target in $maxit iterations")
end

function mu_residuals(pool, H)
    flat = maximum(abs(abs2(v[j]) - BigFloat(1) / 6) for v in pool for j in 1:6)
    mu = maximum(abs(abs2(dot(H[:, k], v)) / 6 - BigFloat(1) / 6) for v in pool for k in 1:6)
    return flat, mu
end

function clique_stats(pool, tol)
    g = SimpleGraph(length(pool))
    for j in eachindex(pool), k in 1:j-1
        abs(dot(pool[j], pool[k])) < tol && add_edge!(g, j, k)
    end
    cliques = maximal_cliques(g)
    six = [sort!(collect(c)) for c in cliques if length(c) == 6]
    worst = isempty(six) ? BigFloat(NaN) :
            maximum(abs(dot(pool[c[j]], pool[c[k]])) for c in six for j in 1:6 for k in 1:j-1)
    return length(six), maximum(length, cliques), worst, six
end

function main()
    # Build the 256-bit reference matrix, then DROP back to Float64 default
    # precision for the HC solve (HC internals slow drastically at 256-bit
    # global precision; discovery is Float64 by design here).
    H = setprecision(BigFloat, 256) do
        dita_pi3_matrix()
    end
    hdef = maximum(abs, H * H' - 6 * I(6))
    @assert hdef < big"1e-70" "256-bit H failed unitarity: $hdef"
    H64 = ComplexF64.(H)
    hgap = maximum(abs, H - Complex{BigFloat}.(H64))
    println("H256 unitarity defect = $hdef;  |H256 - H64| = $hgap"); flush(stdout)
    @assert hgap < 1e-13 "branch mismatch between 256-bit and Float64 construction"

    Random.seed!(20260917)
    println("Float64 polyhedral solve ..."); flush(stdout)
    out = solve_pool_audited(H64)
    a = out.audit
    println("audit: mv=$(a.mixed_volume) attempted=$(a.paths_attempted) sols=$(a.numerical_solutions) cert=$(a.certified_solutions) conj=$(a.conjugate_locus_solutions) valid=$(a.valid_mu_vectors) dedup=$(a.deduplicated_vectors)")
    flush(stdout)
    # Completeness gate = the repo's CORRECTED gate (certify-to-track):
    # every tracked path certified. mv (252) is only a generic upper bound;
    # 240 < 252 here is degeneracy, not lost paths (cf. docs/ENGINEERING.md §4).
    @assert a.numerical_solutions == a.certified_solutions &&
            a.certified_solutions > 0 && a.valid_mu_vectors > 0 "incomplete tracking - no claim"

    # ---- everything below runs at 256-bit global precision ----
    setprecision(BigFloat, 256) do

        pool64 = out.pool
        pool_big = [Complex{BigFloat}.(v) for v in pool64]
        refined = Vector{Vector{Complex{BigFloat}}}()
        maxsteps = 0
        for v in pool64
            w, st = refine(v, H)
            push!(refined, w)
            maxsteps = max(maxsteps, st)
        end
        println("refined $(length(refined)) vectors (max Newton steps $maxsteps)"); flush(stdout)

        minsep = minimum(norm(refined[j] - refined[k]) for j in eachindex(refined) for k in 1:j-1)
        @assert minsep > big"1e-20" "refinement collapsed distinct vectors: $minsep"

        flat0, mu0 = mu_residuals(pool_big, H)
        flat1, mu1 = mu_residuals(refined, H)
        println("residuals after refinement: flat=$flat1  mu=$mu1"); flush(stdout)
        # per-vector worst-offender diagnostics (reported, not just asserted)
        worst_v, worst_flat, worst_mu = 0, BigFloat(0), BigFloat(0)
        for (i, w) in enumerate(refined)
            fi = maximum(abs(abs2(w[j]) - BigFloat(1) / 6) for j in 1:6)
            mi = maximum(abs(abs2(dot(H[:, k], w)) / 6 - BigFloat(1) / 6) for k in 1:6)
            fi > worst_flat && (worst_v, worst_flat = i, fi)
            mi > worst_mu && (worst_mu = mi)
        end
        println("worst per-vector: v#$(worst_v) flat=$(worst_flat) mu=$(worst_mu)"); flush(stdout)
        @assert flat1 < big"1e-60" && mu1 < big"1e-60" "refinement did not reach 1e-60"


        lines = String[]
        push!(lines, "=== lambda=pi/3 precision experiment (Float64 discovery + 256-bit refinement) ===")
        push!(lines, "timestamp=$(now())  Julia=$(VERSION)  seed=20260917  precision=256 bits")
        push!(lines, "method: Float64 polyhedral HC solve (audited pipeline), then 256-bit phase-gauge")
        push!(lines, "Newton refinement against the 256-bit Karlsson matrix. NOT 256-bit path tracking.")
        push!(lines, "H256_unitarity_defect = $hdef")
        push!(lines, "H256_vs_H64_max_diff  = $hgap")
        push!(lines, "pool_audit: mv=$(a.mixed_volume) attempted=$(a.paths_attempted) solutions=$(a.numerical_solutions) certified=$(a.certified_solutions) conjugate=$(a.conjugate_locus_solutions) valid=$(a.valid_mu_vectors) dedup=$(a.deduplicated_vectors)")
        push!(lines, "Float64_pipeline_worst_MU_defect = $(out.worst_defect)")
        push!(lines, "residuals_before_refinement: flat=$(flat0)  mu=$(mu0)")
        push!(lines, "residuals_after_refinement:  flat=$(flat1)  mu=$(mu1)")
        push!(lines, "min_refined_separation = $minsep;  max_Newton_steps = $maxsteps")
        push!(lines, "")
        for tol in (big"1e-8", big"1e-12", big"1e-14", big"1e-20", big"1e-40")
            c0, m0, e0, idx0 = clique_stats(pool_big, tol)
            c1, m1, r1, idx1 = clique_stats(refined, tol)
            push!(lines, @sprintf("tol=%s  six-cliques: float=%d refined=%d  max_clique: %d/%d  worst_clique_edge: %.3e / %.3e  index_sets_identical=%s",
                                   tol, c0, c1, m0, m1, Float64(e0), Float64(r1), idx0 == idx1))
        end
        push!(lines, "")
        push!(lines, "scope: numerical stability of the recovered pool at lambda=pi/3 only.")
        push!(lines, "No W1 emptiness, no exact elimination, no pool-completeness certificate,")
        push!(lines, "no statement about F6 arc boundaries, and no T3 dependence.")
        open(joinpath(ROOT, "results", "pi3_precision_256.txt"), "w") do io
            for ln in lines; println(io, ln); end
        end
        println(); println.(lines); flush(stdout)
        println("done=$(now())"); flush(stdout)
    end
end

main()

