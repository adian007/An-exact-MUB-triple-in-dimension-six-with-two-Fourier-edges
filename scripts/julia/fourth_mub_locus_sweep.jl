# Phase 1.3c: Fourth-MUB sweep on interior locus points (>=20).
# STOP RULE: halt with exit code 2 if found_fourth=true anywhere.
# Usage: julia --project=. scripts/julia/fourth_mub_locus_sweep.jl [--quick]
include(joinpath(@__DIR__, "_paths.jl"))
using Printf, Dates, LinearAlgebra

const DITA_THETA = acos(1 / sqrt(3))
const OUT = joinpath(RESULTS_DIR, "fourth_mub_locus_sweep.txt")
const HP_BITS = 400
const MU_TOL = 1e-12
const ORTHO_TOL = 1e-12
const QUICK = "--quick" in ARGS

"""Interior points on Dita lambda curve and F6 theta=0 lambda curve."""
function locus_points()
    pts = NamedTuple[]
    # Dita: spread lambda along axis where clique-6 persists
    lam_vals = QUICK ? [0.35, 0.4, 0.45] : collect(range(0.25, 0.55; length = 12))
    for lam in lam_vals
        push!(pts, (name = "Dita", theta = DITA_THETA, phi = pi / 4, lambda = lam))
    end
    # F6 theta=0: lambda axis (phi gauge)
    lam_f6 = QUICK ? [0.25, 0.3, 0.35] : collect(range(0.22, 0.38; length = 10))
    for lam in lam_f6
        push!(pts, (name = "F6_theta0", theta = 0.0, phi = 0.5, lambda = lam))
    end
    return pts
end

function sweep_one(pt)
    H = build_karlsson_family(pt.theta, pt.phi, pt.lambda)
    !is_hadamard(H; tol = 1e-8) && return (
        ok = false, max_clique = 0, found_third = false, found_fourth = false,
        hp_ok = false, pool_complete = false,
    )
    pool, _ = generate_candidate_pool_fresh(H; verbose = false)
    pool = deduplicate_pool(pool)
    pc = pool_completeness_report(H; verbose = false)
    ext = check_four_mub_extension(pool, H; mu_tol = MU_TOL, ortho_tol = ORTHO_TOL)
    hp_ok = true
    if ext.found_third && ext.max_clique >= 6
        g = _orthogonality_graph(pool; ortho_tol = ORTHO_TOL)
        cliques = [c for c in maximal_cliques(g) if length(c) >= 6]
        if !isempty(cliques)
            best_c = argmax(c -> length(c), cliques)
            hp = verify_clique_hp([pool[i] for i in best_c]; bits = HP_BITS)
            hp_ok = hp.ortho_max < 1e-8 && hp.mu_norm_max < 1e-8
        end
    end
    per_basis = fourth_mub_per_basis_report(pool, H; ortho_tol = ORTHO_TOL, mu_tol = MU_TOL)
    any_fourth = any(r.found_fourth for r in per_basis)
    return (
        ok = true,
        max_clique = ext.max_clique,
        found_third = ext.found_third,
        found_fourth = ext.found_fourth || any_fourth,
        hp_ok = hp_ok,
        pool_complete = pc.pool_complete_flag,
        n_pool = length(pool),
        n_third_bases = length(per_basis),
    )
end

function main()
    pts = locus_points()
    @assert length(pts) >= (QUICK ? 3 : 20) "Need >=20 interior points (use full run)"

    lines = String[]
    push!(lines, "=== Fourth-MUB locus sweep (Phase 1.3c) ===")
    push!(lines, "Date: $(Dates.now())  n_points=$(length(pts))  HP_bits=$HP_BITS")
    push!(lines, "")

    found_any_fourth = false
    n_third = 0
    for (i, pt) in enumerate(pts)
        r = sweep_one(pt)
        @printf("[%d/%d] %s (%.4f, %.4f, %.4f) clique=%d third=%s fourth=%s hp=%s\n",
            i, length(pts), pt.name, pt.theta, pt.phi, pt.lambda,
            r.max_clique, r.found_third, r.found_fourth, r.hp_ok)
        push!(lines, @sprintf("%s  theta=%.6g phi=%.6g lambda=%.6g  clique=%d third=%s fourth=%s hp=%s pool_complete=%s",
            pt.name, pt.theta, pt.phi, pt.lambda, r.max_clique, r.found_third,
            r.found_fourth, r.hp_ok, r.pool_complete))
        r.found_third && (n_third += 1)
        if r.found_fourth
            found_any_fourth = true
            push!(lines, "  *** STOP: found_fourth=true — escalate to investigate_candidate.jl ***")
        end
    end

    push!(lines, "")
    push!(lines, "Summary: $(length(pts)) points, $n_third with third MUB, found_fourth=$(found_any_fourth)")

    open(OUT, "w") do io
        for ln in lines
            println(io, ln)
        end
    end
    for ln in lines
        println(ln)
    end
    println("Wrote $OUT")

    if found_any_fourth
        println("\n*** PIPELINE STOP: found_fourth=true — run investigate_candidate.jl ***")
        exit(2)
    end
end

main()
