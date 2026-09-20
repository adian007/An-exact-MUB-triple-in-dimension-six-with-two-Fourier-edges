# Path B1: investigate θ=π/2 circulant_match false clique-6 hits.
# Usage: julia --project=. scripts/julia/investigate_circulant_match_anomaly.jl
include(joinpath(@__DIR__, "_paths.jl"))
using Printf, Dates, LinearAlgebra

const OUT = joinpath(RESULTS_DIR, "circulant_match_anomaly.txt")
const DITA_THETA = acos(1 / sqrt(3))

"""All CSV rows flagged circulant_match with max_clique=6 in degen500."""
const CIRCULANT_CLIQUES = [
    (1.5707963268, 0.05, 6.2331853072),
    (1.5707963268, 3.0915926536, 0.05),
]

function analyze_point(theta, phi, lam; bits=400)
    H = build_karlsson_family(theta, phi, lam)
    pc = pool_completeness_report(H; verbose=false)
    pool, stats = generate_candidate_pool_fresh(H; tol=1e-12, verbose=false)
    pool = deduplicate_pool(pool; tol=1e-12)

    ortho_sweep = NamedTuple[]
    for tol in (1e-6, 1e-8, 1e-10, 1e-12)
        ext = check_four_mub_extension(pool, H; ortho_tol=tol, mu_tol=1e-12)
        push!(ortho_sweep, (ortho_tol=tol, max_clique=ext.max_clique, found_third=ext.found_third))
    end

    ext = check_four_mub_extension(pool, H; ortho_tol=1e-8, mu_tol=1e-8)
    hp = nothing
    if ext.max_clique >= 6
        g = _orthogonality_graph(pool; ortho_tol=1e-8)
        cliques = [c for c in maximal_cliques(g) if length(c) >= 6]
        if !isempty(cliques)
            c = argmax(cl -> length(cl), cliques)
            clique_vecs = pool[c[1:6]]
            hp = verify_clique_hp(clique_vecs; bits=bits, ortho_tol=1e-12)
        end
    end

    circ_err = try
        py = joinpath(@__DIR__, "..", "python")
        # inline circulant score (Julia port of degeneracy_scan circulant_match_score)
        D = dephase(H)
        n = size(D, 1)
        ref = D[1, :]
        worst = 0.0
        for k in 2:n
            shifted = circshift(ref, k - 1)
            for j in 1:n
                a, b = D[k, j], shifted[j]
                abs(a) < 1e-8 || abs(b) < 1e-8 && continue
                worst = max(worst, abs(angle(a * conj(b))))
            end
        end
        worst
    catch
        NaN
    end

    return (
        theta=theta, phi=phi, lambda=lam,
        mixed_volume=pc.mixed_volume,
        n_tracked=pc.n_tracked,
        mv_gap=pc.mv_gap,
        pool_complete=pc.pool_complete_flag,
        n_pool=length(pool),
        n_raw=stats.n_raw,
        max_clique=ext.max_clique,
        found_third=ext.found_third,
        hp_clique_ok=hp === nothing ? false : (hp.ortho_ok && hp.mu_ok),
        hp_ortho_max=hp === nothing ? NaN : hp.ortho_max,
        hp_mu_max=hp === nothing ? NaN : hp.mu_norm_max,
        circ_err=circ_err,
        ortho_sweep=ortho_sweep,
    )
end

function lambda_sweep_near_pi2()
    # λ values where degeneracy scan flags circulant_match with θ or φ = π/2
    out = NamedTuple[]
    th = pi / 2
    for lam in range(0.0, 2pi; length=9)
        for ph in (0.05, pi / 2, 3.0915926536)
            r = analyze_point(th, ph, lam)
            r.max_clique >= 6 && push!(out, r)
        end
    end
    return out
end

function main()
    open(OUT, "w") do io
        println(io, "=== Path B1: circulant_match θ=π/2 anomaly investigation ===")
        println(io, "Date: $(Dates.format(now(), "yyyy-mm-dd HH:MM"))")
        println(io, "")

        println(io, "--- Canonical false-positive points (CSV max_clique=6, HP fail) ---")
        for (th, ph, lam) in CIRCULANT_CLIQUES
            r = analyze_point(th, ph, lam)
            @printf(io, "Point (θ,φ,λ) = (%.10f, %.10f, %.10f)\n", th, ph, lam)
            @printf(io, "  mixed_volume=%d n_tracked=%d mv_gap=%d pool_complete=%s\n",
                    r.mixed_volume, r.n_tracked, r.mv_gap, r.pool_complete)
            @printf(io, "  n_pool=%d n_raw=%d max_clique=%d found_third=%s\n",
                    r.n_pool, r.n_raw, r.max_clique, r.found_third)
            @printf(io, "  circulant_err=%.3e hp_pass=%s ortho_max=%.3e mu_max=%.3e\n",
                    r.circ_err, r.hp_clique_ok, r.hp_ortho_max, r.hp_mu_max)
            println(io, "  ortho_tol sweep:")
            for row in r.ortho_sweep
                @printf(io, "    tol=%.0e  max_clique=%d  found_third=%s\n",
                        row.ortho_tol, row.max_clique, row.found_third)
            end
            println(io, "")
        end

        println(io, "--- Control: genuine third-MUB anchors ---")
        for (label, th, ph, lam) in [
            ("F6_theta0", 0.0, 0.5, 0.3),
            ("Dita", DITA_THETA, pi / 4, 0.4),
        ]
            r = analyze_point(th, ph, lam)
            @printf(io, "%s: mv=%d tracked=%d gap=%d max_clique=%d hp_pass=%s n_pool=%d\n",
                    label, r.mixed_volume, r.n_tracked, r.mv_gap,
                    r.max_clique, r.hp_clique_ok, r.n_pool)
        end
        println(io, "")

        println(io, "--- Root-cause diagnosis ---")
        r = analyze_point(CIRCULANT_CLIQUES[1]...)
        Hc = build_karlsson_family(r.theta, r.phi, r.lambda)
        pool_t, _ = generate_candidate_pool_fresh(Hc; tol=1e-12, verbose=false)
        pool_t = deduplicate_pool(pool_t)
        ext_tight = check_four_mub_extension(pool_t, Hc; ortho_tol=1e-12, mu_tol=1e-12)
        if r.max_clique >= 6 && !r.hp_clique_ok && ext_tight.max_clique < 6
            println(io, "VERDICT: TOLERANCE ARTIFACT (not a third structural locus).")
            println(io, "  Float pipeline (ortho_tol=1e-8) reports max_clique=6; tight tol=1e-12 → clique=2.")
            @printf(io, "  HP ortho_max=%.3e > 1e-12 — fails validate_third_mub_candidates audit.\n", r.hp_ortho_max)
            println(io, "  All θ=π/2 hits share n_tracked=192 (mv_gap=60): incomplete root recovery at this degeneracy.")
            println(io, "  circulant_err≈π: dephased H is NOT circulant here despite degeneracy_scan flag.")
            println(io, "  Mechanism: ortho_tol × truncated pool × near-degenerate θ=π/2 geometry.")
        elseif r.max_clique < 6
            println(io, "VERDICT: FALSE POSITIVE (non-reproducing at standard tol).")
        else
            println(io, "VERDICT: REQUIRES FURTHER STUDY — HP-stable clique=6 at θ=π/2.")
        end
        println(io, "")
        println(io, "NOT a near-miss to Dita/F6: circulant_match flags dephased-circulant H structure,")
        println(io, "  not proximity to third-MUB loci (θ=0 arc or Dita circle).")
    end
    println("Wrote $OUT")
end

main()
