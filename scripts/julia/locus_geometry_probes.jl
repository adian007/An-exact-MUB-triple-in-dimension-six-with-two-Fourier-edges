# L5-L6 locus geometry probes: semicontinuity boundaries at Dita and F6 theta=0.
# Usage: julia --project=. scripts/julia/locus_geometry_probes.jl

include(joinpath(@__DIR__, "_paths.jl"))

using Printf
using Dates

const OUT = joinpath(RESULTS_DIR, "locus_geometry_probes.txt")
const DITA_THETA = acos(1 / sqrt(3))
const HP_BITS = 256

function probe_point(theta, phi, lam)
    H = build_karlsson_family(theta, phi, lam)
    pool, = generate_candidate_pool_fresh(H; verbose=false)
    pool = deduplicate_pool(pool)
    ext = check_four_mub_extension(pool, H; mu_tol=1e-10, ortho_tol=1e-8)
    mc = ext.max_clique
    hp_ok = false
    if mc >= 6
        g = _orthogonality_graph(pool; ortho_tol=1e-8)
        cliques = [c for c in maximal_cliques(g) if length(c) >= 6]
        if !isempty(cliques)
            vecs = [pool[i] for i in cliques[1][1:6]]
            hp = verify_clique_hp(vecs; bits=HP_BITS)
            hp_ok = hp.ortho_ok && hp.mu_ok
        end
    end
    return mc, hp_ok
end

function dita_phi_boundary_probes()
    lam_vals = [0.3, 0.4, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
    deltas = [0.0, 0.001]
    results = NamedTuple[]
    for lam in lam_vals, d in deltas
        mc, hp = probe_point(DITA_THETA, pi / 4 + d, lam)
        push!(results, (lam=lam, delta_phi=d, max_clique=mc, hp=hp))
    end
    return results
end

function f6_lambda_boundary_probes()
    lam0 = 0.3
    d_plus = 0.055044803035
    d_minus = 0.0625984193874
    probes = [
        ("interior", lam0),
        ("plus_boundary_in", lam0 + d_plus - 1e-6),
        ("plus_boundary_out", lam0 + d_plus + 1e-6),
        ("minus_boundary_in", lam0 - d_minus + 1e-6),
        ("minus_boundary_out", lam0 - d_minus - 1e-6),
    ]
    results = NamedTuple[]
    for (label, lam) in probes
        mc, hp = probe_point(0.0, 0.5, lam)
        push!(results, (label=label, lambda=lam, max_clique=mc, hp=hp))
    end
    return results
end

function main()
    println("Running Dita phi boundary probes...")
    dita = dita_phi_boundary_probes()
    println("Running F6 lambda boundary probes...")
    f6 = f6_lambda_boundary_probes()

    open(OUT, "w") do io
        println(io, "=== Locus geometry probes (L5-L6) ===")
        println(io, "Date: $(Dates.format(now(), "yyyy-mm-ddTHH:MM:SS"))  HP_bits=$HP_BITS\n")

        println(io, "--- L5: Dita phi-offset ---")
        clique6_at_slice = drop_at_001 = 0
        for r in dita
            @printf(io, "  lam=%.1g dphi=%+.4g clique=%d hp=%s\n",
                    r.lam, r.delta_phi, r.max_clique, r.hp)
            r.delta_phi == 0.0 && r.max_clique >= 6 && (clique6_at_slice += 1)
            abs(r.delta_phi - 0.001) < 1e-9 && r.max_clique <= 2 && (drop_at_001 += 1)
        end
        @printf(io, "  clique>=6 at phi=pi/4: %d/9  drop at dphi=0.001: %d/9\n\n",
                clique6_at_slice, drop_at_001)

        println(io, "--- L6: F6 theta=0 lambda arc ---")
        for r in f6
            @printf(io, "  %-20s lam=%.12g clique=%d hp=%s\n",
                    r.label, r.lambda, r.max_clique, r.hp)
        end

        l5_ok = clique6_at_slice == 9 && drop_at_001 == 9
        interior = findfirst(r -> r.label == "interior", f6)
        l6_ok = interior !== nothing && f6[interior].max_clique >= 6
        println(io, "\nVERDICT L5: $(l5_ok ? "PASS" : "PARTIAL")")
        println(io, "VERDICT L6: $(l6_ok ? "PASS" : "CHECK")")
    end
    println("Wrote $OUT")
end

main()
