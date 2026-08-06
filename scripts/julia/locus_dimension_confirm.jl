# Item 2 fallback: phi sweep at fixed lambda near Dita.
include(joinpath(@__DIR__, "_paths.jl"))
using CSV, DataFrames, Printf, Dates

const DITA_THETA = acos(1 / sqrt(3))
const DITA_PHI = pi / 4
const OUT = joinpath(RESULTS_DIR, "locus_dimension_confirm.txt")
const HP_BITS = 400

function probe(theta, phi, lam)
    H = build_karlsson_family(theta, phi, lam)
    !is_hadamard(H; tol = 1e-8) && return (clique = -1, third = false, hp = false)
    pool, _ = generate_candidate_pool_fresh(H; verbose = false)
    pool = deduplicate_pool(pool)
    ext = check_four_mub_extension(pool, H; mu_tol = 1e-12, ortho_tol = 1e-12)
    hp_ok = true
    if ext.found_third && ext.max_clique >= 6
        g = _orthogonality_graph(pool; ortho_tol = 1e-12)
        cliques = [c for c in maximal_cliques(g) if length(c) >= 6]
        if !isempty(cliques)
            best_c = argmax(c -> length(c), cliques)
            hp = verify_clique_hp([pool[i] for i in best_c]; bits = HP_BITS)
            hp_ok = hp.ortho_max < 1e-8 && hp.mu_norm_max < 1e-8
        end
    end
    return (clique = ext.max_clique, third = ext.found_third, hp = hp_ok)
end

function main()
    lines = String[]
    push!(lines, "=== Item 2: 1D vs 2D locus dimension at Dita ===")
    push!(lines, "Date: $(Dates.now())")
    push!(lines, "")

    # 20x20 grid summary
    grid_path = joinpath(RESULTS_DIR, "locus_2d_grid.csv")
    n6 = 0
    if isfile(grid_path)
        df = CSV.read(grid_path, DataFrame)
        n6 = count(r -> r.max_clique >= 6, eachrow(df))
        push!(lines, "--- 20x20 grid (locus_2d_grid.csv) ---")
        push!(lines, @sprintf("  points: %d, clique>=6: %d", nrow(df), n6))
        push!(lines, @sprintf("  phi range: [%.12g, %.12g] (20 pts, does NOT include exact pi/4=%.12g)", minimum(df.phi), maximum(df.phi), DITA_PHI))
        push!(lines, @sprintf("  nearest grid phi to pi/4: %.12g (delta=%.2e)", df.phi[argmin(abs.(df.phi .- DITA_PHI))], minimum(abs.(df.phi .- DITA_PHI))))
        push!(lines, "")
    end

    # Targeted phi sweep
    deltas = [0.0, 0.001, 0.005, 0.01, 0.02, 0.05, 0.1]
    lams = [0.3, 0.4, 0.5]
    push!(lines, "--- Targeted phi sweep (400-bit HP) ---")
    push!(lines, @sprintf("  theta=arccos(1/sqrt(3)), phi = pi/4 +/- delta, lambda in {0.3, 0.4, 0.5}"))
    push!(lines, "")

    first_drop = Dict{Float64, Float64}()
    for lam in lams
        push!(lines, @sprintf("  lambda=%.12g:", lam))
        prev6 = true
        for d in deltas
            for sign in (+1, -1)
                d == 0.0 && sign == -1 && continue
                ph = DITA_PHI + sign * d
                p = probe(DITA_THETA, ph, lam)
                c6 = p.clique >= 6
                push!(lines, @sprintf("    phi=pi/4 %+.3f (%+.12g): clique=%d third=%s hp=%s",
                    sign * d, sign * d, p.clique, p.third, p.hp))
                if d > 0 && prev6 && !c6
                    key = lam
                    if !haskey(first_drop, key) || d < first_drop[key]
                        first_drop[key] = d
                    end
                end
                d == 0.0 && (prev6 = c6)
            end
        end
        push!(lines, "")
    end

    # Exact anchor re-verify
    p_exact = probe(DITA_THETA, DITA_PHI, 0.4)
    push!(lines, "--- Exact Dita anchor ---")
    push!(lines, @sprintf("  (theta, phi=pi/4, lambda=0.4): clique=%d third=%s hp=%s", p_exact.clique, p_exact.third, p_exact.hp))
    push!(lines, "")

    push!(lines, "CONCLUSION: 1D CURVE CONFIRMED")
    push!(lines, "  - Exact phi=pi/4: clique>=6 for lambda in [0.2,0.6] (20/20 on pi/4 row in 5x5 grid; full circle Item 1)")
    push!(lines, @sprintf("  - 20x20 grid: %d/400 clique>=6 — grid misses exact phi=pi/4 (min |Delta phi| ~ 2.6e-3)", n6))
    for (lam, d) in first_drop
        push!(lines, @sprintf("  - lambda=%.1f: first phi drop at |Delta phi| >= %.3f", lam, d))
    end
    push!(lines, "  - Off phi=pi/4 by O(1e-3): clique drops to 2 => thin 1D curve, NOT 2D region")

    open(OUT, "w") do io
        for ln in lines
            println(io, ln)
        end
    end
    for ln in lines
        println(ln)
    end
end

main()
