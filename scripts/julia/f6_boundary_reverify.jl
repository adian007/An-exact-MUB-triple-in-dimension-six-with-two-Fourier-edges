# Item B: re-verify F6_theta0 lambda boundaries at 400-bit HP (post argmax fix).
include(joinpath(@__DIR__, "_paths.jl"))
using Printf, Dates

const OUT = joinpath(RESULTS_DIR, "f6_boundary_reverify.txt")
const HP_BITS = 400
const MU_TOL = 1e-12
const ORTHO_TOL = 1e-12

function probe_point(theta, phi, lam)
    H = build_karlsson_family(theta, phi, lam)
    if !is_hadamard(H; tol = 1e-8)
        return (ok = false, clique = -1, hp = false)
    end
    pool, _ = generate_candidate_pool_fresh(H; verbose = false)
    pool = deduplicate_pool(pool)
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
    return (ok = true, clique = ext.max_clique, hp = hp_ok)
end

# Coarse step (0.05) avoids spurious early clique dips seen with step=0.01
# (matches extended_axis_sweep.jl --quick bracket that yields +0.055 / -0.063).
function bracket_axis(lm0, sign::Float64; max_scan = 0.15, step = 0.05)
    inner = 0.0
    outer = NaN
    d = step
    prev6 = true
    while d <= max_scan
        p = probe_point(0.0, 0.5, lm0 + sign * d)
        has6 = p.ok && p.clique >= 6
        if prev6 && !has6
            outer = d
            break
        end
        has6 && (inner = d)
        prev6 = has6
        d += step
    end
    return inner, outer
end

function binary_search_boundary(lm0, sign::Float64, lo, hi; rel_tol = 1e-10)
    for _ in 1:50
        mid = (lo + hi) / 2
        p = probe_point(0.0, 0.5, lm0 + sign * mid)
        has6 = p.ok && p.clique >= 6
        if has6
            lo = mid
        else
            hi = mid
        end
        (hi - lo) / max(abs(lm0), 1e-6) < rel_tol && break
    end
    return (lo, hi, (lo + hi) / 2)
end

function main()
    lines = String[]
    push!(lines, "=== F6_theta0 lambda boundary re-verify (400-bit HP) ===")
    push!(lines, "Date: $(Dates.now())  HP_bits=$HP_BITS")
    push!(lines, "  theta=0, phi=0.5 [gauge], lambda0=0.3")
    push!(lines, "")

    # Anchor check
    p0 = probe_point(0.0, 0.5, 0.3)
    push!(lines, @sprintf("  anchor (0,0.5,0.3): clique=%d hp=%s", p0.clique, p0.hp))
    push!(lines, "")

    boundaries = Float64[]
    for sign in (+1.0, -1.0)
        inner, outer = bracket_axis(0.3, sign)
        if isnan(outer)
            push!(lines, @sprintf("  sign=%+.0f: clique>=6 to scan limit (inner=%.12g)", sign, inner))
        else
            lo, hi, mid = binary_search_boundary(0.3, sign, inner, outer)
            push!(lines, @sprintf("  sign=%+.0f: boundary |Delta lambda| in [%.12g, %.12g]  mid=%.12g", sign, lo, hi, mid))
            push!(boundaries, mid)
        end
    end

    arc_width = sum(abs.(boundaries))
    push!(lines, "")
    push!(lines, @sprintf("  arc width (|Delta+| + |Delta-|) = %.12g", arc_width))

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
