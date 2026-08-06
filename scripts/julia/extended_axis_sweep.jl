# Phase 1.4: binary-search clique-6 boundaries + 2D locus grid.
# Usage: julia --project=. scripts/julia/extended_axis_sweep.jl [--quick]
include(joinpath(@__DIR__, "_paths.jl"))
using Printf, Dates, DelimitedFiles

const DITA_THETA = acos(1 / sqrt(3))
const OUT_BOUND = joinpath(RESULTS_DIR, "phase1_boundaries.txt")
const OUT_GRID = joinpath(RESULTS_DIR, "locus_2d_grid.csv")
const HP_BITS = 400
const MU_TOL = 1e-12
const ORTHO_TOL = 1e-12

const QUICK = "--quick" in ARGS

function probe_point(theta, phi, lam; hp_bits = HP_BITS, verify_hp = true)
    H = build_karlsson_family(theta, phi, lam)
    if !is_hadamard(H; tol = 1e-8)
        return (ok = false, clique = -1, third = false, fourth = false, hp = false)
    end
    pool, _ = generate_candidate_pool_fresh(H; verbose = false)
    pool = deduplicate_pool(pool)
    ext = check_four_mub_extension(pool, H; mu_tol = MU_TOL, ortho_tol = ORTHO_TOL)
    hp_ok = true
    if verify_hp && ext.found_third && ext.max_clique >= 6
        g = _orthogonality_graph(pool; ortho_tol = ORTHO_TOL)
        cliques = [c for c in maximal_cliques(g) if length(c) >= 6]
        if !isempty(cliques)
            best_c = argmax(c -> length(c), cliques)
            hp = verify_clique_hp([pool[i] for i in best_c]; bits = hp_bits)
            hp_ok = hp.ortho_max < 1e-8 && hp.mu_norm_max < 1e-8
        end
    end
    return (
        ok = true,
        clique = ext.max_clique,
        third = ext.found_third,
        fourth = ext.found_fourth,
        hp = hp_ok,
    )
end

"""Bracket [inner, outer] where inner has clique>=6 and outer does not."""
function bracket_axis(th0, ph0, lm0, axis::Symbol, sign::Float64; max_scan = 2π)
    inner = 0.0
    outer = NaN
    step = QUICK ? 0.05 : 0.01
    d = step
    prev6 = true
    while d <= max_scan
        th, ph, lam = th0, ph0, lm0
        if axis == :lambda
            lam = lm0 + sign * d
        elseif axis == :phi
            ph = max(1e-8, ph0 + sign * d)
        else
            th = th0 + sign * d
        end
        p = probe_point(th, ph, lam)
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

function binary_search_boundary(th0, ph0, lm0, axis::Symbol, sign::Float64, lo, hi; rel_tol = 1e-10)
    # lo: clique>=6, hi: clique<6
    for _ in 1:(QUICK ? 20 : 50)
        mid = (lo + hi) / 2
        th, ph, lam = th0, ph0, lm0
        if axis == :lambda
            lam = lm0 + sign * mid
        elseif axis == :phi
            ph = max(1e-8, ph0 + sign * mid)
        else
            th = th0 + sign * mid
        end
        p = probe_point(th, ph, lam)
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

function run_boundaries()
    lines = String[]
    push!(lines, "=== Phase 1.4: exact clique-6 boundaries (binary search) ===")
    push!(lines, "Date: $(Dates.now())  HP_bits=$HP_BITS  quick=$QUICK")
    push!(lines, "")

    results = NamedTuple[]

    # Dita lambda axis (both directions)
    push!(lines, "--- Dita: lambda axis (theta=$(DITA_THETA), phi=pi/4, lambda0=0.4) ---")
    for sign in (+1.0, -1.0)
        inner, outer = bracket_axis(DITA_THETA, pi / 4, 0.4, :lambda, sign)
        if isnan(outer)
            push!(lines, @sprintf("  sign=%+.0f: clique>=6 to scan limit (inner=%.4g, no drop found)", sign, inner))
            push!(results, (anchor = "Dita", axis = "lambda", sign = sign, lo = inner, hi = NaN, boundary = NaN))
        else
            lo, hi, mid = binary_search_boundary(DITA_THETA, pi / 4, 0.4, :lambda, sign, inner, outer)
            push!(lines, @sprintf("  sign=%+.0f: boundary |Δλ| in [%.12g, %.12g]  mid=%.12g", sign, lo, hi, mid))
            push!(results, (anchor = "Dita", axis = "lambda", sign = sign, lo = lo, hi = hi, boundary = mid))
        end
    end
    push!(lines, "")

    # F6 theta=0 lambda axis (phi is gauge — skip phi axis)
    push!(lines, "--- F6_theta0: lambda axis (theta=0, phi=0.5 [gauge], lambda0=0.3) ---")
    for sign in (+1.0, -1.0)
        inner, outer = bracket_axis(0.0, 0.5, 0.3, :lambda, sign)
        if isnan(outer)
            push!(lines, @sprintf("  sign=%+.0f: clique>=6 to scan limit (inner=%.4g)", sign, inner))
            push!(results, (anchor = "F6_theta0", axis = "lambda", sign = sign, lo = inner, hi = NaN, boundary = NaN))
        else
            lo, hi, mid = binary_search_boundary(0.0, 0.5, 0.3, :lambda, sign, inner, outer)
            push!(lines, @sprintf("  sign=%+.0f: boundary |Δλ| in [%.12g, %.12g]  mid=%.12g", sign, lo, hi, mid))
            push!(results, (anchor = "F6_theta0", axis = "lambda", sign = sign, lo = lo, hi = hi, boundary = mid))
        end
    end
    push!(lines, "")

    # F6 combo point (0, 0.51, 0.31)
    push!(lines, "--- F6 combo (theta=0, phi=0.51, lambda=0.31) re-verify ---")
    p_combo = probe_point(0.0, 0.51, 0.31)
    p_phi = probe_point(0.0, 0.51, 0.3)
    p_lam = probe_point(0.0, 0.5, 0.31)
    push!(lines, @sprintf("  (0,0.51,0.31): clique=%d hp=%s", p_combo.clique, p_combo.hp))
    push!(lines, @sprintf("  (0,0.51,0.30) phi-only: clique=%d", p_phi.clique))
    push!(lines, @sprintf("  (0,0.50,0.31) lambda-only: clique=%d", p_lam.clique))
    push!(lines, "  (phi redundant at theta=0 => H(0,0.51,0.31)=H(0,0.5,0.31) up to gauge)")
    push!(lines, "")

    open(OUT_BOUND, "w") do io
        for ln in lines
            println(io, ln)
        end
    end
    for ln in lines
        println(ln)
    end
    return results
end

"""Evenly spaced range with exact center pinned (for phi-sensitive locus)."""
function anchored_range(center, halfwidth, n)
    n == 1 && return [center]
    pts = collect(center .+ range(-halfwidth, halfwidth; length = n))
    pts[cld(n, 2)] = center
    return pts
end

function run_2d_grid()
    n = QUICK ? 5 : 20
    ph_center = pi / 4
    lam_center = 0.4
    lams = anchored_range(lam_center, 0.2, n)
    phis = anchored_range(ph_center, 0.05, n)

    rows = String[]
    push!(rows, "theta,phi,lambda,max_clique,found_third,found_fourth,hp_ok,is_hadamard")
    n6 = 0
    total = 0
    for ph in phis
        for lam in lams
            p = probe_point(DITA_THETA, ph, lam; verify_hp = false)
            total += 1
            c6 = p.ok && p.clique >= 6
            c6 && (n6 += 1)
            push!(rows, @sprintf("%.12g,%.12g,%.12g,%d,%s,%s,%s,%s",
                DITA_THETA, ph, lam, p.clique, p.third, p.fourth, p.hp, p.ok))
        end
    end
    open(OUT_GRID, "w") do io
        for r in rows
            println(io, r)
        end
    end
    println("\n=== 2D grid (Dita theta fixed): $n6 / $total clique>=6 ===")
    println("Wrote $OUT_GRID")
    return (n6, total, n)
end

println("Phase 1.4 extended axis sweep (quick=$QUICK)")
if "--grid-only" in ARGS
    run_2d_grid()
else
    run_boundaries()
    run_2d_grid()
end
