# Item 1: Lambda periodicity / arc bounds at Dita anchor.
# Usage: julia --project=. scripts/julia/lambda_periodicity_dita.jl
include(joinpath(@__DIR__, "_paths.jl"))
using Printf, Dates

const DITA_THETA = acos(1 / sqrt(3))
const DITA_PHI = pi / 4
const LAMBDA0 = 0.4
const OUT = joinpath(RESULTS_DIR, "lambda_periodicity_dita.txt")
const HP_BITS = 400
const MU_TOL = 1e-12
const ORTHO_TOL = 1e-12

function probe_lambda(lam; hp_bits = HP_BITS)
    H = build_karlsson_family(DITA_THETA, DITA_PHI, lam)
    if !is_hadamard(H; tol = 1e-8)
        return (hadamard = false, clique = -1, third = false, fourth = false, hp = false)
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
            hp = verify_clique_hp([pool[i] for i in best_c]; bits = hp_bits)
            hp_ok = hp.ortho_max < 1e-8 && hp.mu_norm_max < 1e-8
        end
    end
    return (
        hadamard = true,
        clique = ext.max_clique,
        third = ext.found_third,
        fourth = ext.found_fourth,
        hp = hp_ok,
    )
end

has_clique6(p) = p.hadamard && p.clique >= 6

"""Scan lambda values; return vector of (lam, probe)."""
function scan_lambdas(lams)
    [(lam = lam, p = probe_lambda(lam)) for lam in lams]
end

"""Find first drop along increasing lambda from start."""
function find_drop_increasing(start, step, max_lam)
    prev6 = has_clique6(probe_lambda(start))
    lam = start + step
    while lam <= max_lam
        p = probe_lambda(lam)
        cur6 = has_clique6(p)
        if prev6 && !cur6
            return (drop_lo = lam - step, drop_hi = lam)
        end
        prev6 = cur6
        lam += step
    end
    return (drop_lo = NaN, drop_hi = NaN)
end

"""Find first drop along decreasing lambda from start."""
function find_drop_decreasing(start, step, min_lam)
    prev6 = has_clique6(probe_lambda(start))
    lam = start - step
    while lam >= min_lam
        p = probe_lambda(lam)
        cur6 = has_clique6(p)
        if prev6 && !cur6
            return (drop_lo = lam, drop_hi = lam + step)
        end
        prev6 = cur6
        lam -= step
    end
    return (drop_lo = NaN, drop_hi = NaN)
end

function binary_search_drop(lo6, hi6; rel_tol = 1e-10, max_iter = 60)
    # lo6: clique>=6, hi6: clique<6
    for _ in 1:max_iter
        mid = (lo6 + hi6) / 2
        if has_clique6(probe_lambda(mid))
            lo6 = mid
        else
            hi6 = mid
        end
        (hi6 - lo6) < rel_tol && break
    end
    return (lo6, hi6, (lo6 + hi6) / 2)
end

function main()
    lines = String[]
    push!(lines, "=== Item 1: Lambda periodicity at Dita ===")
    push!(lines, "Date: $(Dates.now())")
    push!(lines, @sprintf("theta=%.16g (= arccos(1/sqrt(3))), phi=pi/4=%.16g, lambda0=%.16g", DITA_THETA, DITA_PHI, LAMBDA0))
    push!(lines, "HP_bits=$HP_BITS")
    push!(lines, "")

    # Spot checks at canonical multiples of pi/2 and 2pi wrap
    spot = [0.0, pi / 2, pi, 3pi / 2, 2pi - 0.01, 2pi, 2pi + 0.01]
    push!(lines, "--- Spot checks (absolute lambda) ---")
    n6_spot = 0
    for lam in spot
        p = probe_lambda(lam)
        c6 = has_clique6(p)
        c6 && (n6_spot += 1)
        push!(lines, @sprintf("  lambda=%.12g: hadamard=%s clique=%d third=%s hp=%s clique6=%s",
            lam, p.hadamard, p.clique, p.third, p.hp, c6))
    end
    push!(lines, @sprintf("  clique>=6 count: %d / %d", n6_spot, length(spot)))
    push!(lines, "")

    # Fine step near 2pi wrap
    push!(lines, "--- Fine sweep [6.15, 6.283] step=0.01 ---")
    fine_hi = collect(range(6.15, 6.283; step = 0.01))
    hi_scan = scan_lambdas(fine_hi)
    n6_hi = count(x -> has_clique6(x.p), hi_scan)
    push!(lines, @sprintf("  points=%d clique>=6=%d", length(hi_scan), n6_hi))
    for x in hi_scan
        if !has_clique6(x.p)
            push!(lines, @sprintf("    first drop in window: lambda=%.12g clique=%d", x.lam, x.p.clique))
            break
        end
    end
    push!(lines, "")

    push!(lines, "--- Fine sweep [-0.05, 0.05] step=0.01 ---")
    fine_lo = collect(range(-0.05, 0.05; step = 0.01))
    lo_scan = scan_lambdas(fine_lo)
    n6_lo = count(x -> has_clique6(x.p), lo_scan)
    push!(lines, @sprintf("  points=%d clique>=6=%d", length(lo_scan), n6_lo))
    for x in lo_scan
        if !has_clique6(x.p)
            push!(lines, @sprintf("    first drop in window: lambda=%.12g clique=%d", x.lam, x.p.clique))
            break
        end
    end
    push!(lines, "")

    # Full period coarse scan (step 0.05)
    push!(lines, "--- Full period scan [0, 2pi] step=0.05 ---")
    full_lams = collect(range(0.0, 2pi; step = 0.05))
    full_scan = scan_lambdas(full_lams)
    n6_full = count(x -> has_clique6(x.p), full_scan)
    drops = [x for x in full_scan if !has_clique6(x.p)]
    push!(lines, @sprintf("  points=%d clique>=6=%d clique<6=%d", length(full_scan), n6_full, length(drops)))
    if !isempty(drops)
        push!(lines, "  drop lambda values (first 10):")
        for x in drops[1:min(10, end)]
            push!(lines, @sprintf("    lambda=%.12g clique=%d hadamard=%s", x.lam, x.p.clique, x.hadamard))
        end
    end
    push!(lines, "")

    # Binary search boundaries if not full circle
    all6 = n6_full == length(full_scan)
    if all6
        push!(lines, "VERDICT: FULL CIRCLE — clique>=6 at all 126 samples on [0, 2pi] (step 0.05)")
        push!(lines, "  Spot checks: $n6_spot/$(length(spot)) clique>=6")
        push!(lines, "  Wrap continuity [6.15,6.283] and [-0.05,0.05]: no drops detected")
    else
        # Find arc bounds via binary search from interior point known to have clique6
        interior = LAMBDA0
        has_clique6(probe_lambda(interior)) || (interior = first(x.lam for x in full_scan if has_clique6(x.p)))

        push!(lines, "--- Binary search arc bounds (400-bit HP) from interior lambda=$interior ---")
        # Search upward to 2pi and beyond
        up_br = find_drop_increasing(interior, 0.01, 2pi + 0.5)
        down_br = find_drop_decreasing(interior, 0.01, -0.5)

        lower = NaN
        upper = NaN
        if !isnan(down_br.drop_lo)
            lo, hi, mid = binary_search_drop(down_br.drop_lo, down_br.drop_hi)
            lower = mid
            push!(lines, @sprintf("  lower boundary: [%.12g, %.12g] mid=%.12g", lo, hi, mid))
        else
            push!(lines, "  lower: no drop found scanning to lambda=-0.5")
        end
        if !isnan(up_br.drop_lo)
            lo, hi, mid = binary_search_drop(up_br.drop_lo, up_br.drop_hi)
            upper = mid
            push!(lines, @sprintf("  upper boundary: [%.12g, %.12g] mid=%.12g", lo, hi, mid))
        else
            push!(lines, "  upper: no drop found scanning to lambda=$(2pi + 0.5)")
        end

        if isnan(lower) && isnan(upper)
            push!(lines, "VERDICT: FULL CIRCLE (no boundary bracketed within extended scan)")
        else
            lo_b = isnan(lower) ? 0.0 : lower
            hi_b = isnan(upper) ? 2pi : upper
            width = hi_b - lo_b
            push!(lines, @sprintf("VERDICT: ARC from [%.12g] to [%.12g], width=%.12g", lo_b, hi_b, width))
        end
    end

    open(OUT, "w") do io
        for ln in lines
            println(io, ln)
        end
    end
    for ln in lines
        println(ln)
    end
    println("\nWrote $OUT")
end

main()
