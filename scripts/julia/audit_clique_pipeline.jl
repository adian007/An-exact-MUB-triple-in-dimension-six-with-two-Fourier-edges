# Phase 1: Audit clique-finding pipeline and null-result plausibility.
# Usage: julia --project=. audit_clique_pipeline.jl

include(joinpath(@__DIR__, "_paths.jl"))

using LinearAlgebra
using Random
using Statistics

const MU_TOLS = (1e-6, 1e-8, 1e-10)
const ORTHO_TOLS = (1e-6, 1e-8, 1e-10, 1e-12)

function check(label, ok, detail = "")
    status = ok ? "PASS" : "FAIL"
    println("  [$status] $label" * (detail == "" ? "" : " — $detail"))
    return ok
end

"""Six columns of a random unitary in C^6 (known ONB / planted third MUB)."""
function random_onb6(rng = Random.default_rng())
    X = randn(rng, ComplexF64, 6, 6)
    Q, _ = qr(X)
    return [Q[:, j] for j in 1:6]
end

"""Plant a 6-clique: ONB columns + decoy vectors with controlled orthogonality."""
function planted_clique_pool(n_decoy = 40; rng = Random.default_rng())
    clique = random_onb6(rng)
    decoys = [normalize(randn(rng, ComplexF64, 6)) for _ in 1:n_decoy]
    return vcat(clique, decoys), clique
end

"""Verify clique vectors are pairwise orthogonal."""
function verify_clique(clique; tol = 1e-12)
    for i in 1:length(clique), j in (i + 1):length(clique)
        d = abs(dot(clique[i], clique[j]))
        d >= tol && return false, d
    end
    return true, 0.0
end

function clique_stats(pool; ortho_tol = 1e-8)
    n = length(pool)
    g = _orthogonality_graph(pool; ortho_tol = ortho_tol)
    cliques = maximal_cliques(g)
    max_c = isempty(cliques) ? 0 : maximum(length(c) for c in cliques)
    dots = Float64[]
    for i in 1:n, j in (i + 1):n
        push!(dots, abs(dot(pool[i], pool[j])))
    end
    return (
        max_clique = max_c,
        n_cliques = length(cliques),
        min_dot = minimum(dots),
        median_dot = median(dots),
        max_dot = maximum(dots),
        n_near_ortho = count(d -> d < ortho_tol, dots),
        ortho_tol = ortho_tol,
    )
end

function test_planted_clique()
    println("\n=== [T1] Planted 6-clique recovery (synthetic, no MU constraints) ===")
    all_ok = true
    rng = MersenneTwister(20260803)
    for n_decoy in (0, 10, 100)
        for ortho_tol in ORTHO_TOLS
            pool, clique = planted_clique_pool(n_decoy; rng)
            ok_ortho, worst = verify_clique(clique)
            ext = check_four_mub_extension(pool, F6; ortho_tol = ortho_tol, mu_tol = 1e-8)
            found = ext.max_clique >= 6
            all_ok &= check(
                "decoys=$n_decoy ortho_tol=$ortho_tol → max_clique≥6",
                found,
                "max_clique=$(ext.max_clique), planted ortho worst=$worst",
            )
        end
    end
    all_ok
end

"""Extract a maximal 6-clique from a real pool (returns indices and vectors)."""
function extract_six_clique(pool; ortho_tol = 1e-8)
    g = _orthogonality_graph(pool; ortho_tol = ortho_tol)
    cliques = [c for c in maximal_cliques(g) if length(c) >= 6]
    isempty(cliques) && return nothing, nothing
    best = argmax(c -> length(c), cliques)
    return best, [pool[i] for i in best]
end

"""Plant a literature-grounded 6-clique from the F6 pool plus decoy vectors."""
function planted_f6_clique_pool(H = F6; n_decoy = 30, ortho_tol = 1e-8,
                                rng = Random.default_rng())
    pool_full, = generate_candidate_pool_fresh(H; verbose = false)
    pool_full = deduplicate_pool(pool_full)
    clique_idx, clique_vecs = extract_six_clique(pool_full; ortho_tol = ortho_tol)
    clique_idx === nothing && return nothing, nothing, nothing
    decoy_idx = setdiff(1:length(pool_full), clique_idx)
    if length(decoy_idx) >= n_decoy
        pick = decoy_idx[randperm(rng, length(decoy_idx))[1:n_decoy]]
    else
        pick = decoy_idx
    end
    planted = vcat(clique_vecs, [pool_full[i] for i in pick])
    return planted, clique_vecs, pool_full
end

function test_planted_f6_clique()
    println("\n=== [T2] Planted MU 6-clique from F6 pool (literature-grounded) ===")
    planted, clique, pool_full = planted_f6_clique_pool(F6)
    if planted === nothing
        println("  [FAIL] Could not extract 6-clique from F6 pool.")
        return false
    end
    ok_ortho, worst = verify_clique(clique)
    all_ok = check("F6 pool contains orthogonal 6-clique", ok_ortho,
                    "worst intra-|dot|=$worst, n_pool=$(length(pool_full))")
    for ortho_tol in ORTHO_TOLS
        for mu_tol in MU_TOLS
            ext = check_four_mub_extension(planted, F6; ortho_tol = ortho_tol, mu_tol = mu_tol)
            all_ok &= check(
                "F6-planted clique recovered (ortho=$ortho_tol, mu=$mu_tol)",
                ext.max_clique >= 6,
                "max_clique=$(ext.max_clique), found_third=$(ext.found_third)",
            )
        end
    end
    all_ok
end

"""Negative control: single planted 6-clique must NOT trigger found_fourth."""
function test_planted_fourth_negative_control()
    println("\n=== [T2b] Planted fourth-MUB negative control ===")
    planted, clique, = planted_f6_clique_pool(F6; n_decoy = 50)
    if planted === nothing
        println("  [SKIP] No F6 clique available for negative control.")
        return true
    end
    ext = check_four_mub_extension(planted, F6)
    ok = ext.found_third && !ext.found_fourth
    check("single 6-clique → found_third=true, found_fourth=false", ok,
          "third=$(ext.found_third) fourth=$(ext.found_fourth) max_clique=$(ext.max_clique)")
end

"""Synthetic positive path for found_fourth logic: two disjoint ONBs in one pool."""
function test_synthetic_two_cliques_negative()
    println("\n=== [T2c] Synthetic two-ONB pool (found_fourth should be false without MU to H) ===")
    rng = MersenneTwister(42)
    onb1 = random_onb6(rng)
    onb2 = random_onb6(rng)
    pool = vcat(onb1, onb2)
    ext = check_four_mub_extension(pool, F6)
    # Two random ONBs are orthogonal internally but not MU to F6 → found_fourth stays false
    check("random two-ONB pool: found_fourth=false (not MU to H)", !ext.found_fourth,
          "max_clique=$(ext.max_clique) fourth=$(ext.found_fourth)")
end

function audit_512_from_csv()
    println("\n=== [T3] 512-point sweep CSV re-check (read-only) ===")
    path = joinpath(RESULTS_DIR, "karlsson_perH_sweep.csv")
    if !isfile(path)
        println("  [SKIP] Missing $path")
        return
    end
    n = 0
    max_cliques = Int[]
    n_third = n_fourth = 0
    n_pools = Set{Int}()
    open(path) do io
        header = readline(io)
        for line in eachline(io)
            isempty(strip(line)) && continue
            parts = split(line, ',')
            length(parts) < 7 && continue
            n += 1
            push!(max_cliques, Int(round(parse(Float64, parts[5]))))
            push!(n_pools, Int(round(parse(Float64, parts[4]))))
            parse(Float64, parts[6]) > 0.5 && (n_third += 1)
            parse(Float64, parts[7]) > 0.5 && (n_fourth += 1)
        end
    end
    println("  rows parsed: $n")
    println("  max_clique range: $(minimum(max_cliques))–$(maximum(max_cliques))")
    println("  found_third=$n_third  found_fourth=$n_fourth")
    println("  n_pool values: $(sort(collect(n_pools)))")
    println("  Grid missed θ=0, Dita (θ≈0.955): sweep used θ,φ≥0.02 only.")
    check("512/512 sampled points have max_clique=2", all(==(2), max_cliques))
    check("512/512 found_third=false", n_third == 0)
    check("512/512 found_fourth=false", n_fourth == 0)
end

function analyze_real_pool(H, label; ortho_tols = ORTHO_TOLS)
    println("\n=== Real pool analysis: $label ===")
    pool, stats = generate_candidate_pool_fresh(H; verbose = false)
    pool = deduplicate_pool(pool)
    println("  n_pool=$(length(pool))  raw=$(stats.n_raw)  verified=$(stats.n_verified)")
    for tol in ortho_tols
        s = clique_stats(pool; ortho_tol = tol)
        println("  ortho_tol=$(tol): max_clique=$(s.max_clique)  " *
                "near_ortho_pairs=$(s.n_near_ortho)/$(length(pool)*(length(pool)-1)÷2)  " *
                "min|dot|=$(s.min_dot)  median=$(s.median_dot)")
    end
    ext = check_four_mub_extension(pool, H)
    println("  check_four_mub: max_clique=$(ext.max_clique) found_third=$(ext.found_third)")
    g = _orthogonality_graph(pool; ortho_tol = 1e-8)
    cliques = maximal_cliques(g)
    if !isempty(cliques)
        c = argmax(c -> length(c), cliques)
        worst = 0.0
        for a in 1:length(c), b in (a + 1):length(c)
            worst = max(worst, abs(dot(pool[c[a]], pool[c[b]])))
        end
        println("  best clique size=$(length(c))  worst intra-|dot|=$worst")
    end
    pool, ext
end

function audit_512_null_result()
    println("\n=== [T4] Spot-check sweep points (varying pool sizes) ===")
    spots = [
        (0.02, 0.02, 0.0),
        (0.906891, 0.906891, 0.785398),
        (1.79362, 1.79362, 1.5708),
        (2.35619, 2.35619, 3.14159),
    ]
    for (theta, phi, lam) in spots
        H = build_karlsson_family(theta, phi, lam)
        analyze_real_pool(H, "θ=$theta φ=$phi λ=$lam")
    end
end

function main()
    println("=== Clique pipeline audit ===")
    ok1 = test_planted_clique()
    ok2 = test_planted_f6_clique()
    test_planted_fourth_negative_control()
    test_synthetic_two_cliques_negative()

    println("\n=== [T5] F6 anchor pool (literature: third MUB triples exist for {I,F6}) ===")
    analyze_real_pool(F6, "H=F6")
    pc = pool_completeness_report(F6; verbose = true)

    println("\n=== [T6] Generic Karlsson point ===")
    analyze_real_pool(build_karlsson_family(0.3, 0.5, 0.2), "generic K6")

    println("\n=== [T7] Dita-family point ===")
    analyze_real_pool(build_karlsson_family(acos(1 / sqrt(3)), pi / 4, 0.4), "Dita")

    audit_512_from_csv()
    audit_512_null_result()

    println("\n=== Verdict ===")
    if ok1 && ok2
        println("  Planted-clique tests PASSED — clique finder works on synthetic and F6 pools.")
        println("  F6 pool completeness: distinct_cert=$(pc.n_distinct_certified) dedup=$(pc.n_dedup_pool) flag=$(pc.pool_complete_flag)")
        println("  If real pools still show max_clique=2 at generic points, that is structural, not a bug.")
    else
        println("  Planted-clique test FAILED — fix clique code before any negative claim.")
        exit(1)
    end
end

main()
