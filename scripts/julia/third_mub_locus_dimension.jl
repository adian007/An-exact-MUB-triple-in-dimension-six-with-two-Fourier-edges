# Fast Task 3 probe: anchors + perturbations + 100-sample local scan (2 Karlsson anchors).
include(joinpath(@__DIR__, "_paths.jl"))
using Printf, Random, Dates

const DITA_THETA = acos(1 / sqrt(3))
const OUT = joinpath(RESULTS_DIR, "third_mub_locus_dimension.txt")
const RADIUS = 1e-4
const N_SAMPLES = 100

function probe(theta, phi, lam)
    H = build_karlsson_family(theta, phi, lam)
    !is_hadamard(H) && return (ok = false, clique = 0, third = false)
    pool, _ = generate_candidate_pool_fresh(H; verbose = false)
    pool = deduplicate_pool(pool)
    ext = check_four_mub_extension(pool, H; mu_tol = 1e-10, ortho_tol = 1e-10)
    return (ok = true, clique = ext.max_clique, third = ext.found_third)
end

function probe_F6()
    pool, _ = generate_candidate_pool_fresh(F6; verbose = false)
    pool = deduplicate_pool(pool)
    ext = check_four_mub_extension(pool, F6; mu_tol = 1e-10, ortho_tol = 1e-10)
    return (ok = true, clique = ext.max_clique, third = ext.found_third)
end

function local_scan(th0, ph0, lm0; n = N_SAMPLES, r = RADIUS)
    rng = MersenneTwister(42)
    cl = Int[]
    for _ in 1:n
        th = th0 + r * randn(rng)
        ph = max(1e-6, ph0 + r * randn(rng))
        lm = lm0 + r * randn(rng)
        try
            p = probe(th, ph, lm)
            p.ok && push!(cl, p.clique)
        catch
        end
    end
    return cl
end

open(OUT, "w") do io
    println(io, "=== Task 3: Third-MUB locus dimension ===")
    println(io, "Date: $(Dates.format(now(), "yyyy-mm-dd"))")
    println(io, "Sampling: radius=$RADIUS, n=$N_SAMPLES per Karlsson anchor\n")

    p6 = probe_F6()
    println(io, "F6 (Fourier): clique=$(p6.clique) third=$(p6.third)")

    for (lab, th, ph, lam) in [
        ("Dita", DITA_THETA, pi / 4, 0.4),
        ("F6_theta0", 0.0, 0.5, 0.3),
    ]
        p = probe(th, ph, lam)
        println(io, "$lab exact: clique=$(p.clique) third=$(p.third)")
    end

    println(io, "\n--- Perturbation tests ---")
    for (lab, th, ph, lam) in [
        ("theta=0 + 1e-4", 1e-4, 0.5, 0.3),
        ("theta=0 + 1e-3", 1e-3, 0.5, 0.3),
        ("theta=0 + 0.01", 0.01, 0.5, 0.3),
        ("Dita + 1e-4 theta", DITA_THETA + 1e-4, pi / 4, 0.4),
        ("Dita + 1e-3 theta", DITA_THETA + 1e-3, pi / 4, 0.4),
        ("Dita + 0.01 theta", DITA_THETA + 0.01, pi / 4, 0.4),
    ]
        p = probe(th, ph, lam)
        println(io, "  $lab: clique=$(p.clique) third=$(p.third) ok=$(p.ok)")
    end

    println(io, "\n--- Local Gaussian sampling ---")
    for (lab, th, ph, lam) in [
        ("Dita", DITA_THETA, pi / 4, 0.4),
        ("F6_theta0", 0.0, 0.5, 0.3),
    ]
        cl = local_scan(th, ph, lam)
        n6 = count(==(6), cl)
        n2 = count(==(2), cl)
        println(io, "  $lab: valid=$(length(cl)) clique6=$n6 ($(round(100*n6/max(1,length(cl)), digits=1))%) clique2=$n2 ($(round(100*n2/max(1,length(cl)), digits=1))%) min=$(isempty(cl) ? 0 : minimum(cl)) max=$(isempty(cl) ? 0 : maximum(cl))")
    end

    println(io, "\n--- Parametric system ---")
    try
        sys = build_parametric_pool_system()
        println(io, "  build_parametric_pool_system: OK ($(length(sys)) vars)")
    catch e
        println(io, "  build_parametric_pool_system: FAILED ($(sprint(showerror, e)))")
    end
    println(io, "  HC monodromy/NID on extended system: not run (API issues on host)")
end
println("Wrote $OUT")
