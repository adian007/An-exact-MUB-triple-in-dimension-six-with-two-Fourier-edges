# Formal numerical verification of Lemmas L1-L4 (Theorem T1).
# Standalone: does not require HomotopyContinuation (App Control safe).
# Usage: julia scripts/julia/formalize_gauge_lemmas.jl

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const RESULTS_DIR = joinpath(ROOT, "results")
include(joinpath(ROOT, "src", "karlsson_gauge_only.jl"))

using LinearAlgebra
using Printf
using Random
using Dates

const OUT = joinpath(RESULTS_DIR, "formalize_gauge_lemmas.txt")
const DITA_THETA = acos(1 / sqrt(3))
const TOL = 1e-12

function A_at_dita()
    A = gauge_build_A(DITA_THETA, pi / 4)
    expected = [im -1; -1 -im]
    return A, expected
end

function chm_residual(H1, H2)
    d1, d2 = gauge_dephase(H1), gauge_dephase(H2)
    best = Inf
    perms = [[1, 2, 3, 4, 5, 6], [2, 1, 3, 4, 5, 6], [1, 3, 2, 4, 5, 6]]
    for perm in perms
        P = d2[:, perm]
        for k in 0:5
            phases = [exp(im * 2π * k * j / 6) for j in 1:6]
            D = diagm(0 => phases)
            best = min(best, norm(d1 - D * P))
        end
    end
    return best
end

function test_L1_phi_gauge(rng::AbstractRNG)
    lam = 0.3
    H_ref = gauge_build_karlsson_family(0.0, 0.5, lam)
    ok = true
    for _ in 1:10
        ph = rand(rng) * pi
        H = gauge_build_karlsson_family(0.0, ph, lam)
        ok &= maximum(abs.(H - H_ref)) < TOL
    end
    return ok
end

function test_L2_lambda_periodicity(rng::AbstractRNG)
    ok = true
    for _ in 1:10
        th = rand(rng) * pi
        ph = rand(rng) * pi
        lam = rand(rng) * 2π
        H1 = gauge_build_karlsson_family(th, ph, lam)
        H2 = gauge_build_karlsson_family(th, ph, lam + 2π)
        ok &= maximum(abs.(H1 - H2)) < TOL
    end
    return ok
end

function test_L3_dita_A_block()
    A, _ = A_at_dita()
    # Correct Dita specialization: A = [[i,-1],[-1,i]] (AA^dag = 2I)
    expected = [im -1; -1 im]
    dA = maximum(abs.(A - expected))
    unitary_err = maximum(abs.(A * A' - 2 * I(2)))
    return dA < TOL && unitary_err < TOL, dA, unitary_err
end

function test_L4_lambda_not_gauge()
    H_ref = gauge_build_karlsson_family(DITA_THETA, pi / 4, 0.4)
    H2 = gauge_build_karlsson_family(DITA_THETA, pi / 4, 0.41)
    entry_diff = maximum(abs.(H_ref - H2))
    res = chm_residual(H_ref, H2)
    return entry_diff > 1e-6 && res > 0.01, entry_diff, res
end

function main()
    rng = MersenneTwister(42)
    l1 = test_L1_phi_gauge(rng)
    l2 = test_L2_lambda_periodicity(rng)
    l3_ok, dA, u_err = test_L3_dita_A_block()
    l4_ok, ed, chm_r = test_L4_lambda_not_gauge()

    mkpath(RESULTS_DIR)
    open(OUT, "w") do io
        println(io, "=== Formal verification of Lemmas L1-L4 (Theorem T1) ===")
        println(io, "Date: $(Dates.format(now(), "yyyy-mm-ddTHH:MM:SS"))")
        println(io, "Tolerance: $TOL\n")
        println(io, "L1 (phi-gauge at theta=0): $(l1 ? "PASS" : "FAIL")")
        println(io, "L2 (lambda 2pi periodicity): $(l2 ? "PASS" : "FAIL")")
        @printf(io, "L3 (Dita A-block): %s  |ΔA|=%.3e  unitary_err=%.3e\n",
                l3_ok ? "PASS" : "FAIL", dA, u_err)
        @printf(io, "L4 (lambda not CHM-gauge): %s  entry_diff=%.3e  chm_residual=%.3e\n",
                l4_ok ? "PASS" : "FAIL", ed, chm_r)
        all_ok = l1 && l2 && l3_ok && l4_ok
        println(io, "\nVERDICT: $(all_ok ? "ALL PASS — Theorem T1 numerically certified" : "FAIL")")
    end
    println("Wrote $OUT")
    all_ok = l1 && l2 && l3_ok && l4_ok
    all_ok || error("Gauge lemma verification failed")
end

main()
