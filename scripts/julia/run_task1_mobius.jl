# Task 1: Mobius degeneracy candidate — BigFloat certify + high-precision third/fourth MUB test.
include(joinpath(@__DIR__, "_paths.jl"))
using Printf

const THETA = 0.9553166181
const PHI = 0.7853981634
const LAM = 0.4
const OUT = joinpath(RESULTS_DIR, "task1_mobius_investigation.txt")

function mobius_z4_dev(theta, phi, lam)
    A = build_A(theta, phi)
    B = -[1 1; 1 -1] - A
    alpha_A, beta_A = A[1, 2]^2, A[1, 1]^2
    alpha_B, beta_B = B[1, 2]^2, B[1, 1]^2
    z1sq = exp(2im * lam)
    z3sq = mobius(z1sq, alpha_A, beta_A)
    z4sq = mobius(z1sq, alpha_B, beta_B)
    num = beta_B - z3sq * conj(alpha_B)
    den = alpha_B - z3sq * conj(beta_B)
    z2sq = num / den
    z4sq_check = mobius(z2sq, alpha_A, beta_A)
    return abs(z4sq - z4sq_check)
end

function hp_fourth_test(theta, phi, lam; bits = 400)
    setprecision(bits)
    H = build_karlsson_family(theta, phi, lam)
    pool, stats = generate_candidate_pool_fresh(H; tol = 1e-10, verbose = false)
    pool = deduplicate_pool(pool; tol = 1e-8)
    ext = check_four_mub_extension(pool, H; ortho_tol = 1e-12, mu_tol = 1e-12)
    hp_detail = nothing
    if ext.found_third
        g = _orthogonality_graph(pool; ortho_tol = 1e-12)
        cliques = [c for c in maximal_cliques(g) if length(c) >= 6]
        best_c = argmax(c -> length(c), cliques)
        clique_vecs = [pool[i] for i in best_c]
        hp_detail = verify_clique_hp(clique_vecs; bits = bits, ortho_tol = 1e-40, mu_tol = 1e-40)
    end
    return (pool = pool, stats = stats, ext = ext, hp_detail = hp_detail, n_pool = length(pool))
end

function main()
    open(OUT, "w") do io
        println(io, "=== Task 1: Mobius degeneracy candidate ===")
        @printf(io, "Exact CSV point: theta=%.10f phi=%.10f lambda=%.1f\n", THETA, PHI, LAM)
        z4 = mobius_z4_dev(THETA, PHI, LAM)
        @printf(io, "z4_dev = %.16e\n", z4)
        dita_theta = acos(1 / sqrt(3))
        @printf(io, "Dita theta (exact) = %.16f, delta = %.3e\n", dita_theta, abs(THETA - dita_theta))

        H = build_karlsson_family(THETA, PHI, LAM)
        @printf(io, "Hadamard defect = %.3e\n\n", norm(H * H' - 6 * I(6)))

        println(io, "[A] Pool completeness (max_precision=128)")
        pc = pool_completeness_report(H; verbose = false)
        @printf(io, "  mv=%d tracked=%d certified=%d distinct=%d verified=%d dedup=%d\n",
                pc.mixed_volume, pc.n_tracked, pc.n_certified, pc.n_distinct_certified,
                pc.n_verified, pc.n_dedup_pool)
        @printf(io, "  pool_complete_flag=%s cert_match=%s\n\n", pc.pool_complete_flag, pc.cert_match)

        println(io, "[B] Re-certify at max_precision=128")
        system = build_numeric_pool_system(H)
        res = solve(system; show_progress = false)
        cert = certify(system, res; show_progress = false, threading = false, max_precision = 128)
        @printf(io, "  n_certified=%d n_tracked=%d\n\n", ncertified(cert), length(solutions(res)))

        println(io, "[C] Third/fourth MUB (standard pool, tol 1e-12)")
        pool, stats = generate_candidate_pool_fresh(H; tol = 1e-10, verbose = false)
        pool = deduplicate_pool(pool; tol = 1e-8)
        ext = check_four_mub_extension(pool, H; ortho_tol = 1e-12, mu_tol = 1e-12)
        @printf(io, "  n_pool=%d max_clique=%d found_third=%s found_fourth=%s\n",
                ext.n_pool, ext.max_clique, ext.found_third, ext.found_fourth)
        @printf(io, "  mu_defect_min=%.3e ortho_defect=%.3e\n\n", ext.mu_defect_min, ext.ortho_defect)

        println(io, "[D] High-precision clique verify (BigFloat bits=400, ~120 decimal digits)")
        hp = hp_fourth_test(THETA, PHI, LAM; bits = 400)
        @printf(io, "  found_third=%s found_fourth=%s max_clique=%d\n",
                hp.ext.found_third, hp.ext.found_fourth, hp.ext.max_clique)
        if hp.hp_detail !== nothing
            @printf(io, "  HP ortho_max=%.3e mu_max=%.3e ortho_ok=%s mu_ok=%s\n",
                    hp.hp_detail.ortho_max, hp.hp_detail.mu_norm_max,
                    hp.hp_detail.ortho_ok, hp.hp_detail.mu_ok)
        end

        println(io, "\n[E] PSLQ / algebraic identification")
        dita = acos(1 / sqrt(3))
        for (val, label) in [(0.0, "0"), (pi/4, "pi/4"), (dita, "arccos(1/sqrt(3))")]
            abs(THETA - val) < 1e-6 && println(io, "  theta ≈ $label")
            abs(PHI - val) < 1e-6 && println(io, "  phi ≈ $label")
        end
        abs(PHI - pi/4) < 1e-6 && println(io, "  phi = pi/4 [Dita slice]")
        abs(THETA - dita) < 1e-6 && println(io, "  theta = arccos(1/sqrt(3)) [Dita]")

        println(io, "\n[F] Per-basis fourth-MUB report")
        for rep in fourth_mub_per_basis_report(pool, H; ortho_tol = 1e-12, mu_tol = 1e-12)
            @printf(io, "  basis %d: n_filtered=%d max_clique=%d found_fourth=%s\n",
                    rep.basis_idx, rep.n_filtered, rep.max_clique_filtered, rep.found_fourth)
        end

        if ext.found_fourth
            println(io, "\n*** COUNTEREXAMPLE: found_fourth=true — STOP ***")
        else
            println(io, "\n=== RESULT: found_fourth=false (no counterexample) ===")
        end
    end
    println("Wrote $OUT")
end

main()
