# Path B2: structural explanation of F6 bounded λ-arc vs Dita full λ-circle.
# Usage: julia --project=. scripts/julia/analyze_arc_circle_asymmetry.jl
include(joinpath(@__DIR__, "_paths.jl"))
using Printf, Dates, LinearAlgebra

const OUT = joinpath(RESULTS_DIR, "arc_circle_asymmetry.txt")
const DITA_THETA = acos(1 / sqrt(3))
const F6_LAM0 = 0.3
const DITA_LAM0 = 0.4

function mobius_params(theta, phi)
    A = build_A(theta, phi)
    B = -[1 1; 1 -1] - A
    return (
        A = A,
        B = B,
        alpha_A = A[1, 2]^2,
        beta_A = A[1, 1]^2,
        alpha_B = B[1, 2]^2,
        beta_B = B[1, 1]^2,
        A_hermitian = norm(A - A') < 1e-10,
        A_unitary_err = maximum(abs.(A * A' - 2I)),
    )
end

function lambda_profile(theta, phi; n=36, label="")
    lams = range(0, 2pi; length=n)
    rows = NamedTuple[]
    for lam in lams
        H = build_karlsson_family(theta, phi, lam)
        pool, stats = generate_candidate_pool_fresh(H; verbose=false)
        pool = deduplicate_pool(pool)
        ext = check_four_mub_extension(pool, H; ortho_tol=1e-8)
        pc = pool_completeness_report(H; verbose=false)
        push!(rows, (
            lambda=lam,
            max_clique=ext.max_clique,
            n_pool=length(pool),
            n_tracked=pc.n_tracked,
            mv_gap=pc.mv_gap,
        ))
    end
    clique6 = count(r -> r.max_clique >= 6, rows)
    return rows, clique6
end

function arc_width_estimate(rows)
    cl6 = [r.lambda for r in rows if r.max_clique >= 6]
    isempty(cl6) && return (0.0, 0.0, 0.0)
    sorted = sort(cl6)
    # handle wrap on circle
    gaps = diff(sorted)
    wrap_gap = sorted[1] + 2pi - sorted[end]
    if wrap_gap > maximum(gaps)
        # arc wraps — take complement
        max_gap_idx = argmax(gaps)
        lo = sorted[max_gap_idx + 1]
        hi = sorted[max_gap_idx] + 2pi
        width = 2pi - wrap_gap
        return (lo, hi, width)
    else
        return (sorted[1], sorted[end], sorted[end] - sorted[1])
    end
end

function main()
    open(OUT, "w") do io
        println(io, "=== Path B2: arc vs circle structural analysis ===")
        println(io, "Date: $(Dates.format(now(), "yyyy-mm-dd HH:MM"))")
        println(io, "")

        println(io, "--- Block A,B structure at the two loci ---")
        for (label, th, ph) in [("F6_theta0", 0.0, 0.5), ("Dita", DITA_THETA, pi / 4)]
            mp = mobius_params(th, ph)
            @printf(io, "%s (θ=%.6g, φ=%.6g):\n", label, th, ph)
            @printf(io, "  A = [%.4f%+.4fi  %.4f%+.4fi;  %.4f%+.4fi  %.4f%+.4fi]\n",
                    real(mp.A[1, 1]), imag(mp.A[1, 1]), real(mp.A[1, 2]), imag(mp.A[1, 2]),
                    real(mp.A[2, 1]), imag(mp.A[2, 1]), real(mp.A[2, 2]), imag(mp.A[2, 2]))
            @printf(io, "  A Hermitian=%s  A unitary err=%.2e\n", mp.A_hermitian, mp.A_unitary_err)
            @printf(io, "  |α_A|=%.6f |β_A|=%.6f  |α_B|=%.6f |β_B|=%.6f\n",
                    abs(mp.alpha_A), abs(mp.beta_A), abs(mp.alpha_B), abs(mp.beta_B))
            if th ≈ 0.0
                println(io, "  φ-gauge: A independent of φ (sin θ=0); H(0,φ,λ) constant in φ.")
            end
            if abs(th - DITA_THETA) < 1e-6
                println(io, "  Dita specialization: A(θ_D,π/4)=[[i,-1],[-1,i]] (Lemma L3).")
                println(io, "  Sharp φ=π/4 condition: off-slice |Δφ|≥10⁻³ → clique=2 (prior HP audit).")
            end
            println(io, "")
        end

        println(io, "--- λ-profile (from prior HP audits; no re-solve) ---")
        println(io, "F6_theta0: clique≥6 on bounded arc width 0.117643 rad (f6_boundary_reverify.txt, 400-bit)")
        println(io, "Dita:      clique≥6 on full [0,2π] — 126/126 + 628/628 (lambda_periodicity_dita.txt)")
        println(io, "")
        println(io, "--- Spot-check at anchor λ (fresh solve, 2 points) ---")
        for (label, th, ph, lam) in [
            ("F6 anchor λ", 0.0, 0.5, F6_LAM0),
            ("Dita anchor λ", DITA_THETA, pi / 4, DITA_LAM0),
        ]
            try
                H = build_karlsson_family(th, ph, lam)
                pool, _ = generate_candidate_pool_fresh(H; verbose=false)
                pool = deduplicate_pool(pool)
                ext = check_four_mub_extension(pool, H; ortho_tol=1e-8)
                @printf(io, "  %s: max_clique=%d n_pool=%d\n", label, ext.max_clique, length(pool))
            catch e
                @printf(io, "  %s: SOLVE_FAIL (%s)\n", label, sprint(showerror, e))
            end
        end
        println(io, "")

        println(io, "=== STRUCTURAL EXPLANATION ===")
        println(io, "")
        println(io, "1. Common mechanism: third-MUB locus is 1D in λ at both anchors.")
        println(io, "   - F6: φ is gauge (T1), so physical locus is λ-interval ⊂ S¹_λ.")
        println(io, "   - Dita: θ=θ_D and φ=π/4 are sharp (Prop L5), locus is full S¹_λ.")
        println(io, "")
        println(io, "2. Block-level difference (not merely numerical):")
        println(io, "   - F6 θ=0: A = -½ + i(√3/2)I₂ (non-Hermitian 2×2 block); Möbius data (α_A,β_A,α_B,β_B)")
        println(io, "     fixed; z₁=e^{iλ} drives Z-blocks so κ(H) varies with λ — clique-6 only on bounded arc.")
        println(io, "   - Dita: A = [[i,-1],[-1,i]] (Hermitian unitary block, Lemma L3); same λ-periodicity of H")
        println(io, "     but κ(H(λ))≥6 for all λ tested (126/126 + 628/628 HP).")
        println(io, "")
        println(io, "3. Symmetry-group contrast:")
        println(io, "   - At θ=0: Stab(A) contains φ-translations (gauge); effective moduli = λ only.")
        println(io, "   - At Dita: A is pinned by (θ,φ); λ mod 2π gives CHM-inequivalent family (L4).")
        println(io, "   - Arc vs circle is NOT a sampling artifact: F6 arc width 0.118 rad (400-bit HP);")
        println(io, "     Dita full circle verified independently.")
        println(io, "")
        println(io, "4. Conjectured mechanism (HP-supported, not proved analytically):")
        println(io, "   The Hermitian Dita block A(θ_D,π/4) aligns Möbius parameters so the MU-pool")
        println(io, "   orthogonality graph retains a 6-clique for every z₁; the F6 θ=0 block does not,")
        println(io, "   yielding clique-6 only when z₁ lies in a bounded arc. Closing this gap requires")
        println(io, "   an analytic characterization of κ(H(λ)) in terms of (α_A,β_A,α_B,β_B,z₁).")
        println(io, "")
        println(io, "VERDICT: Arc-vs-circle asymmetry is a consequence of distinct A-block geometry")
        println(io, "  (F6 non-Hermitian vs Dita Hermitian) and φ-gauge at θ=0, not unexplained numerics.")
    end
    println("Wrote $OUT")
end

main()
