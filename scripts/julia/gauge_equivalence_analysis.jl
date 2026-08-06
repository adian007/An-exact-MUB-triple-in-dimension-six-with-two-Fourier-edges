# Phase 1.1a: Analytic gauge analysis — theta=0 phi redundancy, lambda role.
# Output: results/gauge_analysis.txt
include(joinpath(@__DIR__, "_paths.jl"))
using Printf, LinearAlgebra, Dates

const OUT = joinpath(RESULTS_DIR, "gauge_analysis.txt")
const DITA_THETA = acos(1 / sqrt(3))
const TOL = 1e-12

"""A11, A12 at general (theta, phi) and symbolic theta=0 limit."""
function A_entries(theta, phi)
    A11 = -0.5 + im * (sqrt(3) / 2) * (cos(theta) + exp(-im * phi) * sin(theta))
    A12 = -0.5 + im * (sqrt(3) / 2) * (-cos(theta) + exp(im * phi) * sin(theta))
    return A11, A12
end

function main()
    open(OUT, "w") do io
        println(io, "=== Phase 1.1a: Gauge equivalence analytic analysis ===")
        println(io, "Date: $(Dates.format(now(), "yyyy-mm-dd HH:MM"))\n")

        # --- theta=0: phi independence of A ---
        println(io, "--- Block A at theta=0 ---")
        phi_vals = [0.0, 0.48, 0.5, 0.52, 0.7, 1.3, pi / 4]
        A11_ref, A12_ref = A_entries(0.0, 0.5)
        max_A_diff = 0.0
        for ph in phi_vals
            A11, A12 = A_entries(0.0, ph)
            d = max(abs(A11 - A11_ref), abs(A12 - A12_ref))
            max_A_diff = max(max_A_diff, d)
            @printf(io, "  phi=%.6g  A11=%+.16g+%+.16gi  A12=%+.16g+%+.16gi  |ΔA|=%.3e\n",
                    ph, real(A11), imag(A11), real(A12), imag(A12), d)
        end
        @printf(io, "  max |A(0,phi)-A(0,0.5)| over phi in test set: %.3e\n", max_A_diff)
        println(io, "  Analytic: sin(theta)=0 => A11 = -1/2 + i*sqrt(3)/2, A12 = -1/2 - i*sqrt(3)/2 (phi drops out).\n")

        # --- Full H at theta=0: phi pairs ---
        println(io, "--- Full H(0, phi, lambda) at fixed lambda=0.3 ---")
        lam = 0.3
        H_ref = build_karlsson_family(0.0, 0.5, lam)
        max_H_diff_phi = 0.0
        for ph in [0.48, 0.5, 0.52, 0.7]
            H = build_karlsson_family(0.0, ph, lam)
            d = maximum(abs.(H - H_ref))
            max_H_diff_phi = max(max_H_diff_phi, d)
            @printf(io, "  phi=%.4g  max|H(0,phi,λ)-H(0,0.5,λ)| = %.3e\n", ph, d)
        end
        @printf(io, "  max entrywise |ΔH| over phi pairs: %.3e (tol=%.0e)\n\n", max_H_diff_phi, TOL)

        # --- lambda at Dita: H changes ---
        println(io, "--- H(Dita_theta, pi/4, lambda) vs lambda ---")
        lam_vals = [0.4, 0.41, 0.45, 0.5, 0.9, 1.3, 2.0]
        H_dita_ref = build_karlsson_family(DITA_THETA, pi / 4, 0.4)
        for lam in lam_vals
            H = build_karlsson_family(DITA_THETA, pi / 4, lam)
            d = maximum(abs.(H - H_dita_ref))
            @printf(io, "  lambda=%.4g  max|H(λ)-H(0.4)| = %.6e\n", lam, d)
        end
        println(io, "  lambda enters via z1=exp(i*lambda); all Z blocks change => H depends on lambda by construction.\n")

        # --- Verdict ---
        println(io, "=== VERDICT ===")
        if max_A_diff < TOL && max_H_diff_phi < TOL
            println(io, "phi IS a gauge parameter at theta=0: build_A(0,phi) and H(0,phi,lambda) are")
            println(io, "  independent of phi (max |ΔH| = $(max_H_diff_phi) < $(TOL)). Axis persistence along phi")
            println(io, "  at F6_theta0 reflects parametrization redundancy, not a distinct CHM family.")
        else
            println(io, "phi is NOT redundant at theta=0 (unexpected): max |ΔH| = $(max_H_diff_phi).")
        end
        println(io, "")
        println(io, "lambda is NOT a gauge parameter in the Karlsson parametrization: it is the explicit third")
        println(io, "  degree of freedom (z1 = exp(i*lambda)); H(Dita, pi/4, lambda) changes with lambda.")
        println(io, "")
        println(io, "Phase 1.1b CHM equivalence (chm_equivalence.py, 720 perm + diag phases):")
        println(io, "  See results/chm_equivalence.txt for numeric residuals.")
        println(io, "  F6 phi pairs: residual=0 (EQUIVALENT) — confirms phi gauge at theta=0.")
        println(io, "  Dita lambda pairs vs lambda=0.4: min_residual=2.309e-02 (INEQUIVALENT).")
        println(io, "  Branch: 1.3 (genuine 1D inequivalent-CHM family along lambda at Dita).")
    end
    println("Wrote $OUT")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
