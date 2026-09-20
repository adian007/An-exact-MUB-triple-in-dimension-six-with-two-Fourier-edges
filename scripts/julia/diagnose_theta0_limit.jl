# Empirical check: what does z2² approach as θ → 0?
# At θ=0 exactly, num2 = den2 = 0 identically (indeterminate 0/0).
# The legacy code evaluates it via floating-point noise, coherently giving
# z2² ≈ alpha_A/beta_A. Question: is that the TRUE limit (so resolving the
# indeterminacy deterministically is mathematically correct), or path noise?
include(joinpath(@__DIR__, "_paths.jl"))
include(joinpath(ROOT, "src", "MubSearch.jl"))
using .MubSearch
using Printf

for (ph, lm) in [(0.5, 0.3), (0.5, 0.7), (0.9, 1.1), (0.5, 0.0), (1.3, 2.4)]
    # θ=0 exact values of alpha_A, beta_A (φ-independent block)
    A0 = build_A_karlsson_original(0.0, ph)
    aA0, bA0 = A0[1, 2]^2, A0[1, 1]^2
    ratio = aA0 / bA0
    @printf "phi=%.2f lam=%.2f: alpha_A/beta_A = %.15f%+.15fi  |.|=%.15f\n" ph lm real(ratio) imag(ratio) abs(ratio)
    for th in (1e-3, 1e-5, 1e-7, 1e-9, 1e-11)
        A = build_A_karlsson_original(th, ph)
        B = -F2_HAD - A
        α_A, β_A = A[1, 2]^2, A[1, 1]^2
        α_B, β_B = B[1, 2]^2, B[1, 1]^2
        z1sq = exp(2im * lm)
        z3sq = (α_A * z1sq - β_A) / (conj(β_A) * z1sq - conj(α_A))
        num2 = β_B - z3sq * conj(α_B)
        den2 = α_B - z3sq * conj(β_B)
        z2sq = num2 / den2
        @printf "  theta=%.0e: z2^2 = %+.12f%+.12fi   |z2^2 - alpha/beta| = %.3e  |z2^2|-1 = %.3e\n" th real(z2sq) imag(z2sq) abs(z2sq - ratio) abs(abs(z2sq) - 1)
    end
end
