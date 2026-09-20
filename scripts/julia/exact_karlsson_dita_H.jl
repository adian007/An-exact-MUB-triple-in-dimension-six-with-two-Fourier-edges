# E1: exact Karlsson H on the Dita slice at lambda=0 and lambda=pi/3, in Nemo.
# No Groebner claim. No HomotopyContinuation.
#
# Field: Q(i, √3) ≅ Q(ζ_12), which contains the Dita A-block [[i,-1],[-1,i]]
# (Lemma L3) and z1 = exp(i λ) at these two algebraic angles.
# MU-vector field Q(ζ_24, √5) is larger and is used for B3 / W1, not for H itself.
#
# Usage: julia --project=. --compiled-modules=no scripts/julia/exact_karlsson_dita_H.jl
# Output: results/exact_karlsson_dita_H.txt

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(ROOT, "src", "karlsson_gauge_only.jl"))

using LinearAlgebra
using Printf
using Dates
using Nemo

const RESULTS_DIR = joinpath(ROOT, "results")
const OUT = joinpath(RESULTS_DIR, "exact_karlsson_dita_H.txt")

# Q(ζ_12) = Q(i, √3), degree 4. ζ_12 = exp(2π i / 12) = exp(i π / 6).
const K12, Z12 = cyclotomic_field(12)
const II = Z12^3          # i
const SQRT3 = Z12 + Z12^11  # 2 cos(π/6) = √3

function embed_qi(a::Int, b::Int)
    return K12(a) + K12(b) * II
end

"""Complex conjugation on Q(ζ_12): ζ ↦ ζ^{-1} = ζ^{11}."""
function conj12(x)
    n = degree(parent(x))
    s = zero(parent(x))
    for k in 0:(n - 1)
        ck = coeff(x, k)
        iszero(ck) && continue
        s += parent(x)(ck) * Z12^((11 * k) % 12)
    end
    return s
end

function mobius_exact(z, alpha, beta)
    num = alpha * z - beta
    den = conj12(beta) * z - conj12(alpha)
    return num * inv(den)
end

"""Square roots of zsq in Q(ζ_12) by a small integer search (no factor/Hecke)."""
function principal_sqrt(zsq)
    rts = typeof(zsq)[]
    for a in -3:3, b in -3:3, c in -3:3, d in -3:3
        u = K12(a) + K12(b) * Z12 + K12(c) * Z12^2 + K12(d) * Z12^3
        iszero(u) && continue
        if u^2 == zsq
            any(r -> r == u, rts) || push!(rts, u)
        end
    end
    isempty(rts) && error("sqrt not in Q(ζ_12) with small coeffs; need a field extension")
    function toC(u)
        ζ = cis(π / 6)
        acc = 0.0 + 0im
        for k in 0:(degree(K12) - 1)
            acc += Float64(coeff(u, k)) * ζ^k
        end
        return acc
    end
    target = sqrt(toC(zsq))  # Julia principal branch on the complex image
    sort!(rts; by=u -> abs(toC(u) - target))
    return rts[1], rts, toC.(rts)
end

function build_exact_dita(lambda_tag::String, z1)
    # Exact A at Dita (Lemma L3): [[i,-1],[-1,i]]
    A11 = II
    A12 = -one(K12)
    A21 = -one(K12)
    A22 = II
    # F2 = [[1,1],[1,-1]]
    # B = -F2 - A
    B11 = -one(K12) - A11   # -1 - i
    B12 = -one(K12) - A12   # -1 - (-1) = 0
    B21 = -one(K12) - A21   # 0
    B22 = one(K12) - A22    # 1 - i

    alpha_A = A12^2
    beta_A = A11^2
    alpha_B = B12^2
    beta_B = B11^2

    z1sq = z1^2
    z3sq = mobius_exact(z1sq, alpha_A, beta_A)
    z4sq = mobius_exact(z1sq, alpha_B, beta_B)
    num = beta_B - z3sq * conj12(alpha_B)
    den = alpha_B - z3sq * conj12(beta_B)
    z2sq = num * inv(den)

    z2, z2cands, z2C = principal_sqrt(z2sq)
    z3, z3cands, z3C = principal_sqrt(z3sq)
    z4, z4cands, z4C = principal_sqrt(z4sq)

    Zleft(z) = [one(K12) one(K12); z -z]
    Zright(z) = [one(K12) z; one(K12) -z]
    F2 = [one(K12) one(K12); one(K12) -one(K12)]
    A = [A11 A12; A21 A22]
    B = [B11 B12; B21 B22]
    Z1, Z2 = Zleft(z1), Zleft(z2)
    Z3, Z4 = Zright(z3), Zright(z4)
    half = inv(K12(2))
    scale(s, M) = s .* M
    top = hcat(F2, Z1, Z2)
    mid = hcat(Z3, scale(half, Z3 * A * Z1), scale(half, Z3 * B * Z2))
    bot = hcat(Z4, scale(half, Z4 * B * Z1), scale(half, Z4 * A * Z2))
    H = vcat(top, mid, bot)

    return (
        tag=lambda_tag, z1=z1, z2=z2, z3=z3, z4=z4,
        z2sq=z2sq, z3sq=z3sq, z4sq=z4sq,
        z2_branches=z2cands, z3_branches=z3cands, z4_branches=z4cands,
        z2C=z2C, z3C=z3C, z4C=z4C,
        A=A, B=B, H=H,
        alpha_A=alpha_A, beta_A=beta_A, alpha_B=alpha_B, beta_B=beta_B,
    )
end

function toC12(u)
    ζ = cis(π / 6)
    acc = 0.0 + 0im
    for k in 0:(degree(K12) - 1)
        acc += Float64(coeff(u, k)) * ζ^k
    end
    return acc
end

function H_to_complex(H)
    C = zeros(ComplexF64, 6, 6)
    for i in 1:6, j in 1:6
        C[i, j] = toC12(H[i, j])
    end
    return C
end

function fmt_elt(u)
    parts = String[]
    for k in 0:(degree(K12) - 1)
        ck = coeff(u, k)
        iszero(ck) && continue
        if k == 0
            push!(parts, string(ck))
        elseif k == 1
            push!(parts, isone(ck) ? "ζ12" : "$(ck)*ζ12")
        else
            push!(parts, isone(ck) ? "ζ12^$k" : "$(ck)*ζ12^$k")
        end
    end
    return isempty(parts) ? "0" : join(parts, " + ")
end

function main()
    lines = String[]
    push!(lines, "=== E1: exact Karlsson H at Dita λ=0 and λ=π/3 (Nemo) ===")
    push!(lines, "Date: $(Dates.format(now(), "yyyy-mm-ddTHH:MM:SS"))")
    push!(lines, "Number field: Q(ζ_12) = Q(i, √3), degree $(degree(K12))")
    push!(lines, "i = ζ12^3;  √3 = ζ12 + ζ12^11")
    push!(lines, "Conjugation: ζ12 ↦ ζ12^11 (complex conjugation).")
    push!(lines, "Square-root branch: the root in Q(ζ_12) closest to Julia's")
    push!(lines, "  principal sqrt of the complex image (matches build_karlsson_family).")
    push!(lines, "This is NOT a Groebner certificate. H only — not B3 / W1.")
    push!(lines, "MU vectors at D0 live in the larger field Q(ζ_24, √5); do not")
    push!(lines, "confuse that with the Hadamard field used here.")
    push!(lines, "")

    z1_0 = one(K12)
    z1_pi3 = Z12^2   # exp(i π / 3) = ζ_12^2
    jobs = [("lambda=0", z1_0, 0.0), ("lambda=pi/3", z1_pi3, pi / 3)]

    for (tag, z1, lam) in jobs
        data = build_exact_dita(tag, z1)
        Href = gauge_build_karlsson_family(acos(1 / sqrt(3)), pi / 4, lam)
        # Choose the (z2,z3,z4) branch triple matching the float construction.
        best_err = Inf
        best_H = data.H
        best_z = (data.z2, data.z3, data.z4)
        Zleft(z) = [one(K12) one(K12); z -z]
        Zright(z) = [one(K12) z; one(K12) -z]
        F2 = [one(K12) one(K12); one(K12) -one(K12)]
        half = inv(K12(2))
        scale(s, M) = s .* M
        for e2 in data.z2_branches, e3 in data.z3_branches, e4 in data.z4_branches
            Z1, Z2 = Zleft(data.z1), Zleft(e2)
            Z3, Z4 = Zright(e3), Zright(e4)
            top = hcat(F2, Z1, Z2)
            mid = hcat(Z3, scale(half, Z3 * data.A * Z1), scale(half, Z3 * data.B * Z2))
            bot = hcat(Z4, scale(half, Z4 * data.B * Z1), scale(half, Z4 * data.A * Z2))
            Hex = vcat(top, mid, bot)
            e = maximum(abs.(H_to_complex(Hex) - Href))
            if e < best_err
                best_err = e
                best_H = Hex
                best_z = (e2, e3, e4)
            end
        end
        Hnum = H_to_complex(best_H)
        n = size(Hnum, 1)
        unit = norm(Hnum * Hnum' - n * I(n))
        modu = maximum(abs.(abs.(Hnum) .- 1))
        push!(lines, "--- $tag ---")
        push!(lines, "  z1 = $(fmt_elt(data.z1))")
        push!(lines, "  z2^2 = $(fmt_elt(data.z2sq))")
        push!(lines, "  z3^2 = $(fmt_elt(data.z3sq))")
        push!(lines, "  z4^2 = $(fmt_elt(data.z4sq))")
        push!(lines, "  matched branches: z2=$(fmt_elt(best_z[1]))  z3=$(fmt_elt(best_z[2]))  z4=$(fmt_elt(best_z[3]))")
        push!(lines, "  n_roots in Q(ζ_12): z2=$(length(data.z2_branches)) z3=$(length(data.z3_branches)) z4=$(length(data.z4_branches))")
        push!(lines, @sprintf("  |H_exact - H_float| max = %.3e  (branch match to build_karlsson_family)", best_err))
        push!(lines, @sprintf("  CHM check on complex image: unitary_err=%.3e  modulus_err=%.3e", unit, modu))
        push!(lines, "  H (exact, matched branches, rows):")
        for i in 1:6
            row = join((fmt_elt(best_H[i, j]) for j in 1:6), ",  ")
            push!(lines, "    [$row]")
        end
        push!(lines, "")
    end

    push!(lines, "Field diagnosis vs A1 LLL: repo A1 fitted B3(λ=0.4) in Q(i,√2,√3).")
    push!(lines, "That field contains this H(λ=0) and H(λ=π/3), but not the D0 MU")
    push!(lines, "vectors (need √5 and ζ_24). Do not widen A1 at λ=0.4.")
    push!(lines, "No Groebner run. No claim that dim I = -1 over CC.")

    mkpath(RESULTS_DIR)
    open(OUT, "w") do io
        println.(Ref(io), lines)
    end
    foreach(println, lines)
    println("Wrote $OUT")
end

main()
