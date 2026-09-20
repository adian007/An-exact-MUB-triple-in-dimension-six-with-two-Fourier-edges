# E3: exact n_wit=1 export over a named number field, never CC.
#
# Nemo-only (no Oscar / no relative number_field): an element of Q(ζ_24, √5)
# is a pair (a, b) in Q(ζ_24)×Q(ζ_24) meaning a + b√5.
#
# H  = D_bc (Bengtsson eqs. 79–80), entries in Q(ζ_24) ⊂ Q(ζ_24, √5)
# B3 = F_D  (eqs. 78, 61), entries need √5 via b2 = (1-2i)/√5 = (1-2i)√5 / 5
#
# n_wit=1 only. Never R = CC[...]. No Groebner run.
#
# Usage: julia --project=. --compiled-modules=no scripts/julia/export_w1_exact.jl
# Output: symbolic_export/w1_D0_exact.m2
#         results/export_w1_exact.txt

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const RESULTS_DIR = joinpath(ROOT, "results")
const M2_DIR = joinpath(ROOT, "symbolic_export")

using LinearAlgebra
using Printf
using Dates
using Nemo

const OUT_LOG = joinpath(RESULTS_DIR, "export_w1_exact.txt")
const OUT_M2 = joinpath(M2_DIR, "w1_D0_exact.m2")

const K24, Z24 = cyclotomic_field(24)

"""a + b√5 with a, b ∈ Q(ζ_24)."""
struct Flt
    a::AbsSimpleNumFieldElem
    b::AbsSimpleNumFieldElem
end

Base.zero(::Type{Flt}) = Flt(zero(K24), zero(K24))
Base.one(::Type{Flt}) = Flt(one(K24), zero(K24))
Base.:+(x::Flt, y::Flt) = Flt(x.a + y.a, x.b + y.b)
Base.:-(x::Flt, y::Flt) = Flt(x.a - y.a, x.b - y.b)
Base.:-(x::Flt) = Flt(-x.a, -x.b)
Base.:*(x::Flt, y::Flt) = Flt(x.a * y.a + K24(5) * x.b * y.b, x.a * y.b + x.b * y.a)
function Base.inv(x::Flt)
    # 1/(a+b√5) = (a-b√5)/(a^2-5b^2)
    nrm = x.a * x.a - K24(5) * x.b * x.b
    nrmi = inv(nrm)
    return Flt(x.a * nrmi, -x.b * nrmi)
end
embed(u) = Flt(u, zero(K24))
iiF() = embed(Z24^6)

function conj24(u)
    n = degree(parent(u))
    s = zero(parent(u))
    for k in 0:(n - 1)
        ck = coeff(u, k)
        iszero(ck) && continue
        s += parent(u)(ck) * Z24^((23 * k) % 24)
    end
    return s
end

conjF(x::Flt) = Flt(conj24(x.a), conj24(x.b))  # √5 real

function fmt24(u)
    parts = String[]
    for k in 0:(degree(K24) - 1)
        ck = coeff(u, k)
        iszero(ck) && continue
        num = numerator(ck)
        den = denominator(ck)
        cs = den == 1 ? string(num) : "($(num))/$(den)"
        if k == 0
            push!(parts, cs)
        elseif k == 1
            # M2 variable is named zt: `zeta` is a built-in (Riemann zeta) in M2 >= 1.24
            push!(parts, cs == "1" ? "zt" : "($cs)*zt")
        else
            push!(parts, cs == "1" ? "zt^$k" : "($cs)*zt^$k")
        end
    end
    return isempty(parts) ? "0" : join(parts, " + ")
end

function m2_elt(x::Flt)
    fa, fb = fmt24(x.a), fmt24(x.b)
    if fb == "0"
        return fa == "0" ? "0" : "($fa)"
    elseif fa == "0"
        return "(($fb)*s5)"
    else
        return "(($fa) + ($fb)*s5)"
    end
end

function fourier3_exact()
    ω3 = Z24^8
    F = Matrix{Flt}(undef, 3, 3)
    for j in 0:2, k in 0:2
        F[j + 1, k + 1] = embed(ω3^(j * k))
    end
    return F
end

function adjointF(M)
    R, C = size(M)
    A = Matrix{Flt}(undef, C, R)
    for i in 1:R, j in 1:C
        A[j, i] = conjF(M[i, j])
    end
    return A
end

function exact_D_bc()
    ω = Z24
    o = embed(one(K24))
    C3 = [o embed(ω^6) embed(ω^6); embed(ω^6) o embed(ω^6); embed(ω^6) embed(ω^6) o]
    C4 = [embed(ω^15) embed(ω^3) embed(ω^3)
          embed(ω^3) embed(ω^15) embed(ω^3)
          embed(ω^3) embed(ω^3) embed(ω^15)]
    top = hcat(C3, C4)
    bot = hcat(C4, (-iiF()) * adjointF(C3))  # scalar * matrix: define below
    return vcat(top, bot)
end

Base.:*(s::Flt, M::AbstractMatrix{Flt}) = [s * M[i, j] for i in 1:size(M, 1), j in 1:size(M, 2)]
Base.:-(M::AbstractMatrix{Flt}) = map(-, M)

function exact_F_D()
    # b2 = (1-2i)/√5 = (1-2i)√5 / 5
    one24 = one(K24)
    num = one24 - K24(2) * Z24^6          # 1-2i
    b2 = Flt(zero(K24), num * inv(K24(5))) # 0 + ((1-2i)/5) √5
    Ddiag = [embed(Z24^9) * b2, embed(one24), embed(one24)]
    F3 = fourier3_exact()
    F3D = Matrix{Flt}(undef, 3, 3)
    for i in 1:3, j in 1:3
        F3D[i, j] = F3[i, j] * Ddiag[j]
    end
    top = hcat(F3, F3)
    bot = hcat(F3D, -F3D)
    return vcat(top, bot)
end

function w1_eqs17(H, B3)
    eqs = String[]
    for i in 1:5
        push!(eqs, "z$i * w$i - 1")
    end
    for B in (H, B3)
        for k in 1:6
            lhs = String[]
            rhs = String[]
            push!(lhs, m2_elt(conjF(B[1, k])))
            push!(rhs, m2_elt(B[1, k]))
            for j in 1:5
                cj = m2_elt(conjF(B[j + 1, k]))
                bj = m2_elt(B[j + 1, k])
                cj != "0" && push!(lhs, "$cj*z$j")
                bj != "0" && push!(rhs, "$bj*w$j")
            end
            push!(eqs, "($(join(lhs, " + ")))*(($(join(rhs, " + ")))) - 6")
        end
    end
    return eqs
end

function main()
    lines = String[]
    push!(lines, "=== E3: exact n_wit=1 export at D0, never CC ===")
    push!(lines, "Date: $(Dates.format(now(), "yyyy-mm-ddTHH:MM:SS"))")
    push!(lines, "Field: Q(ζ_24, √5) represented as (a,b) ↦ a+b√5, a,b ∈ Q(ζ_24)")
    push!(lines, "Nemo cyclotomic_field(24) only; no Oscar; no relative number_field.")
    push!(lines, "H  = D_bc (Bengtsson eqs. 79-80)")
    push!(lines, "B3 = F_D  (eqs. 78, 61), b2 = (1-2i)√5 / 5")
    push!(lines, "n_wit=1 only (17 eqs). n_wit=2 not started. Ring is not CC.")
    push!(lines, "No Groebner computation is performed.")
    push!(lines, "")

    H = exact_D_bc()
    B3 = exact_F_D()
    eqs17 = w1_eqs17(H, B3)
    push!(lines, "exported generators = $(length(eqs17)) (5 unimodular + 6 MU-H + 6 MU-B3)")

    mkpath(M2_DIR)
    open(OUT_M2, "w") do io
        println(io, "-- Exact n_wit=1 fourth-MUB witness at D0 / D_bc")
        println(io, "-- Field Q(zeta_24, sqrt(5)); NEVER CC")
        println(io, "-- Generated $(Dates.now()) by export_w1_exact.jl")
        println(io, "-- H = Bengtsson D_bc (eqs. 79-80); B3 = F_D (eqs. 78, 61)")
        println(io, "-- n_wit=1 only. Do not saturate with n_wit=2.")
        println(io, "-- This file is an exact-coefficient export, not a Groebner certificate.")
        println(io, "-- zt stands for zeta_24 (M2 reserves `zeta` for the built-in zeta function)")
        println(io, "kk = toField(QQ[zt, s5] / ideal(zt^8 - zt^4 + 1, s5^2 - 5));")
        println(io, "R = kk[z1, z2, z3, z4, z5, w1, w2, w3, w4, w5];")
        println(io, "witness = ideal(")
        for (i, eq) in enumerate(eqs17)
            sep = i < length(eqs17) ? "," : ""
            println(io, "  $eq$sep")
        end
        println(io, ");")
        println(io, "-- Do not run groebner(witness) in this research pass.")
        println(io, "-- print(\"#gens = \" | toString numgens witness | newline);")
    end
    body = read(OUT_M2, String)
    push!(lines, "Wrote $OUT_M2  ($(length(body)) bytes)")
    push!(lines, "Contains 'R = CC'? $(occursin("R = CC", body))")
    push!(lines, "Contains 'toField'? $(occursin("toField", body))")
    push!(lines, "First generator: $(eqs17[1])")
    push!(lines, "VERDICT: exact n_wit=1 export exists over Q(ζ_24, √5). No solver run.")

    mkpath(RESULTS_DIR)
    open(OUT_LOG, "w") do io
        println.(Ref(io), lines)
    end
    foreach(println, lines)
    println("Wrote $OUT_LOG")
    occursin("R = CC", body) && error("refused: export contained R = CC")
end

main()
