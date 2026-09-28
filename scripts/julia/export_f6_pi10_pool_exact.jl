using LinearAlgebra
using Nemo
using Printf

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const OUT_M2 = joinpath(ROOT, "symbolic_export", "f6_pi10_pool_exact.m2")
const OUT_LOG = joinpath(ROOT, "results", "f6_pi10_pool_exact_export.txt")

const K60, ZETA60 = cyclotomic_field(60)
const OMEGA = ZETA60^20
const PARAMETER = ZETA60^3
const ONE_K = one(K60)

zleft(a) = [ONE_K ONE_K; a -a]
zright(a) = [ONE_K a; ONE_K -a]

function exact_f6_pi10()
    f2 = [ONE_K ONE_K; ONE_K -ONE_K]
    a = [OMEGA OMEGA^2; OMEGA -OMEGA^2]
    b = -f2 - a
    zr = zright(ONE_K)
    zl_lambda = zleft(PARAMETER)
    zl_two = zleft(OMEGA)

    return [
        f2 zl_lambda zl_two
        zr zr * a * zl_lambda / 2 zr * b * zl_two / 2
        zr zr * b * zl_lambda / 2 zr * a * zl_two / 2
    ]
end

function conjugate60(x)
    return sum(parent(x)(coeff(x, k)) * ZETA60^(59k)
               for k in 0:(degree(K60) - 1); init = zero(K60))
end

function field_expression(x)
    terms = String[]
    for k in 0:(degree(K60) - 1)
        c = coeff(x, k)
        iszero(c) && continue
        n, d = numerator(c), denominator(c)
        coefficient = d == 1 ? string(n) : "($n/$d)"
        push!(terms, k == 0 ? coefficient : "($coefficient)*zz^$k")
    end
    return isempty(terms) ? "0" : "(" * join(terms, " + ") * ")"
end

function pool_equations(h)
    eqs = String[]
    for j in 1:5
        push!(eqs, "x$j*y$j-1")
    end

    for k in 1:5
        lhs = [field_expression(conjugate60(h[1, k]))]
        rhs = [field_expression(h[1, k])]
        for j in 1:5
            push!(lhs, "$(field_expression(conjugate60(h[j + 1, k])))*x$j")
            push!(rhs, "$(field_expression(h[j + 1, k]))*y$j")
        end
        push!(eqs, "(" * join(lhs, "+") * ")*(" * join(rhs, "+") * ")-6")
    end
    return eqs
end

function main()
    h = exact_f6_pi10()
    gram = h * map(conjugate60, transpose(h))
    is_hadamard = all(iszero(gram[i, j] - (i == j ? K60(6) : zero(K60)))
                      for i in 1:6, j in 1:6)
    is_hadamard || error("Exact F6(pi/10) Gram matrix is not 6I")
    all(iszero(x * conjugate60(x) - ONE_K) for x in h) ||
        error("Exact F6(pi/10) has a non-unimodular entry")

    eqs = pool_equations(h)
    length(eqs) == 10 || error("Expected 10 equations in the MU-pool ideal")

    mkpath(dirname(OUT_M2))
    mkpath(dirname(OUT_LOG))
    open(OUT_M2, "w") do io
        println(io, "-- Exact F6(theta=0, lambda=pi/10) per-vector MU ideal")
        println(io, "-- Coefficient field Q(zeta_60), no floating-point coefficients")
        println(io, "-- One MU equation is omitted by the unitary sum identity")
        println(io, "kk = toField(QQ[zz]/ideal(zz^16 + zz^14 - zz^10 - zz^8 - zz^6 + zz^2 + 1));")
        println(io, "R = kk[x1,x2,x3,x4,x5,y1,y2,y3,y4,y5];")
        println(io, "I = ideal(")
        for (i, eq) in enumerate(eqs)
            println(io, "  ", eq, i == length(eqs) ? "" : ",")
        end
        println(io, ");")
        println(io, "print(\"Starting exact Groebner calculation\" | newline);")
        println(io, "G = gb I;")
        println(io, "print(\"Exact MU-pool ideal computed; dimension=\" | toString dim I | newline);")
    end

    open(OUT_LOG, "w") do io
        println(io, "Exact F6(theta=0, lambda=pi/10) matrix: 36/36 entries unimodular; H*H^*=6I.")
        println(io, "Coefficient field: Q(zeta_60), degree $(degree(K60)).")
        println(io, "Exact per-vector MU-pool ideal exported: 10 equations, 10 variables.")
        println(io, "Macaulay2 input: $OUT_M2")
        println(io, "No Gröbner result is recorded by the exporter; the M2 script runs gb I.")
    end

    @printf("Exact Hadamard check passed over Q(zeta_60); wrote %s\n", OUT_M2)
end

main()
