# Reconstruct B3(λ=0.4) at the Dita anchor via Nemo.jl LLL integer relations.
#
# For each of the 36 entries, seek z ≈ (1/q) Σ_k a_k β_k with
# β ∈ {1, √2, √3, √6, i, i√2, i√3, i√6} and a_k,q ∈ Z, |a_k|,|q| ≤ bound.
#
# Method `lll` (default): Nemo flint LLL on the standard integer-relation lattice
#   generators e_j ⊕ (⌊γ Re w_j⌉, ⌊γ Im w_j⌉),  w_j = β_j (j≤8), w_9 = −z.
# Method `ls-round`: previous least-squares + local neighborhood (weak baseline).
#
# Usage:
#   julia --project=. scripts/julia/reconstruct_b3_algebraic.jl --use-cache --method lll --basis-degree 8
#   julia --project=. scripts/julia/reconstruct_b3_algebraic.jl --use-cache --method lll --bits 512 --basis-degree 12
#
# Outputs: results/reconstruct_b3_algebraic.txt, .meta.txt

include(joinpath(@__DIR__, "_paths.jl"))
include(joinpath(ROOT, "src", "dita_third_mub_construction.jl"))

using LinearAlgebra
using Printf
using Dates
using Nemo

# -------------------- CLI --------------------
function parse_args(args)
    opts = Dict{String,Any}(
        "lambda" => 0.4,
        "bits" => 512,
        "basis_degree" => 8,
        "tol" => nothing,
        "gauge" => "col1",
        "scale" => "raw",
        "method" => "lll",          # lll | ls-round | both
        "out" => joinpath(RESULTS_DIR, "reconstruct_b3_algebraic.txt"),
        "use_cache" => false,
        "gamma_exp" => nothing,     # override γ = 10^gamma_exp
    )
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--lambda" && i < length(args)
            opts["lambda"] = parse(Float64, args[i+1]); i += 2
        elseif a == "--bits" && i < length(args)
            opts["bits"] = parse(Int, args[i+1]); i += 2
        elseif a == "--basis-degree" && i < length(args)
            opts["basis_degree"] = parse(Int, args[i+1]); i += 2
        elseif a == "--tol" && i < length(args)
            opts["tol"] = parse(Float64, args[i+1]); i += 2
        elseif a == "--gauge" && i < length(args)
            opts["gauge"] = args[i+1]; i += 2
        elseif a == "--scale" && i < length(args)
            opts["scale"] = args[i+1]; i += 2
        elseif a == "--method" && i < length(args)
            opts["method"] = args[i+1]; i += 2
        elseif a == "--gamma-exp" && i < length(args)
            opts["gamma_exp"] = parse(Int, args[i+1]); i += 2
        elseif a == "--use-cache"
            opts["use_cache"] = true; i += 1
        elseif a == "--out" && i < length(args)
            opts["out"] = args[i+1]; i += 2
        elseif a in ("-h", "--help")
            println("""
Usage: reconstruct_b3_algebraic.jl [--lambda 0.4] [--bits 512]
       [--basis-degree 8] [--tol T] [--gauge col1|none]
       [--scale raw|sqrt6] [--method lll|ls-round|both]
       [--gamma-exp N] [--use-cache] [--out PATH]
""")
            exit(0)
        else
            error("Unknown argument: $a")
        end
    end
    if opts["tol"] === nothing
        opts["tol"] = max(1e-12, 10.0^(-min(14, opts["bits"] ÷ 40)))
    end
    return opts
end

# -------------------- Algebraic basis --------------------
function algebraic_basis_labels()
    return [
        "1", "sqrt2", "sqrt3", "sqrt6",
        "i", "i*sqrt2", "i*sqrt3", "i*sqrt6",
    ]
end

function algebraic_basis_values(::Type{T}) where {T<:Real}
    s2 = sqrt(T(2)); s3 = sqrt(T(3)); s6 = sqrt(T(6))
    return Complex{T}[
        1, s2, s3, s6,
        im, im * s2, im * s3, im * s6,
    ]
end

function gauge_fix_columns!(B::AbstractMatrix{<:Complex})
    m, n = size(B)
    for j in 1:n
        pivot = zero(eltype(B))
        for i in 1:m
            if abs(B[i, j]) > 1e-14
                pivot = B[i, j]
                break
            end
        end
        if abs(pivot) > 0
            phase = conj(pivot) / abs(pivot)
            B[:, j] .*= phase
        end
    end
    return B
end

function _format_combo(a::Vector{Int}, q::Int, labels::Vector{String})
    terms = String[]
    for (k, ak) in enumerate(a)
        ak == 0 && continue
        push!(terms, ak == 1 ? labels[k] :
              ak == -1 ? "-$(labels[k])" : "$(ak)*$(labels[k])")
    end
    isempty(terms) && return "0"
    body = join(terms, " + ")
    body = replace(body, "+ -" => "- ")
    return q == 1 ? body : "($body)/$q"
end

# -------------------- Nemo LLL integer relation --------------------
"""
Find integers (a_1..a_n, q) with q ≠ 0, |a_k|,|q| ≤ bound, minimizing
|Σ a_k β_k − q z| via flint LLL (Nemo) on the relation lattice.
"""
function fit_entry_lll(z::Complex; bound::Int=8, bits::Int=512, tol::Float64=1e-12,
                       gamma_exp::Union{Nothing,Int}=nothing)
    labels = algebraic_basis_labels()
    n = length(labels)
    setprecision(BigFloat, bits) do
        β = algebraic_basis_values(BigFloat)
        zb = Complex{BigFloat}(BigFloat(real(z)), BigFloat(imag(z)))
        gexp = gamma_exp === nothing ? 16 : gamma_exp  # Float64-sourced B3: ~10^16 is enough
        γ = BigFloat(10)^gexp

        # Lattice generators as ROWS (Nemo.lll reduces the row lattice).
        # Each of the n+1 rows is:
        #   (a_1,...,a_n,q,  γ(Σ a Reβ − q Re z), γ(Σ a Imβ − q Im z))
        # encoded as the standard basis of Z^{n+1} plus the two weighted relation coords.
        nrows = n + 1
        ncols = n + 3
        Ment = zeros(BigInt, nrows, ncols)
        for i in 1:nrows
            Ment[i, i] = 1
        end
        for i in 1:n
            Ment[i, n + 2] = BigInt(round(γ * real(β[i])))
            Ment[i, n + 3] = BigInt(round(γ * imag(β[i])))
        end
        Ment[n + 1, n + 2] = BigInt(round(γ * (-real(zb))))
        Ment[n + 1, n + 3] = BigInt(round(γ * (-imag(zb))))

        A = matrix(Nemo.ZZ, Ment)
        L = Nemo.lll(A)

        best = nothing
        best_err = Inf
        for i in 1:nrows
            ok_int = true
            coeffs = zeros(Int, n + 1)
            for j in 1:(n + 1)
                try
                    coeffs[j] = Int(L[i, j])
                catch
                    ok_int = false
                    break
                end
            end
            ok_int || continue
            q = coeffs[n + 1]
            a = coeffs[1:n]
            q == 0 && continue
            if q < 0
                q = -q
                a = -a
            end
            maximum(abs, a) > bound && continue
            q > bound && continue
            approx = sum(BigFloat(a[k]) * β[k] for k in 1:n) / BigFloat(q)
            err = Float64(abs(zb - approx))
            if err < best_err
                best_err = err
                best = (
                    ok = err < tol,
                    err = err,
                    coeffs = a,
                    denom = q,
                    expr = _format_combo(a, q, labels),
                    method = "nemo-lll",
                )
            end
        end

        # Also try a few smaller γ in case of overflow / scaling issues
        if best === nothing || !best.ok
            for gexp2 in unique([gexp - 10, gexp - 20, max(10, gexp ÷ 2), 12, 16, 20, 24])
                gexp2 < 8 && continue
                γ2 = BigFloat(10)^gexp2
                Ment2 = zeros(BigInt, nrows, ncols)
                for i in 1:nrows
                    Ment2[i, i] = 1
                end
                for i in 1:n
                    Ment2[i, n + 2] = BigInt(round(γ2 * real(β[i])))
                    Ment2[i, n + 3] = BigInt(round(γ2 * imag(β[i])))
                end
                Ment2[n + 1, n + 2] = BigInt(round(γ2 * (-real(zb))))
                Ment2[n + 1, n + 3] = BigInt(round(γ2 * (-imag(zb))))
                L2 = Nemo.lll(matrix(Nemo.ZZ, Ment2))
                for i in 1:nrows
                    coeffs = [Int(L2[i, j]) for j in 1:(n + 1)]
                    q = coeffs[n + 1]
                    a = coeffs[1:n]
                    q == 0 && continue
                    if q < 0
                        q = -q; a = -a
                    end
                    maximum(abs, a) > bound && continue
                    q > bound && continue
                    approx = sum(BigFloat(a[k]) * β[k] for k in 1:n) / BigFloat(q)
                    err = Float64(abs(zb - approx))
                    if err < best_err
                        best_err = err
                        best = (
                            ok = err < tol,
                            err = err,
                            coeffs = a,
                            denom = q,
                            expr = _format_combo(a, q, labels),
                            method = "nemo-lll",
                        )
                    end
                end
            end
        end
        return best
    end
end

"""Legacy weak baseline: LS + local integer neighborhood (not a certificate of absence)."""
function fit_entry_ls_round(z::Complex; bound::Int=8, tol::Float64=1e-12)
    β = algebraic_basis_values(Float64)
    labels = algebraic_basis_labels()
    n = length(β)
    zr = Float64(real(z)); zi = Float64(imag(z))
    zc = Complex{Float64}(z)
    best = nothing
    best_err = Inf
    M = zeros(Float64, 2, n)
    for k in 1:n
        M[1, k] = real(β[k])
        M[2, k] = imag(β[k])
    end
    nbhd = bound <= 4 ? (-2:2) : (-1:1)
    for q in 1:bound
        target = q * [zr, zi]
        a_ls = M \ target
        a0 = round.(Int, a_ls)
        for δ in Iterators.product(ntuple(_ -> nbhd, n)...)
            a = a0 .+ collect(Int, δ)
            maximum(abs, a) > bound && continue
            approx = sum(a[k] * β[k] for k in 1:n) / q
            err = abs(zc - approx)
            if err < best_err
                best_err = err
                best = (ok = err < tol, err = err, coeffs = a, denom = q,
                        expr = _format_combo(a, q, labels), method = "ls-round")
            end
        end
    end
    return best
end

function fit_entry(z::Complex; bound::Int=8, bits::Int=512, tol::Float64=1e-12,
                   method::String="lll", gamma_exp=nothing)
    if method == "ls-round"
        return fit_entry_ls_round(z; bound=bound, tol=tol)
    elseif method == "both"
        a = fit_entry_lll(z; bound=bound, bits=bits, tol=tol, gamma_exp=gamma_exp)
        b = fit_entry_ls_round(z; bound=bound, tol=tol)
        if a === nothing
            return b
        elseif b === nothing
            return a
        else
            return a.err <= b.err ? a : b
        end
    else
        # default / "lll" / legacy "pslq" alias → Nemo LLL
        return fit_entry_lll(z; bound=bound, bits=bits, tol=tol, gamma_exp=gamma_exp)
    end
end

# -------------------- Main --------------------
function main(args)
    opts = parse_args(args)
    λ = opts["lambda"]
    bits = opts["bits"]
    deg = opts["basis_degree"]
    tol = opts["tol"]
    method = opts["method"]
    gamma_exp = opts["gamma_exp"]

    mkpath(RESULTS_DIR)
    out_path = opts["out"]

    println("="^72)
    println("B3 algebraic reconstruction @ Dita λ=$λ  [Nemo LLL]")
    println("bits=$bits  coeff-bound=$deg  tol=$tol  method=$method")
    println("gauge=$(opts["gauge"])  scale=$(opts["scale"])  gamma_exp=$(gamma_exp)")
    println("="^72)

    t0 = time()
    cache_bin = joinpath(RESULTS_DIR, "reconstruct_b3_B3_lambda$(λ).csv")
    use_cache = opts["use_cache"] && isfile(cache_bin)
    B3 = Matrix{ComplexF64}(undef, 6, 6)
    clique_idx = Int[]
    n_pool = 0
    max_clique = 0

    if use_cache
        println("\n[1] Loading cached B3 from $cache_bin ...")
        rows = readlines(cache_bin)
        for part in split(rows[1], ',')
            k, v = split(part, '=')
            k = strip(k); v = strip(v)
            k == "max_clique" && (max_clique = parse(Int, v))
            k == "n_pool" && (n_pool = parse(Int, v))
        end
        clique_idx = parse.(Int, split(split(rows[2], '=')[2], ';'))
        for i in 1:6
            cols = split(rows[2 + i], ';')
            for j in 1:6
                re, im_ = parse.(Float64, split(cols[j], ','))
                B3[i, j] = ComplexF64(re, im_)
            end
        end
        println("    max_clique=$max_clique  n_pool=$n_pool  indices=$clique_idx")
    else
        println("\n[1] Extracting third MUB via dita_third_mub_at($λ) ...")
        d = dita_third_mub_at(λ; verbose=false)
        B3 = Complex{Float64}.(d.third.B3)
        clique_idx = d.third.indices
        n_pool = d.extension.n_pool
        max_clique = d.extension.max_clique
        println("    max_clique=$max_clique  n_pool=$n_pool")
        open(cache_bin, "w") do io
            println(io, "lambda=$λ,max_clique=$max_clique,n_pool=$n_pool")
            println(io, "indices=", join(clique_idx, ";"))
            for i in 1:6
                println(io, join([string(real(B3[i, j]), ",", imag(B3[i, j])) for j in 1:6], ";"))
            end
        end
        println("    cached B3 → $cache_bin")
    end

    vecs = [B3[:, j] for j in 1:6]
    hp = verify_clique_hp(vecs; bits=min(bits, 400))
    println("    HP verify: ortho_ok=$(hp.ortho_ok) mu_ok=$(hp.mu_ok) ortho_max=$(hp.ortho_max)")

    if opts["gauge"] == "col1"
        gauge_fix_columns!(B3)
        println("    Applied column gauge")
    end
    if opts["scale"] == "sqrt6"
        B3 .*= sqrt(6)
        println("    Scaling by √6")
    end

    labels = algebraic_basis_labels()
    println("\n[2] Fitting 36 entries with method=$method (bound=$deg) ...")
    n_ok = 0
    n_total = 36
    open(out_path, "w") do io
        println(io, "B3 algebraic reconstruction (Nemo LLL)")
        println(io, "timestamp = ", Dates.now())
        println(io, "lambda = ", λ)
        println(io, "bits = ", bits)
        println(io, "basis_degree = ", deg)
        println(io, "tol = ", tol)
        println(io, "gauge = ", opts["gauge"])
        println(io, "scale = ", opts["scale"])
        println(io, "method = ", method)
        println(io, "gamma_exp = ", gamma_exp)
        println(io, "engine = Nemo.lll (flint)")
        println(io, "basis = ", join(labels, ", "))
        println(io, "hp_ortho_ok = ", hp.ortho_ok, " hp_mu_ok = ", hp.mu_ok)
        println(io, "clique_indices = ", clique_idx)
        println(io)
        println(io, "PRIOR_NOTE: earlier ls-round/brute runs are weak baselines only;")
        println(io, "  this file supersedes them for claim-ledger purposes when method=lll.")
        println(io)
        println(io, "entry_ij | ok | err | method | expression | numeric")
        println(io, "-"^72)

        for i in 1:6, j in 1:6
            z = B3[i, j]
            fit = fit_entry(z; bound=deg, bits=bits, tol=tol, method=method, gamma_exp=gamma_exp)
            ok = fit !== nothing && fit.ok
            ok && (n_ok += 1)
            if fit === nothing
                @printf("%d,%d  FAIL  (no candidate within bound)\n", i, j)
                println(io, @sprintf("%d,%d | false | Inf | none | — | %s", i, j, string(z)))
            else
                status = ok ? "OK" : "MISS"
                @printf("%d,%d  %s  err=%.3e  [%s]  %s\n",
                        i, j, status, fit.err, fit.method, fit.expr)
                println(io, @sprintf("%d,%d | %s | %.6e | %s | %s | %s",
                                     i, j, string(ok), fit.err, fit.method, fit.expr, string(z)))
            end
        end

        println(io)
        println(io, "="^40)
        println(io, "FIT_RATE = $n_ok / $n_total")
        println(io, "elapsed_s = ", round(time() - t0; digits=2))
        verdict = n_ok == n_total ? "FULL_ALGEBRAIC_FIT" :
                  n_ok >= 30 ? "NEARLY_FULL_FIT" :
                  n_ok >= 18 ? "PARTIAL_FIT" : "NO_SMALL_BASIS_FIT"
        println(io, "verdict = ", verdict)
        println(io)
        println(io, "SEARCH_QUALITY: Nemo flint LLL integer-relation lattice;")
        println(io, "  coeff bound=$deg, bits=$bits, tol=$tol.")
        if verdict == "NO_SMALL_BASIS_FIT"
            println(io, "  Null result under THIS search — may still miss if true")
            println(io, "  representation needs larger coeffs or a bigger field.")
        end
    end

    elapsed = round(time() - t0; digits=2)
    verdict = n_ok == n_total ? "FULL_ALGEBRAIC_FIT" :
              n_ok >= 30 ? "NEARLY_FULL_FIT" :
              n_ok >= 18 ? "PARTIAL_FIT" : "NO_SMALL_BASIS_FIT"
    println("\nFIT RATE: $n_ok / $n_total   verdict=$verdict   elapsed=$(elapsed)s")
    println("wrote: $out_path")

    open(joinpath(RESULTS_DIR, "reconstruct_b3_algebraic.meta.txt"), "w") do io
        println(io, "fit_ok=$n_ok")
        println(io, "fit_total=$n_total")
        println(io, "basis_degree=$deg")
        println(io, "bits=$bits")
        println(io, "lambda=$λ")
        println(io, "verdict=$verdict")
        println(io, "tol=$tol")
        println(io, "method=$method")
        println(io, "engine=Nemo.lll")
        println(io, "gauge=$(opts["gauge"])")
        println(io, "scale=$(opts["scale"])")
        println(io, "search_quality=nemo_flint_lll")
    end
    return (n_ok=n_ok, n_total=n_total, verdict=verdict)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
