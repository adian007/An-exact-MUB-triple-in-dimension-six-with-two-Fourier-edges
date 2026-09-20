# E0: identify which Karlsson (theta_D, pi/4, lambda) is CHM-equivalent to textbook D0.
# No HomotopyContinuation (Application Control safe).
# Usage: julia --project=. scripts/julia/identify_d0_in_karlsson.jl
#        Prefer the Python sibling if this host cannot load Julia deps.

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(ROOT, "src", "karlsson_gauge_only.jl"))
include(joinpath(ROOT, "src", "brierley_weigert_notes.jl"))

using LinearAlgebra
using Printf
using Dates

const RESULTS_DIR = joinpath(ROOT, "results")
const OUT = joinpath(RESULTS_DIR, "identify_d0_in_karlsson_julia.txt")
const DITA_THETA = acos(1 / sqrt(3))
const DITA_PHI = pi / 4
const LAM_CERT = (0.0, 0.4, pi / 3, 2pi / 3, pi, 4pi / 3, 5pi / 3)
const TOL_EQUIV = 1e-8

function _phase(z)
    abs(z) < 1e-15 && return 1.0 + 0im
    return z / abs(z)
end

"""Column-permutation residual after dephasing (same logic as chm_equivalence.py)."""
function residual_for_perm(H1, H2, perm)
    H2p = H2[:, perm]
    Dr = ones(ComplexF64, 6)
    for i in 1:6
        if abs(H2p[i, 1]) > 1e-15
            Dr[i] = _phase(H1[i, 1] / H2p[i, 1])
        end
    end
    Dc = ones(ComplexF64, 6)
    for j in 2:6
        ratios = ComplexF64[]
        for i in 1:6
            if abs(H2p[i, j]) > 1e-15 && abs(Dr[i]) > 1e-15
                push!(ratios, H1[i, j] / (Dr[i] * H2p[i, j]))
            end
        end
        if !isempty(ratios)
            r0 = ratios[1]
            Dc[j] = _phase(mean(_phase(r / r0) for r in ratios) * r0)
        end
    end
    Hfit = Dr .* H2p .* transpose(Dc)
    return norm(H1 - Hfit)
end

function mean(iter)
    s = 0.0 + 0im
    n = 0
    for x in iter
        s += x
        n += 1
    end
    return s / n
end

function permutations6()
    acc = Vector{Vector{Int}}()
    used = falses(6)
    cur = zeros(Int, 6)
    function rec(k)
        if k > 6
            push!(acc, copy(cur))
            return
        end
        for i in 1:6
            used[i] && continue
            used[i] = true
            cur[k] = i
            rec(k + 1)
            used[i] = false
        end
    end
    rec(1)
    return acc
end

const PERMS6 = permutations6()

function chm_residual(H1, H2)
    best = (residual=Inf, variant="H", perm=Int[])
    variants = (("H", H2), ("H.T", collect(transpose(H2))),
                ("conj(H)", conj.(H2)), ("conj(H).T", conj.(collect(transpose(H2)))))
    D1 = gauge_dephase(H1)
    for (vname, H2v) in variants
        D2 = gauge_dephase(H2v)
        for perm in PERMS6
            res = residual_for_perm(D1, D2, perm)
            if res < best.residual
                best = (residual=res, variant=vname, perm=perm)
            end
        end
    end
    return best
end

function lambda_grid()
    algebraic = [0.0, pi/12, pi/8, pi/6, pi/4, pi/3, pi/2, 2pi/3, 3pi/4, pi,
                 4pi/3, 5pi/3, 2pi]
    dense = range(0, 2pi; length=97)
    vals = sort(collect(Set(round(x; digits=12) for x in Iterators.flatten((algebraic, dense, LAM_CERT)))))
    return vals
end

function fmt_perm(p)
    isempty(p) && return "None"
    return "(" * join(p, " ") * ")"
end

function main()
    D0 = dita_D0()
    lines = String[]
    push!(lines, "=== E0 (Julia): identify textbook Dita D0 inside Karlsson Dita slice ===")
    push!(lines, "Date: $(Dates.format(now(), "yyyy-mm-ddTHH:MM:SS"))")
    push!(lines, "No HomotopyContinuation. CHM residual matches chm_equivalence.py.")
    push!(lines, "")

    d0d = chm_defect(D0)
    dbcd = chm_defect(dita_D_bc())
    push!(lines, @sprintf("D0 CHM ok=%s  unitary=%.3e  modulus=%.3e", d0d.ok, d0d.unitary, d0d.modulus))
    push!(lines, @sprintf("D_bc CHM ok=%s  unitary=%.3e  modulus=%.3e", dbcd.ok, dbcd.unitary, dbcd.modulus))
    rdbc = chm_residual(D0, dita_D_bc())
    push!(lines, @sprintf("D0 vs D_bc: residual=%.6e  variant=%s  perm=%s",
                          rdbc.residual, rdbc.variant, fmt_perm(rdbc.perm)))
    push!(lines, "")
    push!(lines, "--- residual(Karlsson(theta_D, pi/4, lambda), D0) [cert / matches / <1e-3] ---")

    best = (lambda=NaN, residual=Inf, variant="?", perm=Int[])
    matches = NamedTuple[]
    cert_rows = NamedTuple[]
    for lam in lambda_grid()
        H = gauge_build_karlsson_family(DITA_THETA, DITA_PHI, lam)
        r = chm_residual(H, D0)
        rec = (lambda=lam, residual=r.residual, variant=r.variant, perm=r.perm,
               equivalent=r.residual < TOL_EQUIV)
        if rec.residual < best.residual
            best = rec
        end
        in_cert = any(abs(lam - c) < 1e-12 for c in LAM_CERT)
        if rec.equivalent || in_cert || rec.residual < 1e-3
            mark = in_cert ? "  LAMBDA_CERT" : ""
            push!(lines, @sprintf("  lambda=%.12g  residual=%.6e  variant=%s  perm=%s  equiv=%s%s",
                                  lam, rec.residual, rec.variant, fmt_perm(rec.perm), rec.equivalent, mark))
        end
        rec.equivalent && push!(matches, rec)
        in_cert && push!(cert_rows, rec)
    end

    push!(lines, "")
    push!(lines, "=== Best match ===")
    push!(lines, @sprintf("lambda=%.12g  residual=%.6e  variant=%s  perm=%s  equiv=%s",
                          best.lambda, best.residual, best.variant, fmt_perm(best.perm),
                          best.residual < TOL_EQUIV))
    push!(lines, "n_matches: $(length(matches))")
    in_cert_match = any(r.equivalent for r in cert_rows)
    push!(lines, "")
    push!(lines, "=== T3 Lambda_cert vs D0 ===")
    for rec in cert_rows
        push!(lines, @sprintf("  lambda=%.12g  residual=%.6e  equiv=%s",
                              rec.lambda, rec.residual, rec.equivalent))
    end
    push!(lines, "")
    if in_cert_match
        push!(lines, "VERDICT: a Lambda_cert point is CHM-equivalent to textbook D0.")
        push!(lines, "T3 at that point is a corollary of BW PRA 79, 052316 / arXiv:0901.4051.")
    else
        push!(lines, "VERDICT: no Lambda_cert point matched D0 at the stated tolerance.")
    end

    mkpath(RESULTS_DIR)
    open(OUT, "w") do io
        println.(Ref(io), lines)
    end
    foreach(println, lines)
    println("Wrote $OUT")
end

main()
