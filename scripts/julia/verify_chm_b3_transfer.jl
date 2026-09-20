# Confirm CHM-equivalence transfers third-MUB / W1 structure under the B3 map.
#
# For each pair (λ, λ+π) with λ ∈ {0, π/3, 2π/3}:
#   1. Recover CHM map dephase(H(λ+π)) ≈ D_r · dephase(H(λ))[:,P] · D_c
#   2. Take one (or all) size-6 clique B3 at λ
#   3. Form B3' with columns D_r * B3[:,j] (renormalized)
#   4. Check: B3' ONB, MU to I, MU to H(λ+π)
#   5. Certify W1(H(λ+π), B3') empty
#
# Usage: julia --project=. scripts/julia/verify_chm_b3_transfer.jl
# Output: results/verify_chm_b3_transfer.txt

include(joinpath(@__DIR__, "_paths.jl"))
include(joinpath(ROOT, "src", "dita_third_mub_construction.jl"))
include(joinpath(ROOT, "scripts", "julia", "locus_classification.jl"))

using Printf, Dates, LinearAlgebra, Statistics

const OUT = joinpath(RESULTS_DIR, "verify_chm_b3_transfer.txt")
const PAIRS = [(0.0, π), (π / 3, 4π / 3), (2π / 3, 5π / 3)]
const MU_TOL = 1e-8
const ORTHO_TOL = 1e-8

function chm_equivalence_map(H1, H2; tol=CHM_TOL)
    best = (residual=Inf, variant="", perm=Int[], D_r=ComplexF64[], D_c=ComplexF64[],
            equivalent=false)
    D1 = dephase_chm(H1)
    for (vname, H2v) in [("H", H2), ("H.T", Matrix(transpose(H2))),
                          ("conj(H)", conj.(H2)), ("conj(H).T", Matrix(transpose(conj.(H2))))]
        D2 = dephase_chm(H2v)
        for perm in _PERMS6
            H2p = D2[:, perm]
            N = 6
            D_r = ones(ComplexF64, N)
            for i in 1:N
                abs(H2p[i, 1]) > 1e-15 && (D_r[i] = _phase(D1[i, 1] / H2p[i, 1]))
            end
            D_c = ones(ComplexF64, N)
            for j in 2:N
                ratios = ComplexF64[]
                for i in 1:N
                    if abs(H2p[i, j]) > 1e-15 && abs(D_r[i]) > 1e-15
                        push!(ratios, D1[i, j] / (D_r[i] * H2p[i, j]))
                    end
                end
                if !isempty(ratios)
                    r0 = ratios[1]
                    phases = [_phase(r / r0) for r in ratios]
                    D_c[j] = mean(phases) * r0
                    abs(D_c[j]) > 1e-15 && (D_c[j] = _phase(D_c[j]))
                end
            end
            H_fit = D_r .* H2p .* transpose(D_c)
            res = norm(D1 - H_fit)
            if res < best.residual
                best = (residual=res, variant=vname, perm=collect(perm), D_r=D_r, D_c=D_c,
                        equivalent=res < tol)
            end
        end
    end
    return best
end

function apply_variant(H, variant)
    variant == "H" && return Matrix(H)
    variant == "H.T" && return Matrix(transpose(H))
    variant == "conj(H)" && return conj.(H)
    variant == "conj(H).T" && return Matrix(transpose(conj.(H)))
    error("unknown variant $variant")
end

function transform_b3(B3, D_r)
    B = Matrix{ComplexF64}(undef, 6, 6)
    for j in 1:6
        v = D_r .* B3[:, j]
        B[:, j] = v / norm(v)
    end
    return B
end

function mu_max_to_basis(B, H)
    m = 0.0
    for j in 1:6, k in 1:6
        m = max(m, abs(abs(dot(B[:, j], H[:, k] / sqrt(6)))^2 - 1 / 6))
    end
    return m
end

function mu_max_to_I(B)
    m = 0.0
    for j in 1:6, i in 1:6
        m = max(m, abs(abs(B[i, j])^2 - 1 / 6))
    end
    return m
end

function ortho_max(B)
    m = 0.0
    for a in 1:6, b in a:6
        target = a == b ? 1.0 : 0.0
        m = max(m, abs(dot(B[:, a], B[:, b]) - target))
    end
    return m
end

function all_six_cliques(H; ortho_tol=1e-8)
    pool, = generate_candidate_pool_fresh(H; verbose=false)
    pool = deduplicate_pool(pool)
    g = _orthogonality_graph(pool; ortho_tol=ortho_tol)
    uniq = Dict{Vector{Int}, Vector{Int}}()
    for c in maximal_cliques(g)
        length(c) < 6 && continue
        idx = collect(c[1:6])
        uniq[sort(idx)] = idx
    end
    return pool, collect(values(uniq))
end

function main()
    mkpath(RESULTS_DIR)
    lines = String[]
    push!(lines, "=== CHM-equivalence B3 transfer check (λ ↔ λ+π) ===")
    push!(lines, "timestamp = $(Dates.now())")
    push!(lines, "Map: dephase(H_tgt) ≈ D_r · dephase(variant(H_src))[:,P] · D_c")
    push!(lines, "B3' columns := normalize(D_r .* B3[:,j])")
    push!(lines, "")

    all_ok = true
    for (λ_src, λ_tgt) in PAIRS
        H_src = build_karlsson_family(DITA_THETA, DITA_PHI, λ_src)
        H_tgt = build_karlsson_family(DITA_THETA, DITA_PHI, λ_tgt)
        mp = chm_equivalence_map(H_tgt, H_src)  # express tgt in terms of src
        push!(lines, @sprintf("--- λ_src=%.6f → λ_tgt=%.6f ---", λ_src, λ_tgt))
        push!(lines, @sprintf("  CHM map: residual=%.3e equiv=%s variant=%s perm=%s",
                              mp.residual, mp.equivalent, mp.variant, mp.perm))
        if !mp.equivalent
            push!(lines, "  FAIL: not CHM-equivalent")
            all_ok = false
            continue
        end

        # Sanity: reconstruct dephased tgt from src
        H_src_v = apply_variant(H_src, mp.variant)
        D_src = dephase_chm(H_src_v)
        D_tgt = dephase_chm(H_tgt)
        H_fit = mp.D_r .* D_src[:, mp.perm] .* transpose(mp.D_c)
        push!(lines, @sprintf("  dephase fit residual=%.3e", norm(D_tgt - H_fit)))

        pool, cliques = all_six_cliques(H_src)
        push!(lines, @sprintf("  src pool=%d n_6cliques=%d (checking all)", length(pool), length(cliques)))

        n_pass = 0
        for (k, idx) in enumerate(cliques)
            B3 = hcat([pool[i] for i in idx]...)
            B3p = transform_b3(B3, mp.D_r)
            o = ortho_max(B3p)
            mi = mu_max_to_I(B3p)
            mh = mu_max_to_basis(B3p, H_tgt)
            # Also MU of B3p to dephased-equivalent column set of H_tgt
            struct_ok = o < ORTHO_TOL && mi < MU_TOL && mh < MU_TOL

            info = certify_fourth_mub_witness_at_H(H_tgt; n_wit=1,
                clique_indices=idx, B3=B3p, verbose=false)
            empty = !info.skipped_solve && info.n_certified > 0 &&
                    info.n_certified_pass_full_residual == 0
            ok = struct_ok && empty
            n_pass += ok ? 1 : 0
            all_ok &= ok
            push!(lines, @sprintf(
                "  clique %d: ortho=%.2e muI=%.2e muH=%.2e struct=%s W1_empty=%s (cert=%d fullpass=%d)",
                k, o, mi, mh, struct_ok, empty, info.n_certified, info.n_certified_pass_full_residual))
        end
        push!(lines, @sprintf("  pair summary: %d/%d transferred cliques OK", n_pass, length(cliques)))
        push!(lines, "")
    end

    push!(lines, all_ok ?
        "PASS: for all three λ↔λ+π pairs, every transformed B3' is a third MUB for H(λ+π) and W1 is empty." :
        "FAIL: at least one transform/W1 check failed.")
    text = join(lines, "\n") * "\n"
    write(OUT, text)
    println(text)
    println("Wrote ", OUT)
end

main()
