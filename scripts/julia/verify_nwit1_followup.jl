# Follow-up: (a) confirm H(λ) ~ H(λ+π) CHM pattern; (b) certify W1 on ALL
# size-6 cliques at λ=0.4; (c) certify λ=0 with fresh B3.
#
# Usage: julia --project=. scripts/julia/verify_nwit1_followup.jl
# Output: results/verify_nwit1_followup.txt

include(joinpath(@__DIR__, "_paths.jl"))
include(joinpath(ROOT, "src", "dita_third_mub_construction.jl"))
include(joinpath(ROOT, "scripts", "julia", "locus_classification.jl"))

using Printf, Dates, Random

const OUT = joinpath(RESULTS_DIR, "verify_nwit1_followup.txt")

function interpret(info)
    n_pass = info.n_certified_pass_full_residual
    if info.skipped_solve
        return "SKIPPED"
    elseif n_pass > 0
        return "SAT_fullpass_$n_pass"
    elseif info.n_certified > 0 && n_pass == 0
        return "EMPTY_SQUARE_$(info.n_certified)"
    else
        return "OTHER"
    end
end

function all_six_cliques(H; ortho_tol=1e-8)
    pool, = generate_candidate_pool_fresh(H; verbose=false)
    pool = deduplicate_pool(pool)
    g = _orthogonality_graph(pool; ortho_tol=ortho_tol)
    cliques = [collect(c[1:6]) for c in maximal_cliques(g) if length(c) >= 6]
    # Deduplicate as sets
    uniq = Dict{Vector{Int}, Vector{Int}}()
    for c in cliques
        uniq[sort(c)] = c
    end
    return pool, collect(values(uniq))
end

function main()
    mkpath(RESULTS_DIR)
    lines = String[]
    push!(lines, "=== n_wit=1 follow-up (CHM λ+π, all-cliques at 0.4, certify λ=0) ===")
    push!(lines, "timestamp = $(Dates.now())")
    push!(lines, "")

    # (a) CHM λ vs λ+π
    push!(lines, "--- (a) CHM residual H(λ) vs H(λ+π) ---")
    for λ in [0.0, 0.4, π / 3, 2π / 3, π / 2]
        H1 = build_karlsson_family(DITA_THETA, DITA_PHI, λ)
        H2 = build_karlsson_family(DITA_THETA, DITA_PHI, λ + π)
        r = chm_equivalence_residual(H1, H2)
        push!(lines, @sprintf("  λ=%.6f vs λ+π: residual=%.6e equiv=%s variant=%s",
                              λ, r.residual, r.equivalent, r.variant))
    end
    # Independent representatives among the seven-point set
    push!(lines, "Proposed CHM-class representatives: {0, 0.4, π/3, 2π/3}")
    reps = [0.0, 0.4, π / 3, 2π / 3]
    min_r = Inf
    any_eq = false
    for i in 1:length(reps), j in (i + 1):length(reps)
        r = chm_equivalence_residual(
            build_karlsson_family(DITA_THETA, DITA_PHI, reps[i]),
            build_karlsson_family(DITA_THETA, DITA_PHI, reps[j]))
        min_r = min(min_r, r.residual)
        any_eq |= r.equivalent
        push!(lines, @sprintf("  rep %.4f vs %.4f: residual=%.6e equiv=%s",
                              reps[i], reps[j], r.residual, r.equivalent))
    end
    push!(lines, any_eq ? "FAIL: representatives not independent" :
                          @sprintf("PASS: four representatives pairwise inequivalent (min_res=%.3e)", min_r))
    push!(lines, "")

    # (b) All 6-cliques at λ=0.4
    push!(lines, "--- (b) W1 certify on EVERY size-6 clique at λ=0.4 ---")
    H04 = build_karlsson_family(DITA_THETA, DITA_PHI, 0.4)
    pool, cliques = all_six_cliques(H04)
    push!(lines, @sprintf("pool_size=%d  n_distinct_6cliques=%d", length(pool), length(cliques)))
    n_empty = 0
    for (k, idx) in enumerate(cliques)
        B3 = hcat([pool[i] for i in idx]...)
        info = certify_fourth_mub_witness_at_H(H04; n_wit=1,
            clique_indices=idx, B3=B3, verbose=false)
        v = interpret(info)
        n_empty += startswith(v, "EMPTY") ? 1 : 0
        push!(lines, @sprintf("  clique %d %s: mv=%d tracked=%d cert=%d fullpass=%d verdict=%s",
                              k, idx, info.mixed_volume, info.n_tracked, info.n_certified,
                              info.n_certified_pass_full_residual, v))
    end
    push!(lines, n_empty == length(cliques) ?
        "PASS: W1 empty for all $(length(cliques)) distinct 6-cliques at λ=0.4" :
        "FAIL: some clique had non-empty W1")
    push!(lines, "")

    # (c) Certify λ=0
    push!(lines, "--- (c) Certify n_wit=1 at λ=0 (fresh B3) ---")
    H0 = build_karlsson_family(DITA_THETA, DITA_PHI, 0.0)
    info0 = certify_fourth_mub_witness_at_H(H0; n_wit=1, verbose=false)
    v0 = interpret(info0)
    push!(lines, @sprintf("  clique=%s mv=%d tracked=%d cert=%d fullpass=%d verdict=%s found4=%s",
                          info0.clique_indices, info0.mixed_volume, info0.n_tracked,
                          info0.n_certified, info0.n_certified_pass_full_residual, v0,
                          info0.found_fourth_combinatorial))
    push!(lines, startswith(v0, "EMPTY") ? "PASS: λ=0 empty" : "FAIL: λ=0 not empty")

    text = join(lines, "\n") * "\n"
    write(OUT, text)
    println(text)
    println("Wrote ", OUT)
end

main()
