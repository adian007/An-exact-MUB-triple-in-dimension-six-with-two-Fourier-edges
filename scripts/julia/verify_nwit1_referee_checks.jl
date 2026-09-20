# Referee checks for n_wit=1 multi-λ certificates:
#   (1) CHM-inequivalence among the five π/3-multiple Dita H(λ)
#   (2) Clique-selection determinism for B3 used in the witness
#
# Usage:
#   julia --project=. scripts/julia/verify_nwit1_referee_checks.jl
#
# Output: results/verify_nwit1_referee_checks.txt

include(joinpath(@__DIR__, "_paths.jl"))
include(joinpath(ROOT, "src", "dita_third_mub_construction.jl"))
include(joinpath(ROOT, "scripts", "julia", "locus_classification.jl"))  # chm_equivalence_residual

using Printf
using Dates
using Random

const OUT = joinpath(RESULTS_DIR, "verify_nwit1_referee_checks.txt")

# Five π/3-multiples from the multi-λ batch (excluding anchor λ=0.4).
const PI3_LAMBDAS = [π / 3, 2π / 3, π, 4π / 3, 5π / 3]
const ALL_CERT_LAMBDAS = [0.0; 0.4; PI3_LAMBDAS]  # seven-point set once λ=0 certified

function select_b3_clique(H; ortho_tol=1e-8, seed=nothing)
    if seed !== nothing
        Random.seed!(seed)
    end
    pool, = generate_candidate_pool_fresh(H; verbose=false)
    pool = deduplicate_pool(pool)
    g = _orthogonality_graph(pool; ortho_tol=ortho_tol)
    cliques = [c for c in maximal_cliques(g) if length(c) >= 6]
    isempty(cliques) && error("no 6-clique")
    # Same rule as certify_fourth_mub_witness_at_H:
    c = argmax(cl -> length(cl), cliques)
    idx = sort(collect(c[1:6]))  # sort for stable comparison of the set
    n6 = length(cliques)
    return (clique_sorted=idx, n_cliques_ge6=n6, pool_size=length(pool))
end

function main()
    mkpath(RESULTS_DIR)
    lines = String[]
    push!(lines, "=== n_wit=1 referee checks (CHM + clique determinism) ===")
    push!(lines, "timestamp = $(Dates.now())")
    push!(lines, "")

    # ---- (1) CHM inequivalence among five π/3-multiples ----
    push!(lines, "--- (1) CHM equivalence among five π/3-multiple H(λ) ---")
    Hs = Dict{Float64, Matrix{ComplexF64}}()
    for λ in PI3_LAMBDAS
        Hs[λ] = build_karlsson_family(DITA_THETA, DITA_PHI, λ)
    end
    H04 = build_karlsson_family(DITA_THETA, DITA_PHI, 0.4)

    min_res = Inf
    max_res = 0.0
    any_equiv = false
    for i in 1:length(PI3_LAMBDAS)
        for j in (i + 1):length(PI3_LAMBDAS)
            la, lb = PI3_LAMBDAS[i], PI3_LAMBDAS[j]
            r = chm_equivalence_residual(Hs[la], Hs[lb])
            min_res = min(min_res, r.residual)
            max_res = max(max_res, r.residual)
            any_equiv |= r.equivalent
            push!(lines, @sprintf("  λ=%.6f vs λ=%.6f: residual=%.6e  equiv=%s  variant=%s",
                                  la, lb, r.residual, r.equivalent, r.variant))
        end
    end
    # Also vs anchor 0.4
    push!(lines, "  vs anchor λ=0.4:")
    for λ in PI3_LAMBDAS
        r = chm_equivalence_residual(H04, Hs[λ])
        min_res = min(min_res, r.residual)
        max_res = max(max_res, r.residual)
        any_equiv |= r.equivalent
        push!(lines, @sprintf("  λ=0.4 vs λ=%.6f: residual=%.6e  equiv=%s  variant=%s",
                              λ, r.residual, r.equivalent, r.variant))
    end
    push!(lines, @sprintf("CHM summary: min_residual=%.6e max_residual=%.6e any_equivalent=%s",
                          min_res, max_res, any_equiv))
    push!(lines, any_equiv ?
        "FAIL: some certified λ are CHM-equivalent — certificates are not independent." :
        "PASS: all five π/3-multiples are pairwise CHM-inequivalent (and vs λ=0.4).")
    push!(lines, "")

    # ---- (2) Clique-selection determinism ----
    push!(lines, "--- (2) Clique-selection determinism (rule: argmax length among maximal 6-cliques) ---")
    push!(lines, "Re-run select_b3_clique 3× per λ with fixed HC seed path (fresh pool each time).")
    det_ok = true
    for λ in [0.4; PI3_LAMBDAS]
        H = build_karlsson_family(DITA_THETA, DITA_PHI, λ)
        runs = [select_b3_clique(H) for _ in 1:3]
        sets = [r.clique_sorted for r in runs]
        n6s = [r.n_cliques_ge6 for r in runs]
        same = all(s == sets[1] for s in sets)
        det_ok &= same
        push!(lines, @sprintf("  λ=%.6f: n_cliques_ge6=%s  pool=%s",
                              λ, n6s, [r.pool_size for r in runs]))
        push!(lines, @sprintf("    run1=%s", sets[1]))
        push!(lines, @sprintf("    run2=%s", sets[2]))
        push!(lines, @sprintf("    run3=%s", sets[3]))
        push!(lines, same ? "    DETERMINISTIC: same clique set across 3 runs" :
                            "    NONDETERMINISTIC: clique set varies (set-equality of indices)")
    end
    push!(lines, "")
    push!(lines, det_ok ?
        "PASS: clique selection returned the same index set on every re-run tested." :
        "WARN: clique indices vary across runs — certificate is still valid for each chosen B3, but report must state the selection rule and that emptiness is relative to the selected B3 (combinatorial found_fourth already scans all third cliques).")
    push!(lines, "")
    push!(lines, "Selection rule (code): among maximal_cliques of the ortho graph with |C|≥6,")
    push!(lines, "  take argmax(|C|), then columns pool[C[1:6]].")
    push!(lines, "Note: emptiness of W1(H,B3) proves no fourth MUB extending {I,H,B3};")
    push!(lines, "  check_four_mub_extension already requires found_fourth=false over all third 6-cliques")
    push!(lines, "  in the pool (combinatorial), which the multi-λ batch recorded as false.")

    text = join(lines, "\n") * "\n"
    write(OUT, text)
    println(text)
    println("Wrote ", OUT)
end

main()
