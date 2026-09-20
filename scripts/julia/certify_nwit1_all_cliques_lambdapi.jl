# Native all-clique W1 certification at the three λ+π partners:
#   π, 4π/3, 5π/3
#
# If Windows Application Control blocks precompiled .dll files in .julia-depot,
# add: julia --compiled-modules=no --project=. ...
#
# Usage: julia --project=. scripts/julia/certify_nwit1_all_cliques_lambdapi.jl
# Output: results/certify_nwit1_all_cliques_lambdapi.txt

include(joinpath(@__DIR__, "_paths.jl"))
include(joinpath(ROOT, "src", "dita_third_mub_construction.jl"))

using Printf, Dates

const OUT = joinpath(RESULTS_DIR, "certify_nwit1_all_cliques_lambdapi.txt")
const LAMBDAS = [π, 4π / 3, 5π / 3]

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

function interpret(info)
    n_pass = info.n_certified_pass_full_residual
    info.skipped_solve && return "SKIPPED"
    n_pass > 0 && return "SAT_fullpass_$n_pass"
    info.n_certified > 0 && n_pass == 0 && return "EMPTY_SQUARE_$(info.n_certified)"
    return "OTHER"
end

function main()
    mkpath(RESULTS_DIR)
    lines = String[]
    push!(lines, "=== Native all-clique n_wit=1 at λ+π partners ===")
    push!(lines, "timestamp = $(Dates.now())")
    push!(lines, "lambdas = {π, 4π/3, 5π/3}")
    push!(lines, "")

    total_cliques = 0
    total_empty = 0
    all_pass = true

    for λ in LAMBDAS
        H = build_karlsson_family(DITA_THETA, DITA_PHI, λ)
        pool, cliques = all_six_cliques(H)
        push!(lines, @sprintf("--- λ=%.10f  pool=%d  n_6cliques=%d ---",
                              λ, length(pool), length(cliques)))
        n_empty = 0
        for (k, idx) in enumerate(cliques)
            B3 = hcat([pool[i] for i in idx]...)
            info = certify_fourth_mub_witness_at_H(H; n_wit=1,
                clique_indices=idx, B3=B3, verbose=false)
            v = interpret(info)
            empty = startswith(v, "EMPTY")
            n_empty += empty ? 1 : 0
            total_cliques += 1
            total_empty += empty ? 1 : 0
            all_pass &= empty
            push!(lines, @sprintf(
                "  clique %d %s: mv=%d tracked=%d cert=%d fullpass=%d found4=%s verdict=%s",
                k, idx, info.mixed_volume, info.n_tracked, info.n_certified,
                info.n_certified_pass_full_residual, info.found_fourth_combinatorial, v))
        end
        push!(lines, @sprintf("  λ summary: empty=%d/%d", n_empty, length(cliques)))
        push!(lines, "")
    end

    push!(lines, "=== GLOBAL ===")
    push!(lines, @sprintf("total_cliques=%d total_empty=%d all_pass=%s",
                          total_cliques, total_empty, all_pass))
    push!(lines, all_pass ?
        "PASS: W1 empty for every native size-6 clique at π, 4π/3, 5π/3." :
        "FAIL: at least one native (λ, clique) was not empty.")

    text = join(lines, "\n") * "\n"
    write(OUT, text)
    println(text)
    println("Wrote ", OUT)
end

main()
