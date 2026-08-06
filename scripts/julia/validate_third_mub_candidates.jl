# Steps 1–3: Audit third-MUB candidates (max_clique>=6), fourth-MUB per basis, clustering.
#
# Usage:
#   julia --project=. scripts/julia/validate_third_mub_candidates.jl
#   julia --project=. scripts/julia/validate_third_mub_candidates.jl --csv results/special_loci_search.csv
#   julia --project=. scripts/julia/validate_third_mub_candidates.jl --dedupe-ref

include(joinpath(@__DIR__, "_paths.jl"))

using LinearAlgebra
using Printf
using Random
using Statistics

const DEFAULT_CSV = joinpath(RESULTS_DIR, "special_loci_search.csv")
const OUT_DIR = RESULTS_DIR
const DITA_THETA = acos(1 / sqrt(3))

function load_search_csv(path)
    lines = filter(!isempty, strip.(readlines(path)))
    header = split(replace(lines[1], "\ufeff" => ""), ',')
    idx = Dict(h => i for (i, h) in enumerate(header))
    rows = NamedTuple[]
    for line in lines[2:end]
        parts = split(line, ',')
        length(parts) < length(header) && continue
        g(k) = parts[idx[k]]
        push!(rows, (
            name = g("name"),
            theta = parse(Float64, g("theta")),
            phi = parse(Float64, g("phi")),
            lambda = parse(Float64, g("lambda")),
            status = g("status"),
            max_clique = isempty(g("max_clique")) ? 0 : parse(Int, g("max_clique")),
            found_third = g("found_third") == "true",
            found_fourth = g("found_fourth") == "true",
        ))
    end
    return rows
end

function third_mub_candidates(rows; dedupe_ref = false, exclude_deep_ref = false)
    cand = [r for r in rows if r.status == "ok" && r.max_clique >= 6]
    exclude_deep_ref && (cand = [r for r in cand if !occursin("ref_ref", r.name)])
    dedupe_ref || return cand
    # One representative per 0.005 parameter cell (matches refinement scale)
    reps = NamedTuple[]
    seen = Set{Tuple{Int,Int,Int}}()
    for r in cand
        occursin("ref_ref", r.name) && continue
        key = (round(Int, r.theta / 0.005), round(Int, r.phi / 0.005), round(Int, r.lambda / 0.005))
        key in seen && continue
        push!(seen, key)
        push!(reps, r)
    end
    return reps
end

function locus_distances(theta, phi, lam)
    d_theta0 = abs(theta)
    d_dita = sqrt((theta - DITA_THETA)^2 + (phi - pi / 4)^2)
    z4 = try
        A = build_A(theta, phi)
        B = -[1 1; 1 -1] - A
        alpha_A, beta_A = A[1, 2]^2, A[1, 1]^2
        alpha_B, beta_B = B[1, 2]^2, B[1, 1]^2
        z1sq = exp(2im * lam)
        z3sq = mobius(z1sq, alpha_A, beta_A)
        z4sq = mobius(z1sq, alpha_B, beta_B)
        num = beta_B - z3sq * conj(alpha_B)
        den = alpha_B - z3sq * conj(beta_B)
        z2sq = num / den
        abs(z4sq - mobius(z2sq, alpha_A, beta_A))
    catch
        Inf
    end
    dists = Dict("Fourier_theta0" => d_theta0, "Dita" => d_dita, "Mobius_z4_dev" => z4)
    near_name = argmin(k -> dists[k], collect(keys(dists)))
    return dists, near_name, dists[near_name]
end

function exact_param_match(theta, phi, lam; atol = 1e-6)
    relations = String[]
    for (val, label) in [(0.0, "0"), (pi / 6, "pi/6"), (pi / 4, "pi/4"), (pi / 3, "pi/3"),
                         (pi / 2, "pi/2"), (DITA_THETA, "arccos(1/sqrt(3))")]
        abs(theta - val) < atol && push!(relations, "theta=$label")
        abs(phi - val) < atol && push!(relations, "phi=$label")
        abs(lam - val) < atol && push!(relations, "lambda=$label")
    end
    for k in 0:11
        v = k * pi / 6
        abs(lam - v) < atol && push!(relations, "lambda=$(k)*pi/6")
    end
    return relations
end

function hp_verify_third(theta, phi, lam; bits = 256)
    H = build_karlsson_family(theta, phi, lam)
    !is_hadamard(H; tol = 1e-10) && return (ok = false, reason = "not_hadamard")
    pool, = generate_candidate_pool_fresh(H; tol = 1e-10, verbose = false)
    pool = deduplicate_pool(pool; tol = 1e-8)
    ext = check_four_mub_extension(pool, H; ortho_tol = 1e-12, mu_tol = 1e-12)
    ext.max_clique < 6 && return (ok = false, reason = "max_clique=$(ext.max_clique)")
    g = _orthogonality_graph(pool; ortho_tol = 1e-12)
    cliques = [c for c in maximal_cliques(g) if length(c) >= 6]
    isempty(cliques) && return (ok = false, reason = "no 6-clique")
    best_c = argmax(c -> length(c), cliques)
    clique_vecs = [pool[i] for i in best_c]
    hp = verify_clique_hp(clique_vecs; bits = bits)
    mu_H = [verify_mu_to_basis_hp(v, H; bits = bits) for v in clique_vecs]
    return (
        ok = hp.ortho_ok && hp.mu_ok && all(m -> m.ok, mu_H),
        ortho_max = hp.ortho_max,
        mu_norm_max = hp.mu_norm_max,
        float_artifact = hp.ortho_max > 1e-14 || hp.mu_norm_max > 1e-14,
    )
end

function fourth_mub_per_point(theta, phi, lam)
    H = build_karlsson_family(theta, phi, lam)
    pool, = generate_candidate_pool_fresh(H; verbose = false)
    pool = deduplicate_pool(pool)
    reports = fourth_mub_per_basis_report(pool, H; ortho_tol = 1e-10, mu_tol = 1e-10)
    aggregate = check_four_mub_extension(pool, H; ortho_tol = 1e-10, mu_tol = 1e-10)
    return reports, aggregate, length(pool)
end

function cluster_analysis(candidates)
    n = length(candidates)
    n < 2 && return (rank = 0, sv = Float64[], tangent_survival = NamedTuple[], n_points = n)
    X = hcat([c.theta for c in candidates], [c.phi for c in candidates], [c.lambda for c in candidates])'
    Xc = X .- mean(X, dims = 2)
    F = svd(Xc)
    rank = count(s -> s > 1e-6 * F.S[1], F.S)
    rng = MersenneTwister(20260803)
    survival = NamedTuple[]
    for i in 1:min(n, 5)
        c = candidates[i]
        dir = randn(rng, 3); dir ./= norm(dir)
        eps = 1e-4
        th2, ph2, lm2 = c.theta + eps * dir[1], c.phi + eps * dir[2], c.lambda + eps * dir[3]
        try
            H2 = build_karlsson_family(th2, ph2, lm2)
            if is_hadamard(H2)
                pool2, = generate_candidate_pool_fresh(H2; verbose = false)
                pool2 = deduplicate_pool(pool2)
                ext2 = check_four_mub_extension(pool2, H2)
                push!(survival, (idx = i, eps = eps, third_survives = ext2.found_third, max_clique = ext2.max_clique))
            else
                push!(survival, (idx = i, eps = eps, third_survives = false, max_clique = 0))
            end
        catch
            push!(survival, (idx = i, eps = eps, third_survives = false, max_clique = 0))
        end
    end
    return (rank = rank, sv = F.S, tangent_survival = survival, n_points = n)
end

function write_csv(path, header, matrix)
    open(path, "w") do io
        println(io, join(header, ","))
        for row in matrix
            println(io, join(row, ","))
        end
    end
end

function run_validation(csv_path = DEFAULT_CSV; dedupe_ref = false, exclude_deep_ref = false)
    println("=== Step 1: Audit third-MUB candidates ===\n")
    rows = load_search_csv(csv_path)
    cand = third_mub_candidates(rows; dedupe_ref = dedupe_ref, exclude_deep_ref = exclude_deep_ref)
    println("Loaded $(length(rows)) rows; $(length(cand)) third-MUB candidates (dedupe=$dedupe_ref)\n")
    isempty(cand) && error("No max_clique>=6 candidates in $csv_path")

    audit_header = ["index", "name", "theta", "phi", "lambda", "d_theta0", "nearest_locus",
                      "nearest_dist", "hp_survival", "ortho_max_hp", "float_artifact", "exact_confirmed"]
    audit_rows = Vector{Vector{String}}()
    fourth_rows = Vector{Vector{String}}()

    for (i, r) in enumerate(cand)
        th, ph, lm = r.theta, r.phi, r.lambda
        dists, near_name, near_dist = locus_distances(th, ph, lm)
        exact = exact_param_match(th, ph, lm)
        hp = hp_verify_third(th, ph, lm)
        exact_str = isempty(exact) ? "" : join(exact, "; ")
        push!(audit_rows, vcat(string(i), r.name, string(th), string(ph), string(lm),
                               string(dists["Fourier_theta0"]), near_name, string(near_dist),
                               string(get(hp, :ok, false)), string(get(hp, :ortho_max, NaN)),
                               string(get(hp, :float_artifact, false)), exact_str))
        @printf("%2d  %-24s  θ=%.6f φ=%.6f λ=%.6f  near=%-16s d=%.2e  hp=%s  exact=%s\n",
                i, r.name, th, ph, lm, near_name, near_dist, get(hp, :ok, false),
                isempty(exact_str) ? "—" : exact_str)
    end
    out1 = joinpath(OUT_DIR, "third_mub_audit.csv")
    write_csv(out1, audit_header, audit_rows)
    println("\nWrote $out1")

    println("\n=== Step 2: Fourth-MUB per third basis ===\n")
    any_fourth = false
    for (i, r) in enumerate(cand)
        reports, agg, n_pool = fourth_mub_per_point(r.theta, r.phi, r.lambda)
        any_fourth |= agg.found_fourth
        @printf("Candidate %2d (%s): n_pool=%d aggregate_fourth=%s bases=%d\n",
                i, r.name, n_pool, agg.found_fourth, length(reports))
        for rep in reports
            @printf("  basis %d: n_filtered=%d max_clique=%d found_fourth=%s\n",
                    rep.basis_idx, rep.n_filtered, rep.max_clique_filtered, rep.found_fourth)
            push!(fourth_rows, vcat(string(i), r.name, string(r.theta), string(r.phi), string(r.lambda),
                                    string(rep.basis_idx), string(n_pool), string(rep.n_filtered),
                                    string(rep.max_clique_filtered), string(rep.found_fourth),
                                    string(agg.found_fourth)))
        end
    end
    if !isempty(fourth_rows)
        out2 = joinpath(OUT_DIR, "fourth_mub_per_basis.csv")
        write_csv(out2, ["candidate_idx", "name", "theta", "phi", "lambda", "basis_idx",
                         "n_pool", "n_filtered", "max_clique_filtered", "found_fourth",
                         "aggregate_found_fourth"], fourth_rows)
        println("\nWrote $out2")
    end

    println("\n=== Step 3: Clustering of third-MUB points ===\n")
    cl = cluster_analysis(cand)
    @printf("Parameter-space SVD rank (of %d points): %d\n", cl.n_points, cl.rank)
    @printf("Singular values: %s\n", join(string.(round.(cl.sv; digits = 6)), ", "))
    for s in cl.tangent_survival
        @printf("  perturb idx=%d eps=%.0e third_survives=%s max_clique=%d\n",
                s.idx, s.eps, s.third_survives, s.max_clique)
    end

    summary_path = joinpath(OUT_DIR, "validation_summary.txt")
    open(summary_path, "w") do io
        println(io, "n_csv_rows=$(length(rows))")
        println(io, "n_candidates=$(length(cand))")
        println(io, "cluster_rank=$(cl.rank)")
        println(io, "cluster_sv=$(join(cl.sv, ","))")
        n_hp = count(row -> row[9] == "true", audit_rows)
        println(io, "hp_survived=$n_hp")
        println(io, "any_fourth=$any_fourth")
        n_surv = count(s -> s.third_survives, cl.tangent_survival)
        println(io, "tangent_note=cluster rank $(cl.rank); tangent test $(n_surv)/$(length(cl.tangent_survival)) preserve third-MUB")
    end
    println("\nWrote $summary_path")
end

function main()
    dedupe = "--dedupe-ref" in ARGS
    exclude_deep = "--exclude-deep-ref" in ARGS
    csv_path = DEFAULT_CSV
    args = copy(ARGS)
    if "--csv" in args
        i = findfirst(==("--csv"), args)
        i !== nothing && length(args) >= i + 1 && (csv_path = args[i + 1])
    end
    run_validation(csv_path; dedupe_ref = dedupe, exclude_deep_ref = exclude_deep)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
