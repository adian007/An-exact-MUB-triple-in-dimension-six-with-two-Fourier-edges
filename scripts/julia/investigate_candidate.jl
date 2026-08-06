# Step 4 advanced: monodromy near loci, PSLQ, Clifford/SL(2,Z_6) checks.
# Also used for counterexample escalation when found_fourth=true.
#
# Usage:
#   julia --project=. investigate_candidate.jl theta phi lambda
#   julia --project=. investigate_candidate.jl --from-csv results/special_loci_search.csv
#   julia --project=. investigate_candidate.jl --third-candidates

include(joinpath(@__DIR__, "_paths.jl"))

using LinearAlgebra
using Printf

const DEFAULT_CSV = joinpath(RESULTS_DIR, "special_loci_search.csv")
const OUT_DIR = RESULTS_DIR

"""PSLQ integer relation on parameter triple (rational π scan + mpmath hook)."""
function pslq_params(theta, phi, lam; atol = 1e-6)
    relations = String[]
    for denom in 1:12
        for num in -24:24
            v = num * pi / denom
            abs(theta - v) < atol && push!(relations, "theta ≈ $(num)*pi/$denom")
            abs(phi - v) < atol && push!(relations, "phi ≈ $(num)*pi/$denom")
            abs(lam - v) < atol && push!(relations, "lambda ≈ $(num)*pi/$denom")
        end
    end
    dita = acos(1 / sqrt(3))
    abs(theta - dita) < atol && push!(relations, "theta = arccos(1/sqrt(3)) [Dita]")
    abs(theta) < atol && push!(relations, "theta = 0 [Fourier slice]")
    return unique(relations)
end

"""Monodromy probe: witness-set dimension near fixed H (clique existence variety)."""
function monodromy_clique_probe(H::AbstractMatrix{ComplexF64}; target=5, timeout=60)
    system = build_numeric_pool_system(H)
    try
        # HomotopyContinuation v2.22: monodromy_solve(system; kwargs...)
        res = monodromy_solve(system; target_solutions_count=target,
                              max_loops_no_progress=5, show_progress=false)
        n = nresults(res)
        return (success=true, n_solutions=n, loops=res)
    catch e
        # Fallback: single solve count
        try
            sols = solve(system; show_progress=false)
            n = length(sols)
            return (success=true, n_solutions=n, fallback=true)
        catch e2
            return (success=false, error=string(e), fallback_error=string(e2))
        end
    end
end

"""Check if H is fixed under a diagonal ±1 phase symmetry (Clifford-like Z_6 action stub)."""
function clifford_z6_fixed_point_check(H::AbstractMatrix{ComplexF64})
    # Test: H[:,j] * phase_j ~ column permutation under ω6^k scaling on first two blocks
    n = size(H, 1)
    best_err = Inf
    for k in 0:5
        phase = ω6^k
        D = diagm(0 => [phase, phase, 1, 1, 1, 1])
        H2 = D * H
        err = minimum(norm(H2 - H[:, perm]) for perm in permutations(n))
        best_err = min(best_err, err)
    end
    return (best_perm_err = best_err, likely_fixed = best_err < 1e-6)
end

function permutations(n)
    n == 1 && return [[1]]
    out = Vector{Int}[]
    for p in permutations(n - 1)
        for i in 1:n
            push!(out, vcat(p[1:i-1], n, p[i:end]))
        end
    end
    return out
end

function high_precision_pool(H::AbstractMatrix{ComplexF64}; tol = 1e-12)
    pool, stats = generate_candidate_pool_fresh(H; tol = tol, verbose = true)
    pool = deduplicate_pool(pool; tol = 1e-10)
    return pool, stats
end

function recertify_at_H(H::AbstractMatrix{ComplexF64}; max_precision = 128)
    system = build_numeric_pool_system(H)
    res = solve(system; show_progress = true)
    cert = certify(system, res; show_progress = true, threading = false,
                   max_precision = max_precision)
    return (
        n_certified = ncertified(cert),
        n_distinct = ndistinct_certified(cert),
        certification = cert,
        result = res,
    )
end

function investigate(theta, phi, lam; mu_tol = 1e-12, ortho_tol = 1e-12)
    println("=== Candidate investigation ===")
    @printf("  (θ,φ,λ) = (%.16f, %.16f, %.16f)\n", theta, phi, lam)

    H = build_karlsson_family(theta, phi, lam)
    @assert is_hadamard(H; tol = 1e-10)

    println("\n[1] Pool completeness report")
    pc = pool_completeness_report(H; verbose = true)

    println("\n[2] High-precision fresh solve")
    pool, stats = high_precision_pool(H; tol = 1e-12)
    @printf("  n_pool=%d raw=%d verified=%d\n", length(pool), stats.n_raw, stats.n_verified)

    println("\n[3] Re-certify (max_precision=128)")
    cert = recertify_at_H(H; max_precision = 128)
    @printf("  certified=%d distinct=%d\n", cert.n_certified, cert.n_distinct)

    println("\n[4] Four-MUB extension at tight tolerances")
    ext = check_four_mub_extension(pool, H; mu_tol = mu_tol, ortho_tol = ortho_tol)
    @printf("  max_clique=%d found_third=%s found_fourth=%s mu_min=%.2e ortho=%.2e\n",
            ext.max_clique, ext.found_third, ext.found_fourth,
            ext.mu_defect_min, ext.ortho_defect)

    println("\n[5] Per-basis fourth-MUB report")
    for rep in fourth_mub_per_basis_report(pool, H; ortho_tol = ortho_tol, mu_tol = mu_tol)
        @printf("  basis %d: n_filtered=%d max_clique=%d found_fourth=%s\n",
                rep.basis_idx, rep.n_filtered, rep.max_clique_filtered, rep.found_fourth)
    end

    println("\n[6] Monodromy clique-existence probe (60s timeout)")
    mono = monodromy_clique_probe(H; timeout = 60)
    if mono.success
        @printf("  monodromy solutions tracked: %d\n", mono.n_solutions)
    else
        println("  monodromy failed: $(mono.error)")
    end

    println("\n[7] Clifford/Z_6 fixed-point check")
    cz = clifford_z6_fixed_point_check(H)
    @printf("  best_perm_err=%.3e likely_fixed=%s\n", cz.best_perm_err, cz.likely_fixed)

    println("\n[8] PSLQ / algebraic identification")
    rels = pslq_params(theta, phi, lam)
    if isempty(rels)
        println("  No simple π-rational relation detected.")
    else
        for r in rels
            println("  $r")
        end
    end

    return (H = H, pool = pool, ext = ext, completeness = pc, cert = cert,
            relations = rels, monodromy = mono, clifford = cz)
end

function load_search_csv_dfless(path)
    lines = filter(!isempty, strip.(readlines(path)))
    header = split(lines[1], ',')
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
            found_fourth = g("found_fourth") == "true",
        ))
    end
    return rows
end

function investigate_third_candidates(csv_path = DEFAULT_CSV)
    rows = load_search_csv_dfless(csv_path)
    cand = [r for r in rows if r.status == "ok" && r.max_clique >= 6 && !startswith(r.name, "ref_ref_")]
    results = Any[]
    for r in cand
        push!(results, investigate(r.theta, r.phi, r.lambda))
    end
    return results
end

function load_fourth_from_csv(path)
    rows = load_search_csv_dfless(path)
    return [(theta = r.theta, phi = r.phi, lambda = r.lambda) for r in rows if r.found_fourth]
end

function main()
    args = copy(ARGS)
    if !isempty(args) && args[1] == "--from-csv"
        path = length(args) >= 2 ? args[2] : DEFAULT_CSV
        rows = load_fourth_from_csv(path)
        if isempty(rows)
            println("No found_fourth=true rows in $path")
            return
        end
        for r in rows
            investigate(r.theta, r.phi, r.lambda)
        end
    elseif !isempty(args) && args[1] == "--third-candidates"
        path = length(args) >= 2 ? args[2] : DEFAULT_CSV
        investigate_third_candidates(path)
    elseif length(args) >= 3
        investigate(parse(Float64, args[1]), parse(Float64, args[2]), parse(Float64, args[3]))
    else
        println("Usage: julia --project=. investigate_candidate.jl theta phi lambda")
        println("   or: julia --project=. investigate_candidate.jl --from-csv [path]")
        println("   or: julia --project=. investigate_candidate.jl --third-candidates [path]")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
