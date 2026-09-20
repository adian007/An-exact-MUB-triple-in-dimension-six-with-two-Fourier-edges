# ============================================================================
# Pool.jl — Parts V, VI, VII of the methodology hardening plan.
#
#   * `PoolAudit`: every stage of the MU-pool pipeline reports its own
#     count, and every count difference between stages is explained.
#     Never report one number for six different things:
#       paths attempted → numerical solutions → certified solutions →
#       conjugate-locus (physical) solutions → MU-valid vectors →
#       deduplicated vectors.
#   * Deduplication audit: the pool count is reported under a sweep of
#     clustering tolerances (1e-6 … 1e-12 and arbitrary precision) so any
#     tolerance-dependent classification is visible instead of hidden.
#   * Global-phase equivalence is handled by a mathematically clean
#     projective (Fubini–Study) distance, not by entrywise rounding.
# ============================================================================

# ---------------------------------------------------------------------------
# Part VI — structured pool audit
# ---------------------------------------------------------------------------

"""
    PoolAudit

Complete bookkeeping of one per-H MU-pool computation. Field semantics:

  * `mixed_volume`        — BKK bound for the square pool system (paths the
                            polyhedral homotopy attempts).
  * `paths_attempted`     — start paths launched (= mixed volume here).
  * `numerical_solutions` — finite solutions returned by the tracker.
  * `failed_paths`        = paths_attempted − numerical_solutions. At
                            degenerate Karlsson points part of the difference
                            is genuine (roots escaping to infinity relative to
                            the generic fiber); it is NEVER assumed harmless —
                            an independent re-solve at a different seed must
                            agree before a negative result is trusted.
  * `certified_solutions` — solutions certified by interval arithmetic.
  * `uncertified_solutions` — returned but not certified (must be 0 for any
                            pool used in a negative claim).
  * `conjugate_locus_solutions` — solutions with w = conj(z) (gauge z0=1):
                            exactly the solutions that correspond to actual
                            complex unit vectors.
  * `conjugate_rejections` — solutions dropped for being off the conjugate
                            locus (they satisfy the algebraic equations but
                            are not vectors; standard, documented behavior).
  * `valid_mu_vectors`    — conjugate-locus solutions passing the independent
                            (non-solver) MU verification |v_i|²=1/6 and
                            |<h_j/√6, v>|²=1/6 for all i, j.
  * `mu_rejections`       — conjugate-locus solutions failing verification.
  * `deduplicated_vectors`— distinct physical vectors up to global phase.
  * `duplicate_count`     — merges performed by deduplication.
"""
struct PoolAudit
    mixed_volume::Int
    paths_attempted::Int
    numerical_solutions::Int
    certified_solutions::Int
    conjugate_locus_solutions::Int
    valid_mu_vectors::Int
    deduplicated_vectors::Int
    failed_paths::Int
    uncertified_solutions::Int
    conjugate_rejections::Int
    mu_rejections::Int
    duplicate_count::Int
    explanations::Vector{String}
end

function audit_explanations(a::PoolAudit)
    lines = String[]
    a.failed_paths > 0 && push!(lines,
        "failed_paths=$(a.failed_paths): paths diverged or failed to converge; at " *
        "degenerate points this can be genuine (roots at infinity). Completeness of the " *
        "physical pool requires an independent re-solve agreement, not this count alone.")
    a.uncertified_solutions > 0 && push!(lines,
        "uncertified_solutions=$(a.uncertified_solutions): returned solutions without " *
        "interval certification — pools with this > 0 may not support negative claims.")
    a.conjugate_rejections > 0 && push!(lines,
        "conjugate_rejections=$(a.conjugate_rejections): algebraic solutions with " *
        "w ≠ conj(z) (no physical vector); expected and documented.")
    a.mu_rejections > 0 && push!(lines,
        "mu_rejections=$(a.mu_rejections): conjugate-locus solutions failing the " *
        "independent MU verification — investigate if nonzero.")
    a.duplicate_count > 0 && push!(lines,
        "duplicate_count=$(a.duplicate_count): merges under projective clustering " *
        "(same vector up to global phase); equivalence relation documented in " *
        "`cluster_pool_hp`.")
    isempty(lines) && push!(lines, "no count gaps beyond the expected conjugate-locus filter")
    return lines
end

function Base.show(io::IO, a::PoolAudit)
    print(io, "PoolAudit(mv=$(a.mixed_volume), attempted=$(a.paths_attempted), " *
          "solutions=$(a.numerical_solutions), certified=$(a.certified_solutions), " *
          "conjugate_locus=$(a.conjugate_locus_solutions), valid_mu=$(a.valid_mu_vectors), " *
          "dedup=$(a.deduplicated_vectors))")
end

# ---------------------------------------------------------------------------
# Parts V + VII — solve with full bookkeeping; projective dedup machinery
# ---------------------------------------------------------------------------

"""
    projective_distance(u, v)

Fubini–Study distance between vectors, invariant under global phase AND
scaling: d = arccos(|<u,v>| / (||u|| ||v||)). The natural equivalence
metric for pool vectors (a vector and any global-phase multiple are the
same quantum state).
"""
function projective_distance(u::AbstractVector{ComplexF64}, v::AbstractVector{ComplexF64})
    c = abs(dot(u, v)) / (norm(u) * norm(v))
    c > 1.0 && (c = 1.0)
    return acos(c)
end

"""
    cluster_pool_hp(pool; bits=128, tol=1e-8)

Single-linkage clustering of pool vectors under the projective metric at
`bits`-bit precision. Returns cluster representatives (highest-precision
canonical phase per cluster: first nonzero component rotated to 1·e^{i0},
using the largest-magnitude component as anchor for numerical stability)
and the merge table. Vectors with projective distance < tol are identified
— this is the explicit equivalence relation required for deduplication.
"""
function cluster_pool_hp(pool::Vector{Vector{ComplexF64}}; bits::Int = 128, tol::Real = 1e-8)
    n = length(pool)
    n == 0 && return Vector{Vector{ComplexF64}}(), Int[]
    parent = collect(1:n)
    find(x) = (x == parent[x] ? x : (parent[x] = find(parent[x])))
    setprecision(bits) do
        pb = [[Complex{BigFloat}(x) for x in v] for v in pool]
        norms = [norm(v) for v in pb]
        for i in 1:n, j in (i + 1):n
            c = abs(dot(pb[i], pb[j])) / (norms[i] * norms[j])
            c > 1 && (c = BigFloat(1))
            d = acos(c)
            if d < BigFloat(tol)
                ri, rj = find(i), find(j)
                ri != rj && (parent[minmax(ri, rj)[2]] = minmax(ri, rj)[1])
            end
        end
    end
    clusters = Dict{Int,Vector{Int}}()
    for i in 1:n
        push!(get!(clusters, find(i), Int[]), i)
    end
    reps = Vector{Vector{ComplexF64}}()
    assignment = Vector{Int}(undef, n)
    for (k, (_, members)) in enumerate(sort(collect(clusters); by = x -> x[1]))
        for m in members
            assignment[m] = k
        end
        push!(reps, _phase_canonical(pool[members[1]]))
    end
    return reps, assignment
end

"""Phase-canonical form: rotate so the largest-magnitude component is real
positive. Well-defined whenever the vector is nonzero (all pool vectors
have all components of modulus 1/√6, so the anchor is unambiguous)."""
function _phase_canonical(v::Vector{ComplexF64})
    w = copy(v)
    _, idx = findmax(abs.(w))
    w ./= w[idx] / abs(w[idx])
    return w
end

"""
    dedup_tolerance_diagnostics(pool; tols=[1e-6,1e-8,1e-10,1e-12], bits=128)

Part VII deliverable: how the distinct-vector count depends on the
clustering tolerance. A pool whose count MOVES between 1e-8 and 1e-10 is
tolerance-fragile and must be flagged (same lesson as the θ=π/2 artifact).
"""
function dedup_tolerance_diagnostics(pool::Vector{Vector{ComplexF64}};
                                     tols = [1e-6, 1e-8, 1e-10, 1e-12], bits::Int = 128)
    rows = Vector{NamedTuple}()
    for tol in tols
        reps, _ = cluster_pool_hp(pool; bits = bits, tol = tol)
        push!(rows, (tol = tol, n_clusters = length(reps), n_input = length(pool)))
    end
    return rows
end

"""
    solve_pool_audited(H; tol=1e-8, hp_verify_bits=0, reseed_crosscheck=false)

Run the complete per-H pool computation (Part V) with full Part VI
bookkeeping. The MU verification of every conjugate-locus solution is
independent of the solver (direct evaluation of |v_i|² and |<h_j,v>|²).
`hp_verify_bits > 0` additionally re-verifies all accepted vectors at
BigFloat precision. `reseed_crosscheck=true` runs a second independent
polyhedral solve (different randomization) and reports whether the
solution-count and physical-pool counts agree — the completeness
cross-check required before negative claims at degenerate points.
"""
function solve_pool_audited(H::AbstractMatrix{ComplexF64};
                            tol::Real = 1e-8,
                            hp_verify_bits::Int = 0,
                            reseed_crosscheck::Bool = false)
    system = build_numeric_pool_system(H)
    mv = mixed_volume(system)
    res = solve(system)
    raw = solutions(res)

    cert = certify(system, res; threading = false)
    n_cert = ncertified(cert)
    n_raw = length(raw)

    conjugate_ok = Bool[]
    for sol in raw
        ok = all(abs(sol[i + 5] - conj(sol[i])) < tol for i in 1:5)
        push!(conjugate_ok, ok)
    end

    valid_vectors = Vector{Vector{ComplexF64}}()
    n_conj = count(conjugate_ok)
    n_mu_rej = 0
    worst_defect = 0.0
    for (k, sol) in enumerate(raw)
        conjugate_ok[k] || continue
        v = _solution_to_v(sol)
        ok, defect = verify_candidate(v, H; tol = tol)
        worst_defect = max(worst_defect, defect)
        if ok
            push!(valid_vectors, v)
        else
            n_mu_rej += 1
        end
    end

    reps, assignment = cluster_pool_hp(valid_vectors; bits = 128, tol = 1e-8)
    n_dupes = length(valid_vectors) - length(reps)

    _audit_fields = (Int(mv), Int(mv), n_raw, Int(n_cert), n_conj, length(valid_vectors),
                     length(reps), Int(mv) - n_raw, n_raw - Int(n_cert),
                     n_raw - n_conj, n_mu_rej, n_dupes)
    audit = PoolAudit(_audit_fields..., String[])
    audit = PoolAudit(audit.mixed_volume, audit.paths_attempted, audit.numerical_solutions,
                      audit.certified_solutions, audit.conjugate_locus_solutions,
                      audit.valid_mu_vectors, audit.deduplicated_vectors,
                      audit.failed_paths, audit.uncertified_solutions,
                      audit.conjugate_rejections, audit.mu_rejections,
                      audit.duplicate_count, audit_explanations(audit))

    hp_ok = true
    if hp_verify_bits > 0
        for v in reps
            hp = verify_candidate_hp(v, H; bits = hp_verify_bits)
            hp_ok &= hp.ok
        end
    end

    crosscheck = missing
    if reseed_crosscheck
        # Force an independent randomization: re-seed the global RNG so the
        # polyhedral start system differs from run 1.
        Random.seed!(Random.default_rng(), rand(1:1_000_000_000))
        system2 = build_numeric_pool_system(H)
        res2 = solve(system2)   # HC draws fresh randomization per call
        raw2 = solutions(res2)
        pool2 = Vector{Vector{ComplexF64}}()
        for sol in raw2
            all(abs(sol[i + 5] - conj(sol[i])) < tol for i in 1:5) || continue
            v = _solution_to_v(sol)
            ok, _ = verify_candidate(v, H; tol = tol)
            ok && push!(pool2, v)
        end
        reps2, _ = cluster_pool_hp(pool2; bits = 128, tol = 1e-8)
        crosscheck = (n_solutions_run1 = n_raw, n_solutions_run2 = length(raw2),
                      n_physical_run1 = length(reps), n_physical_run2 = length(reps2),
                      agree = (length(raw2) == n_raw && length(reps2) == length(reps)))
    end

    return (pool = reps, audit = audit, worst_defect = worst_defect,
            hp_verified = hp_verify_bits > 0 ? hp_ok : missing,
            crosscheck = crosscheck)
end

"""BigFloat MU verification of a full vector (independent of the solver)."""
function verify_candidate_hp(v::Vector{ComplexF64}, H::AbstractMatrix{ComplexF64};
                             bits::Int = 256, tol::Real = 1e-12)
    out = nothing
    setprecision(bits) do
        vb = [Complex{BigFloat}(x) for x in v]
        defect = maximum(abs(abs2(vb[i]) - BigFloat(1) / 6) for i in 1:6)
        for k in 1:6
            c = [Complex{BigFloat}(H[i, k]) for i in 1:6]
            c ./= sqrt(BigFloat(6))
            defect = max(defect, abs(abs2(dot(c, vb)) - BigFloat(1) / 6))
        end
        out = (defect = Float64(defect), ok = defect < BigFloat(tol))
    end
    return out
end
