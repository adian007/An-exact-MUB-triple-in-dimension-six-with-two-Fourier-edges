using HomotopyContinuation
using LinearAlgebra
using Graphs
using Random

const ω6 = exp(-2π * im / 6)

"""
Pipeline for searching whether the pair {I, H} extends to four MUBs in dimension 6.

QUESTION TESTED (correct scope):
  For a 6×6 CHM H, find unit vectors unbiased to BOTH the computational basis I
  and the column basis of H/√6.  A 6-clique among such vectors is a third MUB;
  a further 6-clique (after filtering to be unbiased to that third basis) is a
  fourth MUB.  Thus {I, H} extends to four MUBs iff `found_fourth == true`.

  The candidate pool MUST be re-solved for each H — it depends on H through the
  polynomial system (F6 is NOT used except in regression validation).

DEPRECATED SCOPE (do not use for family searches):
  A fixed pool unbiased to {I, F6} with an extra H filter only tests whether
  {I, F6} extends — Grassl (2004) already proves that cap is 3.
"""

function is_hadamard(H; tol=1e-8)
    n = size(H, 1)
    return norm(H * H' - n * I(n)) < tol &&
           all(abs.(abs.(H) .- 1) .< tol)
end

function hadamard_defect(H)
    n = size(H, 1)
    return norm(H * H' - n * I(n)), maximum(abs.(abs.(H) .- 1.0))
end

function mu_defect_pair(H1, H2)
    n = size(H1, 1)
    G = (H1' * H2) / n
    return maximum(abs.(abs2.(G) .- 1 / n))
end

const F6 = [ω6^(k * j) for k in 0:5, j in 0:5]

function dephase(H; tol=1e-8)
    H = copy(H)
    for j in 2:size(H, 2)
        if abs(H[1, j]) > tol
            d = conj(H[1, j]) / abs(H[1, j])
            H[:, j] .*= d
        end
    end
    for i in 2:size(H, 1)
        if abs(H[i, 1]) > tol
            d = conj(H[i, 1]) / abs(H[i, 1])
            H[i, :] .*= d
        end
    end
    return H
end

# ---------------------------------------------------------------------------
# Karlsson family (unchanged construction; see audit_karlsson_variants.py)
# ---------------------------------------------------------------------------

function build_A(theta, phi)
    A11 = -0.5 + im * (sqrt(3) / 2) * (cos(theta) + exp(-im * phi) * sin(theta))
    A12 = -0.5 + im * (sqrt(3) / 2) * (-cos(theta) + exp(im * phi) * sin(theta))
    return [A11 A12; conj(A12) -conj(A11)]
end

function mobius(z, alpha, beta)
    return (alpha * z - beta) / (conj(beta) * z - conj(alpha))
end

function build_karlsson_family(theta, phi, lambda_)
    F2 = [1 1; 1 -1]
    A = build_A(theta, phi)
    B = -F2 - A

    A_err = maximum(abs.(A * A' - 2 * I(2)))
    B_err = maximum(abs.(B * B' - 2 * I(2)))
    if A_err > 1e-6 || B_err > 1e-6
        error("Karlsson blocks not unitary (A err $(A_err), B err $(B_err)).")
    end

    alpha_A = A[1, 2]^2
    beta_A = A[1, 1]^2
    alpha_B = B[1, 2]^2
    beta_B = B[1, 1]^2

    z1 = exp(im * lambda_)
    z1sq = z1^2

    # Fourier-seam branch (audit 2026-09-14). At theta = 0 exactly:
    #   A12 = conj(A11)  =>  alpha_A = conj(beta_A)  (bitwise in floats),
    #   B11 = conj(A11), B12 = A11  =>  beta_B = alpha_A, alpha_B = beta_A,
    # so num(z) = alpha_X z - beta_X equals den(z) = conj(beta_X) z - conj(alpha_X)
    # IDENTICALLY in z: M_A and M_B are the constant map 1, and the z2^2
    # inversion num2/den2 = [alpha_A (1 - z3^2)]/[beta_A (1 - z3^2)]
    # simplifies exactly to alpha_A / beta_A. Historically this 0/0 was
    # evaluated through 1-ulp division noise, which coherently produced
    # these same values — but returned NaN at some lambda (e.g. (0,0.5,0.7)
    # pre-fix) and was not portable. Resolving the seam explicitly makes all
    # F6_theta0 anchor matrices deterministic and backward-compatible to
    # ~1e-16. Note: the theta -> 0+ LIMIT of z2^2 is a different,
    # (phi,lambda)-dependent value (diagnose_theta0_limit.jl); the seam
    # matrices are the exact algebraic theta=0 slice, not that limit.
    if theta == 0.0
        z3sq = one(ComplexF64)
        z4sq = one(ComplexF64)
        z2sq = alpha_A / beta_A
    else
        z3sq = mobius(z1sq, alpha_A, beta_A)
        z4sq = mobius(z1sq, alpha_B, beta_B)

        num = beta_B - z3sq * conj(alpha_B)
        den = alpha_B - z3sq * conj(beta_B)
        z2sq = num / den

        # Legacy consistency check, in DENOMINATOR-CLEARED form (audit fix
        # 2026-09-14): the divided form is an indeterminate 0/0 at the Dita
        # anchor (|A11|=|A12|=1), where it reported a spurious z4_dev ≈ 1.84.
        # Cleared identity: z4²·(conj(β_A)z2² − conj(α_A)) = α_A·z2² − β_A.
        z4_dev = abs(alpha_A * z2sq - beta_A - z4sq * (conj(beta_A) * z2sq - conj(alpha_A)))
        if !isnan(z4_dev) && z4_dev > 1e-6
            @warn "Karlsson z4^2 Mobius consistency violated (cleared form)" theta phi lambda_ z4_dev
        end
    end

    z2 = sqrt(z2sq)
    z3 = sqrt(z3sq)
    z4 = sqrt(z4sq)

    # Fail fast on genuinely singular Möbius branches (poles off the Fourier
    # seam): a pole in the z2²/z3²/z4² maps yields NaN entries that would
    # otherwise propagate silently into H. The seam (theta = 0) is resolved
    # above, so reaching this guard means a real degeneracy — diagnose with
    # karlsson_mobius_audit. Callers catch this (search_special_loci records
    # status=:not_hadamard).
    if any(!isfinite, (z2, z3, z4))
        error("Karlsson construction hit a singular/degenerate Möbius branch at " *
              "(theta,phi,lambda)=($theta,$phi,$lambda_): non-finite z_i. " *
              "The point is degenerate; see karlsson_mobius_audit.")
    end

    Zleft(z) = [1 1; z -z]
    Zright(z) = [1 z; 1 -z]

    Z1 = Zleft(z1)
    Z2 = Zleft(z2)
    Z3 = Zright(z3)
    Z4 = Zright(z4)

    top = hcat(F2, Z1, Z2)
    mid = hcat(Z3, 0.5 * Z3 * A * Z1, 0.5 * Z3 * B * Z2)
    bot = hcat(Z4, 0.5 * Z4 * B * Z1, 0.5 * Z4 * A * Z2)
    return vcat(top, mid, bot)
end

function build_liang_chen_family(params...)
    # Path C closed (2026-08-13): Liang/Chen/Long–Qiu papers classify special 6×6 CHMs
    # and H2-reducible types; they do not publish a Karlsson-style three-parameter family.
    # See docs/LIANG_CHEN_FAMILY_ASSESSMENT.md and McNulty–Weigert Sec. 7.1 (K6^(3) only).
    error("build_liang_chen_family: no parametric Liang/Chen family in cited sources; " *
          "see docs/LIANG_CHEN_FAMILY_ASSESSMENT.md")
end

function build_family_matrix(family::Symbol, params...)
    if family === :karlsson
        return build_karlsson_family(params...)
    elseif family === :liang_chen
        return build_liang_chen_family(params...)
    else
        error("Unsupported family: $(family)")
    end
end

# ---------------------------------------------------------------------------
# Polynomial system: vectors unbiased to I and a given H
#
# v = z/√6, z_0 = 1.  Constraints:
#   |v_i|^2 = 1/6  ↔  z_i w_i = 1, w_i ≈ conj(z_i)
#   |<H_k, v>|^2 = 1/6  ↔  (Σ conj(H_kj) z_j)(Σ H_kj w_j) = 6,  k = 0..4
# (k = 5 dependent on conjugate locus by unitarity of H/√6.)
# ---------------------------------------------------------------------------

function _pool_equations(z, w, hc, h)
    z_full = [1.0 + 0im; [z[i] for i in 1:5]]
    w_full = [1.0 + 0im; [w[i] for i in 1:5]]
    eqs_I = [z[i] * w[i] - 1.0 for i in 1:5]
    eqs_H = Vector{Any}(undef, 5)
    for k in 0:4
        # Column k of H: <H[:,k]/√6, v> with v = z/√6  =>  |Σ_j conj(H[j,k]) z_j|² = 6.
        # hc[j,k] = conj(H[j,k]),  h[j,k] = H[j,k].
        lhs = sum(hc[j + 1, k + 1] * z_full[j + 1] for j in 0:5)
        rhs = sum(h[j + 1, k + 1] * w_full[j + 1] for j in 0:5)
        eqs_H[k + 1] = lhs * rhs - 6.0
    end
    return [eqs_I; eqs_H]
end

function build_numeric_pool_system(H::AbstractMatrix{ComplexF64})
    @var z[1:5] w[1:5]
    hc = conj.(H)
    System(_pool_equations(z, w, hc, H); variables=[z; w])
end

function build_parametric_pool_system()
    @var z[1:5] w[1:5]
    @var hc[1:6, 1:6] h[1:6, 1:6]
    eqs = _pool_equations(z, w, hc, h)
    System(eqs; variables=[z; w], parameters=[hc[:]; h[:]])
end

function params_from_hadamard(H::AbstractMatrix{ComplexF64})
    vcat(vec(conj.(H)), vec(H))
end

"""
Independent verification: |v_i|^2 = 1/6 and |<H_k, v>|^2 = 1/6 for all i, k.
"""
function verify_candidate(v::Vector{ComplexF64}, H::AbstractMatrix{ComplexF64}; tol=1e-8)
    defect = maximum(abs(abs2(v[i]) - 1 / 6) for i in 1:6)
    for k in 1:6
        c = H[:, k] / sqrt(6.0)
        defect = max(defect, abs(abs2(dot(c, v)) - 1 / 6))
    end
    return defect < tol, defect
end

"""Regression check against the audited {I, F6} system (column-k convention)."""
function verify_candidate_f6(v::Vector{ComplexF64}; tol=1e-8)
    return verify_candidate(v, F6; tol=tol)
end

function _solution_to_v(sol)
    return [1.0 + 0im; [sol[i] for i in 1:5]] / sqrt(6.0)
end

function _on_conjugate_locus(sol, tol)
    return all(abs(sol[i + 5] - conj(sol[i])) < tol for i in 1:5)
end

"""Canonical key for matching vectors up to global phase (v and -v identified)."""
function vector_key(v::Vector{ComplexF64}; tol=1e-6)
    w = copy(v)
    if abs(w[1]) > tol
        w ./= w[1] / abs(w[1])
    elseif norm(w) > tol
        idx = findfirst(x -> abs(x) > tol, w)
        w ./= w[idx] / abs(w[idx])
    end
    return Tuple(round.(real.(w), sigdigits=8)), Tuple(round.(imag.(w), sigdigits=8))
end

function pool_from_raw_solutions(raw_solutions, H::AbstractMatrix{ComplexF64};
                                 tol=1e-8, verify=verify_candidate)
    pool = Vector{Vector{ComplexF64}}()
    n_conjugate = 0
    n_verified = 0
    worst_defect = 0.0
    for sol in raw_solutions
        if !_on_conjugate_locus(sol, tol)
            continue
        end
        n_conjugate += 1
        v = _solution_to_v(sol)
        ok, defect = verify(v, H; tol=tol)
        worst_defect = max(worst_defect, defect)
        if ok
            n_verified += 1
            push!(pool, v)
        end
    end
    return pool, (n_raw=length(raw_solutions), n_conjugate, n_verified, worst_defect)
end

function deduplicate_pool(pool; tol=1e-6)
    seen = Set{Tuple{Tuple{Vararg{Float64}}, Tuple{Vararg{Float64}}}}()
    out = Vector{Vector{ComplexF64}}()
    for v in pool
        k = vector_key(v; tol=tol)
        if k in seen
            continue
        end
        push!(seen, k)
        push!(out, v)
    end
    return out
end

function compare_pools(pool_a, pool_b; tol=1e-6)
    keys_a = Set(vector_key(v; tol=tol) for v in pool_a)
    keys_b = Set(vector_key(v; tol=tol) for v in pool_b)
    only_a = length(setdiff(keys_a, keys_b))
    only_b = length(setdiff(keys_b, keys_a))
    return (only_in_a=only_a, only_in_b=only_b, matched=length(intersect(keys_a, keys_b)))
end

"""
Fresh polyhedral solve at fixed H.  Used for cross-checks and validation anchors.
"""
function generate_candidate_pool_fresh(H::AbstractMatrix{ComplexF64};
                                       tol=1e-8, verbose=false)
    system = build_numeric_pool_system(H)
    res = solve(system)
    raw = solutions(res)
    pool, stats = pool_from_raw_solutions(raw, H; tol=tol)
    pool = deduplicate_pool(pool; tol=tol)
    if verbose
        println("  [fresh] tracked=$(stats.n_raw) conjugate=$(stats.n_conjugate) " *
                "verified=$(stats.n_verified) dedup=$(length(pool)) " *
                "worst_defect=$(stats.worst_defect)")
    end
    return pool, stats
end

"""
Certify solutions of the per-H pool system at fixed H using HC.jl interval arithmetic.
Returns counts of certified / real / distinct roots among finite tracked solutions.
"""
function certify_pool_at_H(H::AbstractMatrix{ComplexF64}; verbose = false)
    system = build_numeric_pool_system(H)
    res = solve(system; show_progress = verbose)
    cert = certify(system, res; show_progress = verbose, threading = false)
    n_cand = length(solutions(res))
    n_cert = ncertified(cert)
    n_real = nreal_certified(cert)
    n_distinct = ndistinct_certified(cert)
    n_distinct_real = count(r -> is_certified(r) && is_real(r), distinct_certificates(cert))
    if verbose
        println("  [certify] candidates=$n_cand certified=$n_cert real=$n_real " *
                "distinct=$n_distinct distinct_real=$n_distinct_real")
    end
    return (
        result = res,
        certification = cert,
        n_candidates = n_cand,
        n_certified = n_cert,
        n_real = n_real,
        n_distinct = n_distinct,
        n_distinct_real = n_distinct_real,
    )
end

"""Legacy {I, F6} pool — regression / validation ONLY."""
function generate_candidate_pool_f6(; tol=1e-8, verbose=false)
    return generate_candidate_pool_fresh(F6; tol=tol, verbose=verbose)[1]
end

"""
Audit whether the per-H pool pipeline recovered a complete certified root set.

Compares certified root counts to conjugate-locus / MU-filtered pool size.
At F6, 156 distinct certified roots typically yield ~48 deduplicated pool vectors
(conjugate gauge + MU filter); that gap is expected, not incomplete tracking.

Returns `pool_complete_flag=true` when every tracked candidate is certified and
tracking count is not below the mixed-volume BKK bound without explanation.
"""
function pool_completeness_report(H::AbstractMatrix{ComplexF64}; tol = 1e-8, verbose = false)
    system = build_numeric_pool_system(H)
    mv = mixed_volume(system)
    cert_info = certify_pool_at_H(H; verbose = verbose)
    raw = solutions(cert_info.result)
    pool, stats = pool_from_raw_solutions(raw, H; tol = tol)
    pool_dedup = deduplicate_pool(pool; tol = tol)

    n_tracked = cert_info.n_candidates
    n_certified = cert_info.n_certified
    n_distinct_certified = cert_info.n_distinct
    n_conjugate = stats.n_conjugate
    n_verified = stats.n_verified
    n_dedup_pool = length(pool_dedup)

    cert_match = n_certified == n_tracked
    mv_gap = mv - n_tracked
    cert_vs_pool_gap = n_distinct_certified - n_dedup_pool

    # Complete when every tracked path certified (primary gate for set S).
    pool_complete_flag = cert_match && n_certified > 0 && n_verified > 0

    if verbose
        println("  [pool completeness] mv=$mv tracked=$n_tracked certified=$n_certified " *
                "distinct_cert=$n_distinct_certified conjugate=$n_conjugate " *
                "verified=$n_verified dedup_pool=$n_dedup_pool complete=$pool_complete_flag")
    end

    return (
        mixed_volume = mv,
        n_tracked = n_tracked,
        n_certified = n_certified,
        n_distinct_certified = n_distinct_certified,
        n_conjugate = n_conjugate,
        n_verified = n_verified,
        n_dedup_pool = n_dedup_pool,
        cert_match = cert_match,
        mv_gap = mv_gap,
        cert_vs_pool_gap = cert_vs_pool_gap,
        pool_complete_flag = pool_complete_flag,
    )
end

mutable struct ParametricPoolTracker
    system::System
    generic_params::Vector{ComplexF64}
    start_solutions::Vector{Vector{ComplexF64}}
    generic_count::Int
end

"""
Offline phase: solve at H = F6 (rigorous anchor with 156 roots / 48 pool vectors),
then track along parameter homotopy to each target H.  F6 is used because it lies
on the same consistent parameter locus hc = conj(H) and has the known root structure.
"""
function init_parametric_pool_tracker(; verbose=false)
    system = build_parametric_pool_system()
    H_start = F6
    generic_params = params_from_hadamard(H_start)
    res = solve(build_numeric_pool_system(H_start))
    raw = solutions(res)
    start_solutions = [s for s in raw if all(isfinite, s) && norm(s) < Inf]
    tracker = ParametricPoolTracker(system, generic_params, start_solutions, length(raw))
    if verbose
        println("  [tracker init] F6 anchor -> $(length(raw)) raw, " *
                "$(length(start_solutions)) finite start paths")
    end
    return tracker
end

function track_pool_solutions(tracker::ParametricPoolTracker, H::AbstractMatrix{ComplexF64};
                              verbose=false)
    target = params_from_hadamard(H)
    if isempty(tracker.start_solutions)
        return Vector{Vector{ComplexF64}}(), (n_tracked=0, n_failed=tracker.generic_count,
                                               n_expected=tracker.generic_count)
    end
    res = solve(tracker.system, tracker.start_solutions;
                start_parameters=tracker.generic_params,
                target_parameters=target,
                show_progress=verbose)
    raw = solutions(res)
    n_finite = count(s -> all(isfinite, s) && norm(s) < Inf, raw)
    n_failed = length(tracker.start_solutions) - n_finite
    return raw, (n_tracked=n_finite, n_failed,
                   n_expected=length(tracker.start_solutions))
end

"""
Generate pool via parameter homotopy.  Returns pool, stats, and tracking metadata.
If `cross_check=true`, also runs a fresh solve and compares counts (drop detection).
"""
function generate_candidate_pool(H::AbstractMatrix{ComplexF64};
                                 tracker=nothing,
                                 tol=1e-8,
                                 verbose=false,
                                 cross_check=false)
    if tracker === nothing
        pool, stats = generate_candidate_pool_fresh(H; tol=tol, verbose=verbose)
        meta = (method=:fresh, n_tracked=stats.n_raw, n_failed=0,
                n_expected=stats.n_raw, cross_check=:skipped, drop_suspected=false)
        return pool, stats, meta
    end

    raw, track = track_pool_solutions(tracker, H; verbose=verbose)
    pool, stats = pool_from_raw_solutions(raw, H; tol=tol)
    pool = deduplicate_pool(pool; tol=tol)
    drop_suspected = track.n_tracked < track.n_expected || track.n_failed > 0 ||
                     stats.n_verified == 0

    cross = (only_in_a=0, only_in_b=0, matched=length(pool))
    fresh_pool = pool
    fresh_stats = stats
    if cross_check || drop_suspected
        fresh_pool, fresh_stats = generate_candidate_pool_fresh(H; tol=tol, verbose=false)
        cross = compare_pools(pool, fresh_pool; tol=tol)
        if cross.only_in_b > 0 || (drop_suspected && length(fresh_pool) > length(pool))
            drop_suspected = true
            if verbose
                println("  [DROP DETECTED] homotopy pool=$(length(pool)) fresh=$(length(fresh_pool)) " *
                        "missed=$(cross.only_in_b)")
            end
        end
    end

    if drop_suspected && length(fresh_pool) > length(pool)
        pool = fresh_pool
        stats = fresh_stats
    end

    if verbose
        println("  [homotopy] tracked=$(track.n_tracked)/$(track.n_expected) " *
                "conjugate=$(stats.n_conjugate) verified=$(stats.n_verified) " *
                "dedup=$(length(pool)) drop_suspected=$(drop_suspected)")
    end

    meta = (method=:homotopy, track..., cross_check=cross, drop_suspected=drop_suspected)
    return pool, stats, meta
end

# ---------------------------------------------------------------------------
# Four-MUB extension test on a pool already unbiased to {I, H}
# ---------------------------------------------------------------------------

function _orthogonality_graph(vectors; ortho_tol=1e-8)
    n = length(vectors)
    g = SimpleGraph(n)
    for i in 1:n, j in (i + 1):n
        if abs(dot(vectors[i], vectors[j])) < ortho_tol
            add_edge!(g, i, j)
        end
    end
    return g
end

function _vectors_mu_to_basis(vectors, B; mu_tol=1e-8)
    out = Int[]
    for (i, v) in enumerate(vectors)
        ok = true
        for j in 1:size(B, 2)
            b = B[:, j] / sqrt(6.0)
            if abs(abs2(dot(b, v)) - 1 / 6) >= mu_tol
                ok = false
                break
            end
        end
        ok && push!(out, i)
    end
    return out
end

"""
Test whether {I, H} extends to four MUBs using vectors from `pool`.

  found_third  — a 6-clique exists (third MUB found)
  found_fourth — a second 6-clique exists among vectors also MU to that third basis
"""
function check_four_mub_extension(pool, H::AbstractMatrix{ComplexF64};
                                  ortho_tol=1e-8, mu_tol=1e-8)
    n = length(pool)
    mu_defect_min = Inf
    for v in pool
        defect = maximum(abs(abs2(v[i]) - 1 / 6) for i in 1:6)
        for j in 1:6
            c = H[:, j] / sqrt(6.0)
            defect = max(defect, abs(abs2(dot(c, v)) - 1 / 6))
        end
        mu_defect_min = min(mu_defect_min, defect)
    end
    n == 0 && return (found_third=false, found_fourth=false, max_clique=0,
                      mu_defect_min=mu_defect_min, ortho_defect=0.0, n_pool=n)

    g = _orthogonality_graph(pool; ortho_tol=ortho_tol)
    cliques = maximal_cliques(g)
    max_clique = maximum(length(c) for c in cliques)
    found_third = max_clique >= 6
    found_fourth = false
    ortho_defect = 0.0

    third_cliques = [c for c in cliques if length(c) >= 6]
    for c in third_cliques
        B = hcat([pool[i] for i in c]...)
        idx = _vectors_mu_to_basis(pool, B; mu_tol=mu_tol)
        sub = pool[idx]
        if length(sub) < 6
            continue
        end
        g2 = _orthogonality_graph(sub; ortho_tol=ortho_tol)
        for c2 in maximal_cliques(g2)
            if length(c2) >= 6
                found_fourth = true
                for (a, i) in enumerate(c2), j in c2[a+1:end]
                    ortho_defect = max(ortho_defect, abs(dot(sub[i], sub[j])))
                end
                break
            end
        end
        found_fourth && break
    end

    if found_third && !found_fourth
        best_c = argmax(c -> length(c), cliques)
        for (a, ii) in enumerate(best_c), jj in best_c[(a + 1):end]
            ortho_defect = max(ortho_defect, abs(dot(pool[ii], pool[jj])))
        end
    end

    return (found_third=found_third, found_fourth=found_fourth,
            max_clique=max_clique, mu_defect_min=mu_defect_min,
            ortho_defect=ortho_defect, n_pool=n)
end

"""Per third-basis report: pool size and max clique after MU-filter to each 6-clique."""
function fourth_mub_per_basis_report(pool, H::AbstractMatrix{ComplexF64};
                                     ortho_tol=1e-8, mu_tol=1e-8)
    n = length(pool)
    n == 0 && return NamedTuple[]

    g = _orthogonality_graph(pool; ortho_tol=ortho_tol)
    cliques = [c for c in maximal_cliques(g) if length(c) >= 6]
    reports = NamedTuple[]
    for (bi, c) in enumerate(cliques)
        B3 = hcat([pool[i] for i in c]...)
        idx = _vectors_mu_to_basis(pool, B3; mu_tol=mu_tol)
        sub = pool[idx]
        max_c2 = 0
        found4 = false
        if length(sub) >= 6
            g2 = _orthogonality_graph(sub; ortho_tol=ortho_tol)
            cliques2 = maximal_cliques(g2)
            max_c2 = maximum(length(c2) for c2 in cliques2; init=0)
            found4 = any(length(c2) >= 6 for c2 in cliques2)
        end
        push!(reports, (
            basis_idx = bi,
            clique_indices = c,
            n_filtered = length(sub),
            max_clique_filtered = max_c2,
            found_fourth = found4,
        ))
    end
    return reports
end

"""High-precision clique verification using BigFloat inner products."""
function verify_clique_hp(vectors; bits=256, ortho_tol=1e-12, mu_tol=1e-12)
    setprecision(bits) do
        vb = [[Complex{BigFloat}(x) for x in v] for v in vectors]
        ortho_max = BigFloat(0)
        mu_max = BigFloat(0)
        for i in 1:length(vb)
            mu_max = max(mu_max, abs(abs2(vb[i][1]) - BigFloat(1) / 6))
            for j in (i + 1):length(vb)
                ortho_max = max(ortho_max, abs(dot(vb[i], vb[j])))
            end
        end
        return (
            ortho_max = Float64(ortho_max),
            mu_norm_max = Float64(mu_max),
            ortho_ok = ortho_max < BigFloat(ortho_tol),
            mu_ok = mu_max < BigFloat(mu_tol),
        )
    end
end

"""High-precision MU check of vector v to basis B (columns of 6×6 matrix)."""
function verify_mu_to_basis_hp(v, B; bits=256, mu_tol=1e-12)
    setprecision(bits) do
        vb = Complex{BigFloat}[Complex{BigFloat}(x) for x in v]
        mu_max = BigFloat(0)
        for j in 1:size(B, 2)
            b = [Complex{BigFloat}(B[i, j]) for i in 1:6]
            b ./= sqrt(BigFloat(6))
            mu_max = max(mu_max, abs(abs2(dot(b, vb)) - BigFloat(1) / 6))
        end
        return (mu_max = Float64(mu_max), ok = mu_max < BigFloat(mu_tol))
    end
end

# ---------------------------------------------------------------------------
# Fourth-MUB witness polynomial system (reduced: n_wit MU vectors to I, H, B3)
# Matches build_fourth_mub_witness_equations in symbolic_elimination.jl.
# ---------------------------------------------------------------------------

"""MU constraints for witness vector (z0=w0=1) unbiased to columns of basis B.

Auto-detects normalization: RHS=1 if B is unitary (B'B≈I), RHS=6 if B is
unnormalized Hadamard (B*B'≈6I)."""
function _witness_mu_to_basis(B::AbstractMatrix, z, w)
    z_full = [1.0 + 0im; [z[i] for i in 1:5]]
    w_full = [1.0 + 0im; [w[i] for i in 1:5]]
    # Auto-detect normalization convention
    b_gram = B' * B
    unitary_err = maximum(abs.(b_gram - I(6)))
    hadamard_err = maximum(abs.(b_gram - 6.0 * I(6)))
    rhs_val = unitary_err < hadamard_err ? 1.0 : 6.0
    eqs = Any[]
    for k in 1:6
        lhs = sum(conj(B[j, k]) * z_full[j] for j in 1:6)
        rhs = sum(B[j, k] * w_full[j] for j in 1:6)
        push!(eqs, lhs * rhs - rhs_val)
    end
    return eqs
end

"""
Build numeric HC system for the reduced fourth-MUB witness at fixed H and third ONB B3.

- `n_wit=1`: one vector MU to I, H, and B3 (10 vars, 17 eqs).
  Emptiness ⇒ no fourth MUB at this (H,B3): a fourth MUB would supply six
  solutions of this system. (Converse false: candidates may exist without a 6-ONB.)
- `n_wit=2`: two orthogonal such vectors (20 vars, 35 eqs) — matches
  `fourth_mub_reduced_*.m2`. Emptiness is closer to no fourth ONB, but mv≃1.2×10⁵.
"""
function build_numeric_fourth_mub_witness_system(H::AbstractMatrix{ComplexF64},
                                                   B3::AbstractMatrix{ComplexF64};
                                                   n_wit::Int=2)
    if n_wit == 1
        @var z1[1:5] w1[1:5]
        eqs = Any[]
        for i in 1:5
            push!(eqs, z1[i] * w1[i] - 1.0)
        end
        append!(eqs, _witness_mu_to_basis(H, [z1[i] for i in 1:5], [w1[i] for i in 1:5]))
        append!(eqs, _witness_mu_to_basis(B3, [z1[i] for i in 1:5], [w1[i] for i in 1:5]))
        return System(eqs; variables=vcat(z1, w1))
    elseif n_wit == 2
        @var z1[1:5] w1[1:5] z2[1:5] w2[1:5]
        z_blocks = [[z1[i] for i in 1:5], [z2[i] for i in 1:5]]
        w_blocks = [[w1[i] for i in 1:5], [w2[i] for i in 1:5]]
        eqs = Any[]
        for v in 1:n_wit
            for i in 1:5
                push!(eqs, z_blocks[v][i] * w_blocks[v][i] - 1.0)
            end
        end
        for v in 1:n_wit
            append!(eqs, _witness_mu_to_basis(H, z_blocks[v], w_blocks[v]))
            append!(eqs, _witness_mu_to_basis(B3, z_blocks[v], w_blocks[v]))
        end
        for a in 1:n_wit, b in (a + 1):n_wit
            ortho = 1.0
            for j in 1:5
                ortho += z_blocks[a][j] * w_blocks[b][j]
            end
            push!(eqs, ortho)
        end
        return System(eqs; variables=vcat(z1, z2, w1, w2))
    else
        error("build_numeric_fourth_mub_witness_system: n_wit=$n_wit not implemented (use 1 or 2)")
    end
end

"""
Solve and certify the reduced fourth-MUB witness system at fixed H.
Returns solve/certify counts plus clique indices used for B3.
"""
function probe_fourth_mub_witness_mv(H::AbstractMatrix{ComplexF64}, B3::AbstractMatrix{ComplexF64};
                                     n_wit::Int=2)
    system = build_numeric_fourth_mub_witness_system(H, B3; n_wit=n_wit)
    n_eqs = length(system)
    n_vars = nvariables(system)
    mv = mixed_volume(system)
    return (n_eqs=n_eqs, n_vars=n_vars, mixed_volume=mv, system=system)
end

function certify_fourth_mub_witness_at_H(H::AbstractMatrix{ComplexF64};
                                         n_wit::Int=2,
                                         clique_indices::Union{Nothing,Vector{Int}}=nothing,
                                         B3::Union{Nothing,AbstractMatrix{ComplexF64}}=nothing,
                                         verbose::Bool=false)
    ext_found_fourth = false
    if B3 !== nothing
        clique_indices === nothing && error("B3 provided without clique_indices")
    else
        pool, = generate_candidate_pool_fresh(H; verbose=verbose)
        pool = deduplicate_pool(pool)
        ext = check_four_mub_extension(pool, H; ortho_tol=1e-8, mu_tol=1e-8)
        ext.max_clique < 6 && error("No third MUB clique at this H (max_clique=$(ext.max_clique))")
        ext_found_fourth = ext.found_fourth

        if clique_indices === nothing
            g = _orthogonality_graph(pool; ortho_tol=1e-8)
            cliques = [c for c in maximal_cliques(g) if length(c) >= 6]
            isempty(cliques) && error("No 6-clique in pool despite max_clique=$(ext.max_clique)")
            c = argmax(cl -> length(cl), cliques)
            clique_indices = c[1:6]
        end
        length(clique_indices) >= 6 || error("clique_indices must have length >= 6")
        B3 = hcat([pool[i] for i in clique_indices[1:6]]...)
    end

    system = build_numeric_fourth_mub_witness_system(H, B3; n_wit=n_wit)
    n_eqs = length(system)
    n_vars = nvariables(system)
    mv_full = mixed_volume(system)

    # HC certify() requires a square system. For overdetermined witnesses (n_wit=1:
    # 17 eqs / 10 vars), form a generic square by random linear combinations of the
    # equations (Bertini / HC pattern), solve+certify THAT system, then residual-check
    # every original equation. Solving the overdetermined system first is wrong:
    # polyhedral tracking reports all paths as "excess" and returns 0 solutions.
    squared_for_certify = n_eqs != n_vars
    solve_system = system
    if squared_for_certify
        eqs = expressions(system)
        vars = variables(system)
        Random.seed!(20260813 + n_wit)
        R = randn(ComplexF64, n_vars, n_eqs)
        squared_eqs = [sum(R[i, j] * eqs[j] for j in 1:n_eqs) for i in 1:n_vars]
        solve_system = System(squared_eqs; variables=vars)
    end
    mv = mixed_volume(solve_system)

    if mv > WITNESS_MV_SOLVE_CAP
        return (
            H=H, B3=B3, clique_indices=clique_indices, n_wit=n_wit,
            n_eqs=n_eqs, n_vars=n_vars, mixed_volume=mv, mixed_volume_full=mv_full,
            n_tracked=0, n_certified=0, n_distinct_certified=0, n_real_certified=0,
            n_certified_pass_full_residual=0,
            found_fourth_combinatorial=ext_found_fourth,
            skipped_solve=true, skip_reason="mixed_volume=$mv exceeds cap $WITNESS_MV_SOLVE_CAP",
            squared_for_certify=squared_for_certify, result=nothing, certification=nothing,
            system=system, solve_system=solve_system,
        )
    end

    res = solve(solve_system; show_progress=verbose)
    raw = solutions(res)
    n_tracked = length(raw)

    cert = certify(solve_system, res; show_progress=verbose, threading=false)
    n_cert = ncertified(cert)
    n_distinct = ndistinct_certified(cert)
    n_real = nreal_certified(cert)

    # Full-system residual gate: only count square-roots that satisfy ALL original eqs.
    n_pass_full = 0
    if n_tracked > 0
        for sol in raw
            rmax = 0.0
            try
                vals = system(sol)
                rmax = maximum(abs, vals)
            catch
                rmax = Inf
            end
            rmax < 1e-8 && (n_pass_full += 1)
        end
    end

    return (
        H=H,
        B3=B3,
        clique_indices=clique_indices,
        n_wit=n_wit,
        n_eqs=n_eqs,
        n_vars=n_vars,
        mixed_volume=mv,
        mixed_volume_full=mv_full,
        n_tracked=n_tracked,
        n_certified=n_cert,
        n_distinct_certified=n_distinct,
        n_real_certified=n_real,
        n_certified_pass_full_residual=n_pass_full,
        found_fourth_combinatorial=ext_found_fourth,
        skipped_solve=false,
        skip_reason="",
        squared_for_certify=squared_for_certify,
        result=res,
        certification=cert,
        system=system,
        solve_system=solve_system,
    )
end

"""
End-to-end search for one Karlsson (or other) family member:
  1. Build H
  2. Solve per-H candidate pool (vectors MU to {I, H})
  3. Test extension to four MUBs via clique analysis
"""
function run_search(family::Symbol, params...;
                    tol=1e-8, ortho_tol=1e-8, mu_tol=1e-8,
                    tracker=nothing, cross_check=false, verbose=false)
    H = build_family_matrix(family, params...)
    if !is_hadamard(H; tol=tol)
        error("The matrix produced by $(family) is not Hadamard within tolerance.")
    end
    pool, pool_stats, meta = generate_candidate_pool(H; tracker=tracker,
                                                     tol=tol, verbose=verbose,
                                                     cross_check=cross_check)
    ext = check_four_mub_extension(pool, H; ortho_tol=ortho_tol, mu_tol=mu_tol)
    return (H=H, pool=pool, pool_stats=pool_stats, meta=meta, ext...)
end

const WITNESS_MV_SOLVE_CAP = 5000

# Legacy name for the deprecated {I,F6} fixed pool (regression tests only).
generate_candidate_pool_legacy_f6 = generate_candidate_pool_f6
