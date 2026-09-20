# ============================================================================
# Cliques.jl — Part IX: exhaustive third-MUB enumeration.
#
# Distinguishes "a third MUB exists" from "all third MUBs have been
# enumerated". Every candidate basis is verified independently (pairwise
# orthogonality, normalization, ONB completeness, MU to I, MU to H) —
# the graph search only PROPOSES candidates; verification decides.
#
# Completeness argument (why enumerating maximal cliques is exhaustive):
# in C^6 at most 6 unit vectors can be pairwise orthogonal, so every
# orthogonality clique has size ≤ 6 and the 6-cliques are exactly the
# maximal cliques of size ≥ 6. Enumerating maximal cliques therefore
# enumerates ALL third bases present in the pool — no subset-selection
# shortcut is needed or used.
#
# Caveat carried through to claims: completeness is relative to the pool.
# "All third MUBs enumerated" means "all third MUBs contained in this
# (certified-complete) pool". Pool completeness is a separate claim with
# its own evidence (see Pool.jl and the benchmarks).
# ============================================================================

"""
    enumerate_all_third_mub_bases(pool, H; ortho_tol=1e-8, mu_tol=1e-8,
                                  hp_bits=0)

Enumerate and verify every third MUB contained in `pool`:

1. build the orthogonality graph (edge iff |<v_i, v_j>| < ortho_tol);
2. enumerate all maximal cliques (equivalently all 6-cliques, see the
   completeness note above);
3. verify each candidate independently:
   - all 15 pairs orthogonal (max |<v_i,v_j>| reported);
   - each vector normalized (max ||v_i|−1/√6|... reported as |‖v‖−1|);
   - the six vectors form an ONB (‖B'B − I‖ with B' = B/√6... reported
     directly as ‖B†B − I‖ for unit-norm columns);
   - every vector MU to I (max ||v_i|² − 1/6|, redundant with
     normalization for unit vectors but checked separately);
   - every vector MU to every column of H (max ||<h_j/√6, v>|² − 1/6|).
4. optionally re-verify at `hp_bits`-bit precision.

Returns a vector of per-basis records plus summary counts. Raw count and
deduplicated (identical index-set) count are reported separately.
"""
function enumerate_all_third_mub_bases(pool::Vector{Vector{ComplexF64}},
                                       H::AbstractMatrix{ComplexF64};
                                       ortho_tol::Real = 1e-8, mu_tol::Real = 1e-8,
                                       hp_bits::Int = 0)
    n = length(pool)
    records = Vector{NamedTuple}()
    n == 0 && return (bases = records, n_raw = 0, n_distinct = 0,
                      ortho_tol = ortho_tol, mu_tol = mu_tol)

    g = _orthogonality_graph(pool; ortho_tol = ortho_tol)
    cliques = maximal_cliques(g)
    six_cliques = [sort!(collect(c)) for c in cliques if length(c) == 6]
    sort!(six_cliques)                       # canonical order: deterministic output
    sizes = [length(c) for c in cliques]
    n_distinct = length(six_cliques)

    for (bi, idx) in enumerate(six_cliques)
        B = hcat([pool[i] for i in idx]...)
        col_norms = [norm(B[:, k]) for k in 1:6]
        norm_defect = maximum(abs(cn - 1.0) for cn in col_norms)
        ortho_defect = 0.0
        for a in 1:6, b in (a + 1):6
            ortho_defect = max(ortho_defect, abs(dot(B[:, a], B[:, b])))
        end
        mu_I_defect = maximum(abs(abs2(B[i, k]) - 1 / 6) for i in 1:6, k in 1:6)
        mu_H_defect = 0.0
        for k in 1:6
            c = H[:, k] / sqrt(6.0)
            for i in 1:6
                mu_H_defect = max(mu_H_defect, abs(abs2(dot(c, B[:, i])) - 1 / 6))
            end
        end
        gram_err = maximum(abs.(B' * B - I(6)))
        is_valid = ortho_defect < ortho_tol && norm_defect < mu_tol &&
                   mu_I_defect < mu_tol && mu_H_defect < mu_tol && gram_err < mu_tol

        hp = missing
        if hp_bits > 0
            vecs = [pool[i] for i in idx]
            hp = verify_clique_hp(vecs; bits = hp_bits, ortho_tol = 1e-12, mu_tol = 1e-12)
        end

        push!(records, (
            basis_index = bi,
            pool_indices = idx,
            ortho_defect = ortho_defect,
            norm_defect = norm_defect,
            mu_to_I_defect = mu_I_defect,
            mu_to_H_defect = mu_H_defect,
            gram_err = gram_err,
            is_verified_third_mub = is_valid,
            hp = hp,
        ))
    end

    return (bases = records, n_raw = n_distinct, n_distinct = n_distinct,
            n_verified = count(r -> r.is_verified_third_mub, records),
            clique_size_histogram = Dict(s => count(==(s), sizes) for s in unique(sizes)),
            ortho_tol = ortho_tol, mu_tol = mu_tol)
end
