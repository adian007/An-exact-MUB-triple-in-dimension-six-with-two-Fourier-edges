# Explicit third-MUB construction on the Dita lambda-circle.
# Route A: extract certified 6-clique from MU-pool; verify parametric persistence.
#
# Brierley-Weigert (Phys. Rev. A 79, 062316, 2009) established third MUBs at
# Dita-type Karlsson parameters. This module operationalizes a constructive
# certificate: at each lambda, the pool orthogonality graph yields a 6-clique
# forming a third ONB unbiased to {I, H(lambda)}.

include(joinpath(@__DIR__, "mub_zauner_6d_liang_chen.jl"))

const DITA_THETA = acos(1 / sqrt(3))
const DITA_PHI = pi / 4

"""Return first 6-clique indices and third ONB matrix (6x6 columns) from pool."""
function extract_third_mub_basis(pool, H::AbstractMatrix{ComplexF64};
                                 ortho_tol=1e-8, mu_tol=1e-8)
    g = _orthogonality_graph(pool; ortho_tol=ortho_tol)
    cliques = [c for c in maximal_cliques(g) if length(c) >= 6]
    isempty(cliques) && return (indices=Int[], B3=Matrix{ComplexF64}(undef, 6, 0))
    c = cliques[1][1:6]
    B3 = hcat([pool[i] for i in c]...)
    # Verify MU to I and H
    for j in 1:6
        v = B3[:, j]
        for k in 1:6
            bI = zeros(ComplexF64, 6); bI[k] = 1.0
            abs(abs2(dot(bI, v)) - 1/6) >= mu_tol && error("Third basis not MU to I")
            cH = H[:, k] / sqrt(6.0)
            abs(abs2(dot(cH, v)) - 1/6) >= mu_tol && error("Third basis not MU to H")
        end
    end
    return (indices=c, B3=B3)
end

"""Build H and third ONB at Dita slice for given lambda."""
function dita_third_mub_at(lambda_; verbose=false)
    H = build_karlsson_family(DITA_THETA, DITA_PHI, lambda_)
    pool, stats = generate_candidate_pool_fresh(H; verbose=verbose)
    pool = deduplicate_pool(pool)
    ext = check_four_mub_extension(pool, H)
    ext.max_clique < 6 && error("No third MUB at lambda=$lambda_ (max_clique=$(ext.max_clique))")
    third = extract_third_mub_basis(pool, H)
    return (H=H, pool=pool, stats=stats, extension=ext, third=third)
end

"""Canonical third basis at anchor lambda=0.4 (reference for deformation)."""
function dita_reference_third(lambda_ref=0.4)
    return dita_third_mub_at(lambda_ref)
end

"""Verify third MUB exists and is HP-orthonormal at a list of lambda values."""
function verify_dita_lambda_circle(lambdas; hp_bits=256)
    results = NamedTuple[]
    for lam in lambdas
        try
            d = dita_third_mub_at(lam; verbose=false)
            vecs = [d.third.B3[:, j] for j in 1:6]
            hp = verify_clique_hp(vecs; bits=hp_bits)
            push!(results, (
                lambda=lam, ok=true, max_clique=d.extension.max_clique,
                hp_ortho=hp.ortho_ok, hp_mu=hp.mu_ok,
                n_pool=d.extension.n_pool,
            ))
        catch e
            push!(results, (lambda=lam, ok=false, error=string(e)))
        end
    end
    return results
end

export DITA_THETA, DITA_PHI, extract_third_mub_basis, dita_third_mub_at,
       dita_reference_third, verify_dita_lambda_circle
