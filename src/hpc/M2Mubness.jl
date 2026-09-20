# ============================================================================
# M2Mubness.jl — MODULE 2: objective, gradients, family parametrizations.
#
# MATHEMATICAL CONVENTION (corrects the task spec's target value):
#   W = Σ_{k<m} Σ_{i,j} |<b_{k,i} | b_{m,j}>|^4
# For any two orthonormal bases in C^6, Σ_ij |<u_i|v_j>|^4 >= 1 (Cauchy–
# Schwarz on the 36 entries of G = B_k' B_m with Σ|G_ij|^2 = 6), with
# equality iff the bases are mutually unbiased. A set of 4 MUBs therefore
# has W = C(4,2) x 1 = 6 EXACTLY, and W is to be MINIMIZED, not maximized:
#   W ∈ [6, 36];  W = 6  ⟺  {B_0..B_3} are 4 MUBs.
# The spec's "maximize to 48 / 48.7982" is inconsistent with its own
# formula; we implement the formula. (Raynal–Lü–Englert's near-MUB analysis
# for d=6 studies exactly this quantity's per-pair minimum being >1 for
# some near-miss triples — our trap detector reports the per-pair profile.)
#
# All kernels are allocation-light 6x6 complex operations; the parallelism
# lives in MODULE 3 (over seeds). GFLOP counting is exposed for MODULE 5.
# ============================================================================

const FLOPS_PER_MATMUL_6 = 2 * 6^3          # complex mul! ≈ 2·216 real-ish ops scale

"""Σ_ij |(B_k' B_m)_ij|^4 — the per-pair MUBness distance (≥1, =1 iff MU)."""
@inline function pair_mubness(Bk::AbstractMatrix{ComplexF64}, Bm::AbstractMatrix{ComplexF64},
                              G::AbstractMatrix{ComplexF64} = zeros(ComplexF64, 6, 6))
    mul!(G, adjoint(Bk), Bm)
    s = 0.0
    @inbounds for ij in eachindex(G)
        a = abs2(G[ij])
        s += a * a
    end
    return s
end

"""
Total W over a vector of bases, plus per-pair values. Allocates nothing
beyond the scratch G when provided.
"""
function total_W(Bs::Vector{Matrix{ComplexF64}};
                 pairs = [(k, m) for k in 1:length(Bs) for m in (k + 1):length(Bs)],
                 G::AbstractMatrix{ComplexF64} = zeros(ComplexF64, 6, 6))
    pp = Vector{Float64}(undef, length(pairs))
    for (c, (k, m)) in enumerate(pairs)
        pp[c] = pair_mubness(Bs[k], Bs[m], G)
    end
    return (W = sum(pp), pairs = pairs, pair_values = pp)
end

# ---------------------------------------------------------------------------
# Analytic Euclidean gradients (validated against finite differences in
# test/runtests.jl — Enzyme/Zygote are unnecessary for 6x6 complex chains
# and slower than the closed form).
#
#   W_km = Σ_ij |g_ij|^4,  g = B_k' B_m.
#   dW = 4 Re tr( M' dg ),  M = |g|^2 ∘ g  (elementwise).
#   dg = dB_k' B_m + B_k' dB_m  ⟹
#   ∇_{B_m} W_km = 4 B_k M,      ∇_{B_k} W_km = 4 M B_m'.
# ---------------------------------------------------------------------------

@inline function _grad_pair!(Gk, Gm, Bk, Bm, G)
    mul!(G, adjoint(Bk), Bm)
    M = similar(G)
    @inbounds for ij in eachindex(G)
        M[ij] = abs2(G[ij]) * G[ij]
    end
    Gm .+= 4 .* (Bk * M)
    Gk .+= 4 .* (M * adjoint(Bm))
    return nothing
end

"""Accumulate Euclidean ∇W w.r.t. every basis into `grads` (preallocated)."""
function grad_total_W!(grads::Vector{Matrix{ComplexF64}},
                       Bs::Vector{Matrix{ComplexF64}},
                       G::AbstractMatrix{ComplexF64} = zeros(ComplexF64, 6, 6))
    n = length(Bs)
    for g in grads
        fill!(g, 0.0 + 0.0im)
    end
    for k in 1:n, m in (k + 1):n
        _grad_pair!(grads[k], grads[m], Bs[k], Bs[m], G)
    end
    return grads
end

# ---------------------------------------------------------------------------
# Riemannian machinery on U(6)^3: tangent projection + Gram–Schmidt
# retraction (modified Gram–Schmidt with a fixed phase convention, as in
# the spec; Q is unitary and Q → closest unitary in practice for small
# steps).
# ---------------------------------------------------------------------------

"""Project Euclidean gradient to the unitary-group tangent at B:
T = G − B·sym(B'G), sym(X) = (X + X')/2."""
function tangent_project!(T, G, B)
    # T = G - B*(B'G + (B'G)')/2
    BG = adjoint(B) * G
    S = (BG + adjoint(BG)) ./ 2
    mul!(T, B, S)
    @inbounds for ij in eachindex(T)
        T[ij] = G[ij] - T[ij]
    end
    return T
end

"""Modified Gram–Schmidt retraction onto U(6): returns unitary Q with the
first nonzero entry of each column made real positive (gauge)."""
function unitary_retract!(B)
    n = size(B, 1)
    for j in 1:n
        for i in 1:j - 1
            r = dot(view(B, :, i), view(B, :, j))
            if r !== 0
                @inbounds for k in 1:n
                    B[k, j] -= r * B[k, i]
                end
            end
        end
        nrm = norm(view(B, :, j))
        nrm > 0 || return false
        p = B[argmax(abs.(view(B, :, j))), j]
        ph = conj(p) / abs(p)
        @inbounds for k in 1:n
            B[k, j] = B[k, j] * ph / nrm
        end
    end
    return true
end

# ---------------------------------------------------------------------------
# Family parametrizations for the CHM-constrained search. Every builder
# returns a 6x6 matrix H with |H_ij| = 1, dephased gauge H[1,:] = H[:,1] = 1;
# bases are B = H/√6. Builders use ONLY repo-verified transcriptions:
#   * Karlsson K6(θ,φ,λ)   — the audited, seam-resolved MubSearch builder
#   * Diţă D(x)            — brierley_weigert_notes.jl (Bengtsson eq. (11))
# (The task spec's F6(a,b) is NOT implemented: its phase pattern could not
# be transcribed reliably from memory; add from Tadej–Życzkowski Table 1
# when the reference is at hand. Do not guess Hadamard patterns.)
# ---------------------------------------------------------------------------

"""One-parameter Diţă family D(x), x ∈ [0,1) (z = e^{2πix}); verified
transcription from src/brierley_weigert_notes.jl."""
function build_dita_x(x::Float64)
    return dita_D(x)
end

"""Karlsson K6(θ,φ,λ) — re-export of the audited, seam-resolved builder
from MubSearch (validated)."""
function build_karlsson3(theta::Float64, phi::Float64, lam::Float64)
    return build_karlsson_family(theta, phi, lam)
end

"""Assemble the 4-basis stack {I, H1, H2, H3}/√6 from a flat parameter
vector:
  mode :triple_karlsson — p = [θ1 φ1 λ1 θ2 φ2 λ2 θ3 φ3 λ3]  (9 params)
  mode :karlsson_dita   — p = [θ1 φ1 λ1 x2 θ3 φ3 λ3]        (7 params)
"""
function family_bases(p::Vector{Float64}; mode::Symbol = :triple_karlsson)
    s6 = sqrt(6.0)
    Hs = Matrix{ComplexF64}[]
    if mode === :triple_karlsson
        length(p) == 9 || error("triple_karlsson needs 9 parameters")
        for c in 1:3
            push!(Hs, build_karlsson3(p[3c - 2], p[3c - 1], p[3c]) ./ s6)
        end
    elseif mode === :karlsson_dita
        length(p) == 7 || error("karlsson_dita needs 7 parameters")
        push!(Hs, build_karlsson3(p[1], p[2], p[3]) ./ s6)
        push!(Hs, build_dita_x(p[4]) ./ s6)
        push!(Hs, build_karlsson3(p[5], p[6], p[7]) ./ s6)
    else
        error("unknown family mode $mode")
    end
    return [Matrix{ComplexF64}(I, 6, 6), Hs...]
end

"""Box bounds for the family multi-start (θ, φ, λ periodic; θ interior to
avoid the Fourier seam and poles)."""
function family_bounds(mode::Symbol = :triple_karlsson)
    kbox = [(1e-2, π - 1e-2), (0.0, 2π), (0.0, 2π)]
    if mode === :triple_karlsson
        return reduce(vcat, kbox for _ in 1:3)
    else
        return vcat(kbox..., [(0.0, 1.0)], kbox...)
    end
end

"""W as a function of the family parameters (multi-start objective).
Invalid points (poles, NaN) return the penalty 72 (= 2 × W_max)."""
function family_W(p::Vector{Float64}; mode::Symbol = :triple_karlsson)
    try
        Bs = family_bases(p; mode = mode)
        w = total_W(Bs).W
        return isnan(w) ? 72.0 : w
    catch e
        return 72.0
    end
end

"""Random dephased unimodular 6x6 matrix (for manifold multi-start)."""
function random_flat_matrix!(rng, H)
    for ij in eachindex(H)
        H[ij] = exp(2im * π * rand(rng))
    end
    @inbounds for j in 2:6
        H[:, j] ./= H[1, j]
    end
    @inbounds for i in 2:6
        H[i, :] ./= H[i, 1]
    end
    return H
end
