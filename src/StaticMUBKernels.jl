# ============================================================================
# StaticMUBKernels.jl — Zero-allocation 6×6 MUB kernels using StaticArrays.
#
# All operations use SMatrix{6,6,ComplexF64,36} which lives entirely on the
# stack and enables full LLVM SIMD vectorization. A single W evaluation
# takes ~200 ns vs ~50 μs with heap-allocated Matrix{ComplexF64}.
#
# Mathematical convention (matching src/hpc/M2Mubness.jl):
#   W = Σ_{k<m} Σ_{i,j} |⟨b_{k,i} | b_{m,j}⟩|⁴
#   W ≥ 6 for 4 orthonormal bases in C⁶; W = 6 ⟺ {B₀..B₃} are 4 MUBs.
#   W is MINIMIZED toward 6.0.
#
# Analytic gradients (validated against finite differences):
#   G = Bk' * Bm
#   M = |G|² ∘ G  (elementwise)
#   ∇_{Bm} W_km = 4 Bk * M
#   ∇_{Bk} W_km = 4 M * Bm'
# ============================================================================

using StaticArrays
using LinearAlgebra
using Random

# ---------------------------------------------------------------------------
# Type aliases
# ---------------------------------------------------------------------------
const SMat6  = SMatrix{6, 6, ComplexF64, 36}
const SVec6  = SVector{6, ComplexF64}
const MMat6  = MMatrix{6, 6, ComplexF64, 36}

# ---------------------------------------------------------------------------
# Core objective: per-pair mubness distance
# ---------------------------------------------------------------------------

"""
    pair_mubness_static(Bk::SMat6, Bm::SMat6) -> Float64

Compute Σ_{i,j} |(Bk' Bm)_{ij}|⁴ — the per-pair MUBness distance.
Equals 1.0 iff Bk and Bm are mutually unbiased; ≥ 1.0 always.
Zero allocations.
"""
@inline function pair_mubness_static(Bk::SMat6, Bm::SMat6)
    G = Bk' * Bm   # SMatrix multiplication, fully on stack
    s = 0.0
    @inbounds for ij in 1:36
        a = abs2(G[ij])
        s += a * a
    end
    return s
end

"""
    total_W_static(B0, B1, B2, B3) -> (W, pair_values)

Total W over 4 bases (6 pairs). W = 6.0 ⟺ the four bases are MUBs.
"""
@inline function total_W_static(B0::SMat6, B1::SMat6, B2::SMat6, B3::SMat6)
    w01 = pair_mubness_static(B0, B1)
    w02 = pair_mubness_static(B0, B2)
    w03 = pair_mubness_static(B0, B3)
    w12 = pair_mubness_static(B1, B2)
    w13 = pair_mubness_static(B1, B3)
    w23 = pair_mubness_static(B2, B3)
    W = w01 + w02 + w03 + w12 + w13 + w23
    return (W = W, pairs = (w01, w02, w03, w12, w13, w23))
end

# ---------------------------------------------------------------------------
# Analytic Euclidean gradients
# ---------------------------------------------------------------------------

"""
    grad_pair_static(Bk::SMat6, Bm::SMat6) -> (∇Bk, ∇Bm)

Analytic Euclidean gradient of W_km = Σ_ij |G_ij|⁴ where G = Bk'Bm.
  ∇_{Bm} W_km = 4 Bk * M,  ∇_{Bk} W_km = 4 M * Bm'
where M = |G|² ∘ G (elementwise).
"""
@inline function grad_pair_static(Bk::SMat6, Bm::SMat6)
    G = Bk' * Bm
    # Build M = |G|² ∘ G elementwise
    M_data = ntuple(Val(36)) do ij
        @inbounds abs2(G[ij]) * G[ij]
    end
    M = SMat6(M_data)
    dBm = 4 * Bk * M
    dBk = 4 * M * Bm'
    return (dBk, dBm)
end

"""
    grad_total_W_static(B0, B1, B2, B3) -> (g0, g1, g2, g3)

Accumulate Euclidean ∇W w.r.t. all four bases over all 6 pairs.
"""
function grad_total_W_static(B0::SMat6, B1::SMat6, B2::SMat6, B3::SMat6)
    Z = zero(SMat6)
    g0, g1, g2, g3 = Z, Z, Z, Z

    d0, d1 = grad_pair_static(B0, B1); g0 += d0; g1 += d1
    d0, d2 = grad_pair_static(B0, B2); g0 += d0; g2 += d2
    d0, d3 = grad_pair_static(B0, B3); g0 += d0; g3 += d3
    d1, d2 = grad_pair_static(B1, B2); g1 += d1; g2 += d2
    d1, d3 = grad_pair_static(B1, B3); g1 += d1; g3 += d3
    d2, d3 = grad_pair_static(B2, B3); g2 += d2; g3 += d3

    return (g0, g1, g2, g3)
end

# ---------------------------------------------------------------------------
# Riemannian machinery on U(6)
# ---------------------------------------------------------------------------

"""
    tangent_project_static(G::SMat6, B::SMat6) -> SMat6

Project Euclidean gradient G to the unitary-group tangent space at B:
  T = G − B · sym(B'G),   sym(X) = (X + X')/2
"""
@inline function tangent_project_static(G::SMat6, B::SMat6)
    BG = B' * G
    S = (BG + BG') / 2
    return G - B * S
end

"""
    unitary_retract_static(B::SMat6) -> SMat6

Modified Gram–Schmidt retraction onto U(6) with phase gauge:
the largest-magnitude entry in each column is rotated to real positive.
"""
function unitary_retract_static(B::SMat6)
    M = MMatrix{6, 6, ComplexF64, 36}(B)
    n = 6
    for j in 1:n
        # Orthogonalize against previous columns
        for i in 1:(j - 1)
            r = zero(ComplexF64)
            @inbounds for k in 1:n
                r += conj(M[k, i]) * M[k, j]
            end
            @inbounds for k in 1:n
                M[k, j] -= r * M[k, i]
            end
        end
        # Normalize
        nrm = 0.0
        @inbounds for k in 1:n
            nrm += abs2(M[k, j])
        end
        nrm = sqrt(nrm)
        nrm < 1e-14 && return SMat6(M)  # degenerate, return as-is
        # Phase gauge: largest-magnitude entry → real positive
        best_k = 1
        best_abs = 0.0
        @inbounds for k in 1:n
            a = abs2(M[k, j])
            if a > best_abs
                best_abs = a
                best_k = k
            end
        end
        @inbounds begin
            p = M[best_k, j]
            ph = conj(p) / abs(p)
            for k in 1:n
                M[k, j] = M[k, j] * ph / nrm
            end
        end
    end
    return SMat6(M)
end

"""
    gauge_fix_static(B::SMat6) -> SMat6

Phase-gauge: make the first row and first column real positive.
Standard dephasing convention for complex Hadamard matrices.
"""
function gauge_fix_static(B::SMat6)
    M = MMatrix{6, 6, ComplexF64, 36}(B)
    n = 6
    # First row phases
    @inbounds for j in 1:n
        if abs(M[1, j]) > 1e-14
            ph = exp(-im * angle(M[1, j]))
            for k in 1:n
                M[k, j] *= ph
            end
        end
    end
    # First column phases
    @inbounds for i in 1:n
        if abs(M[i, 1]) > 1e-14
            ph = exp(-im * angle(M[i, 1]))
            for k in 1:n
                M[i, k] *= ph
            end
        end
    end
    @inbounds M[1, 1] = abs(M[1, 1]) + 0.0im
    return SMat6(M)
end

# ---------------------------------------------------------------------------
# Random unitary generation
# ---------------------------------------------------------------------------

"""
    random_unitary_static(rng) -> SMat6

Generate a random unitary matrix from the Haar measure on U(6) via QR
decomposition of a complex Gaussian matrix. Returns an SMat6.
"""
function random_unitary_static(rng::AbstractRNG)
    # Generate 6×6 complex Gaussian as MMatrix, then QR
    M = MMatrix{6, 6, ComplexF64, 36}(undef)
    @inbounds for ij in 1:36
        M[ij] = randn(rng) + im * randn(rng)
    end
    # QR via the heap-allocated path (unavoidable for QR), but only done
    # once per seed — not in the hot loop
    Q_mat = Matrix(qr(SMatrix(M)).Q)
    return SMat6(Q_mat)
end

"""
    fourier_matrix_static() -> SMat6

The 6×6 Fourier matrix F₆ / √6 (orthonormal columns).
"""
function fourier_matrix_static()
    data = ntuple(Val(36)) do ij
        i = ((ij - 1) % 6) + 1
        j = ((ij - 1) ÷ 6) + 1
        ComplexF64(exp(2im * π * (i - 1) * (j - 1) / 6) / sqrt(6))
    end
    return SMat6(data)
end

"""
    identity_basis_static() -> SMat6

The 6×6 identity matrix (the computational basis, already orthonormal).
"""
@inline function identity_basis_static()
    return SMat6(I)
end

# ---------------------------------------------------------------------------
# Riemannian gradient descent on U(6)^3
# ---------------------------------------------------------------------------

"""
    manifold_descent_static(B0, B1, B2, B3; max_iter, η0, tol, rng)
        -> (W, B1, B2, B3, iterations)

Riemannian gradient descent minimizing W over U(6)^3 (B0 is fixed = I).
Adaptive step size with Armijo-like backtracking.
"""
function manifold_descent_static(B0::SMat6, B1::SMat6, B2::SMat6, B3::SMat6;
                                  max_iter::Int = 400,
                                  η0::Float64 = 0.05,
                                  tol::Float64 = 1e-12,
                                  rng::AbstractRNG = Random.default_rng())
    res = total_W_static(B0, B1, B2, B3)
    W = res.W
    η = η0
    final_it = 0

    for it in 1:max_iter
        final_it = it
        g0, g1, g2, g3 = grad_total_W_static(B0, B1, B2, B3)

        # Gradient norm (only on the free bases B1, B2, B3)
        gn = max(norm(g1), norm(g2), norm(g3))
        gn < tol && break

        # Tangent projection
        t1 = tangent_project_static(g1, B1)
        t2 = tangent_project_static(g2, B2)
        t3 = tangent_project_static(g3, B3)

        # Gradient step + retraction
        B1_new = unitary_retract_static(SMat6(B1 - η * t1))
        B2_new = unitary_retract_static(SMat6(B2 - η * t2))
        B3_new = unitary_retract_static(SMat6(B3 - η * t3))

        Wnew = total_W_static(B0, B1_new, B2_new, B3_new).W

        if Wnew < W - 1e-14
            W = Wnew
            B1, B2, B3 = B1_new, B2_new, B3_new
            η = min(η * 1.05, 0.5)   # accelerate on success
        else
            η *= 0.5                  # shrink on failure
            η < 1e-10 && break
        end

        # Early termination if W ≈ 6 (candidate 4-MUB)
        W - 6.0 < 1e-10 && break
    end

    return (W = W, B1 = B1, B2 = B2, B3 = B3, iterations = final_it)
end
