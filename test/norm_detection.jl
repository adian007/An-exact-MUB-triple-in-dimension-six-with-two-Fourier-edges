using Test
using LinearAlgebra
using Random

function detect_normalization(B3::AbstractMatrix; tolerance::Real=1e-10)
    size(B3) == (6, 6) || throw(ArgumentError("B3 must be a 6×6 matrix"))
    isfinite(tolerance) && tolerance >= 0 || throw(ArgumentError("tolerance must be finite and nonnegative"))
    all(isfinite, B3) || throw(ArgumentError("B3 must contain only finite values"))

    gram = B3' * B3
    unitary_err = maximum(abs.(gram - I(6)))
    hadamard_err = maximum(abs.(gram - 6.0 * I(6)))
    if unitary_err <= tolerance && hadamard_err > tolerance
        rhs_val = 1.0
        scale = "unitary"
    elseif hadamard_err <= tolerance && unitary_err > tolerance
        rhs_val = 6.0
        scale = "unnormalized"
    else
        rhs_val = nothing
        scale = "ambiguous"
    end
    return (unitary_err, hadamard_err, rhs_val, scale, gram)
end

# TEST A — Unit-normalized B3 (RHS → 1.0)
@testset "TEST A: Unit-normalized B3" begin
    # A1: B3 = I → B3'B3 = I
    # unitary_err = max|I-I|=0, hadamard_err = max|I-6I|=5
    # Gap = 5.0 → UNAMBIGUOUS
    B3_a1 = Matrix{ComplexF64}(I, 6, 6)
    (uerr, herr, rhs, sc, gram) = detect_normalization(B3_a1)
    @test uerr == 0.0; @test herr == 5.0; @test rhs == 1.0; @test sc == "unitary"
    @test gram ≈ I(6); @test (herr - uerr) == 5.0
    println("[A1] gap=$(herr-uerr) → rhs=$rhs ($sc) ✓")

    # A2: B3 = random unitary Q (QR), B3'B3 ≈ I
    Random.seed!(42)
    Q, _ = qr(randn(ComplexF64, 6, 6))
    B3_a2 = Q
    (uerr, herr, rhs, sc, gram) = detect_normalization(B3_a2)
    @test uerr < 1e-14; @test herr ≈ 5.0 atol=1e-14
    @test rhs == 1.0; @test sc == "unitary"
    @test (herr - uerr) > 4.999
    println("[A2] gap=$(herr-uerr) → rhs=$rhs ($sc) ✓")

    # A3: B3 = permutation matrix (unitary), B3'B3 = I
    P = [0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 1 0 0;
         0 0 0 0 1 0; 0 0 0 0 0 1; 1 0 0 0 0 0] |> x->x+0im
    B3_a3 = P
    (uerr, herr, rhs, sc, gram) = detect_normalization(B3_a3)
    @test uerr == 0.0; @test herr == 5.0; @test rhs == 1.0; @test sc == "unitary"
    println("[A3] gap=$(herr-uerr) → rhs=$rhs ($sc) ✓")
end

# TEST B — Unnormalized B3 (RHS → 6.0)
@testset "TEST B: Unnormalized B3" begin
    # B1: B3 = sqrt(6)*I → B3'B3 = 6*I
    # unitary_err = |6-1|=5, hadamard_err = |6-6|=0
    # Gap = -5.0 → UNAMBIGUOUS
    B3_b1 = sqrt(6.0) * Matrix{ComplexF64}(I, 6, 6)
    (uerr, herr, rhs, sc, gram) = detect_normalization(B3_b1)
    @test uerr == 5.0; @test herr == 0.0; @test rhs == 6.0; @test sc == "unnormalized"
    @test gram ≈ 6.0 * I(6); @test (herr - uerr) == -5.0
    println("[B1] gap=$(herr-uerr) → rhs=$rhs ($sc) ✓")

    # B2: B3 = sqrt(6)*Q for random unitary Q
    Random.seed!(123)
    Q, _ = qr(randn(ComplexF64, 6, 6))
    B3_b2 = sqrt(6.0) * Q
    (uerr, herr, rhs, sc, gram) = detect_normalization(B3_b2)
    @test uerr ≈ 5.0 atol=1e-14; @test herr < 1e-14
    @test rhs == 6.0; @test sc == "unnormalized"
    @test (herr - uerr) < -4.999
    println("[B2] gap=$(herr-uerr) → rhs=$rhs ($sc) ✓")

    # B3: B3 = sqrt(6)*P (permutation matrix)
    P = [0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 1 0 0;
         0 0 0 0 1 0; 0 0 0 0 0 1; 1 0 0 0 0 0] |> x->x+0im
    B3_b3 = sqrt(6.0) * P
    (uerr, herr, rhs, sc, gram) = detect_normalization(B3_b3)
    @test uerr == 5.0; @test herr == 0.0; @test rhs == 6.0; @test sc == "unnormalized"
    println("[B3] gap=$(herr-uerr) → rhs=$rhs ($sc) ✓")
end


# TEST C — Ambiguous/edge cases
@testset "TEST C: Edge cases" begin
    # C1: B3 = I + 1e-12*E (tiny perturbation)
    # B3'B3 ≈ I + 2e-12*E → unitary_err ≈ O(1e-12), hadamard_err ≈ 5
    Random.seed!(7)
    E = randn(ComplexF64, 6, 6)
    B3_c1 = Matrix{ComplexF64}(I, 6, 6) + 1e-12 * E
    (uerr, herr, rhs, sc, gram) = detect_normalization(B3_c1)
    gap = herr - uerr
    @test uerr < 1e-11; @test herr ≈ 5.0 atol=1e-11
    @test rhs == 1.0; @test sc == "unitary"
    @test gap > 4.999
    println("[C1] B3=I+1e-12·E: gap=$gap → rhs=$rhs ($sc) ✓")

    # C2: B3 = (1+ε)*I, ε=1e-12
    eps2 = 1e-12
    B3_c2 = (1.0 + eps2) * Matrix{ComplexF64}(I, 6, 6)
    (uerr, herr, rhs, sc, gram) = detect_normalization(B3_c2)
    gap = herr - uerr
    @test uerr ≈ 2*eps2 atol=5e-16 rtol=0
    @test herr ≈ 5.0-2*eps2 atol=5e-16 rtol=0
    @test rhs == 1.0; @test sc == "unitary"; @test gap > 4.999
    println("[C2] B3=(1+ε)·I: gap=$gap → rhs=$rhs ($sc) ✓")

    # C3: B3 = sqrt(3.5)*I — EXACTLY equidistant from I and 6I
    # B3'B3 = 3.5I
    # unitary_err=2.5, hadamard_err=2.5, gap=0 → TIE
    # ⚠️ Code: 2.5 < 2.5 → FALSE → ELSE: rhs=6.0 (silent!)
    B3_c3 = sqrt(3.5) * Matrix{ComplexF64}(I, 6, 6)
    (uerr, herr, rhs, sc, gram) = detect_normalization(B3_c3)
    gap = herr - uerr
    @test uerr ≈ 2.5 atol=1e-14 rtol=0
    @test herr ≈ 2.5 atol=1e-14 rtol=0
    @test abs(gap) < 1e-14
    @test isnothing(rhs); @test sc == "ambiguous"
    println("[C3] B3=√3.5·I: gap=$gap → rhs=$rhs ($sc)")

    # C4: B3 = sqrt(3.5+ε)*I, ε=1e-12 (just above midpoint)
    # d=3.5+ε → unitary_err=2.5+ε, hadamard_err=2.5-ε
    # Gap=-2ε≈-2e-12 → hadamard barely wins → rhs=6.0
    # The gap is below the classification tolerance, so neither scale is accepted.
    eps4 = 1e-12
    B3_c4 = sqrt(3.5 + eps4) * Matrix{ComplexF64}(I, 6, 6)
    (uerr, herr, rhs, sc, gram) = detect_normalization(B3_c4)
    gap = herr - uerr
    @test uerr ≈ 2.5+eps4 atol=2e-15 rtol=0
    @test herr ≈ 2.5-eps4 atol=2e-15 rtol=0
    @test gap ≈ -2*eps4 atol=2e-15 rtol=0
    @test isnothing(rhs); @test sc == "ambiguous"
    println("[C4] B3=√(3.5+ε)·I: gap=$gap → rhs=$rhs ($sc)")

    # C5: B3 = sqrt(3.5-ε)*I, ε=1e-12 (just below midpoint)
    # d=3.5-ε → unitary_err=2.5-ε, hadamard_err=2.5+ε
    # Gap=+2ε≈+2e-12 → unitary barely wins → rhs=1.0
    # Same midpoint; the candidate is not close to either normalization target.
    B3_c5 = sqrt(3.5 - eps4) * Matrix{ComplexF64}(I, 6, 6)
    (uerr, herr, rhs, sc, gram) = detect_normalization(B3_c5)
    gap = herr - uerr
    @test uerr ≈ 2.5-eps4 atol=2e-15 rtol=0
    @test herr ≈ 2.5+eps4 atol=2e-15 rtol=0
    @test gap ≈ 2*eps4 atol=2e-15 rtol=0
    @test isnothing(rhs); @test sc == "ambiguous"
    println("[C5] B3=√(3.5-ε)·I: gap=$gap → rhs=$rhs ($sc)")


# C6: Mixed normalization — 4 unitary + 2 unnormalized columns
# B3 = diag(1,1,1,1,sqrt(6),sqrt(6))
# B3'B3 = diag(1,1,1,1,6,6)
# unitary_err=max(0,0,0,0,5,5)=5, hadamard_err=max(5,5,5,5,0,0)=5
# Gap=0 → EXACT TIE (even with 4/6 unitary!)
B3_c6a = Diagonal([1.0+0im, 1.0+0im, 1.0+0im, 1.0+0im, sqrt(6)+0im, sqrt(6)+0im])
(uerr, herr, rhs, sc, gram) = detect_normalization(B3_c6a)
gap = herr - uerr
@test uerr == 5.0; @test herr == 5.0; @test gap == 0.0
@test isnothing(rhs); @test sc == "ambiguous"
println("[C6a] Mixed(4U+2H): gap=$gap → rhs=$rhs ($sc)")

# C6b: Perturbed mixed — one unnormalized column slightly larger
# diag(1,1,1,1,6.01,6): unitary_err=5.01, hadamard_err=5
# The unnormalized target is closer, but this mixed Gram matrix is ambiguous.
B3_c6b = Diagonal([1.0+0im, 1.0+0im, 1.0+0im, 1.0+0im, sqrt(6.01)+0im, sqrt(6)+0im])
(uerr, herr, rhs, sc, gram) = detect_normalization(B3_c6b)
gap = herr - uerr
@test uerr ≈ 5.01 atol=1e-12 rtol=0
@test herr ≈ 5.0 atol=1e-12 rtol=0
@test gap ≈ -0.01 atol=1e-12 rtol=0
@test isnothing(rhs); @test sc == "ambiguous"
println("[C6b] Mixed(4U+2H,6.01): gap=$gap → rhs=$rhs ($sc)")

# C7: Mixed — 5 unitary + 1 unnormalized
# B3'B3 = diag(1,1,1,1,1,6)
# unitary_err=max(0,0,0,0,0,5)=5, hadamard_err=max(5,5,5,5,5,0)=5
# Gap=0 → TIE → unnormalized (even 5/6 unitary can't win!)
B3_c7 = Diagonal([1.0+0im, 1.0+0im, 1.0+0im, 1.0+0im, 1.0+0im, sqrt(6)+0im])
(uerr, herr, rhs, sc, gram) = detect_normalization(B3_c7)
gap = herr - uerr
@test uerr == 5.0; @test herr == 5.0; @test gap == 0.0
@test isnothing(rhs); @test sc == "ambiguous"
println("[C7] Mixed(5U+1H): gap=$gap → rhs=$rhs ($sc)")

# C8: Moderate perturbation near midpoint (gap < 1e-6)
B3_c8 = sqrt(3.5 + 1e-8) * Matrix{ComplexF64}(I, 6, 6)
(uerr, herr, rhs, sc, gram) = detect_normalization(B3_c8)
gap = herr - uerr
@test abs(gap) < 1e-6
@test isnothing(rhs); @test sc == "ambiguous"
println("[C8] B3=√(3.5+1e-8)·I: gap=$gap → rhs=$rhs ($sc)")

# C9: Random non-orthogonal matrix (heuristic pick)
Random.seed!(99)
B3_c9 = randn(ComplexF64, 6, 6)
(uerr, herr, rhs, sc, gram) = detect_normalization(B3_c9)
gap = herr - uerr
@test !isnan(uerr) && !isinf(uerr)
@test isnothing(rhs); @test sc == "ambiguous"
println("[C9] B3=random: gap=$gap → rhs=$rhs ($sc)")
end

@testset "Input validation" begin
    @test_throws ArgumentError detect_normalization(zeros(ComplexF64, 5, 6))
    @test_throws ArgumentError detect_normalization(Matrix{ComplexF64}(I, 6, 6); tolerance=-1)
    nonfinite = Matrix{ComplexF64}(I, 6, 6)
    nonfinite[1, 1] = ComplexF64(NaN, 0)
    @test_throws ArgumentError detect_normalization(nonfinite)
end

# Margin safety analysis
@testset "Margin safety analysis" begin
    unitary = detect_normalization(Matrix{ComplexF64}(I, 6, 6))
    unnormalized = detect_normalization(sqrt(6.0) * Matrix{ComplexF64}(I, 6, 6))
    midpoint = detect_normalization(sqrt(3.5) * Matrix{ComplexF64}(I, 6, 6))
    near_midpoint = detect_normalization(sqrt(3.5 + 1e-8) * Matrix{ComplexF64}(I, 6, 6))

    @test unitary[4] == "unitary"
    @test unnormalized[4] == "unnormalized"
    @test midpoint[4] == "ambiguous"
    @test near_midpoint[4] == "ambiguous"
end

# Analytical notes (static analysis, not executable)
#=
KEY ANALYTICAL RESULTS
==========================================

1. MIDPOINT: For diagonal entry d of B3'B3, the decision boundary
   is d = 3.5 where |d-1| = |d-6| = 2.5.

2. OFF-DIAGONAL ENTRIES: The target matrices have zero off-diagonal entries,
    so off-diagonal Gram entries contribute equally to both errors. However,
    the maximum can be dominated by those entries, so they cannot be ignored.

3. MAX-ABS DOMINANCE: Because metric = maximum(abs.()), ONE entry
   with large deviation overrides ALL others. This means:
    - A single large residual can dominate either maximum.
    - Mixed-scale matrices can tie or be far from both accepted targets.

4. AMBIGUOUS INPUTS: Ties, mixed-scale inputs, and inputs not close to either
    target return no RHS and the `ambiguous` label.

5. NEAR-MIDPOINT INPUTS: B3 = sqrt(3.5±ε)·I gives |gap| ≈ 2ε.
    For ε=1e-12, this is larger than Float64 machine epsilon but below the
    classification tolerance; neither normalization target is accepted.

6. SAFE FOR NORMAL CASES: For genuinely unitary B3'B3≈I or
   unnormalized B3'B3≈6I, |gap|≈5.0 — enormous safety margin.

7. RANDOM MATRICES: A random matrix is not labeled as a normalization unless
    its Gram matrix is within tolerance of one of the two target matrices.

8. TOLERANCE: The default acceptance tolerance is `1e-10` in the maximum
    absolute Gram-entry residual and can be overridden by the caller.
=#


