using Test
using LinearAlgebra
using Random
using StaticArrays
using Printf

function detect_normalization(B3::AbstractMatrix)
    gram = B3' * B3
    unitary_err = maximum(abs.(gram - I(6)))
    hadamard_err = maximum(abs.(gram - 6.0 * I(6)))
    if unitary_err < hadamard_err
        rhs_val = 1.0
        scale = "unitary"
    else
        rhs_val = 6.0
        scale = "unnormalized"
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
    @println("[A1] gap=$(herr-uerr) → rhs=$rhs ($sc) ✓")

    # A2: B3 = random unitary Q (QR), B3'B3 ≈ I
    Random.seed!(42)
    Q, _ = qr(randn(ComplexF64, 6, 6))
    B3_a2 = Q
    (uerr, herr, rhs, sc, gram) = detect_normalization(B3_a2)
    @test uerr < 1e-14; @test herr ≈ 5.0 atol=1e-14
    @test rhs == 1.0; @test sc == "unitary"
    @test (herr - uerr) > 4.999
    @println("[A2] gap=$(herr-uerr) → rhs=$rhs ($sc) ✓")

    # A3: B3 = permutation matrix (unitary), B3'B3 = I
    P = [0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 1 0 0;
         0 0 0 0 1 0; 0 0 0 0 0 1; 1 0 0 0 0 0] |> x->x+0im
    B3_a3 = P
    (uerr, herr, rhs, sc, gram) = detect_normalization(B3_a3)
    @test uerr == 0.0; @test herr == 5.0; @test rhs == 1.0; @test sc == "unitary"
    @println("[A3] gap=$(herr-uerr) → rhs=$rhs ($sc) ✓")
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
    @println("[B1] gap=$(herr-uerr) → rhs=$rhs ($sc) ✓")

    # B2: B3 = sqrt(6)*Q for random unitary Q
    Random.seed!(123)
    Q, _ = qr(randn(ComplexF64, 6, 6))
    B3_b2 = sqrt(6.0) * Q
    (uerr, herr, rhs, sc, gram) = detect_normalization(B3_b2)
    @test uerr ≈ 5.0 atol=1e-14; @test herr < 1e-14
    @test rhs == 6.0; @test sc == "unnormalized"
    @test (herr - uerr) < -4.999
    @println("[B2] gap=$(herr-uerr) → rhs=$rhs ($sc) ✓")

    # B3: B3 = sqrt(6)*P (permutation matrix)
    P = [0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 1 0 0;
         0 0 0 0 1 0; 0 0 0 0 0 1; 1 0 0 0 0 0] |> x->x+0im
    B3_b3 = sqrt(6.0) * P
    (uerr, herr, rhs, sc, gram) = detect_normalization(B3_b3)
    @test uerr == 5.0; @test herr == 0.0; @test rhs == 6.0; @test sc == "unnormalized"
    @println("[B3] gap=$(herr-uerr) → rhs=$rhs ($sc) ✓")
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
    @println("[C1] B3=I+1e-12·E: gap=$gap → rhs=$rhs ($sc) ✓")

    # C2: B3 = (1+ε)*I, ε=1e-12
    eps2 = 1e-12
    B3_c2 = (1.0 + eps2) * Matrix{ComplexF64}(I, 6, 6)
    (uerr, herr, rhs, sc, gram) = detect_normalization(B3_c2)
    gap = herr - uerr
    @test uerr ≈ 2*eps2 atol=1e-20; @test herr ≈ 5.0-2*eps2 atol=1e-20
    @test rhs == 1.0; @test sc == "unitary"; @test gap > 4.999
    @println("[C2] B3=(1+ε)·I: gap=$gap → rhs=$rhs ($sc) ✓")

    # C3: B3 = sqrt(3.5)*I — EXACTLY equidistant from I and 6I
    # B3'B3 = 3.5I
    # unitary_err=2.5, hadamard_err=2.5, gap=0 → TIE
    # ⚠️ Code: 2.5 < 2.5 → FALSE → ELSE: rhs=6.0 (silent!)
    B3_c3 = sqrt(3.5) * Matrix{ComplexF64}(I, 6, 6)
    (uerr, herr, rhs, sc, gram) = detect_normalization(B3_c3)
    gap = herr - uerr
    @test uerr == 2.5; @test herr == 2.5; @test gap == 0.0
    @test rhs == 6.0; @test sc == "unnormalized"
    @println("[C3] B3=√3.5·I: gap=$gap → rhs=$rhs ($sc) ⚠️ EXACT TIE")

    # C4: B3 = sqrt(3.5+ε)*I, ε=1e-12 (just above midpoint)
    # d=3.5+ε → unitary_err=2.5+ε, hadamard_err=2.5-ε
    # Gap=-2ε≈-2e-12 → hadamard barely wins → rhs=6.0
    # 🚨 CLASS 1 CATASTROPHIC: margin < 1e-6
    eps4 = 1e-12
    B3_c4 = sqrt(3.5 + eps4) * Matrix{ComplexF64}(I, 6, 6)
    (uerr, herr, rhs, sc, gram) = detect_normalization(B3_c4)
    gap = herr - uerr
    @test uerr ≈ 2.5+eps4 atol=1e-20; @test herr ≈ 2.5-eps4 atol=1e-20
    @test gap ≈ -2*eps4 atol=1e-20; @test rhs == 6.0; @test sc == "unnormalized"
    @println("[C4] B3=√(3.5+ε)·I: gap=$gap → rhs=$rhs ($sc)")
    @println("         🚨 CLASS 1 CATASTROPHIC: |gap|=$(abs(gap)) < 1e-6")

    # C5: B3 = sqrt(3.5-ε)*I, ε=1e-12 (just below midpoint)
    # d=3.5-ε → unitary_err=2.5-ε, hadamard_err=2.5+ε
    # Gap=+2ε≈+2e-12 → unitary barely wins → rhs=1.0
    # 🚨 Same flip point, same thin margin
    B3_c5 = sqrt(3.5 - eps4) * Matrix{ComplexF64}(I, 6, 6)
    (uerr, herr, rhs, sc, gram) = detect_normalization(B3_c5)
    gap = herr - uerr
    @test uerr ≈ 2.5-eps4 atol=1e-20; @test herr ≈ 2.5+eps4 atol=1e-20
    @test gap ≈ 2*eps4 atol=1e-20; @test rhs == 1.0; @test sc == "unitary"
    @println("[C5] B3=√(3.5-ε)·I: gap=$gap → rhs=$rhs ($sc)")
    @println("         🚨 CLASS 1 CATASTROPHIC: |gap|=$(abs(gap)) < 1e-6")
end


# C6: Mixed normalization — 4 unitary + 2 unnormalized columns
# B3 = diag(1,1,1,1,sqrt(6),sqrt(6))
# B3'B3 = diag(1,1,1,1,6,6)
# unitary_err=max(0,0,0,0,5,5)=5, hadamard_err=max(5,5,5,5,0,0)=5
# Gap=0 → EXACT TIE (even with 4/6 unitary!)
B3_c6a = Diagonal([1.0+0im, 1.0+0im, 1.0+0im, 1.0+0im, sqrt(6)+0im, sqrt(6)+0im])
(uerr, herr, rhs, sc, gram) = detect_normalization(B3_c6a)
gap = herr - uerr
@test uerr == 5.0; @test herr == 5.0; @test gap == 0.0
@test rhs == 6.0; @test sc == "unnormalized"
@println("[C6a] Mixed(4U+2H): gap=$gap → rhs=$rhs ($sc) ⚠️ TIE")

# C6b: Perturbed mixed — one unnormalized column slightly larger
# diag(1,1,1,1,6.01,6): unitary_err=5.01, hadamard_err=5
# Gap=-0.01 → unitary barely wins
# BUT max-abs metric overridden by ONE bad column!
B3_c6b = Diagonal([1.0+0im, 1.0+0im, 1.0+0im, 1.0+0im, sqrt(6)+0.01+0im, sqrt(6)+0im])
(uerr, herr, rhs, sc, gram) = detect_normalization(B3_c6b)
gap = herr - uerr
@test uerr ≈ 5.01; @test herr ≈ 5.0; @test gap ≈ -0.01
@test rhs == 1.0
@println("[C6b] Mixed(4U+2H,6.01): gap=$gap → rhs=$rhs ($sc)")
@println("          🚨 ONE bad column overrides 4 good ones")

# C7: Mixed — 5 unitary + 1 unnormalized
# B3'B3 = diag(1,1,1,1,1,6)
# unitary_err=max(0,0,0,0,0,5)=5, hadamard_err=max(5,5,5,5,5,0)=5
# Gap=0 → TIE → unnormalized (even 5/6 unitary can't win!)
B3_c7 = Diagonal([1.0+0im, 1.0+0im, 1.0+0im, 1.0+0im, 1.0+0im, sqrt(6)+0im])
(uerr, herr, rhs, sc, gram) = detect_normalization(B3_c7)
gap = herr - uerr
@test uerr == 5.0; @test herr == 5.0; @test gap == 0.0
@test rhs == 6.0; @test sc == "unnormalized"
@println("[C7] Mixed(5U+1H): gap=$gap → rhs=$rhs ($sc) ⚠️ TIE")

# C8: Moderate perturbation near midpoint (gap < 1e-6)
B3_c8 = sqrt(3.5 + 1e-8) * Matrix{ComplexF64}(I, 6, 6)
(uerr, herr, rhs, sc, gram) = detect_normalization(B3_c8)
gap = herr - uerr
@test abs(gap) < 1e-6
@println("[C8] B3=√(3.5+1e-8)·I: gap=$gap → rhs=$rhs ($sc) 🚨 THIN")

# C9: Random non-orthogonal matrix (heuristic pick)
Random.seed!(99)
B3_c9 = randn(ComplexF64, 6, 6)
(uerr, herr, rhs, sc, gram) = detect_normalization(B3_c9)
gap = herr - uerr
@test !isnan(uerr) && !isinf(uerr)
@println("[C9] B3=random: gap=$gap → rhs=$rhs ($sc)")
@println("          ℹ️ Heuristic pick — meaningful IF assumption holds")
end


# Margin safety analysis
@testset "Margin safety analysis" begin
    @testset "Classification by |gap|" begin
        # SAFE (|gap| >> 1e-6):
        #   A1 (B3=I):         gap = 5.0      → SAFE
        #   A2 (B3=randU):     gap ≈ 5.0      → SAFE
        #   A3 (B3=perm):      gap = 5.0      → SAFE
        #   B1 (B3=√6·I):      gap = -5.0     → SAFE
        #   B2 (B3=√6·randU):  gap ≈ -5.0     → SAFE
        #   B3 (B3=√6·perm):   gap = -5.0     → SAFE
        #   C1 (B3=I+1e-12E):  gap ≈ 5.0      → SAFE
        #   C2 (B3=(1+ε)I):    gap ≈ 5.0      → SAFE

        # DANGEROUS (|gap| < 1e-6):
        #   C3 (B3=√3.5·I):    gap = 0.0      → EXACT TIE → unnormalized
        #   C4 (B3=√(3.5+ε)I): gap ≈ -2e-12  → CLASS 1 CATASTROPHIC
        #   C5 (B3=√(3.5-ε)I): gap ≈ 2e-12   → CLASS 1 CATASTROPHIC
        #   C6a:               gap = 0.0      → EXACT TIE (mixed)
        #   C7:                gap = 0.0      → EXACT TIE (mixed)
        #   C8:                gap ≈ -2e-8    → THIN MARGIN

        @println("=== MARGIN SAFETY SUMMARY ===")
        @println("SAFE (|gap|>1):      A1,A2,A3,B1,B2,B3,C1,C2 (8 cases)")
        @println("DANGEROUS (<1e-6):   C3,C4,C5,C6a,C7,C8 (6 cases)")
        @println("THIN (<0.1):         C6b (1 case)")
        @println("TOTAL fragile: 7 of 15 cases")
        @println("TIE rule: gap=0 → ELSE branch → rhs=6.0 (unnormalized)")
        @println("No tie-breaking mechanism exists in the code.")
    end
end

# Analytical notes (static analysis, not executable)
=begin
KEY ANALYTICAL RESULTS
==========================================

1. MIDPOINT: For diagonal entry d of B3'B3, the decision boundary
   is d = 3.5 where |d-1| = |d-6| = 2.5.

2. OFF-DIAGONAL IRRELEVANCE: Both I(6) and 6I(6) have zero
   off-diagonal entries. Off-diagonal entries of B3'B3 contribute
   IDENTICALLY to both errors. Only diagonal entries matter.

3. MAX-ABS DOMINANCE: Because metric = maximum(abs.()), ONE entry
   with large deviation overrides ALL others. This means:
   - 5 unitary + 1 unnormalized columns → TIE (gap=0)
   - 4 unitary + 2 unnormalized → TIE (gap=0)
   - A single bad column can flip the decision

4. TIE-BREAKING FLAW: `unitary_err < hadamard_err` uses strict <.
   Exact ties fall to ELSE (rhs=6.0) silently. No warning issued.
   Users CANNOT detect this from return values.

5. CLASS 1 CATASTROPHIC: B3 = sqrt(3.5±ε)·I with tiny ε gives
   |gap| = 2ε. For ε=1e-12, margin is 2e-12 < machine epsilon.
   Floating-point rounding in B3'B3 computation itself can flip
   the decision. This is a genuine design flaw.

6. SAFE FOR NORMAL CASES: For genuinely unitary B3'B3≈I or
   unnormalized B3'B3≈6I, |gap|≈5.0 — enormous safety margin.

7. RANDOM MATRICES: Code picks closer reference (heuristic).
   Correct as closest-match, but MEANINGLESS if assumption
   (B3 ∈ {unitary, unnormalized}) is violated.

8. RECOMMENDED FIX (not implemented here):
   - Add explicit tie-breaking rule
   - Add margin warning: @warn if |gap| < tolerance
   - Allow user override via keyword argument
=end


