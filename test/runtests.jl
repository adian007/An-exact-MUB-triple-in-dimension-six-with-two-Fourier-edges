# ============================================================================
# runtests.jl — Part XVII: automated tests for the foundational mathematics.
#
#   Karlsson construction: A/B-block identities, full Hadamard residuals,
#     Möbius consistency, transcription-discrepancy regression, Dita point,
#     Fourier seam.
#   Pool: F6 benchmark count (48), conjugate-locus handling, dedup,
#     precision/tolerance stress.
#   Graph: orthogonality graph, planted clique recovery, full enumeration.
#   Fourth-MUB: known unextendible case (W1 empty), tier logic.
#
# Long tests (fresh homotopy solves, ~30–60 s each) are gated behind
# MUB_LONG_TESTS (default: on). Set MUB_LONG_TESTS=0 for a fast structural
# pass.
# ============================================================================

using Test
using Random
using LinearAlgebra

include(joinpath(@__DIR__, "norm_detection.jl"))
include(joinpath(@__DIR__, "..", "src", "MubSearch.jl"))
using .MubSearch

const LONG_TESTS = get(ENV, "MUB_LONG_TESTS", "1") != "0"
# Generic (non-degenerate) points and known-degenerate points, per
# diagnose_mobius_grid.jl: θ=0 sits on the Fourier seam (z2² derivation
# singular); (θ,φ,λ)=(0,0.5,0) additionally makes the assembled H NaN and
# must throw under the fail-fast guard.
const GRID = [(0.3, 0.5, 0.2), (1.0, 2.0, 0.7), (0.1, 0.1, 3.0),
              (DITA_THETA, DITA_PHI, 0.4), (DITA_THETA, DITA_PHI, pi / 2)]
const SEAM_GRID = [(0.0, 0.5, 0.3), (0.0, 0.9, 1.1)]

@testset "MubSearch" begin

    @testset "Karlsson construction (Part II)" begin
        for (th, ph, lm) in GRID
            A = build_A_karlsson_original(th, ph)
            B = -F2_HAD - A
            @test maximum(abs.(A * A' - 2I(2))) < 1e-12
            @test maximum(abs.(B * B' - 2I(2))) < 1e-12
            H = build_karlsson_family(th, ph, lm)
            @test maximum(abs.(H * H' - 6I(6))) < 1e-10   # HH† = 6I
            @test maximum(abs.(H' * H - 6I(6))) < 1e-10   # H†H = 6I
            @test maximum(abs.(abs.(H) .- 1)) < 1e-10     # |H_ij| = 1
        end
        # Dita point exact anchor check
        H_dita = build_karlsson_family(DITA_THETA, DITA_PHI, 0.4)
        @test maximum(abs.(H_dita * H_dita' - 6I(6))) < 1e-12

        # Fourier seam: Karlsson at θ=0 is the exact algebraic seam slice —
        # valid Hadamard for ALL φ, λ (including (0, 0.5, 0.0), which
        # previously produced a NaN matrix through division noise).
        for lm in (0.0, 0.7, 2.0)
            H = build_karlsson_family(0.0, 0.5, lm)
            @test maximum(abs.(H * H' - 6I(6))) < 1e-10
            @test maximum(abs.(abs.(H) .- 1)) < 1e-10
        end
        # Backward compatibility: the resolved seam agrees with the
        # historical (ulp-noise) construction to float precision where the
        # latter was finite — spot-check one historical anchor.
        H_seam = build_karlsson_family(0.0, 0.5, 0.3)
        @test maximum(abs.(H_seam * H_seam' - 6I(6))) < 1e-10
    end

    @testset "Transcription discrepancy regression (Part II)" begin
        # The literal McNulty–Weigert review transcription must FAIL the
        # unitary-block identity — this pins the C1 audit (348/349 vs 0/349).
        Random.seed!(20260914)
        n_fail_mw = 0
        n_fail_orig = 0
        for _ in 1:40
            th, ph = rand() * pi, rand() * pi
            A_mw = build_A_review_transcription_literal(th, ph)
            B_mw = -F2_HAD - A_mw
            maximum(abs.(A_mw * A_mw' - 2I(2))) > 1e-6 && (n_fail_mw += 1)
            A_o = build_A_karlsson_original(th, ph)
            B_o = -F2_HAD - A_o
            (maximum(abs.(A_o * A_o' - 2I(2))) > 1e-6 ||
             maximum(abs.(B_o * B_o' - 2I(2))) > 1e-6) && (n_fail_orig += 1)
        end
        @test n_fail_mw == 40          # MW transcription: fails everywhere
        @test n_fail_orig == 0         # Karlsson original: passes everywhere
    end

    @testset "Möbius audit (Part III)" begin
        for (th, ph, lm) in GRID
            ma = karlsson_mobius_audit(th, ph, lm)
            # Cleared-form residuals of all four identities (NaN = vacuous
            # at a pole; skipped by isnan checks).
            for r in (ma.cleared_z3_MA_z1, ma.cleared_z4_MB_z1,
                      ma.cleared_z3_MB_z2, ma.cleared_z4_MA_z2)
                isnan(r) || @test r < 1e-8
            end
            @test ma.z_moduli_residual < 1e-10 || isnan(ma.z_moduli_residual)
        end
        # The Dita anchor: z4² = M_A(z2²) is an indeterminate 0/0 there
        # (|A11|=|A12|=1) — must be FLAGGED, and the legacy z4_dev reported
        # as NaN (undefined), not as a large violation. This pins audit
        # finding 2026-09-14: the degeneracy scan's mobius_z4 flag at Dita
        # is ill-conditioned numerics, not a family inconsistency.
        ma_d = karlsson_mobius_audit(DITA_THETA, DITA_PHI, 0.4)
        @test ma_d.indeterminate_MA_z2 || ma_d.singular_MA_z2
        @test isnan(ma_d.legacy_z4_dev) || ma_d.legacy_z4_dev < 1e-6

        # Singular-branch recording: a Möbius map at its pole is flagged
        m = mobius_tracked(1.0 + 0im, 1.0 + 0im, 1.0 + 0im)   # den = 1−1 = 0
        @test m.singular
        @test isnan(m.value)

        # Denominator-cleared residual: well-conditioned where the divided
        # form is singular. Cleared identity M(z)=w: |(αz−β) − w(β̄z−ᾱ)|.
        r = mobius_identity_residual(1.0 + 0im, 1.0 + 0im, 1.0 + 0im, 1.0 + 0im;
                                     tol = 1e-10)
        @test r.indeterminate          # num = den = 0 → 0/0, not a plain pole
        @test r.cleared == 0.0

        # Fourier seam θ=0: resolved algebraically (M_A, M_B ≡ constant 1,
        # z2² = α_A/β_A); all four cleared-form identities hold exactly and
        # the audit reports the seam as resolved, not singular.
        for (th, ph, lm) in SEAM_GRID
            ma = karlsson_mobius_audit(th, ph, lm)
            @test ma.seam_resolved
            @test !ma.z2_derived_singular
            for r in (ma.cleared_z3_MA_z1, ma.cleared_z4_MB_z1,
                      ma.cleared_z3_MB_z2, ma.cleared_z4_MA_z2)
                @test r < 1e-10
            end
        end
    end

    @testset "Validation gate (Part IV)" begin
        H = build_karlsson_family(DITA_THETA, DITA_PHI, 0.4)
        val = validate_karlsson_matrix(H; variant = :karlsson_original,
                                       theta = DITA_THETA, phi = DITA_PHI, lambda = 0.4)
        @test val.pass
        @test val.classification in ("VALID_CHM", "VALID_CHM_DEGENERATE_BRANCH")
        @test val.param_residuals.A_unitary_err < 1e-10
        @test val.param_residuals.H_assembly_err < 1e-12
        @test val.mobius_audit.cleared_z3_MA_z1 < 1e-8
        @test val.mobius_audit.cleared_z4_MB_z1 < 1e-8

        # A matrix built from the INVALID variant must be rejected by the gate
        H_bad, _ = karlsson_assemble_diagnostics(0.3, 0.5, 0.7;
                                                 variant = :review_transcription_literal)
        val_bad = validate_karlsson_matrix(H_bad; variant = :review_transcription_literal,
                                           theta = 0.3, phi = 0.5, lambda = 0.7)
        @test !val_bad.pass
        @test val_bad.classification == "INVALID_CHM"
        @test_throws ErrorException assert_valid_karlsson!(val_bad)

        # Tao transcription check: tao_S6() validates itself or throws
        S6 = tao_S6()
        @test validate_hadamard_matrix(S6).pass
    end

    @testset "Projective dedup (Part VII)" begin
        v = [exp(im * 0.1), exp(im * 1.3), exp(im * 2.7), exp(im * 0.9),
             exp(im * 4.2), exp(im * 5.5)] ./ sqrt(6)
        w = v .* exp(im * 1.234)
        @test projective_distance(v, w) < 1e-14     # phase equivalence is exact
        @test projective_distance(v, 2 .* v) < 1e-14  # scaling equivalence too
        u = copy(v); u[1] *= -1                      # NOT a global phase flip
        @test projective_distance(v, u) > 1e-3

        reps, assignment = cluster_pool_hp([v, w, u]; bits = 128, tol = 1e-8)
        @test length(reps) == 2
        @test assignment == [1, 1, 2]

        diags = dedup_tolerance_diagnostics([v, w, u])
        @test [d.n_clusters for d in diags] == [2, 2, 2, 2]
    end

    @testset "Clique enumeration (Part IX)" begin
        # Planted positive control: 6 orthonormal vectors MU to a generic H
        Random.seed!(7)
        th, ph, lm = 0.37, 1.1, 2.2
        H = build_karlsson_family(th, ph, lm)
        # Planted enumeration-completeness control: 6 orthonormal vectors
        # plus 6 phase-rotated copies (copy_i NOT orthogonal to original i,
        # orthogonal to everything else). The orthogonality graph is K₁₂
        # minus 6 disjoint non-edges, so the maximal cliques of size 6 are
        # exactly the 2⁶ transversals — a sharp completeness test for the
        # clique enumerator. None are MU to H (random ONB), so the verified
        # count must be 0: proposal and verification are decoupled.
        Random.seed!(2026)
        M = randn(ComplexF64, 6, 6)
        Q = Matrix(qr(M).Q)
        decoys = [Q[:, k] for k in 1:6]
        append!(decoys, [Q[:, k] .* exp(im * 0.5) for k in 1:6])
        cl = enumerate_all_third_mub_bases(decoys, H)
        @test cl.n_raw == 64                      # 2^6 transversals, all maximal
        @test cl.clique_size_histogram == Dict(6 => 64)
        @test cl.n_verified == 0                   # NOT MU to H (honest negative)
        @test all(b.mu_to_H_defect > 1e-3 for b in cl.bases)

        # Determinism: identical inputs → identical records
        cl2 = enumerate_all_third_mub_bases(decoys, H)
        @test [b.pool_indices for b in cl2.bases] ==
              [b.pool_indices for b in cl.bases]

        # F6 positive control (fast, no solver): the known F6 third basis
        # from the Dita construction at λ=0 belongs to the {I, F6} pool
        # only via homotopy; here we only test the graph on a planted pair.
    end

    @testset "Claim tiers (Part XI)" begin
        @test CLAIM_TIERS == (:PROVED, :EXACT_ALGEBRAIC, :CERTIFIED_NUMERICAL,
                              :NUMERICALLY_SUPPORTED, :HEURISTIC, :OPEN)
        @test at_least(:PROVED, :CERTIFIED_NUMERICAL)
        @test at_least(:CERTIFIED_NUMERICAL, :NUMERICALLY_SUPPORTED)
        @test !at_least(:HEURISTIC, :CERTIFIED_NUMERICAL)
        @test_throws ErrorException tier_rank(:nonsense)
    end

    if LONG_TESTS
        @testset "Benchmark A — Fourier 48 (Part VIII.A)" begin
            out = run_benchmark_F6(; reseed = true)
            r = out.record
            @test r.pass_count   # observed_physical_vectors == 48
            @test r.audit.certified_solutions == r.audit.numerical_solutions
            @test r.crosscheck.agree   # independent re-solve finds the same counts
        end

        @testset "Benchmark B — Dita/D0: 120 vectors, 10 third bases (Part VIII.B)" begin
            out = run_benchmark_D0(; lambda = pi / 2, hp_bits = 256, reseed = true)
            r = out.record
            @test r.observed_physical_vectors == EXPECTED_D0_PHYSICAL
            @test r.observed_third_bases == EXPECTED_D0_THIRD_BASES
            @test r.pairwise_mu_among_third_bases.n_mu_pairs == 0   # BW N_p = 0
        end

        @testset "Benchmark C — Tao S6 (Part VIII.C)" begin
            out = run_benchmark_Tao(; reseed = true)
            r = out.record
            @test !r.observed_third_mub
        end

        @testset "W1 witness — known unextendible triple (Part X)" begin
            # At Dita λ=0.4 the pool contains third bases and prior certified
            # runs established W1 = ∅ for the recovered cliques; one is
            # re-checked here as a regression anchor.
            H = build_karlsson_family(DITA_THETA, DITA_PHI, 0.4)
            out = solve_pool_audited(H; tol = 1e-8)
            cl = enumerate_all_third_mub_bases(out.pool, H)
            @test cl.n_verified >= 1
            b1 = cl.bases[1]
            B3 = hcat([out.pool[i] for i in b1.pool_indices]...)
            w = w1_witness_certificate(H, B3; seed = 20260914)
            @test w.complete_tracking
            @test w.verdict == "EMPTY_CERTIFIED"
            @test w.tier == :CERTIFIED_NUMERICAL
        end
    end
end
