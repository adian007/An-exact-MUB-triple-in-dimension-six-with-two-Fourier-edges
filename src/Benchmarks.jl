# ============================================================================
# Benchmarks.jl — Part VIII: hard regression benchmarks.
#
#   A. Fourier benchmark: {I, F6} must yield exactly 48 physical MU vectors.
#   B. Dita/D0 benchmark: at the Karlsson point CHM-equivalent to the
#      Brierley–Weigert matrix D0 (λ = π/2 and 3π/2 on the Dita slice),
#      the pool must yield 120 physical MU vectors forming 10 third bases
#      (BW 2009, PRA 79, 052316: N_v = 120, N_t = 10, N_p = 0).
#      This is the correctness gate for the whole research program.
#   C. Tao benchmark: the isolated order-6 CHM S6 of Terence Tao (entries
#      3rd roots of unity; explicit transcription from Wuttig–Tindall 2026,
#      arXiv:2608.18053, Definition 5 / Eq. (10)).
#
# Every benchmark returns a structured record with expected vs observed
# values, the full PoolAudit, the clique enumeration, and a claim tier.
# A benchmark FAILS loudly; results are never adjusted to pass.
# ============================================================================

"""
    tao_S6()

Terence Tao's isolated 6×6 complex Hadamard matrix, transcribed from
Wuttig–Tindall (arXiv:2608.18053), Definition 5, Eq. (10), with
ω = exp(2πi/3):

    [1  1   1   1   1   1 ]
    [1  1   ω   ω   ω²  ω²]
    [1  ω   1   ω²  ω²  ω ]
    [1  ω   ω²  1   ω   ω²]
    [1  ω²  ω²  ω   1   ω ]
    [1  ω²  ω   ω²  ω   1 ]

The first call validates the transcription (HH† = 6I etc.) and throws if
the transcription is wrong — the validation IS the transcription check.
"""
function tao_S6()
    ω = exp(2π * im / 3)
    w2 = ω^2
    S6 = ComplexF64.([
        1 1 1 1 1 1
        1 1 ω ω w2 w2
        1 ω 1 w2 w2 ω
        1 ω w2 1 ω w2
        1 w2 w2 ω 1 ω
        1 w2 ω w2 ω 1
    ])
    val = validate_hadamard_matrix(S6)
    val.pass || error("Tao S6 transcription FAILED validation: $(val)")
    return S6
end

const DITA_THETA = acos(1 / sqrt(3))
const DITA_PHI = pi / 4

# Expected values (literature pins — do not modify to make tests pass):
const EXPECTED_F6_PHYSICAL = 48                     # known MU-vector count for {I, F6}
const EXPECTED_D0_PHYSICAL = 120                    # BW 2009: N_v = 120
const EXPECTED_D0_THIRD_BASES = 10                  # BW 2009: N_t = 10

"""
    run_benchmark_F6(; tol=1e-8, hp_bits=0, reseed=true)

Benchmark A. Returns (record, pool, audit, cliques).
"""
function run_benchmark_F6(; tol::Real = 1e-8, hp_bits::Int = 0, reseed::Bool = true)
    val = validate_hadamard_matrix(F6)
    val.pass || error("F6 failed intrinsic CHM validation: $(val)")
    out = solve_pool_audited(F6; tol = tol, hp_verify_bits = hp_bits,
                             reseed_crosscheck = reseed)
    cl = enumerate_all_third_mub_bases(out.pool, F6; ortho_tol = tol, mu_tol = tol,
                                       hp_bits = 0)
    obs = out.audit.deduplicated_vectors
    record = (
        benchmark = "A_fourier",
        target = "F6",
        expected_physical_vectors = EXPECTED_F6_PHYSICAL,
        observed_physical_vectors = obs,
        pass_count = obs == EXPECTED_F6_PHYSICAL,
        n_third_bases = cl.n_verified,
        audit = out.audit,
        crosscheck = out.crosscheck,
        hp_verified = out.hp_verified,
        tier = obs == EXPECTED_F6_PHYSICAL ? :CERTIFIED_NUMERICAL : :OPEN,
    )
    return (record = record, pool = out.pool, cliques = cl)
end

"""
    run_benchmark_D0(; lambda=π/2, tol=1e-8, hp_bits=256, reseed=true,
                       certify_w1=false)

Benchmark B at the Dita-slice point H(θ_D, π/4, λ) that is CHM-equivalent
to D0 for λ ∈ {π/2, 3π/2} (E0, numerical residual ~1e-13). Expected:
120 physical vectors, 10 third bases, no two third bases pairwise MU
(N_p = 0 ⇒ no fourth MUB). `certify_w1=true` additionally runs the
one-vector witness on every enumerated third basis.
"""
function run_benchmark_D0(; lambda::Real = pi / 2, tol::Real = 1e-8,
                          hp_bits::Int = 256, reseed::Bool = true,
                          certify_w1::Bool = false, w1_seed::Integer = 20260914)
    H = build_karlsson_family(DITA_THETA, DITA_PHI, lambda)
    val = validate_karlsson_matrix(H; variant = :karlsson_original,
                                   theta = DITA_THETA, phi = DITA_PHI, lambda = lambda,
                                   hp_bits = 0)
    assert_valid_karlsson!(val)

    out = solve_pool_audited(H; tol = tol, hp_verify_bits = hp_bits,
                             reseed_crosscheck = reseed)
    cl = enumerate_all_third_mub_bases(out.pool, H; ortho_tol = tol, mu_tol = tol,
                                       hp_bits = hp_bits)

    obs_v = out.audit.deduplicated_vectors
    obs_b = cl.n_verified
    count_pass = (obs_v == EXPECTED_D0_PHYSICAL && obs_b == EXPECTED_D0_THIRD_BASES)

    w1 = missing
    if certify_w1 && count_pass
        w1 = [w1_witness_certificate(H, hcat([out.pool[i] for i in b.pool_indices]...);
                                     seed = w1_seed + bi)
              for (bi, b) in enumerate(cl.bases)]
    end

    diagnosis = _diagnose_benchmark_b(out, cl, obs_v, obs_b)

    record = (
        benchmark = "B_dita_D0",
        target = "H(θ_D, π/4, λ=$(lambda))",
        chm_equivalent_to = "D0 (BW 2009) at λ ∈ {π/2, 3π/2}",
        expected_physical_vectors = EXPECTED_D0_PHYSICAL,
        observed_physical_vectors = obs_v,
        expected_third_bases = EXPECTED_D0_THIRD_BASES,
        observed_third_bases = obs_b,
        pass_count = count_pass,
        pairwise_mu_among_third_bases = _pairwise_mu_among_bases(out.pool, cl; tol = tol),
        audit = out.audit,
        crosscheck = out.crosscheck,
        hp_verified = out.hp_verified,
        diagnosis = diagnosis,
        w1 = w1,
        tier = count_pass ? :CERTIFIED_NUMERICAL : :OPEN,
    )
    return (record = record, pool = out.pool, cliques = cl, validation = val)
end

function _pairwise_mu_among_bases(pool, cl; tol::Real = 1e-8)
    # BW's N_p = 0: no two vectors from distinct third bases are MU to each
    # other, equivalently no pair of third bases is MU. Test all basis pairs.
    n_pairs_mu = 0
    worst = 0.0
    nb = length(cl.bases)
    for a in 1:nb, b in (a + 1):nb
        ia, ib = cl.bases[a].pool_indices, cl.bases[b].pool_indices
        all_mu = all(minimum(abs(abs2(dot(pool[i], pool[j])) - 1 / 6) for j in ib) < tol
                     for i in ia)
        pair_worst = maximum(minimum(abs(abs2(dot(pool[i], pool[j])) - 1 / 6)
                                     for j in ib) for i in ia)
        worst = max(worst, pair_worst)
        all_mu && (n_pairs_mu += 1)
    end
    return (n_mu_pairs = n_pairs_mu, n_basis_pairs = nb * (nb - 1) ÷ 2, worst_defect = worst)
end

function _diagnose_benchmark_b(out, cl, obs_v, obs_b)
    steps = String[]
    a = out.audit
    if obs_v != EXPECTED_D0_PHYSICAL
        push!(steps, "pool count mismatch: $(obs_v) ≠ 120")
        a.failed_paths > 0 && push!(steps,
            "failed_paths=$(a.failed_paths): tracking lost paths at this (degenerate) point")
        out.crosscheck !== missing &&
            push!(steps, "reseed crosscheck: run2 found $(out.crosscheck.n_solutions_run2) " *
                 "solutions / $(out.crosscheck.n_physical_run2) physical (agree=$(out.crosscheck.agree))")
    end
    if obs_b != EXPECTED_D0_THIRD_BASES
        push!(steps, "third-basis count mismatch: $(obs_b) ≠ 10")
        push!(steps, "clique histogram: $(cl.clique_size_histogram)")
    end
    isempty(steps) && push!(steps, "counts match BW 2009 (N_v=120, N_t=10)")
    return steps
end

"""
    run_benchmark_Tao(; tol=1e-8, hp_bits=0, reseed=true)

Benchmark C: the isolated Tao matrix S6. Expectation per the research
plan: the recovered pool does NOT contain a third MUB. If a 6-clique IS
found, that is a FINDING — it must be reported as-is and escalated to
high-precision verification, never suppressed.
"""
function run_benchmark_Tao(; tol::Real = 1e-8, hp_bits::Int = 0, reseed::Bool = true)
    S6 = tao_S6()
    out = solve_pool_audited(S6; tol = tol, hp_verify_bits = hp_bits,
                             reseed_crosscheck = reseed)
    cl = enumerate_all_third_mub_bases(out.pool, S6; ortho_tol = tol, mu_tol = tol,
                                       hp_bits = 0)
    record = (
        benchmark = "C_tao_isolated",
        target = "S6 (Tao, transcription via Wuttig–Tindall 2026 Eq. 10)",
        expected_third_mub = false,
        observed_third_mub = cl.n_verified > 0,
        observed_physical_vectors = out.audit.deduplicated_vectors,
        n_third_bases = cl.n_verified,
        audit = out.audit,
        crosscheck = out.crosscheck,
        hp_verified = out.hp_verified,
        tier = cl.n_verified == 0 ? :CERTIFIED_NUMERICAL :
               (cl.n_verified > 0 ? :NUMERICALLY_SUPPORTED : :OPEN),
        note = cl.n_verified > 0 ?
            "FINDING: third MUB present in the Tao pool — escalate to HP verification " *
            "and W1 certification; do not suppress." : "",
    )
    return (record = record, pool = out.pool, cliques = cl)
end
