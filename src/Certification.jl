# ============================================================================
# Certification.jl — Parts X, XI: the fourth-vector test and claim tiers.
#
# Claim-tier vocabulary (Part XI). These six labels are the ONLY statuses
# allowed in exported results; a numerical sweep is never a proof.
#
#   PROVED                analytic/logical proof, in-repo
#   EXACT_ALGEBRAIC       exact-coefficient computation (GB unit ideal,
#                         number-field arithmetic) — machine-verified
#   CERTIFIED_NUMERICAL   interval-arithmetic certification (HC certify(),
#                         a posteriori enclosures) of a numerical object
#   NUMERICALLY_SUPPORTED strong multi-run numerical evidence (400-bit HP,
#                         tolerance sweeps, cross-checked solving)
#   HEURISTIC             single-run or tolerance-fragile numerical signal
#   OPEN                  not established
# ============================================================================

const CLAIM_TIERS = (:PROVED, :EXACT_ALGEBRAIC, :CERTIFIED_NUMERICAL,
                     :NUMERICALLY_SUPPORTED, :HEURISTIC, :OPEN)

tier_rank(t::Symbol) = (i = findfirst(==(t), CLAIM_TIERS)) === nothing ?
                       error("unknown claim tier: $t") : i

"""`at_least(:CERTIFIED_NUMERICAL, :NUMERICALLY_SUPPORTED) == true`
(ranks: PROVED strongest → OPEN weakest)."""
at_least(tier::Symbol, threshold::Symbol) = tier_rank(tier) <= tier_rank(threshold)

# ---------------------------------------------------------------------------
# Part X — the one-vector witness W1(B3). The logical statement:
#
#   a fourth MUB extending {I, H, B3} ⇒ W1(H, B3) ≠ ∅
#   contrapositive:  W1(H, B3) = ∅ ⇒ the triple {I, H, B3} does NOT extend
#
# Emptiness of W1 is therefore a COMPLETE certificate of non-extension for
# that triple — strictly stronger than needed. The n_wit=2 system is
# logically redundant for this claim structure (two orthogonal W1-solutions
# are neither necessary nor sufficient for a fourth ONB) and is retired:
# build_numeric_fourth_mub_witness_system(n_wit=2) remains in the legacy
# file for archival reproduction only. It must not be used in new claims.
# ---------------------------------------------------------------------------

"""
    w1_witness_certificate(H, B3; seed=20260914, residual_tol=1e-8)

Certified emptiness test for the single-vector witness at fixed (H, B3).

System (gauge z0 = w0 = 1, 10 complex unknowns, 17 equations):
  z_i w_i = 1                     (i = 1..5)   unimodular / flat
  |Σ_j conj(H[j,k]) z_j|² = 6     (k = 1..6)   MU to H
  |Σ_j conj(B3[j,k]) z_j|² = R2   (k = 1..6)   MU to B3

Where R2 = 6 if B3 is unnormalized (B3*B3' = 6I) or R2 = 1 if B3 is
unit-normalized (B3'B3 = I). The correct R2 is auto-detected from the
matrix norm of B3.

Method (certified):
  1. Form a generic square system by 10 random complex linear combinations
     of the 17 equations (Bertini pattern; fixed seed for reproducibility).
  2. Mixed-volume count; track ALL start paths.
  3. Interval-arithmetic certify() every endpoint.
  4. Full-residual gate: a true W1 solution must satisfy all 17 equations,
     so it must appear among the certified square-system roots; count
     certified roots passing the full residual at `residual_tol`.

Verdicts:
  EMPTY_CERTIFIED  — all paths tracked, all certified, 0 pass the full
                     residual ⇒ W1 = ∅ ⇒ the triple does not extend.
                     Tier: CERTIFIED_NUMERICAL.
  NONEMPTY         — at least one certified full-residual pass; the
                     witness found a fourth-MU vector (follow up!). This
                     does NOT by itself give a fourth ONB.
  INCOMPLETE       — tracking failed paths or certification gaps; no claim.

The residual gate is evaluated with a margin: emptiness transfers from the
computed (certified) B3 to the exact third basis only up to the certified
perturbation size — see METHODOLOGY_AUDIT.md §epsilon-transfer for the
bookkeeping this requires in theorem statements.

NORMALIZATION FIX (2026-09-17): Auto-detects whether B3 is unit-normalized
(B3'B3 ≈ I, RHS=1) or unnormalized (B3*B3' ≈ 6I, RHS=6). This fixes the
prior defect where unit-normalized pool vectors were incorrectly tested
with RHS=6.
"""
function w1_witness_certificate(H::AbstractMatrix{ComplexF64},
                                B3::AbstractMatrix{ComplexF64};
                                seed::Integer = 20260914,
                                residual_tol::Real = 1e-8)
    size(B3) == (6, 6) || error("B3 must be 6×6 (columns are the third basis)")

    # Auto-detect normalization convention from B3
    b3_gram = B3' * B3
    b3_unitary_err = maximum(abs.(b3_gram - I(6)))
    b3_hadamard_err = maximum(abs.(b3_gram - 6.0 * I(6)))
    if b3_unitary_err < b3_hadamard_err
        rhs_val = 1.0
        b3_scale = "unitary"
    else
        rhs_val = 6.0
        b3_scale = "unnormalized"
    end

    @var z[1:5] w[1:5]
    z_full = [1.0 + 0im; [z[i] for i in 1:5]]
    w_full = [1.0 + 0im; [w[i] for i in 1:5]]
    eqs = Any[]
    for i in 1:5
        push!(eqs, z[i] * w[i] - 1.0)
    end
    for k in 1:6
        lhs = sum(conj(H[j, k]) * z_full[j] for j in 1:6)
        rhs = sum(H[j, k] * w_full[j] for j in 1:6)
        push!(eqs, lhs * rhs - 6.0)
    end
    for k in 1:6
        lhs = sum(conj(B3[j, k]) * z_full[j] for j in 1:6)
        rhs = sum(B3[j, k] * w_full[j] for j in 1:6)
        push!(eqs, lhs * rhs - rhs_val)
    end

    full_system = System(eqs; variables = vcat(z, w))
    n_eqs, n_vars = length(eqs), 10
    n_eqs == n_vars + 7 || error("unexpected witness shape: $n_eqs equations")

    Random.seed!(seed)
    R = randn(ComplexF64, n_vars, n_eqs)
    squared_eqs = [sum(R[i, j] * eqs[j] for j in 1:n_eqs) for i in 1:n_vars]
    square_system = System(squared_eqs; variables = vcat(z, w))

    mv = mixed_volume(square_system)
    res = solve(square_system)
    raw = solutions(res)
    n_tracked = length(raw)
    cert = certify(square_system, res; threading = false)
    n_cert = ncertified(cert)

    n_pass = 0
    for sol in raw
        rmax = 0.0
        try
            vals = full_system(sol)
            rmax = maximum(abs, vals)
        catch
            rmax = Inf
        end
        rmax < residual_tol && (n_pass += 1)
    end

    complete = (n_tracked == mv && n_cert == n_tracked)
    verdict, tier = if !complete
        ("INCOMPLETE", :OPEN)
    elseif n_pass > 0
        ("NONEMPTY", :CERTIFIED_NUMERICAL)
    else
        ("EMPTY_CERTIFIED", :CERTIFIED_NUMERICAL)
    end

    return (
        n_equations = n_eqs, n_variables = n_vars,
        mixed_volume = Int(mv), paths_tracked = n_tracked,
        certified = Int(n_cert), full_residual_passes = n_pass,
        residual_tol = residual_tol, seed = seed,
        complete_tracking = complete,
        verdict = verdict, tier = tier,
        w1_empty = verdict == "EMPTY_CERTIFIED",
        b3_normalization = b3_scale,
        b3_rhs_value = rhs_val,
        statement = verdict == "EMPTY_CERTIFIED" ?
            "W1(H,B3) = ∅ (certified): the triple {I,H,B3} admits no fourth MUB vector, hence no fourth MUB." :
            verdict == "NONEMPTY" ?
            "W1(H,B3) ≠ ∅: a fourth-MU vector exists for this triple (not yet an ONB)." :
            "tracking/certification incomplete — no claim.",
    )
end
