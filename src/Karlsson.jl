# ============================================================================
# Karlsson.jl — Parts II, III, IV of the methodology hardening plan.
#
#   * Two independent transcriptions of the Karlsson A-block
#     (`karlsson_original` vs `review_transcription_literal`), with the
#     regression verdict fixed by audit (C1: original passes, literal
#     review transcription fails).
#   * Tracked Möbius map: all four identities of the construction verified
#     independently, singular branches recorded, never silently passed.
#   * `validate_karlsson_matrix`: single authoritative gate returning a
#     structured report. No pool solver may run on an invalid matrix.
#
# Provenance of the two transcriptions:
#   karlsson_original: B. Karlsson, "Three-parameter complex Hadamard
#     matrices of order 6", LAA 434 (2011), arXiv:1003.4177. A-block
#     A = [[A11, A12], [conj(A12), -conj(A11)]] with
#     A11 = -1/2 + (i√3/2)(cosθ + e^{-iφ} sinθ),
#     A12 = -1/2 + (i√3/2)(-cosθ + e^{+iφ} sinθ).
#   review_transcription_literal: McNulty–Weigert review
#     (arXiv:2410.23997) Eqs. (7.3)–(7.5) as printed:
#     A11 = -1/2 + (i√3/2)(cosθ + e^{-iφ} sinθ), A12 unchanged,
#     A = [[A11, A12], [A12, -A11]]  (conjugates dropped).
#     This variant FAILS AA† = 2I on a full (θ,φ) grid — kept only as a
#     regression witness for the transcription discrepancy. Do not use it
#     for any MUB computation.
# ============================================================================

const F2_HAD = [1 1; 1 -1]

"""Karlsson's original A-block (the audited, load-bearing transcription)."""
function build_A_karlsson_original(theta::Real, phi::Real)
    A11 = -0.5 + im * (sqrt(3) / 2) * (cos(theta) + exp(-im * phi) * sin(theta))
    A12 = -0.5 + im * (sqrt(3) / 2) * (-cos(theta) + exp(im * phi) * sin(theta))
    return [A11 A12; conj(A12) -conj(A11)]
end

"""McNulty–Weigert review transcription, exactly as printed (Eqs. 7.3–7.5).

Known to violate the unitary-block identity; kept for the regression test
that pins down the transcription discrepancy. NEVER use for MUB work."""
function build_A_review_transcription_literal(theta::Real, phi::Real)
    A11 = -0.5 + im * (sqrt(3) / 2) * (cos(theta) + exp(-im * phi) * sin(theta))
    A12 = -0.5 + im * (sqrt(3) / 2) * (-cos(theta) + exp(im * phi) * sin(theta))
    return [A11 A12; A12 -A11]
end

const KARLSSON_VARIANTS = (:karlsson_original, :review_transcription_literal)

function build_A_variant(variant::Symbol, theta::Real, phi::Real)
    variant === :karlsson_original && return build_A_karlsson_original(theta, phi)
    variant === :review_transcription_literal &&
        return build_A_review_transcription_literal(theta, phi)
    error("Unknown Karlsson variant: $variant (use one of $KARLSSON_VARIANTS)")
end

"""
    mobius_tracked(z, alpha, beta; tol=1e-12)

Möbius map M(z) = (αz − β)/(conj(β)z − conj(α)) with the singular branch
recorded instead of silently passed through. `singular=true` when the
denominator magnitude drops below `tol`; the value is then NaN and the
caller must record the event in provenance.
"""
function mobius_tracked(z::Number, alpha::Number, beta::Number; tol::Real = 1e-12)
    den = conj(beta) * z - conj(alpha)
    den_abs = abs(den)
    singular = den_abs < tol
    value = singular ? ComplexF64(NaN, NaN) : ComplexF64((alpha * z - beta) / den)
    return (value = value, singular = singular, den_abs = den_abs)
end

"""
    mobius_identity_residual(z_target, z_source, alpha, beta)

Residual of the identity M(z_source) = z_target evaluated in
DENOMINATOR-CLEARED (polynomial) form:

    | (alpha·z_source − beta) − z_target·(conj(beta)·z_source − conj(alpha)) |

Unlike the divided form, this is well-conditioned at Möbius poles and is
the primary audit quantity. Also reports the divided-form residual (NaN
when the map is singular at z_source) and the indeterminacy flag
(num ≈ 0 AND den ≈ 0 simultaneously — the 0/0 case; at the Dita anchor
|A11|=|A12|=1 makes M_A(z2²) exactly this).
"""
function mobius_identity_residual(z_target::Number, z_source::Number,
                                  alpha::Number, beta::Number; tol::Real = 1e-10)
    num = alpha * z_source - beta
    den = conj(beta) * z_source - conj(alpha)
    cleared = abs(num - z_target * den)
    indeterminate = abs(num) < tol && abs(den) < tol
    singular = abs(den) < tol && !indeterminate
    divided = (singular || indeterminate || !isfinite(z_source)) ?
              ComplexF64(NaN, NaN) : num / den
    return (cleared = cleared, divided = isnan(divided) ? NaN : abs(divided - z_target),
            singular = singular, indeterminate = indeterminate,
            num_abs = abs(num), den_abs = abs(den))
end

"""
    karlsson_mobius_audit(theta, phi, lambda_; variant, tol)

Verify the four Möbius identities of the Karlsson construction
(= Wuttig–Tindall 2026, Supp. Eq. (S.1.128)):

    z3² = M_A(z1²) = M_B(z2²),    z4² = M_B(z1²) = M_A(z2²),

independently of how `build_karlsson_family` internally derives z2, z3, z4.
All residuals are computed in denominator-cleared form (primary,
well-conditioned at poles) with the divided form reported alongside.

Audit finding (2026-09-14): at the Dita anchor the identity z4² = M_A(z2²)
is an indeterminate 0/0 (|A11|=|A12|=1 there), which is why the legacy
`z4_dev` diagnostic jumps from 1.3e-5 (rounded θ, degeneracy scan) to 1.84
(full precision) — it is not a family inconsistency. Singular and
indeterminate branches are flagged, never silently passed.
"""
function karlsson_mobius_audit(theta::Real, phi::Real, lambda_::Real;
                               variant::Symbol = :karlsson_original,
                               tol::Real = 1e-12)
    A = build_A_variant(variant, theta, phi)
    B = -F2_HAD - A
    alpha_A, beta_A = A[1, 2]^2, A[1, 1]^2
    alpha_B, beta_B = B[1, 2]^2, B[1, 1]^2

    z1sq = exp(2im * lambda_)
    seam = (theta == 0.0)

    if seam
        # Fourier seam (see build_karlsson_family): M_A, M_B are the constant
        # map 1 and z2² = alpha_A/beta_A is the exact algebraic resolution.
        z3sq = one(ComplexF64)
        z4sq = one(ComplexF64)
        z2sq_derived = alpha_A / beta_A
        z2_derived_singular = false
        seam_resolved = true
    else
        m_A1 = mobius_tracked(z1sq, alpha_A, beta_A; tol = tol)   # z3² = M_A(z1²)
        m_B1 = mobius_tracked(z1sq, alpha_B, beta_B; tol = tol)   # z4² = M_B(z1²)
        z3sq = m_A1.value
        z4sq = m_B1.value

        # Independent derivation of z2² from the identity z3² = M_B(z2²):
        #   z2² = (β_B − z3² conj(α_B)) / (α_B − z3² conj(β_B))
        num2 = beta_B - z3sq * conj(alpha_B)
        den2 = alpha_B - z3sq * conj(beta_B)
        z2_derived_singular = !isfinite(num2 / den2) || abs(den2) < tol
        z2sq_derived = z2_derived_singular ? ComplexF64(NaN, NaN) : num2 / den2
        seam_resolved = false
    end

    # All four identities, cleared form (primary) + divided form (diagnostic)
    r31 = mobius_identity_residual(z3sq, z1sq, alpha_A, beta_A; tol = 1e-10)
    r41 = mobius_identity_residual(z4sq, z1sq, alpha_B, beta_B; tol = 1e-10)
    r32 = isnan(z2sq_derived) ? missing :
          mobius_identity_residual(z3sq, z2sq_derived, alpha_B, beta_B; tol = 1e-10)
    r42 = isnan(z2sq_derived) ? missing :
          mobius_identity_residual(z4sq, z2sq_derived, alpha_A, beta_A; tol = 1e-10)

    finite_sq = filter(isfinite, (z2sq_derived, z3sq, z4sq))
    z_mod_res = isempty(finite_sq) ? NaN :
                maximum(abs(abs(z2) - 1) for z2 in finite_sq)

    return (
        variant = variant,
        theta = theta, phi = phi, lambda = lambda_,
        z2sq = z2sq_derived, z3sq = z3sq, z4sq = z4sq,
        id_z3_MA_z1 = r31,      # z3² = M_A(z1²)
        id_z4_MB_z1 = r41,      # z4² = M_B(z1²)
        id_z3_MB_z2 = r32,      # z3² = M_B(z2²)
        id_z4_MA_z2 = r42,      # z4² = M_A(z2²)  (indeterminate at the Dita anchor)
        # Convenience accessors (cleared-form residuals; NaN/missing when undefined):
        cleared_z3_MA_z1 = r31.cleared,
        cleared_z4_MB_z1 = r41.cleared,
        cleared_z3_MB_z2 = r32 === missing ? NaN : r32.cleared,
        cleared_z4_MA_z2 = r42 === missing ? NaN : r42.cleared,
        singular_MA_z1 = r31.singular || !isfinite(z3sq),
        singular_MB_z1 = r41.singular || !isfinite(z4sq),
        indeterminate_MA_z2 = r42 !== missing && r42.indeterminate,
        indeterminate_MB_z2 = r32 !== missing && r32.indeterminate,
        singular_MA_z2 = r42 !== missing && r42.singular,
        seam_resolved = seam_resolved,
        z2_derived_singular = z2_derived_singular,
        z_moduli_residual = z_mod_res,
        # Legacy-compat divided-form z4_dev (NaN = undefined at a pole/indeterminate):
        legacy_z4_dev = r42 === missing ? NaN :
                        (r42.singular || r42.indeterminate ? NaN : r42.divided),
    )
end

"""
    karlsson_assemble_diagnostics(theta, phi, lambda_; variant)

Assemble the full 6×6 Karlsson matrix from a given A-variant WITHOUT the
hard unitarity error of `build_karlsson_family`, so that invalid variants
can be audited instead of crashing. Returns (H, diagnostics NamedTuple).
"""
function karlsson_assemble_diagnostics(theta::Real, phi::Real, lambda_::Real;
                                       variant::Symbol = :karlsson_original)
    A = build_A_variant(variant, theta, phi)
    B = -F2_HAD - A
    a_err = maximum(abs.(A * A' - 2 * I(2)))
    b_err = maximum(abs.(B * B' - 2 * I(2)))

    ma = karlsson_mobius_audit(theta, phi, lambda_; variant = variant)
    z2sq, z3sq, z4sq = ma.z2sq, ma.z3sq, ma.z4sq

    Zleft(z) = [1 1; z -z]
    Zright(z) = [1 z; 1 -z]
    z1 = exp(im * lambda_)
    Z1 = Zleft(z1)
    Z2 = Zleft(sqrt(z2sq))
    Z3 = Zright(sqrt(z3sq))
    Z4 = Zright(sqrt(z4sq))

    top = hcat(F2_HAD, Z1, Z2)
    mid = hcat(Z3, 0.5 * Z3 * A * Z1, 0.5 * Z3 * B * Z2)
    bot = hcat(Z4, 0.5 * Z4 * B * Z1, 0.5 * Z4 * A * Z2)
    H = vcat(top, mid, bot)

    u_err = maximum(abs.(H * H' - 6 * I(6)))
    m_err = maximum(abs.(abs.(H) .- 1))
    diag = (
        variant = variant,
        theta = theta, phi = phi, lambda = lambda_,
        A_unitary_err = a_err, B_unitary_err = b_err,
        H_unitary_err = u_err, H_unimodular_err = m_err,
        mobius = ma,
        degenerate_family = ma.singular_MA_z1 || ma.singular_MB_z1 ||
                            ma.singular_MA_z2 || ma.z2_derived_singular ||
                            ma.indeterminate_MA_z2 || ma.indeterminate_MB_z2,
    )
    return H, diag
end

# ---------------------------------------------------------------------------
# Part IV — the single authoritative validation gate.
# ---------------------------------------------------------------------------

"""
    validate_hadamard_matrix(H; tol_unimodular=1e-10, tol_unitary=1e-10)

General CHM gate for an arbitrary 6×6 matrix (used for Tao's S₆ and for
cross-checking Karlsson output). Returns a structured report; `pass`
requires BOTH HH† = 6I and H†H = 6I and unimodular entries.
"""
function validate_hadamard_matrix(H::AbstractMatrix{ComplexF64};
                                  tol_unimodular::Real = 1e-10,
                                  tol_unitary::Real = 1e-10)
    n = size(H, 1)
    n == 6 || error("expect a 6×6 matrix")
    hh = maximum(abs.(H * H' - n * I(6)))
    hth = maximum(abs.(H' * H - n * I(6)))
    mod_err = maximum(abs.(abs.(H) .- 1))
    pass = hh < tol_unitary && hth < tol_unitary && mod_err < tol_unimodular
    return (
        n = n,
        hh_norm = hh,               # max |HH† − 6I|
        hth_norm = hth,             # max |H†H − 6I|
        unimodularity_residual = mod_err,
        tol_unimodular = tol_unimodular,
        tol_unitary = tol_unitary,
        pass = pass,
        classification = pass ? "VALID_CHM" : "INVALID_CHM",
    )
end

"""
    validate_karlsson_matrix(H; variant, theta, phi, lambda, <tols>)

The single authoritative gate (Part IV). Checks the assembled matrix H
against its claimed parameters and variant, at Float64 and (optionally)
BigFloat precision. Returns a structured report; `pass=true` is REQUIRED
before any MUB pool computation on H.

When `theta/phi/lambda` are omitted, only the intrinsic CHM identities are
checked (parameter traceability is then `missing`).
"""
function validate_karlsson_matrix(H::AbstractMatrix{ComplexF64};
                                  variant::Union{Nothing,Symbol} = nothing,
                                  theta::Union{Nothing,Real} = nothing,
                                  phi::Union{Nothing,Real} = nothing,
                                  lambda::Union{Nothing,Real} = nothing,
                                  tol_unimodular::Real = 1e-10,
                                  tol_unitary::Real = 1e-10,
                                  tol_mobius::Real = 1e-8,
                                  hp_bits::Int = 0)
    base = validate_hadamard_matrix(H; tol_unimodular = tol_unimodular,
                                    tol_unitary = tol_unitary)

    param_residuals = missing
    mobius = missing
    degenerate_family = false
    if variant !== nothing && theta !== nothing && phi !== nothing && lambda !== nothing
        H2, diag = karlsson_assemble_diagnostics(theta, phi, lambda; variant = variant)
        assembly_err = maximum(abs.(H2 - H))
        mobius = diag.mobius
        degenerate_family = diag.degenerate_family
        param_residuals = (
            A_unitary_err = diag.A_unitary_err,
            B_unitary_err = diag.B_unitary_err,
            H_assembly_err = assembly_err,
            mobius_cleared_max = max(isnan(mobius.cleared_z3_MA_z1) ? 0.0 : mobius.cleared_z3_MA_z1,
                                     isnan(mobius.cleared_z4_MB_z1) ? 0.0 : mobius.cleared_z4_MB_z1,
                                     isnan(mobius.cleared_z3_MB_z2) ? 0.0 : mobius.cleared_z3_MB_z2,
                                     isnan(mobius.cleared_z4_MA_z2) ? 0.0 : mobius.cleared_z4_MA_z2),
            z_moduli_residual = mobius.z_moduli_residual,
            legacy_z4_dev = mobius.legacy_z4_dev,
        )
        # Möbius identities in cleared form are the consistency check; at
        # poles/indeterminacies the residuals are NaN (vacuous) and the
        # matrix is judged by its direct CHM identities alone. Branch
        # degeneracy (Fourier seam θ=0, Dita anchor 0/0) is recorded as a
        # STATUS, never silently passed, and never fails a valid matrix.
        mobius_ok = (isnan(mobius.cleared_z3_MA_z1) ? true : mobius.cleared_z3_MA_z1 < tol_mobius) &&
                    (isnan(mobius.cleared_z4_MB_z1) ? true : mobius.cleared_z4_MB_z1 < tol_mobius) &&
                    (isnan(mobius.cleared_z3_MB_z2) ? true : mobius.cleared_z3_MB_z2 < tol_mobius) &&
                    (isnan(mobius.cleared_z4_MA_z2) ? true : mobius.cleared_z4_MA_z2 < tol_mobius)
    else
        mobius_ok = true
        assembly_err = missing
    end

    hp = missing
    if hp_bits > 0
        setprecision(hp_bits) do
            Hb = [Complex{BigFloat}(x) for x in H]
            n = BigFloat(6)
            hh = maximum(abs.(Hb * Hb' - n * I(6)))
            hth = maximum(abs.(Hb' * Hb - n * I(6)))
            mod_err = maximum(abs.(abs.(Hb) .- 1))
            hp = (bits = hp_bits,
                  hh_norm = Float64(hh), hth_norm = Float64(hth),
                  unimodularity_residual = Float64(mod_err))
        end
    end

    core_pass = base.pass && mobius_ok
    return (
        variant = variant,
        theta = theta, phi = phi, lambda = lambda,
        hh_norm = base.hh_norm,
        hth_norm = base.hth_norm,
        unimodularity_residual = base.unimodularity_residual,
        param_residuals = param_residuals,
        mobius_audit = mobius,
        degenerate_family = degenerate_family,
        high_precision = hp,
        tol_unimodular = tol_unimodular,
        tol_unitary = tol_unitary,
        tol_mobius = tol_mobius,
        pass = core_pass,
        classification = core_pass ? (degenerate_family ? "VALID_CHM_DEGENERATE_BRANCH"
                                                       : "VALID_CHM")
                                   : "INVALID_CHM",
    )
end

"""Gate helper: refuse to run pool solvers on an unvalidated matrix."""
function assert_valid_karlsson!(val)
    val.pass || error("Karlsson validation FAILED (classification=$(val.classification)); " *
                      "refusing to run MUB pool solvers. Report: $(val)")
    return val
end
