# ============================================================================
# MubSearch.jl — the audited, module-scoped entry point for the MUB search
# pipeline (Part XVI restructuring).
#
# Layering:
#   1. mub_zauner_6d_liang_chen.jl — legacy core (single source of truth for
#      the pool system, per-H solving, orthogonality graph, HP verifiers).
#      Kept byte-identical so all existing scripts keep working when they
#      include the file directly.
#   2. Karlsson.jl      — Parts II–IV: transcription variants, tracked
#                         Möbius map, authoritative validation gate.
#   3. Pool.jl          — Parts V–VII: PoolAudit bookkeeping, projective
#                         dedup, tolerance diagnostics.
#   4. Cliques.jl       — Part IX: exhaustive third-MUB enumeration with
#                         independent per-basis verification.
#   5. Certification.jl — Parts X–XI: one-vector witness W1, claim tiers.
#   6. Benchmarks.jl    — Part VIII: Fourier / D0 / Tao regression gates.
#   7. Provenance.jl    — Part XV: metadata + result writer.
#
# Usage:
#   include("src/MubSearch.jl"); using .MubSearch
# ============================================================================
module MubSearch

using HomotopyContinuation
using LinearAlgebra
using Graphs
using Random
using Printf
using Dates

include("mub_zauner_6d_liang_chen.jl")   # legacy core
include("Karlsson.jl")
include("Pool.jl")
include("Cliques.jl")
include("Certification.jl")
include("Benchmarks.jl")
include("Provenance.jl")

# — Karlsson construction & validation (Parts II–IV) —
export KARLSSON_VARIANTS, F2_HAD,
       build_A_karlsson_original, build_A_review_transcription_literal,
       build_A_variant, mobius_tracked, mobius_identity_residual,
       karlsson_mobius_audit, karlsson_assemble_diagnostics,
       validate_hadamard_matrix, validate_karlsson_matrix, assert_valid_karlsson!

# — Pool audit & dedup (Parts V–VII) —
export PoolAudit, audit_explanations, projective_distance,
       cluster_pool_hp, dedup_tolerance_diagnostics,
       solve_pool_audited, verify_candidate_hp

# — Cliques (Part IX) —
export enumerate_all_third_mub_bases

# — Certification (Parts X–XI) —
export CLAIM_TIERS, tier_rank, at_least, w1_witness_certificate

# — Benchmarks (Part VIII) —
export tao_S6, run_benchmark_F6, run_benchmark_D0, run_benchmark_Tao,
       DITA_THETA, DITA_PHI, EXPECTED_F6_PHYSICAL, EXPECTED_D0_PHYSICAL,
       EXPECTED_D0_THIRD_BASES

# — Provenance (Part XV) —
export provenance_record, write_result

# — Re-exports from the legacy core used by new scripts —
export F6, build_karlsson_family, build_numeric_pool_system,
       check_four_mub_extension, verify_clique_hp, _orthogonality_graph

end # module
