# ============================================================================
# run_benchmarks.jl — Part VIII production run.
#
# Runs the three hard regression benchmarks and writes structured results
# with full provenance:
#
#   A. Fourier {I, F6}: expect 48 physical MU vectors.
#   B. Dita/D0 at λ=π/2 (CHM-equivalent to BW's D0): expect 120 vectors,
#      10 third bases, N_p = 0; W1 certified empty on every third basis.
#      Also λ=0.4 (72-vector constellation) for contrast.
#   C. Tao S6 (isolated): expect no third MUB.
#
# Output: results/benchmarks/benchmarks.json (+ .md summary).
# Usage: julia --project=. scripts/julia/run_benchmarks.jl [--quick]
#   --quick skips W1 certification inside benchmark B.
# ============================================================================

include(joinpath(@__DIR__, "_paths.jl"))
include(joinpath(ROOT, "src", "MubSearch.jl"))
using .MubSearch

const OUT_DIR = joinpath(RESULTS_DIR, "benchmarks")
mkpath(OUT_DIR)
const QUICK = "--quick" in ARGS

prov = provenance_record(
    seed = 20260914,
    parameters = Dict("benchmarks" => ["A_fourier", "B_dita_D0(lam=pi/2)",
                                       "B_dita_D0(lam=0.4)", "C_tao"],
                      "dita_theta" => DITA_THETA, "dita_phi" => DITA_PHI),
    solver = Dict("pool" => "HC polyhedral + certify() interval arithmetic",
                  "witness" => "W1 n_wit=1, Bertini square + certify()",
                  "dedup" => "projective (Fubini-Study) clustering @128 bit"),
    tolerances = Dict("solve_verify_tol" => 1e-8, "cluster_tol" => 1e-8,
                      "w1_residual_tol" => 1e-8),
    extra = Dict("runner" => "scripts/julia/run_benchmarks.jl",
                 "quick_mode" => QUICK),
)

results = Dict{String,Any}()

println("=== Benchmark A: Fourier {I, F6} ===")
A = run_benchmark_F6(; hp_bits = 256, reseed = true)
results["A_fourier"] = A.record
println(A.record.audit)
println("  observed=", A.record.observed_physical_vectors,
        " expected=", A.record.expected_physical_vectors,
        " PASS=", A.record.pass_count,
        " crosscheck_agree=", A.record.crosscheck.agree)

println("=== Benchmark B: Dita λ=π/2 (D0-equivalent) ===")
B = run_benchmark_D0(; lambda = pi / 2, hp_bits = 256, reseed = true,
                     certify_w1 = !QUICK)
results["B_dita_D0_pi2"] = B.record
println(B.record.audit)
println("  vectors: observed=", B.record.observed_physical_vectors,
        " expected=", B.record.expected_physical_vectors)
println("  third bases: observed=", B.record.observed_third_bases,
        " expected=", B.record.expected_third_bases)
println("  N_p (MU basis pairs): ", B.record.pairwise_mu_among_third_bases)
println("  diagnosis: ", B.record.diagnosis)
if !QUICK
    for w in B.record.w1
        println("  W1 clique ", w.n_equations > 0 ? "" : "", "...",
                " verdict=", w.verdict, " tier=", w.tier)
    end
end

println("=== Benchmark B2: Dita λ=0.4 (72-vector constellation) ===")
B2 = run_benchmark_D0(; lambda = 0.4, hp_bits = 256, reseed = true,
                      certify_w1 = false)
results["B_dita_lam04"] = B2.record
println(B2.record.audit)
println("  vectors: observed=", B2.record.observed_physical_vectors,
        " third bases: ", B2.record.observed_third_bases)

println("=== Benchmark C: Tao S6 (isolated) ===")
C = run_benchmark_Tao(; hp_bits = 256, reseed = true)
results["C_tao"] = C.record
println(C.record.audit)
println("  third MUB present: ", C.record.observed_third_mub,
        " (expected false)   vectors=", C.record.observed_physical_vectors)

# Dedup tolerance diagnostics on every pool (Part VII deliverable)
println("=== Dedup tolerance diagnostics ===")
for (name, pool) in [("A_fourier", A.pool), ("B_dita_D0_pi2", B.pool),
                     ("B_dita_lam04", B2.pool), ("C_tao", C.pool)]
    diags = dedup_tolerance_diagnostics(pool)
    results["dedup_diagnostics_" * name] = diags
    println("  ", name, ": ", [(d.tol, d.n_clusters) for d in diags])
end

# Verdict block
results["verdicts"] = Dict(
    "A_fourier" => A.record.pass_count ? "PASS" : "FAIL",
    "B_dita_D0_pi2_counts" =>
        B.record.observed_physical_vectors == EXPECTED_D0_PHYSICAL &&
        B.record.observed_third_bases == EXPECTED_D0_THIRD_BASES ? "PASS" : "MISMATCH",
    "B_dita_D0_pi2_w1" => QUICK ? "SKIPPED" :
        all(w.verdict == "EMPTY_CERTIFIED" for w in B.record.w1) ? "ALL_EMPTY" : "SEE_RECORD",
    "C_tao" => !C.record.observed_third_mub ? "PASS" : "FINDING_SEE_RECORD",
)

path = write_result(joinpath(OUT_DIR, "benchmarks.json"), results; provenance = prov)
println("Wrote ", path)

# Human-readable summary
open(joinpath(OUT_DIR, "benchmarks.md"), "w") do io
    println(io, "# Regression benchmarks — run $(prov["timestamp_utc"])")
    println(io)
    println(io, "| Benchmark | Expected | Observed | Verdict |")
    println(io, "|-----------|----------|----------|---------|")
    println(io, "| A: F6 physical MU vectors | 48 | $(A.record.observed_physical_vectors) | $(results["verdicts"]["A_fourier"]) |")
    println(io, "| B: D0 vectors (λ=π/2) | 120 | $(B.record.observed_physical_vectors) | $(results["verdicts"]["B_dita_D0_pi2_counts"]) |")
    println(io, "| B: D0 third bases (λ=π/2) | 10 | $(B.record.observed_third_bases) | $(results["verdicts"]["B_dita_D0_pi2_counts"]) |")
    println(io, "| B2: Dita λ=0.4 vectors | (72 historical) | $(B2.record.observed_physical_vectors) | — |")
    println(io, "| C: Tao S6 third MUB | none | $(C.record.observed_third_mub ? "FOUND" : "none") | $(results["verdicts"]["C_tao"]) |")
    println(io)
    println(io, "Provenance: commit `$(prov["git_commit"])`, Julia $(prov["julia_version"]), seed $(prov["random_seed"]).")
end
println("Wrote ", joinpath(OUT_DIR, "benchmarks.md"))
