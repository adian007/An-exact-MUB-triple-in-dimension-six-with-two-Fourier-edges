include(joinpath(@__DIR__, "_paths.jl"))
include(joinpath(ROOT, "src", "MubSearch.jl"))

using HomotopyContinuation
using LinearAlgebra
using Dates
using Printf
using Random
using CSV
using DataFrames

const LAMBDA_STAR = 0.1114802243779665542913031975274172717685818097497045
const DELTA_LAMBDA = 2e-6
const CAMPAIGN_SEED = 20260928
const OUTPUT_DIR = joinpath(RESULTS_DIR, "campaigns", "i3_singular_locus")
const DITA_Z = exp(im * LAMBDA_STAR)

function build_dita_H(z::ComplexF64)
    return ComplexF64[
        1 1 1 1 1 1
        1 -1 z -z im -im
        1 -im im im -im -1
        1 im -z z -1 -im
        1 inv(z) -im -1 -inv(z) im
        1 -inv(z) -1 -im inv(z) im
    ]
end

function build_symbolic_pool_equations()
    @var x[1:5] y[1:5] z t u[1:10]
    H = [
        1 1 1 1 1 1
        1 -1 z -z im -im
        1 -im im im -im -1
        1 im -z z -1 -im
        1 t -im -1 -t im
        1 -t -1 -im t im
    ]
    Hbar = [
        1 1 1 1 1 1
        1 -1 t -t -im im
        1 im -im -im im -1
        1 -im -t t -1 im
        1 z im -1 -z -im
        1 -z -1 im z -im
    ]

    xfull = [1; x]
    yfull = [1; y]
    eqs_inverse = [x[j] * y[j] - 1 for j in 1:5]
    eqs_mu = [
        sum(Hbar[j, k] * xfull[j] for j in 1:6) *
        sum(H[j, k] * yfull[j] for j in 1:6) - 6 for k in 1:5
    ]
    pool_equations = [eqs_inverse; eqs_mu]
    pool_variables = [x; y]
    jac = differentiate(pool_equations, pool_variables)
    null_equations = [
        sum(jac[row, col] * u[col] for col in 1:10) for row in 1:10
    ]
    augmented_equations = [
        pool_equations;
        null_equations;
        z * t - 1;
        sum(u[j]^2 for j in 1:10) - 1
    ]
    variables = [x; y; z; t; u]
    return System(augmented_equations; variables), variables
end

function json_real_pool(pool)
    [[[real(value), imag(value)] for value in vector] for vector in pool]
end

function physical_pool_at(lambda::Real; verbose::Bool=true)
    H = build_dita_H(ComplexF64(cis(lambda)))
    defect = maximum(abs.(H * H' - 6I(6)))
    defect < 1e-10 || error("H_D column convention failed CHM check: $defect")
    Random.seed!(CAMPAIGN_SEED)
    pool, stats = generate_candidate_pool_fresh(H; verbose=verbose)
    return H, pool, stats
end

function run_transition_counts()
    mkpath(OUTPUT_DIR)
    samples = [
        ("below", LAMBDA_STAR - DELTA_LAMBDA),
        ("critical", LAMBDA_STAR),
        ("above", LAMBDA_STAR + DELTA_LAMBDA),
    ]
    results = Dict{String,Any}()
    for (label, lambda) in samples
        println("=== I3 pool at $label λ=$(repr(lambda)) ===")
        H, pool, stats = physical_pool_at(lambda)
        result = Dict{String,Any}(
            "lambda" => lambda,
            "hadamard_defect" => maximum(abs.(H * H' - 6I(6))),
            "raw_solutions" => stats.n_raw,
            "conjugate_locus_solutions" => stats.n_conjugate,
            "verified_mu_vectors" => stats.n_verified,
            "distinct_physical_vectors" => length(pool),
            "expected_regular_physical_vectors" => label == "below" ? 120 : 72,
            "expected_total_physical_vectors_including_singular" =>
                label == "critical" ? 96 : label == "below" ? 120 : 72,
            "expected_singular_roots" => label == "critical" ? 24 : 0,
            "pool_vectors_re_im" => json_real_pool(pool),
        )
        result["regular_count_matches_expected"] =
            result["distinct_physical_vectors"] ==
            result["expected_regular_physical_vectors"]
        results[label] = result
        @printf("distinct physical vectors: %d (expected %d regular roots)\n",
                result["distinct_physical_vectors"],
                result["expected_regular_physical_vectors"])
    end
    path = joinpath(OUTPUT_DIR, "transition_counts.json")
    MubSearch.write_result(path, results;
                 provenance=MubSearch.provenance_record(
                     seed=CAMPAIGN_SEED,
                     parameters=Dict("lambda_star" => LAMBDA_STAR,
                                     "delta_lambda" => DELTA_LAMBDA,
                                     "hadamard_convention" => "H_D columns; unnormalized"),
                     solver=Dict("method" => "fresh HomotopyContinuation solve"),
                     tolerances=Dict("pool_filter" => 1e-8),
                     extra=Dict("completion_status" =>
                         all(v["regular_count_matches_expected"] for v in values(results)) ?
                         "completed" : "count_mismatch")))
    return results, path
end

function run_augmented_solve(; certify_roots::Bool=false)
    mkpath(OUTPUT_DIR)
    system, variables = build_symbolic_pool_equations()
    println("Augmented I3 singular system: $(length(variables)) variables, " *
            "$(length(system)) equations")
    println("Computing mixed volume...")
    mv = mixed_volume(system)
    println("Augmented-system mixed volume: $mv")
    println("Solving Jacobian-null system...")
    Random.seed!(CAMPAIGN_SEED)
    result = solve(system; show_progress=true)
    raw = solutions(result)
    println("Finite augmented solutions: $(length(raw))")

    records = Vector{Dict{String,Any}}()
    for (solution_index, solution) in enumerate(raw)
        all(isfinite, solution) || continue
        vals = Dict(variables[i] => solution[i] for i in eachindex(variables))
        xvals = ComplexF64[vals[variables[j]] for j in 1:5]
        yvals = ComplexF64[vals[variables[j + 5]] for j in 1:5]
        zval = ComplexF64(vals[variables[11]])
        tval = ComplexF64(vals[variables[12]])
        uvals = ComplexF64[vals[variables[j + 12]] for j in 1:10]
        residual = maximum(abs.(evaluate(system, solution)))
        physical_defect = maximum((
            abs(abs(zval) - 1),
            abs(tval - conj(zval)),
            maximum(abs(abs(value) - 1) for value in xvals),
            maximum(abs(yvals[j] - conj(xvals[j])) for j in 1:5),
            abs(sum(value^2 for value in uvals) - 1),
        ))
        push!(records, Dict{String,Any}(
            "solution_index" => solution_index,
            "z_re_im" => [real(zval), imag(zval)],
            "t_re_im" => [real(tval), imag(tval)],
            "lambda" => angle(zval),
            "x_re_im" => [[real(v), imag(v)] for v in xvals],
            "y_re_im" => [[real(v), imag(v)] for v in yvals],
            "u_re_im" => [[real(v), imag(v)] for v in uvals],
            "system_residual" => residual,
            "physical_defect" => physical_defect,
            "physical_candidate" => physical_defect < 1e-6,
        ))
    end

    certified_count = missing
    uncertified_count = missing
    if certify_roots
        println("Certifying augmented-system solutions...")
        certificates = certify(system, result; show_progress=true, threading=false)
        certified_count = ncertified(certificates)
        all_certificates = HomotopyContinuation.certificates(certificates)
        uncertified_count = length(all_certificates) - certified_count
        for record in records
            index = record["solution_index"]
            record["certified"] = is_certified(all_certificates[index])
        end
    end
    physical_records = filter(record -> record["physical_candidate"], records)
    sort!(physical_records; by=record -> record["system_residual"])
    unique_records = Dict{String,Any}[]
    for record in physical_records
        zvalue = complex(record["z_re_im"]...)
        xvalue = ComplexF64[complex(pair...) for pair in record["x_re_im"]]
        duplicate = any(unique_records) do existing
            abs(complex(existing["z_re_im"]...) - zvalue) < 1e-6 &&
            maximum(abs.(ComplexF64[complex(pair...) for pair in existing["x_re_im"]] .-
                         xvalue)) < 1e-6
        end
        duplicate || push!(unique_records, record)
    end
    physical_records = unique_records
    target_z = DITA_Z
    sort!(physical_records; by=record ->
        abs(complex(record["z_re_im"]...) - target_z))
    for record in physical_records
        record["distance_to_target_z"] =
            abs(complex(record["z_re_im"]...) - target_z)
    end
    output = Dict{String,Any}(
        "mixed_volume" => mv,
        "finite_augmented_solutions" => length(raw),
        "physical_candidates" => length(physical_records),
        "certified_solutions" => certified_count,
        "uncertified_solutions" => uncertified_count,
        "certified_physical_candidates" =>
            certify_roots ? count(record -> get(record, "certified", false),
                                  physical_records) : missing,
        "root_isolation_status" =>
            certify_roots ? "certified augmented-system roots" :
                            "uncertified numerical roots",
        "normalization" => "sum(u_j^2)=1",
        "kernel_variables" => "Jacobian columns are x1..x5,y1..y5",
        "mu_columns" => "H_D columns 1..5; sixth MU equation is Parseval-dependent",
        "physical_roots" => physical_records,
        "all_finite_roots" => records,
    )
    path = joinpath(OUTPUT_DIR, "augmented_singular_roots.json")
    MubSearch.write_result(path, output;
                 provenance=MubSearch.provenance_record(
                     seed=CAMPAIGN_SEED,
                     parameters=Dict("lambda_star_target" => LAMBDA_STAR),
                     solver=Dict("method" => "Jacobian-null augmented polynomial system",
                                 "mixed_volume" => mv),
                     tolerances=Dict("physical_filter" => 1e-6),
                     extra=Dict("completion_status" =>
                         certify_roots ? "certified" : "uncertified_numerical",
                                "target_root_isolation" => "physical candidates are numerically deduplicated")))
    return output, path
end

function run_clique_exports()
    mkpath(OUTPUT_DIR)
    results = Dict{String,Any}()
    for (label, lambda) in [
        ("lambda_0_4", 0.4),
        ("lambda_pi_over_3", pi / 3),
        ("lambda_2pi_over_3", 2pi / 3),
    ]
        println("=== Verified I3 cliques at $label ===")
        H = build_dita_H(ComplexF64(cis(lambda)))
        Random.seed!(CAMPAIGN_SEED)
        pool, stats = generate_candidate_pool_fresh(H; verbose=true)
        enumeration = MubSearch.enumerate_all_third_mub_bases(
            pool, H; ortho_tol=1e-8, mu_tol=1e-8, hp_bits=128
        )
        verified = NamedTuple[]
        for record in enumeration.bases
            (record.is_verified_third_mub && record.hp.ortho_ok && record.hp.mu_ok) ||
                continue
            mu_to_H = [
                MubSearch.verify_mu_to_basis_hp(pool[index], H; bits=128)
                for index in record.pool_indices
            ]
            all(check -> check.ok, mu_to_H) || continue
            push!(verified, (record=record,
                             mu_to_H_max=maximum(check.mu_max for check in mu_to_H)))
        end
        cliques = [
            [[[real(value), imag(value)] for value in pool[index]]
             for index in item.record.pool_indices]
            for item in verified
        ]
        results[label] = Dict{String,Any}(
            "lambda" => lambda,
            "pool_count" => length(pool),
            "raw_solutions" => stats.n_raw,
            "conjugate_locus_solutions" => stats.n_conjugate,
            "verified_mu_vectors" => stats.n_verified,
            "clique_count" => length(verified),
            "max_hp_mu_to_H_defect" =>
                isempty(verified) ? 0.0 : maximum(item.mu_to_H_max for item in verified),
            "verified_cliques_re_im" => cliques,
        )
        @printf("verified pool=%d, verified third bases=%d\n",
                length(pool), length(verified))
    end
    path = joinpath(OUTPUT_DIR, "third_mub_cliques.json")
    MubSearch.write_result(path, results;
                 provenance=MubSearch.provenance_record(
                     seed=CAMPAIGN_SEED,
                     parameters=Dict("lambda_values" => [0.4, pi / 3, 2pi / 3],
                                     "hadamard_convention" => "H_D columns; unnormalized"),
                     solver=Dict("method" => "fresh pool solve and 128-bit clique verification"),
                     tolerances=Dict("pool_filter" => 1e-8,
                                     "clique_ortho" => 1e-8,
                                     "clique_mu" => 1e-8),
                     extra=Dict("completion_status" => "completed",
                                "claim_scope" => "verified cliques in recovered pools")))
    return results, path
end

function run_saved_root_certification()
    mkpath(OUTPUT_DIR)
    input_path = joinpath(OUTPUT_DIR, "roots_to_certify.csv")
    isfile(input_path) || error(
        "Missing $input_path; generate it with " *
        "`python scripts/python/analyze_i3_singular.py prepare-certification`."
    )
    frame = CSV.read(input_path, DataFrame)
    system, _ = build_symbolic_pool_equations()
    root_solutions = Vector{Vector{ComplexF64}}(undef, nrow(frame))
    for row_index in 1:nrow(frame)
        row = frame[row_index, :]
        values = ComplexF64[]
        for prefix in ("x", "y")
            for variable_index in 1:5
                push!(values, complex(row[Symbol("$(prefix)$(variable_index)r")],
                                      row[Symbol("$(prefix)$(variable_index)i")]))
            end
        end
        for prefix in ("z", "t")
            push!(values, complex(row[Symbol("$(prefix)r")], row[Symbol("$(prefix)i")]))
        end
        for variable_index in 1:10
            push!(values, complex(row[Symbol("u$(variable_index)r")],
                                  row[Symbol("u$(variable_index)i")]))
        end
        root_solutions[row_index] = values
    end
    println("Certifying $(length(root_solutions)) saved augmented-system roots...")
    result = certify(system, root_solutions; max_precision=512,
                     show_progress=true, threading=false)
    certs = HomotopyContinuation.certificates(result)
    length(certs) == length(root_solutions) ||
        error("Certification returned $(length(certs)) results for $(length(root_solutions)) inputs.")
    roots = [
        Dict{String,Any}(
            "root_index" => index,
            "certified" => is_certified(certs[index]),
            "approximation" => [[real(value), imag(value)]
                                for value in root_solutions[index]],
        )
        for index in eachindex(certs)
    ]
    output = Dict{String,Any}(
        "input_root_count" => length(root_solutions),
        "certified_root_count" => ncertified(result),
        "uncertified_root_count" => length(certs) - ncertified(result),
        "roots" => roots,
    )
    path = joinpath(OUTPUT_DIR, "root_certification.json")
    MubSearch.write_result(path, output;
                 provenance=MubSearch.provenance_record(
                     seed=CAMPAIGN_SEED,
                     parameters=Dict("max_precision" => 512),
                     solver=Dict("method" => "HomotopyContinuation interval certification",
                                 "input" => "roots_to_certify.csv"),
                     tolerances=Dict("physical_root_filter" => 1e-6),
                     extra=Dict("completion_status" => "completed")))
    println("Certified $(ncertified(result))/$(length(certs)) roots; wrote $path")
    return output, path
end

function main(args)
    mode = isempty(args) ? "counts" : first(args)
    certify_roots = "--certify" in args
    if mode == "counts"
        run_transition_counts()
    elseif mode == "singular"
        run_augmented_solve(; certify_roots)
    elseif mode == "cliques"
        run_clique_exports()
    elseif mode == "certify-saved"
        run_saved_root_certification()
    elseif mode == "all"
        run_transition_counts()
        run_augmented_solve(; certify_roots)
        run_clique_exports()
        run_saved_root_certification()
    elseif mode == "inspect"
        system, variables = build_symbolic_pool_equations()
        println("Augmented I3 singular system: $(length(variables)) variables, " *
                "$(length(system)) equations")
    else
        error("Usage: julia --project=. scripts/julia/i3_singular_locus.jl " *
              "[counts|singular|cliques|certify-saved|all] [--certify]")
    end
end

main(ARGS)
