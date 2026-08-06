# Task 1: Re-certify the 20-point quick-run anchor set (no new search points).
include(joinpath(@__DIR__, "_paths.jl"))
using Printf

function anchor_points()
    dita_theta = acos(1 / sqrt(3))
    anchors = NamedTuple[]
    push!(anchors, (:name => "F6_cyclic", :theta => 0.0, :phi => 0.0, :lambda => 0.0))
    push!(anchors, (:name => "F6_theta0", :theta => 0.0, :phi => 0.5, :lambda => 0.3))
    push!(anchors, (:name => "Dita", :theta => dita_theta, :phi => pi / 4, :lambda => 0.4))
    push!(anchors, (:name => "Mobius_audit_fail",
                     :theta => 0.9553166181245092, :phi => 0.7853981633974483, :lambda => 0.4))
    for (name, th) in [("theta_0", 0.0), ("theta_pi2", pi / 2), ("theta_pi3", pi / 3),
                       ("theta_pi4", pi / 4), ("theta_pi6", pi / 6)]
        push!(anchors, (:name => name, :theta => th, :phi => pi / 4, :lambda => pi / 6))
    end
    for k in 0:11
        push!(anchors, (:name => "lambda_kpi6_k=$k", :theta => 0.3, :phi => 0.5, :lambda => k * pi / 6))
    end
    return anchors
end

function main()
    seen = Set{Tuple{Float64,Float64,Float64}}()
    rows = NamedTuple[]
    for pt in anchor_points()
        key = (round(pt.theta, digits = 10), round(pt.phi, digits = 10), round(pt.lambda, digits = 10))
        key in seen && continue
        push!(seen, key)

        H = build_karlsson_family(pt.theta, pt.phi, pt.lambda)
        hadamard = is_hadamard(H)
        if !hadamard
            push!(rows, (
                name = pt.name, theta = pt.theta, phi = pt.phi, lambda = pt.lambda,
                status = :not_hadamard, pool_complete = missing, n_tracked = missing,
                n_certified = missing, n_verified = missing, max_clique = missing,
                found_third = missing, found_fourth = missing, cert_error = "",
            ))
            continue
        end

        pc = pool_completeness_report(H; verbose = false)
        pool, _ = generate_candidate_pool_fresh(H; verbose = false)
        pool = deduplicate_pool(pool)
        ext = check_four_mub_extension(pool, H)
        push!(rows, (
            name = pt.name, theta = pt.theta, phi = pt.phi, lambda = pt.lambda,
            status = :ok, pool_complete = pc.pool_complete_flag,
            n_tracked = pc.n_tracked, n_certified = pc.n_certified,
            n_verified = pc.n_verified, max_clique = ext.max_clique,
            found_third = ext.found_third, found_fourth = ext.found_fourth, cert_error = "",
        ))
    end

    ok = filter(r -> r.status == :ok, rows)
    complete = filter(r -> r.pool_complete === true, ok)
    incomplete = filter(r -> r.pool_complete !== true, ok)

    println("=== Task 1 quick-run anchor verification ($(length(rows)) unique points) ===")
    @printf("%-18s %12s %12s %12s %6s %6s %6s %6s %s\n",
            "name", "theta", "phi", "lambda", "nt", "nc", "nv", "cliq", "complete")
    for r in rows
        if r.status == :not_hadamard
            @printf("%-18s %12.10f %12.10f %12.10f  — not_hadamard\n",
                    r.name, r.theta, r.phi, r.lambda)
        else
            @printf("%-18s %12.10f %12.10f %12.10f %6d %6d %6d %6d %6s %s\n",
                    r.name, r.theta, r.phi, r.lambda,
                    r.n_tracked, r.n_certified, r.n_verified, r.max_clique,
                    r.pool_complete, r.found_third ? "third" : "")
        end
    end
    println("\nSummary: $(length(ok)) Hadamard-ok, $(length(complete)) pool_complete=true, $(length(incomplete)) incomplete")
    for r in incomplete
        println("  INCOMPLETE: $(r.name) tracked=$(r.n_tracked) certified=$(r.n_certified) verified=$(r.n_verified)")
    end
end

main()
