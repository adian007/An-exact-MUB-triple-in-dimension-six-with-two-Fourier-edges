# Task 2: Pre/post argmax fix audit for all third-MUB-positive points.
include(joinpath(@__DIR__, "_paths.jl"))
using Printf

const DITA_THETA = acos(1 / sqrt(3))

"""Simulate pre-fix argmax bug: cliques[argmax(c->length(c), cliques)]."""
function check_four_mub_prefixed(pool, H; ortho_tol = 1e-8, mu_tol = 1e-8)
    n = length(pool)
    n == 0 && return (found_third = false, found_fourth = false, max_clique = 0, bug_hit = false)
    g = _orthogonality_graph(pool; ortho_tol = ortho_tol)
    cliques = maximal_cliques(g)
    max_clique = maximum(length(c) for c in cliques)
    found_third = max_clique >= 6
    found_fourth = false
    bug_hit = false
    third_cliques = [c for c in cliques if length(c) >= 6]
    for c in third_cliques
        B = hcat([pool[i] for i in c]...)
        idx = _vectors_mu_to_basis(pool, B; mu_tol = mu_tol)
        sub = pool[idx]
        length(sub) < 6 && continue
        g2 = _orthogonality_graph(sub; ortho_tol = ortho_tol)
        for c2 in maximal_cliques(g2)
            if length(c2) >= 6
                found_fourth = true
                break
            end
        end
        found_fourth && break
    end
    if found_third && !found_fourth
        try
            _ = cliques[argmax(c -> length(c), cliques)]  # Julia 1.12 bug path
        catch
            bug_hit = true
        end
    end
    return (found_third = found_third, found_fourth = found_fourth,
            max_clique = max_clique, bug_hit = bug_hit)
end

function third_mub_points()
    pts = NamedTuple[]
    push!(pts, (:label => "F6 (Fourier)", :theta => missing, :phi => missing, :lambda => missing, :use_F6 => true))
    push!(pts, (:label => "F6_theta0 (θ=0 subfamily)", :theta => 0.0, :phi => 0.5, :lambda => 0.3, :use_F6 => false))
    push!(pts, (:label => "Dita anchor", :theta => DITA_THETA, :phi => pi / 4, :lambda => 0.4, :use_F6 => false))
    # θ=0 subfamily slice from symbolic_elimination (30 pts) — sample representative + CSV anchors
    for (th, ph, lam, lab) in [
        (0.0, pi / 4, pi / 6, "theta_0 slice (π/4,π/6)"),
        (0.0, 0.5, 0.3, "F6_theta0 duplicate"),
    ]
        push!(pts, (:label => lab, :theta => th, :phi => ph, :lambda => lam, :use_F6 => false))
    end
    return pts
end

function eval_point(pt)
    H = pt.use_F6 ? F6 : build_karlsson_family(pt.theta, pt.phi, pt.lambda)
    pool, _ = try
        generate_candidate_pool_fresh(H; verbose = false)
    catch e
        @warn "pool generation failed" pt.label exception = e
        return (post = (found_third = false, found_fourth = false, max_clique = 0),
                pre = (found_third = false, found_fourth = false, max_clique = 0, bug_hit = false),
                n_pool = 0, error = string(e))
    end
    pool = deduplicate_pool(pool)
    post = check_four_mub_extension(pool, H)
    pre = check_four_mub_prefixed(pool, H)
    return (post = post, pre = pre, n_pool = length(pool), error = "")
end

function main()
    src_mtime = stat(joinpath(ROOT, "src", "mub_zauner_6d_liang_chen.jl")).mtime
    quick_log = joinpath(RESULTS_DIR, "special_loci_run_log.txt")
    quick_mtime = isfile(quick_log) ? stat(quick_log).mtime : missing
    out = joinpath(RESULTS_DIR, "task2_argmax_audit.txt")

    open(out, "w") do io
        println(io, "=== Task 2: argmax fix audit ===")
        println(io, "  src/mub_zauner_6d_liang_chen.jl mtime: $src_mtime")
        println(io, "  special_loci_run_log.txt mtime:    $quick_mtime")
        println(io, "  Quick run predates fix: $(quick_mtime !== missing && quick_mtime < src_mtime)")
        println(io, "  check_four_mub_extension: NO try/catch — argmax bug only affects ortho_defect branch")
        println(io, "  found_fourth set BEFORE buggy argmax block → not silently false-negative\n")

        @printf(io, "%-28s %6s %6s %6s %6s %6s %8s\n",
                "point", "cliq", "3rd", "4th", "bug?", "pool", "timing")
        for pt in third_mub_points()
            r = eval_point(pt)
            timing = (quick_mtime !== missing && quick_mtime < src_mtime) ? "pre-fix" : "post-fix"
            @printf(io, "%-28s %6d %6s %6s %6s %6d %8s\n",
                    pt.label, r.post.max_clique,
                    string(r.post.found_third), string(r.post.found_fourth),
                    get(r, :pre, (bug_hit = false,)).bug_hit ? "CRASH" : "ok",
                    r.n_pool, timing)
        end

        csv_path = joinpath(RESULTS_DIR, "special_loci_search_backup848.csv")
        if isfile(csv_path)
            println(io, "\n--- backup848 found_third=true unique coords (quick-run anchors) ---")
            seen = Set{Tuple{Float64,Float64,Float64}}()
            open(csv_path) do csv_io
                header = split(strip(readline(csv_io)), ',')
                for line in eachline(csv_io)
                    parts = split(line, ',')
                    length(parts) < length(header) && continue
                    row = Dict(zip(header, parts))
                    row["found_third"] == "true" || continue
                    (row["name"] in ("F6_theta0", "Dita")) || continue
                    th, ph, lam = parse.(Float64, (row["theta"], row["phi"], row["lambda"]))
                    key = (round(th, digits = 10), round(ph, digits = 10), round(lam, digits = 10))
                    key in seen && continue
                    push!(seen, key)
                    @printf(io, "  %s (%.10f, %.10f, %.10f): csv_4th=%s timing=pre-fix\n",
                            row["name"], th, ph, lam, row["found_fourth"])
                end
            end
        end
    end
    println("Wrote $out")
end

main()
