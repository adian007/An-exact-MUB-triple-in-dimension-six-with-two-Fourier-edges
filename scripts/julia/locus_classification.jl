# Phase 2 Track D: cluster HP-verified clique-6 loci, CHM equivalence, classification JSON.
#
# Usage:
#   julia --project=. scripts/julia/locus_classification.jl
#   julia --project=. scripts/julia/locus_classification.jl --csv results/special_loci_search.csv

include(joinpath(@__DIR__, "_paths.jl"))

using Dates
using LinearAlgebra
using Printf
using Statistics

const DEFAULT_SEARCH_CSV = joinpath(RESULTS_DIR, "special_loci_search.csv")
const DEFAULT_AUDIT_CSV = joinpath(RESULTS_DIR, "third_mub_audit.csv")
const OUT_JSON = joinpath(RESULTS_DIR, "locus_classification.json")
const OUT_TXT = joinpath(RESULTS_DIR, "locus_classification.txt")
const DITA_THETA = acos(1 / sqrt(3))
const CHM_TOL = 1e-10

function load_search_csv(path)
    lines = filter(!isempty, strip.(readlines(path)))
    header = split(replace(lines[1], "\ufeff" => ""), ',')
    idx = Dict(h => i for (i, h) in enumerate(header))
    rows = NamedTuple[]
    for line in lines[2:end]
        parts = split(line, ',')
        length(parts) < length(header) && continue
        g(k) = parts[idx[k]]
        push!(rows, (
            name = g("name"),
            theta = parse(Float64, g("theta")),
            phi = parse(Float64, g("phi")),
            lambda = parse(Float64, g("lambda")),
            status = g("status"),
            max_clique = isempty(g("max_clique")) ? 0 : parse(Int, g("max_clique")),
        ))
    end
    return rows
end

function load_audit_csv(path)
    lines = filter(!isempty, strip.(readlines(path)))
    header = split(replace(lines[1], "\ufeff" => ""), ',')
    idx = Dict(h => i for (i, h) in enumerate(header))
    rows = NamedTuple[]
    for line in lines[2:end]
        parts = split(line, ',')
        length(parts) < length(header) && continue
        g(k) = parts[idx[k]]
        push!(rows, (
            name = g("name"),
            theta = parse(Float64, g("theta")),
            phi = parse(Float64, g("phi")),
            lambda = parse(Float64, g("lambda")),
            hp_survival = g("hp_survival"),
            nearest_locus = g("nearest_locus"),
        ))
    end
    return rows
end

function hp_verified_clique6(search_rows, audit_rows)
    audit_keys = Set{Tuple{String, Float64, Float64, Float64}}()
    for r in audit_rows
        r.hp_survival == "true" || continue
        push!(audit_keys, (r.name, r.theta, r.phi, r.lambda))
    end
    out = NamedTuple[]
    for r in search_rows
        r.status != "ok" && continue
        r.max_clique < 6 && continue
        key = (r.name, r.theta, r.phi, r.lambda)
        key in audit_keys || continue
        push!(out, r)
    end
    return out
end

function dephase_chm(H; tol = 1e-8)
    H = copy(H)
    for j in 2:size(H, 2)
        abs(H[1, j]) > tol && (H[:, j] .*= conj(H[1, j]) / abs(H[1, j]))
    end
    for i in 2:size(H, 1)
        abs(H[i, 1]) > tol && (H[i, :] .*= conj(H[i, 1]) / abs(H[i, 1]))
    end
    return H
end

function _phase(z)
    abs(z) < 1e-15 && return 1.0 + 0im
    return z / abs(z)
end

function chm_residual_for_perm(H1, H2, perm)
    H2p = H2[:, perm]
    N = 6
    D_r = ones(ComplexF64, N)
    for i in 1:N
        abs(H2p[i, 1]) > 1e-15 && (D_r[i] = _phase(H1[i, 1] / H2p[i, 1]))
    end
    D_c = ones(ComplexF64, N)
    for j in 2:N
        ratios = ComplexF64[]
        for i in 1:N
            if abs(H2p[i, j]) > 1e-15 && abs(D_r[i]) > 1e-15
                push!(ratios, H1[i, j] / (D_r[i] * H2p[i, j]))
            end
        end
        if !isempty(ratios)
            r0 = ratios[1]
            phases = [_phase(r / r0) for r in ratios]
            D_c[j] = mean(phases) * r0
            abs(D_c[j]) > 1e-15 && (D_c[j] = _phase(D_c[j]))
        end
    end
    H_fit = D_r .* H2p .* D_c'
    res = norm(H1 - H_fit)
    return isfinite(res) ? res : Inf
end

const _PERMS6 = begin
    function perms(n)
        n == 1 && return [[1]]
        out = Vector{Int}[]
        for p in perms(n - 1), i in 1:n
            q = copy(p); insert!(q, i, n); push!(out, q)
        end
        return out
    end
    perms(6)
end

function chm_equivalence_residual(H1, H2; tol = CHM_TOL)
    best = (residual = Inf, variant = "", equivalent = false)
    D1 = dephase_chm(H1)
    for (vname, H2v) in [("H", H2), ("H.T", transpose(H2)),
                          ("conj(H)", conj.(H2)), ("conj(H).T", transpose(conj.(H2)))]
        D2 = dephase_chm(H2v)
        for perm in _PERMS6
            res = chm_residual_for_perm(D1, D2, perm)
            res < best.residual && (best = (residual = res, variant = vname, equivalent = res < tol))
        end
    end
    return best
end

function locus_label(theta, phi, lam; atol = 0.02)
    abs(theta) < atol && return "F6_theta0"
    abs(theta - DITA_THETA) < atol && abs(phi - pi / 4) < atol && return "Dita"
    abs(theta - pi / 2) < atol && return "circulant_match"
    return "other"
end

function cluster_hp_points(points)
    clusters = Dict{String, Vector{NamedTuple}}()
    for p in points
        push!(get!(clusters, locus_label(p.theta, p.phi, p.lambda), NamedTuple[]), p)
    end
    return clusters
end

function cluster_representative(pts)
    # Prefer primary anchor (non-ref) for valid Karlsson build at theta=0
    primaries = sort([p for p in pts if !occursin("ref_", p.name)]; by = p -> p.name)
    anchor = !isempty(primaries) ? primaries[1] : pts[1]
    return (
        theta = anchor.theta,
        phi = anchor.phi,
        lambda = anchor.lambda,
        name = anchor.name,
    )
end

function topology_for_label(lbl, pts)
    lbl == "Dita" && return "1D circle in lambda (Dita, phi=pi/4)"
    lbl == "F6_theta0" && begin
        lams = extrema(p.lambda for p in pts)
        return @sprintf("1D bounded arc in lambda at theta=0 [%.4f, %.4f]", lams[1], lams[2])
    end
    lbl == "circulant_match" && return "failed HP at theta=pi/2 (heuristic artifact)"
    return "unknown"
end

function chm_class_for(lbl, residual_f6, residual_dita)
    lbl == "F6_theta0" && return residual_f6 < CHM_TOL ? "F6_family_phi_gauge" : "F6_related"
    lbl == "Dita" && return "Dita_lambda_family_CHM_inequivalent"
    return "unclassified"
end

function json_escape(s)
    replace(s, "\\" => "\\\\", "\"" => "\\\"")
end

function write_json(path, entries)
    open(path, "w") do io
        println(io, "[")
        for (i, e) in enumerate(entries)
            sep = i < length(entries) ? "," : ""
            println(io, "  {")
            println(io, "    \"name\": \"$(json_escape(e.name))\",")
            @printf(io, "    \"coords\": {\"theta\": %.17g, \"phi\": %.17g, \"lambda\": %.17g},\n",
                    e.coords_theta, e.coords_phi, e.coords_lambda)
            println(io, "    \"topology\": \"$(json_escape(e.topology))\",")
            println(io, "    \"chm_class\": \"$(json_escape(e.chm_class))\",")
            println(io, "    \"hp_count\": $(e.hp_count)")
            print(io, "  }$sep\n")
        end
        println(io, "]")
    end
end

function run_classification(search_csv = DEFAULT_SEARCH_CSV, audit_csv = DEFAULT_AUDIT_CSV)
    println("=== Phase 2: Track D locus classification ===\n")
    search_rows = load_search_csv(search_csv)
    audit_rows = load_audit_csv(audit_csv)
    hp_pts = hp_verified_clique6(search_rows, audit_rows)
    println("HP-verified clique-6: $(length(hp_pts)) (from audit + search CSV)\n")
    isempty(hp_pts) && error("No HP-verified clique-6 points")

    clusters = cluster_hp_points(hp_pts)
    H_f6 = build_karlsson_family(0.0, 0.5, 0.3)
    H_dita = build_karlsson_family(DITA_THETA, pi / 4, 0.4)

    entries = NamedTuple[]
    chm_lines = String[]
    for (lbl, pts) in sort(collect(clusters); by = x -> -length(x[2]))
        rep = cluster_representative(pts)
        H_rep = build_karlsson_family(rep.theta, rep.phi, rep.lambda)
        r_f6 = chm_equivalence_residual(H_rep, H_f6).residual
        r_di = chm_equivalence_residual(H_rep, H_dita).residual
        @printf("Cluster %-18s n=%3d  rep=(%.6f, %.6f, %.6f)  CHM vs F6=%.2e vs Dita=%.2e\n",
                lbl, length(pts), rep.theta, rep.phi, rep.lambda, r_f6, r_di)
        push!(entries, (
            name = lbl,
            coords_theta = rep.theta, coords_phi = rep.phi, coords_lambda = rep.lambda,
            topology = topology_for_label(lbl, pts),
            chm_class = chm_class_for(lbl, r_f6, r_di),
            hp_count = length(pts),
        ))
    end

    # CHM between cluster representatives
    labels = [e.name for e in entries]
    for i in 1:length(labels), j in (i + 1):length(labels)
        Hi = build_karlsson_family(entries[i].coords_theta, entries[i].coords_phi, entries[i].coords_lambda)
        Hj = build_karlsson_family(entries[j].coords_theta, entries[j].coords_phi, entries[j].coords_lambda)
        r = chm_equivalence_residual(Hi, Hj)
        ln = @sprintf("CHM %s vs %s: residual=%.6e equiv=%s", labels[i], labels[j], r.residual, r.equivalent)
        push!(chm_lines, ln)
        println(ln)
    end

    write_json(OUT_JSON, entries)

    open(OUT_TXT, "w") do io
        println(io, "=== Track D locus classification (Phase 2) ===")
        println(io, "Date: $(Dates.now())")
        println(io, "Search CSV: $search_csv")
        println(io, "Audit CSV: $audit_csv")
        println(io, "HP-verified clique-6: $(length(hp_pts))")
        println(io, "")
        for e in entries
            println(io, "--- $(e.name) ---")
            @printf(io, "  coords: theta=%.17g phi=%.17g lambda=%.17g\n",
                    e.coords_theta, e.coords_phi, e.coords_lambda)
            println(io, "  topology: $(e.topology)")
            println(io, "  chm_class: $(e.chm_class)")
            println(io, "  hp_count: $(e.hp_count)")
            println(io, "")
        end
        println(io, "=== CHM between cluster representatives ===")
        for ln in chm_lines
            println(io, ln)
        end
        println(io, "")
        r_cross = chm_equivalence_residual(H_f6, H_dita)
        @printf(io, "F6 vs Dita cross-locus: residual=%.6e equiv=%s (distinct components)\n",
                r_cross.residual, r_cross.equivalent)
        println(io, "JSON: $OUT_JSON")
    end
    println("\nWrote $OUT_JSON")
    println("Wrote $OUT_TXT")
    return entries
end

function main()
    search_csv = DEFAULT_SEARCH_CSV
    if "--csv" in ARGS
        i = findfirst(==("--csv"), ARGS)
        i !== nothing && length(ARGS) >= i + 1 && (search_csv = ARGS[i + 1])
    end
    run_classification(search_csv)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
