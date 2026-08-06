# Top-50 z4_dev degeneracy batch (supplementary locus stress test).
# Usage: julia --project=. scripts/julia/run_top50_degen.jl

include(joinpath(@__DIR__, "_paths.jl"))

using Printf
using Dates

const DEGEN_CSV = joinpath(RESULTS_DIR, "degeneracy_candidates.csv")
const OUT = joinpath(RESULTS_DIR, "top50_degen_stress.txt")

function load_top50_degen()
    isfile(DEGEN_CSV) || error("Missing $DEGEN_CSV")
    lines = filter(!isempty, strip.(readlines(DEGEN_CSV)))
    header = lines[1]
    rows = lines[2:end]
    # Sort by z4_dev descending if column present
    parsed = Vector{NamedTuple}()
    for line in rows
        parts = split(line, ',')
        length(parts) < 4 && continue
        th = parse(Float64, parts[1])
        ph = parse(Float64, parts[2])
        lam = parse(Float64, parts[3])
        z4 = length(parts) >= 4 ? tryparse(Float64, parts[4]) : 0.0
        push!(parsed, (theta=th, phi=ph, lambda=lam, z4_dev=something(z4, 0.0)))
    end
    sort!(parsed, by=r -> -r.z4_dev)
    return parsed[1:min(50, length(parsed))]
end

function main()
    pts = load_top50_degen()
    n_clique6 = 0
    n_third = 0
    n_fourth = 0

    open(OUT, "w") do io
        println(io, "=== Top-50 z4_dev degeneracy stress test ===")
        println(io, "Date: $(Dates.format(now(), "yyyy-mm-ddTHH:MM:SS"))\n")
        for (i, p) in enumerate(pts)
            try
                H = build_karlsson_family(p.theta, p.phi, p.lambda)
                !is_hadamard(H) && continue
                pool, = generate_candidate_pool_fresh(H; verbose=false)
                pool = deduplicate_pool(pool)
                rep = pool_completeness_report(H; pool=pool)
                ext = check_four_mub_extension(pool, H)
                ext.max_clique >= 6 && (n_clique6 += 1)
                ext.found_third && (n_third += 1)
                ext.found_fourth && (n_fourth += 1)
                @printf(io, "%3d z4=%.2e th=%.4g ph=%.4g lam=%.4g clique=%d third=%s fourth=%s complete=%s\n",
                        i, p.z4_dev, p.theta, p.phi, p.lambda,
                        ext.max_clique, ext.found_third, ext.found_fourth, rep.pool_complete_flag)
            catch e
                @printf(io, "%3d FAIL th=%.4g ph=%.4g lam=%.4g err=%s\n",
                        i, p.theta, p.phi, p.lambda, e)
            end
        end
        println(io, "\nSUMMARY: clique6=$n_clique6 third=$n_third fourth=$n_fourth / $(length(pts))")
    end
    println("Wrote $OUT")
end

main()
