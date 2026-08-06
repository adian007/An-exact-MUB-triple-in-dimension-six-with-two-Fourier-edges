# Track A: dense Dita lambda-circle fourth-MUB test (step 0.01 on full 2π circle).
#
# Usage:
#   julia --project=. scripts/julia/dita_lambda_fourth_dense.jl --max-points 10
#   julia --project=. scripts/julia/dita_lambda_fourth_dense.jl --step 0.01
#   julia --project=. scripts/julia/dita_lambda_fourth_dense.jl --step 0.01 --resume

include(joinpath(@__DIR__, "_paths.jl"))
using Printf, Dates, CSV, DataFrames

const DITA_THETA = acos(1 / sqrt(3))
const DITA_PHI = pi / 4
const OUT_CSV = joinpath(RESULTS_DIR, "dita_lambda_fourth_dense.csv")
const OUT_LOG = joinpath(RESULTS_DIR, "dita_lambda_fourth_dense.log")
const MU_TOLS = (1e-6, 1e-8, 1e-10)

function parse_args()
    step = 0.01
    max_pts = typemax(Int)
    resume = false
    for (i, a) in enumerate(ARGS)
        a == "--step" && i < length(ARGS) && (step = parse(Float64, ARGS[i + 1]))
        a == "--max-points" && i < length(ARGS) && (max_pts = parse(Int, ARGS[i + 1]))
        a == "--resume" && (resume = true)
    end
    return (; step = step, max_points = max_pts, resume = resume)
end

function lambda_grid(step)
    n = max(1, round(Int, 2 * pi / step))
    return [k * step for k in 0:(n - 1)]
end

function probe(lam)
    H = build_karlsson_family(DITA_THETA, DITA_PHI, lam)
    !is_hadamard(H; tol = 1e-8) &&
        return (status = :not_hadamard, max_clique = -1, found_third = false, found_fourth = false,
                found_fourth_1e6 = false, found_fourth_1e8 = false, found_fourth_1e10 = false, n_pool = 0)
    pool, _ = generate_candidate_pool_fresh(H; verbose = false)
    pool = deduplicate_pool(pool)
    ext = check_four_mub_extension(pool, H)
    by_tol = Dict{Float64, NamedTuple}()
    for mu_tol in MU_TOLS
        by_tol[mu_tol] = check_four_mub_extension(pool, H; mu_tol = mu_tol)
    end
    return (
        status = :ok,
        max_clique = ext.max_clique,
        found_third = ext.found_third,
        found_fourth = ext.found_fourth,
        found_fourth_1e6 = by_tol[1e-6].found_fourth,
        found_fourth_1e8 = by_tol[1e-8].found_fourth,
        found_fourth_1e10 = by_tol[1e-10].found_fourth,
        n_pool = length(pool),
    )
end

function load_done()
    seen = Set{Float64}()
    !isfile(OUT_CSV) && return seen
    df = CSV.read(OUT_CSV, DataFrame)
    for r in eachrow(df)
        push!(seen, round(r.lambda, digits = 12))
    end
    return seen
end

function main()
    cfg = parse_args()
    lams = lambda_grid(cfg.step)
    length(lams) > cfg.max_points && (lams = lams[1:cfg.max_points])
    done = cfg.resume ? load_done() : Set{Float64}()
    cfg.resume && !isempty(done) && println("Resume: $(length(done)) lambdas already in CSV")

    mkpath(dirname(OUT_CSV))
    io_csv = open(OUT_CSV, cfg.resume && isfile(OUT_CSV) ? "a" : "w")
    if !(cfg.resume && isfile(OUT_CSV))
        println(io_csv, "lambda,status,max_clique,found_third,found_fourth,found_fourth_1e6,found_fourth_1e8,found_fourth_1e10,n_pool")
    end

    lines = String[]
    push!(lines, "=== Dita lambda-circle dense fourth-MUB test ===")
    push!(lines, "Date: $(Dates.now())")
    push!(lines, @sprintf("  theta=arccos(1/sqrt(3)), phi=pi/4, step=%.6g, n_lambda=%d", cfg.step, length(lams)))
    push!(lines, "")

    n_ok = n_third = n_fourth = n_c6 = 0
    fourth_hits = Float64[]

    for lam in lams
        key = round(lam, digits = 12)
        key in done && continue
        p = probe(lam)
        if p.status == :ok
            n_ok += 1
            p.found_third && (n_third += 1)
            p.found_fourth && (n_fourth += 1; push!(fourth_hits, lam))
            p.max_clique >= 6 && (n_c6 += 1)
        end
        println(io_csv, @sprintf("%.17g,%s,%d,%s,%s,%s,%s,%s,%d",
            lam, p.status, p.max_clique, p.found_third, p.found_fourth,
            p.found_fourth_1e6, p.found_fourth_1e8, p.found_fourth_1e10, p.n_pool))
        flush(io_csv)
        push!(lines, @sprintf("  lambda=%.6g: clique=%d third=%s fourth=%s",
            lam, p.max_clique, p.found_third, p.found_fourth))
        p.found_fourth && push!(lines, "  *** COUNTEREXAMPLE at lambda=$lam ***")
    end
    close(io_csv)

    push!(lines, "")
    push!(lines, @sprintf("SUMMARY: probed=%d ok=%d clique>=6=%d third=%d fourth=%d",
        length(lams) - length(done), n_ok, n_c6, n_third, n_fourth))
    if !isempty(fourth_hits)
        push!(lines, "FOURTH MUB HITS: $(fourth_hits)")
    end

    open(OUT_LOG, cfg.resume && isfile(OUT_LOG) ? "a" : "w") do io
        for ln in lines
            println(io, ln)
        end
    end
    for ln in lines
        println(ln)
    end
    println("\nWrote $OUT_CSV and $OUT_LOG")
end

main()
