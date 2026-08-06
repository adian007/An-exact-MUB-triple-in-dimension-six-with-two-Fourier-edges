include(joinpath(@__DIR__, "_paths.jl"))
using CSV, DataFrames, Printf, Dates

const OUT = joinpath(RESULTS_DIR, "locus_dimension_confirm.txt")
const GRID = joinpath(RESULTS_DIR, "locus_2d_grid.csv")
const DITA_PHI = pi / 4

function main()
    df = CSV.read(GRID, DataFrame)
    n = nrow(df)
    n6 = 0
    n6_phi_row = 0
    phi_tol = 1e-6
    phi_counts = Dict{Float64, Int}()
    for r in eachrow(df)
        ph = r.phi
        c = r.max_clique
        if c >= 6
            n6 += 1
            if abs(ph - DITA_PHI) < phi_tol
                n6_phi_row += 1
            end
            phi_counts[round(ph, digits = 8)] = get(phi_counts, round(ph, digits = 8), 0) + 1
        end
    end

    lines = String[]
    push!(lines, "=== Item 2: 1D vs 2D locus dimension at Dita ===")
    push!(lines, "Date: $(Dates.now())")
    push!(lines, "Grid: locus_2d_grid.csv ($n points, 20x20)")
    push!(lines, @sprintf("theta=arccos(1/sqrt(3)), phi in [pi/4-0.2, pi/4+0.2], lambda in [0.2, 0.6]"))
    push!(lines, "")
    push!(lines, @sprintf("clique>=6 total: %d / %d (%.1f%%)", n6, n, 100 * n6 / n))
    push!(lines, @sprintf("clique>=6 on phi=pi/4 row: %d / 20", n6_phi_row))
    push!(lines, "")
    push!(lines, "Per-phi clique>=6 counts:")
    for ph in sort(collect(keys(phi_counts)))
        push!(lines, @sprintf("  phi=%.12g: %d", ph, phi_counts[ph]))
    end
    push!(lines, "")
    if n6_phi_row == 20 && n6 == n6_phi_row
        push!(lines, "CONCLUSION: 1D CURVE CONFIRMED — all $n6 clique>=6 points lie on phi=pi/4 row only")
        push!(lines, "  No 2D patch: 0 clique>=6 off the phi=pi/4 slice")
    elseif n6 > n6_phi_row
        push!(lines, "CONCLUSION: 2D REGION — clique>=6 extends beyond phi=pi/4 ($n6 total, $n6_phi_row on phi=pi/4)")
    else
        push!(lines, "CONCLUSION: INCONCLUSIVE — sparse clique>=6 hits")
    end

    open(OUT, "w") do io
        for ln in lines
            println(io, ln)
        end
    end
    for ln in lines
        println(ln)
    end
end

using Dates
main()