# Finalize degen1845 CSV after regen (archive corrupted, promote regen, verify).
# Usage:
#   julia --project=. scripts/julia/finalize_degen1845_csv.jl
include(joinpath(@__DIR__, "_paths.jl"))
using Dates, Printf

const REGEN = joinpath(RESULTS_DIR, "special_loci_degen1845_regen.csv")
const TARGET = joinpath(RESULTS_DIR, "special_loci_degen1845.csv")
const CORRUPT = joinpath(RESULTS_DIR, "special_loci_degen1845_CORRUPTED592.csv")
const META = joinpath(RESULTS_DIR, "special_loci_degen1845.meta.txt")

function count_rows(path)
    n = 0
    open(path) do io
        readline(io)
        for _ in eachline(io)
            n += 1
        end
    end
    return n
end

function main()
    !isfile(REGEN) && error("Missing $REGEN — run search_special_loci.jl --degen-cap 1845 --no-refine first")
    n = count_rows(REGEN)
    n < 1800 && @warn "Regen has only $n rows (target ~1863); promote anyway?"

    if isfile(CORRUPT) && !isfile(TARGET)
        cp(CORRUPT, TARGET)
        println("Archived corrupt copy as $TARGET (from CORRUPTED592)")
    elseif isfile(TARGET)
        backup = TARGET * ".pre_regen_$(Dates.format(now(), "yyyymmdd_HHMMSS")).bak"
        cp(TARGET, backup)
        println("Backed up existing $TARGET → $backup")
    end
    cp(REGEN, TARGET; force=true)
    println("Promoted $REGEN → $TARGET ($n data rows)")

    open(META, "a") do io
        println(io, "")
        println(io, "# Regeneration $(Dates.format(now(), "yyyy-mm-dd HH:MM"))")
        println(io, "csv_status=REGENERATED from special_loci_degen1845_regen.csv")
        println(io, "csv_rows=$n")
    end
    println("Appended status to $META")
    println("Next: python scripts/python/finish_project_audit.py")
end

main()
