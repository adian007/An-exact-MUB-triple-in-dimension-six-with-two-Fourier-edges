# Attempt Oscar.jl install and Groebner probe (optional).
# Usage: julia --project=. scripts/julia/install_oscar_probe.jl

using Pkg

const LOG = joinpath(@__DIR__, "..", "..", "results", "track_c_elimination", "oscar_install.log")
mkpath(dirname(LOG))

function main()
    open(LOG, "w") do io
        println(io, "=== Oscar.jl install probe ===")
        println(io, "Date: $(Dates.now())\n")
        try
            @eval using Oscar
            println(io, "Oscar already installed: $(Oscar.versioninfo())")
            println("Oscar available")
        catch
            println(io, "Oscar not installed; attempting Pkg.add (may take 30-60 min)...")
            try
                Pkg.add("Oscar")
                @eval using Oscar
                println(io, "Install OK: $(Oscar.versioninfo())")
                println("Oscar installed successfully")
            catch e
                println(io, "Install FAILED: $e")
                println(io, "Use Macaulay2 Docker: scripts/docker/run_m2.ps1 pool_Dita.m2")
                println("Oscar install failed — use M2 Docker fallback")
            end
        end
    end
end

using Dates
main()
