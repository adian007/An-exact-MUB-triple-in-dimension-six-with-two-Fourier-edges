println("=== Julia env check ===")
println("julia version: ", VERSION)
println("project: ", Base.active_project())
try
    using Pkg
    Pkg.status()
catch e
    println("Pkg.status failed: ", e)
end

println("\n=== Checking packages ===")
for pkg in ["Nemo", "AbstractAlgebra", "FLINT", "Hecke", "LinearAlgebra", "Combinatorics"]
    try
        @eval using $(Symbol(pkg))
        println(pkg, " — OK")
    catch e
        println(pkg, " — MISSING: ", typeof(e))
    end
end
