# Certification smoke test at F6.
include(joinpath(@__DIR__, "_paths.jl"))
c = certify_pool_at_H(F6; verbose = true)
println(c)
