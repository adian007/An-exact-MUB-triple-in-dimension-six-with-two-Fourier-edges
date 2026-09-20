include(joinpath(@__DIR__, "_paths.jl"))
using LinearAlgebra

H = build_karlsson_family(acos(1 / sqrt(3)), pi / 4, 0.4)
pool, = generate_candidate_pool_fresh(H; verbose=false)
pool = deduplicate_pool(pool)
g = _orthogonality_graph(pool; ortho_tol=1e-8)
cliques = [c for c in maximal_cliques(g) if length(c) >= 6]
c = argmax(cl -> length(cl), cliques)
B3 = hcat([pool[i] for i in c]...)

ω = exp(-2π * im / 6)
roots = [ω^k for k in 0:5]

function root_dev(z)
    return minimum(abs.(z .- roots))
end

println("clique: ", c)
println("max |H| dev from 6th roots: ", maximum(root_dev(x) for x in H))
println("max |B3| dev from 6th roots: ", maximum(root_dev(x) for x in B3))
println("1/sqrt(6) ≈ ", 1 / sqrt(6))
println("unique rounded |H|: ", sort(unique(round.(abs.(H[:]), sigdigits=8))))
println("unique rounded |B3|: ", sort(unique(round.(abs.(B3[:]), sigdigits=8))))

# Check if entries match simple algebraic templates
templates = Dict(
    "0" => 0.0 + 0im,
    "1" => 1.0 + 0im,
    "-1" => -1.0 + 0im,
    "i" => im,
    "-i" => -im,
    "1/sqrt(6)" => 1 / sqrt(6),
    "i/sqrt(6)" => im / sqrt(6),
    "sqrt3/2" => sqrt(3) / 2,
    "cos(pi/4)" => cos(pi / 4),
    "sin(pi/4)" => sin(pi / 4),
)
for (name, M) in [("H", H), ("B3", B3)]
    println("\n--- $name ---")
    for z in M[:]
        best = ("", Inf)
        for (label, t) in templates
            d = abs(z - t)
            if d < best[2]
                best = (label, d)
            end
            d = abs(z - im * t)
            if d < best[2]
                best = ("i*($label)", d)
            end
        end
        for k in 0:5
            d = abs(z - roots[k + 1])
            if d < best[2]
                best = ("ω^$k", d)
            end
        end
        if best[2] > 1e-8
            println("  non-template: ", z, " best=", best)
        end
    end
end
