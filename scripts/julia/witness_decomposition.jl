# Task 3: Numerical irreducible decomposition at F6 and Dita anchors.
# Uses HomotopyContinuation.jl on the per-H MU pool polynomial system.
#
# Usage: julia --project=. witness_decomposition.jl

include(joinpath(@__DIR__, "_paths.jl"))
using Printf
import HomotopyContinuation: degree, dim, codim, witness_set,
    numerical_irreducible_decomposition, witness_sets, is_irreducible,
    monodromy_solve, nvariables

const OUT = joinpath(RESULTS_DIR, "witness_decomposition.txt")
const DITA_THETA = acos(1 / sqrt(3))

function analyze_H(name, H::AbstractMatrix{ComplexF64})
    println("\n=== $name ===")
    @assert is_hadamard(H; tol = 1e-8) "H is not Hadamard"
    system = build_numeric_pool_system(H)
    mv = mixed_volume(system)
    nvars = nvariables(system)
    neqs = length(system)
    println("  System: $neqs equations, $nvars variables, mixed_volume=$mv")

    res = solve(system; show_progress = false)
    n_sols = length(solutions(res))
    cert = certify(system, res; show_progress = false, threading = false, max_precision = 64)
    n_cert = ncertified(cert)
    println("  Tracked solutions: $n_sols, certified: $n_cert")

    wit_result = nothing
    for attempt in [(:dim, 0), (:codim, nvars)]
        try
            W = witness_set(system; attempt[1] => attempt[2], show_progress = false)
            wit_result = (
                dim = dim(W), codim = codim(W), degree = degree(W),
                n_points = length(solutions(W)), method = attempt,
            )
            @printf("  witness_set(%s=%d): dim=%d codim=%d degree=%d\n",
                    attempt[1], attempt[2], wit_result.dim, wit_result.codim, wit_result.degree)
            break
        catch e
            println("  witness_set($(attempt[1])=$(attempt[2])) failed: $(sprint(showerror, e))")
        end
    end

    nid_result = nothing
    try
        println("  Running numerical_irreducible_decomposition...")
        N = numerical_irreducible_decomposition(system; show_progress = false)
        Ws = witness_sets(N)
        components = NamedTuple[]
        for (d, ws) in Ws
            for W in ws
                ir = is_irreducible(W)
                deg = degree(W)
                pts = solutions(W)
                c = (dim = d, degree = deg, irreducible = ir, n_points = length(pts))
                push!(components, c)
                @printf("    dim=%d degree=%d irreducible=%s n_points=%d\n", d, deg, ir, length(pts))
            end
        end
        nid_result = components
    catch e
        println("  NID failed: $(sprint(showerror, e))")
        try
            println("  Fallback: monodromy_solve (target=$mv)...")
            mres = monodromy_solve(system; target_solutions_count = mv, show_progress = false)
            ns = length(solutions(mres))
            println("  monodromy tracked: $ns solutions")
            nid_result = [(monodromy_solutions = ns, target = mv, complete = ns >= mv)]
        catch e2
            println("  monodromy failed: $(sprint(showerror, e2))")
        end
    end

    reps = [s[1:min(3, length(s))] for s in solutions(res)[1:min(3, n_sols)]]
    zero_dim = wit_result !== nothing ? wit_result.dim == 0 : true

    return (
        name = name,
        n_eqs = neqs, n_vars = nvars, mixed_volume = mv,
        n_tracked = n_sols, n_certified = n_cert,
        witness = wit_result,
        nid = nid_result,
        representatives = reps,
        zero_dim = zero_dim,
    )
end

function main()
    anchors = [
        ("F6 (Fourier matrix)", F6),
        ("Dita Karlsson", build_karlsson_family(DITA_THETA, pi / 4, 0.4)),
    ]

    results = [analyze_H(name, H) for (name, H) in anchors]

    open(OUT, "w") do io
        println(io, "=== Task 3: Witness decomposition at F6 and Dita ===")
        for r in results
            println(io, "\n--- $(r.name) ---")
            @printf(io, "  eqs=%d vars=%d mv=%d tracked=%d certified=%d\n",
                    r.n_eqs, r.n_vars, r.mixed_volume, r.n_tracked, r.n_certified)
            if r.witness !== nothing
                @printf(io, "  witness: dim=%d codim=%d degree=%d zero-dim=%s\n",
                        r.witness.dim, r.witness.codim, r.witness.degree, r.zero_dim)
            end
            if r.nid !== nothing
                println(io, "  NID / monodromy components:")
                for c in r.nid
                    println(io, "    $c")
                end
            end
            if r.zero_dim && r.witness !== nothing && r.witness.degree <= 300
                println(io, "  RECOMMEND: Groebner elimination feasible (zero-dim, degree=$(r.witness.degree))")
            elseif r.zero_dim
                println(io, "  RECOMMEND: zero-dim but degree=$(r.witness.degree) — Groebner costly")
            else
                println(io, "  RECOMMEND: deprioritize symbolic Groebner (positive-dim or witness_set failed)")
            end
            if !isempty(r.representatives)
                println(io, "  Representative points (first 3 coords of first solutions):")
                for (i, rep) in enumerate(r.representatives)
                    println(io, "    sol[$i]: $rep")
                end
            end
        end
    end
    println("\nWrote $OUT")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
