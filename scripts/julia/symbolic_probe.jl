# Symbolic probe: mixed volume of per-H pool system (root count upper bound at generic H).
# Full clique-existence elimination is intractable in this script; we certify the
# MU-pool subsystem structure instead.
#
# Usage: julia --project=. symbolic_probe.jl

include(joinpath(@__DIR__, "_paths.jl"))

function mixed_volume_at_H(H::AbstractMatrix{ComplexF64})
    system = build_numeric_pool_system(H)
    return mixed_volume(system)
end

function main()
    println("=== Symbolic / combinatorial probes ===\n")

    println("[S1] Mixed volume of per-H pool system (BKK bound on root count)")
    for (name, H) in [
        ("F6", F6),
        ("generic", build_karlsson_family(0.3, 0.5, 0.2)),
        ("Dita", build_karlsson_family(acos(1 / sqrt(3)), pi / 4, 0.4)),
        ("theta=0", build_karlsson_family(0.0, 0.5, 0.3)),
    ]
        mv = mixed_volume_at_H(H)
        println("  $name: mixed_volume = $mv")
    end

    println("\n[S2] Parametric system structure")
    sys = build_parametric_pool_system()
    println("  variables: ", nvariables(sys))
    println("  equations: ", length(sys))
    println("  parameters: ", length(parameters(sys)))

    println("\n[S3] Clique-existence as algebraic variety — NOT attempted here")
    println("  Encoding 'exists 6-clique among MU vectors' requires ~60+ witness variables")
    println("  (6 vectors × 10 poly vars) plus orthogonality and MU constraints.")
    println("  Groebner elimination over (theta,phi,lambda) is beyond this pipeline.")
    println("  Negative clique results remain numerical/heuristic unless elimination is done.")

    println("\n[S4] What WOULD constitute a proof")
    println("  (a) Certified root count = pool size at all parameters in a region, AND")
    println("  (b) Symbolic proof that orthogonality graph clique number < 6 on that region, OR")
    println("  (c) An explicit fourth-MUB counterexample with certified arithmetic.")

    println("\n[S5] Tractable subsystems (see symbolic_elimination.jl)")
    println("  theta=0 slice: numerical scan only — run symbolic_elimination.jl")
    println("  Local fourth-MUB witness at F6/Dita: ~60 witness vars — Macaulay2 export fallback")
    println("  Mobius degeneracy locus: exported to symbolic_export/*.m2")
    summary_path = joinpath(RESULTS_DIR, "symbolic_elimination_summary.txt")
    if isfile(summary_path)
        println("\n  Latest elimination summary:")
        for line in eachline(summary_path)
            println("    $line")
        end
    end
end

main()
