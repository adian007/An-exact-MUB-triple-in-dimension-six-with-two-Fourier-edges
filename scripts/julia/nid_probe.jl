# Part A3: HomotopyContinuation NID / witness_set probe (Phase 2.3 fixed API).
include(joinpath(@__DIR__, "_paths.jl"))
using Pkg, Printf, Dates
import HomotopyContinuation

const OUT = joinpath(RESULTS_DIR, "nid_probe.txt")
const DITA_THETA = acos(1 / sqrt(3))

function hc_version(io)
    println(io, "=== HomotopyContinuation environment ===")
    Pkg.status("HomotopyContinuation"; io = io)
    println(io, "\nLoaded module: HomotopyContinuation")
    try
        println(io, "  pkgversion: ", Pkg.pkgversion(HomotopyContinuation))
    catch
        println(io, "  (pkgversion unavailable)")
    end
    println(io, "Julia: ", VERSION)
end

function try_witness(name, system; io)
    println(io, "\n--- witness_set: $name ---")
    nvars = HomotopyContinuation.nvariables(system)
    neqs = length(system)
    mv = mixed_volume(system)
    println(io, "  eqs=$neqs vars=$nvars mv=$mv")
    # Square 10-eq/10-var: use default corank (dim=0) — do NOT pass codim=nvars (invalid).
    for attempt in [(:default, nothing), (:dim, 0)]
        try
            W = if attempt[1] == :default
                HomotopyContinuation.witness_set(system; show_progress = false)
            else
                HomotopyContinuation.witness_set(system; dim = attempt[2], show_progress = false)
            end
            d = HomotopyContinuation.dim(W)
            c = HomotopyContinuation.codim(W)
            deg = HomotopyContinuation.degree(W)
            ns = length(HomotopyContinuation.solutions(W))
            label = attempt[1] == :default ? "witness_set()" : "witness_set(dim=$(attempt[2]))"
            @printf(io, "  SUCCESS %s: dim=%d codim=%d degree=%d n_solutions=%d\n",
                    label, d, c, deg, ns)
            return (ok = true, dim = d, codim = c, degree = deg)
        catch e
            label = attempt[1] == :default ? "witness_set()" : "witness_set(dim=$(attempt[2]))"
            println(io, "  FAIL $label:")
            showerror(io, e)
            println(io)
        end
    end
    return (ok = false,)
end

function try_nid(name, system; io)
    println(io, "\n--- numerical_irreducible_decomposition: $name ---")
    try
        N = HomotopyContinuation.numerical_irreducible_decomposition(system; show_progress = false)
        Ws = HomotopyContinuation.witness_sets(N)
        dims_deg = NamedTuple[]
        for (d, wlist) in Ws
            for W in wlist
                ir = HomotopyContinuation.is_irreducible(W)
                deg = HomotopyContinuation.degree(W)
                @printf(io, "  component: dim=%d degree=%d irreducible=%s\n", d, deg, ir)
                push!(dims_deg, (dim = d, degree = deg))
            end
        end
        return (ok = true, components = dims_deg)
    catch e
        println(io, "  FAIL NID:")
        showerror(io, e)
        println(io)
        return (ok = false,)
    end
end

function nid_at_locus(name, theta, phi, lam; io)
    println(io, "\n--- Locus NID: $name (θ=$(theta), φ=$(phi), λ=$(lam)) ---")
    H = build_karlsson_family(theta, phi, lam)
    !is_hadamard(H; tol = 1e-8) && (println(io, "  not Hadamard — skip"); return)
    sys = build_numeric_pool_system(H)
    wit = try_witness("$name pool", sys; io)
    nid = try_nid("$name pool", sys; io)
    return wit, nid
end

function main()
    open(OUT, "w") do io
        println(io, "=== Part A3: NID / witness_set probe (Phase 2.3) ===")
        println(io, "Date: $(Dates.format(now(), "yyyy-mm-dd HH:MM"))\n")
        hc_version(io)

        H_dita = build_karlsson_family(DITA_THETA, pi / 4, 0.4)
        sys_dita = build_numeric_pool_system(H_dita)
        wit_d = try_witness("Dita fixed-H pool (λ=0.4)", sys_dita; io)
        nid_d = try_nid("Dita fixed-H pool (λ=0.4)", sys_dita; io)

        # Multi-point NID along Dita lambda slice (Phase 2.3)
        println(io, "\n=== NID along Dita lambda slice (fixed θ,φ) ===")
        for lam in [0.4, 0.5, 0.6, 0.8, 1.0]
            nid_at_locus("Dita_λ=$lam", DITA_THETA, pi / 4, lam; io)
        end

        println(io, "\n=== NID at F6_theta0 (λ=0.3) ===")
        nid_at_locus("F6_theta0", 0.0, 0.5, 0.3; io)

        println(io, "\n--- Parametric pool system (hc,h parameters) ---")
        try
            psys = build_parametric_pool_system()
            println(io, "  vars: ", HomotopyContinuation.variables(psys))
            println(io, "  params: ", length(HomotopyContinuation.parameters(psys)), " parameter slots")
            try_witness("parametric pool (hc,h)", psys; io)
        catch e
            println(io, "  build_parametric_pool_system failed:")
            showerror(io, e)
            println(io)
        end

        println(io, "\n=== Summary ===")
        if wit_d.ok
            @printf(io, "Fixed-H Dita witness: dim=%d degree=%d (0-dim isolated roots at fixed H)\n",
                    wit_d.dim, wit_d.degree)
        else
            println(io, "Fixed-H witness_set: use NID (witness_set API edge case on square systems)")
        end
        if nid_d.ok
            nd = nid_d.components
            all_zero = all(c.dim == 0 for c in nd)
            all_deg1 = all(c.degree == 1 for c in nd)
            @printf(io, "NID: %d components, all dim=0: %s, all degree=1: %s\n",
                    length(nd), all_zero, all_deg1)
        end
        println(io, "Parameter-space (θ,φ,λ) third-MUB locus NID: NOT OBTAINED (combinatorial clique constraint)")
    end
    println("Wrote $OUT")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
