# Part A: Multi-radius isotropic + directional sampling for third-MUB locus dimension.
# Usage: julia --project=. scripts/julia/third_mub_locus_sampling.jl
include(joinpath(@__DIR__, "_paths.jl"))
using Printf, Random, LinearAlgebra, Dates

const DITA_THETA = acos(1 / sqrt(3))
const OUT = joinpath(RESULTS_DIR, "third_mub_locus_sampling.txt")
const RADII = [1e-6, 1e-5, 1e-4, 1e-3, 1e-2]
const N_ISO = 100
const MAGS = [1e-6, 1e-5, 1e-4, 1e-3, 1e-2]
const MU_TOL = 1e-10
const ORTHO_TOL = 1e-10

const ANCHORS = [
    ("Dita", DITA_THETA, pi / 4, 0.4),
    ("F6_theta0", 0.0, 0.5, 0.3),
]

function probe_clique(theta, phi, lam; hp_bits = 256)
    H = build_karlsson_family(theta, phi, lam)
    !is_hadamard(H; tol = 1e-8) && return (ok = false, clique = 0, third = false, hp = false)
    pool, _ = generate_candidate_pool_fresh(H; verbose = false)
    pool = deduplicate_pool(pool; tol = 1e-10)
    ext = check_four_mub_extension(pool, H; mu_tol = MU_TOL, ortho_tol = ORTHO_TOL)
    hp_ok = true
    if ext.found_third && ext.max_clique >= 6
        g = _orthogonality_graph(pool; ortho_tol = ORTHO_TOL)
        cliques = [c for c in maximal_cliques(g) if length(c) >= 6]
        if !isempty(cliques)
            best_c = argmax(c -> length(c), cliques)
            hp = verify_clique_hp([pool[i] for i in best_c]; bits = hp_bits)
            hp_ok = hp.ortho_max < 1e-8 && hp.mu_norm_max < 1e-8
        end
    end
    return (ok = true, clique = ext.max_clique, third = ext.found_third, hp = hp_ok)
end

function count_cliques(cliques)
    n6 = count(>=(6), cliques)
    n2 = count(==(2), cliques)
    n_other = length(cliques) - n6 - n2
    return (n6, n2, n_other)
end

function isotropic_scan(th0, ph0, lm0; r, n = N_ISO, rng = MersenneTwister(42))
    cliques = Int[]
    hp_flags = Bool[]
    for _ in 1:n
        th = th0 + r * randn(rng)
        ph = max(1e-8, ph0 + r * randn(rng))
        lm = lm0 + r * randn(rng)
        try
            p = probe_clique(th, ph, lm)
            p.ok && (push!(cliques, p.clique); push!(hp_flags, p.hp))
        catch e
            @warn "probe failed" exception = e
        end
    end
    n6, n2, no = count_cliques(cliques)
    return (valid = length(cliques), n6, n2, no, cliques, hp_flags)
end

"""Smallest |delta| in MAGS where clique drops below 6 (both +/- if applicable)."""
function axis_sweep(th0, ph0, lm0, axis::Symbol)
    results = NamedTuple[]
    max6 = 0.0
    first_drop = Inf
    for sign in (-1.0, 1.0)
        prev6 = true  # exact anchor has clique 6
        for m in MAGS
            d = sign * m
            th, ph, lam = th0, ph0, lm0
            if axis == :theta
                th = th0 + d
            elseif axis == :phi
                ph = max(1e-8, ph0 + d)
            else
                lam = lm0 + d
            end
            p = probe_clique(th, ph, lam)
            cl = p.ok ? p.clique : 0
            cl >= 6 && (max6 = max(max6, abs(d)))
            if prev6 && cl < 6 && first_drop == Inf
                first_drop = abs(d)
            end
            prev6 = cl >= 6
            push!(results, (axis = axis, sign = sign, delta = d, clique = cl, third = p.third))
        end
    end
    return results, max6, first_drop
end

function pair_sweep(th0, ph0, lm0, pair::Tuple{Symbol, Symbol})
    results = NamedTuple[]
    for m in MAGS
        for s1 in (-1.0, 1.0), s2 in (-1.0, 1.0)
            th, ph, lam = th0, ph0, lm0
            d1, d2 = s1 * m, s2 * m
            a, b = pair
            if a == :theta
                th += d1
            elseif a == :phi
                ph = max(1e-8, ph + d1)
            else
                lam += d1
            end
            if b == :theta
                th += d2
            elseif b == :phi
                ph = max(1e-8, ph + d2)
            else
                lam += d2
            end
            p = probe_clique(th, ph, lam)
            push!(results, (pair = pair, m = m, clique = p.ok ? p.clique : 0, third = p.third))
        end
    end
    persist = maximum(r.m for r in results if r.clique >= 6; init = 0.0)
    return results, persist
end

function random_direction_walk(th0, ph0, lm0, dir::Vector{Float64}; seed = 0)
    dir = dir / norm(dir)
    results = NamedTuple[]
    persist = 0.0
    for m in MAGS
        th = th0 + m * dir[1]
        ph = max(1e-8, ph0 + m * dir[2])
        lam = lm0 + m * dir[3]
        p = probe_clique(th, ph, lam)
        cl = p.ok ? p.clique : 0
        cl >= 6 && (persist = m)
        push!(results, (seed = seed, m = m, clique = cl, third = p.third, dir = dir))
        # also negative direction
        th2 = th0 - m * dir[1]
        ph2 = max(1e-8, ph0 - m * dir[2])
        lam2 = lm0 - m * dir[3]
        p2 = probe_clique(th2, ph2, lam2)
        cl2 = p2.ok ? p2.clique : 0
        cl2 >= 6 && (persist = max(persist, m))
        push!(results, (seed = seed, m = -m, clique = cl2, third = p2.third, dir = dir))
    end
    return results, persist
end

function main()
    rng = MersenneTwister(12345)
    priority_hits = NamedTuple[]

    open(OUT, "w") do io
        println(io, "=== Part A: Third-MUB locus local sampling ===")
        println(io, "Date: $(Dates.format(now(), "yyyy-mm-dd HH:MM"))")
        println(io, "Anchors: Dita, F6_theta0")
        println(io, "Isotropic: $N_ISO Gaussian samples per radius in $RADII")
        println(io, "Tolerances: mu=$MU_TOL ortho=$ORTHO_TOL HP bits=256 on any clique>=6\n")

        # Exact anchors
        println(io, "--- Exact anchor probes ---")
        for (name, th, ph, lam) in ANCHORS
            p = probe_clique(th, ph, lam; hp_bits = 400)
            println(io, "  $name: clique=$(p.clique) third=$(p.third) hp=$(p.hp)")
        end

        # A1 Multi-radius
        println(io, "\n=== A1: Multi-radius isotropic Gaussian ===")
        println(io, "anchor\tradius\tvalid\tn_clique6\tn_clique2\tn_other")
        iso_table = NamedTuple[]
        for (name, th, ph, lam) in ANCHORS
            for r in RADII
                sc = isotropic_scan(th, ph, lam; r = r, rng = MersenneTwister(hash((name, r)) % typemax(Int)))
                @printf(io, "%s\t%.0e\t%d\t%d\t%d\t%d\n",
                        name, r, sc.valid, sc.n6, sc.n2, sc.no)
                push!(iso_table, (anchor = name, r = r, sc...))
                if sc.n6 > 0
                    for (i, c) in enumerate(sc.cliques)
                        c >= 6 && push!(priority_hits, (anchor = name, kind = "isotropic", r = r, clique = c, hp = sc.hp_flags[i]))
                    end
                end
            end
        end

        # A2 Directional
        println(io, "\n=== A2: Directional sweeps ===")
        println(io, "\n--- Pure axis (smallest drop below 6) ---")
        println(io, "anchor\taxis\tmax|delta|_clique6\tfirst_drop_mag")
        for (name, th, ph, lam) in ANCHORS
            for ax in (:theta, :phi, :lambda)
                _, persist, first_drop = axis_sweep(th, ph, lam, ax)
                fd = first_drop == Inf ? ">$(maximum(MAGS))" : @sprintf("%.0e", first_drop)
                @printf(io, "%s\t%s\t%.0e\t%s\n", name, ax, persist, fd)
                persist >= 1e-4 && persist > 0 && push!(priority_hits,
                    (anchor = name, kind = "axis_$ax", r = persist, clique = 6, hp = true))
            end
        end

        println(io, "\n--- Pairwise combinations ---")
        println(io, "anchor\tpair\tmax_m_clique6")
        pairs = [(:theta, :phi), (:theta, :lambda), (:phi, :lambda)]
        for (name, th, ph, lam) in ANCHORS
            for pr in pairs
                _, persist = pair_sweep(th, ph, lam, pr)
                @printf(io, "%s\t%s+%s\t%.0e\n", name, pr[1], pr[2], persist)
                persist >= 1e-4 && persist > 0 && push!(priority_hits,
                    (anchor = name, kind = "pair_$pr", r = persist, clique = 6, hp = true))
            end
        end

        println(io, "\n--- 3 fixed random directions ---")
        dirs = [
            LinearAlgebra.normalize([1.0, 0.3, -0.2]),
            LinearAlgebra.normalize([0.1, 1.0, 0.4]),
            LinearAlgebra.normalize([0.5, -0.3, 1.0]),
        ]
        println(io, "anchor\tseed\tmax_m_clique6\tdir")
        for (name, th, ph, lam) in ANCHORS
            for (si, dir) in enumerate(dirs)
                _, persist = random_direction_walk(th, ph, lam, copy(dir); seed = si)
                @printf(io, "%s\t%d\t%.0e\t[%.3f,%.3f,%.3f]\n", name, si, persist, dir...)
                persist >= 1e-4 && persist > 0 && push!(priority_hits,
                    (anchor = name, kind = "rand_dir_$si", r = persist, clique = 6, hp = true))
            end
        end

        println(io, "\n=== PRIORITY HITS (clique>=6 beyond anchor) ===")
        if isempty(priority_hits)
            println(io, "NONE")
        else
            for h in priority_hits
                println(io, h)
            end
        end
    end
    println("Wrote $OUT")
    if !isempty(priority_hits)
        println("WARNING: $(length(priority_hits)) priority clique>=6 hits — review before paper update")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
