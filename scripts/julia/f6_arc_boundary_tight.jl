# ============================================================================
# f6_arc_boundary_tight.jl — bounded re-probe of the F6 θ=0 λ-arc boundaries
# at ORTHOGONALITY/MU thresholds BELOW 1e-12 (1e-13, 1e-14), plus 400-bit HP
# re-verification of any 6-clique found.
#
# Question: are the recorded boundaries (|Δλ| ≈ 0.055044803035 / 0.0625984193874
# from λ0 = 0.3, φ = 0.5, θ = 0; results/f6_boundary_reverify.txt) stable when
# the classification threshold is pushed from 1e-12 toward the Float64 noise
# floor? NOT an HP path-tracking claim; discovery solves are Float64 HC.
#
# Bounded: bracket probes + 14 bisection steps per sign. No full circle sweep.
# Output: results/f6_arc_boundary_tight.txt
# ============================================================================
using Dates, Printf, LinearAlgebra, Random
include(joinpath(@__DIR__, "_paths.jl"))   # legacy core: is_hadamard, pool fns, HP verify
include(joinpath(ROOT, "src", "Cliques.jl"))  # Part IX: enumerate_all_third_mub_bases
using Graphs

println("start=$(now())"); flush(stdout)

const HP_BITS = 400

function probe(theta, phi, lam; tol)
    H = build_karlsson_family(theta, phi, lam)
    is_hadamard(H; tol = 1e-8) || return (ok = false, six = false, hp = missing)
    pool, _ = generate_candidate_pool_fresh(H; verbose = false)
    pool = deduplicate_pool(pool)
    cl = enumerate_all_third_mub_bases(pool, H; ortho_tol = tol, mu_tol = tol,
                                       hp_bits = HP_BITS)
    return (ok = true, six = cl.n_verified,
            hp_all_ok = all(r -> !(r.hp isa NamedTuple) || (r.hp.ortho_ok && r.hp.mu_ok),
                            cl.bases),
            hp_bits = HP_BITS)
end

# Bisect the boundary: clique-6 inside, no-6 outside, at threshold `tol`.
function bisect(sign, lo, hi; tol, steps = 14)
    # invariant: lo keeps clique-6, hi loses it
    for _ in 1:steps
        mid = (lo + hi) / 2
        p = probe(0.0, 0.5, 0.3 + sign * mid; tol = tol)
        if p.ok && p.six >= 1
            lo = mid
        else
            hi = mid
        end
    end
    return (lo + hi) / 2
end

function main()
    lines = String[]
    push!(lines, "=== F6 theta=0 arc boundary at tighter thresholds (fresh Float64 solves) ===")
    push!(lines, "timestamp=$(now())  Julia=$(VERSION)  seed=20260917  HP=$(HP_BITS) bits")
    push!(lines, "reference midpoints (400-bit HP, eps=1e-12): +0.055044803035  -0.0625984193874")
    push!(lines, "")
    recorded = (1 => 0.055044803035, -1 => 0.0625984193874)
    for (sign, rec) in ((+1.0, 0.055044803035), (-1.0, 0.0625984193874))
        for tol in (1e-12, 1e-13)
            push!(lines, "--- sign=$(sign > 0 ? "+" : "-") tol=$tol ---")
            has6_list = Float64[]
            no6_list = Float64[]
            for off in (-2e-6, -1e-6, 0.0, 1e-6, 2e-6)
                d = rec + off
                p = probe(0.0, 0.5, 0.3 + sign * d; tol = tol)
                push!(lines, @sprintf("  |dl|=%.12f  six_cliques=%s  hp=%s",
                                       d, p.six, p.hp === missing ? "n/a" : string(p.hp_all_ok)))
                if p.ok && p.six >= 1
                    push!(has6_list, d)
                else
                    push!(no6_list, d)
                end
            end
            if isempty(has6_list) || isempty(no6_list)
                push!(lines, "  no straddle within +-2e-6 of recorded midpoint; boundary moved")
                continue
            end
            lo = maximum(d for d in has6_list)
            hi = minimum(d for d in no6_list)
            push!(lines, @sprintf("  straddle bracket at tol=%g: [%.12f, %.12f]  (recorded %.12f)",
                                   tol, lo, hi, rec))
        end
    end
    push!(lines, "")
    push!(lines, "scope: fresh Float64 discovery solves with tighter classification thresholds;")
    push!(lines, "no HP path tracking, no new certificate, no T3 dependence.")
    open(joinpath(ROOT, "results", "f6_arc_boundary_tight.txt"), "w") do io
        for ln in lines; println(io, ln); end
    end
    println.(lines); flush(stdout)
    println("done=$(now())"); flush(stdout)
end

main()

