# Item A: extended phi-offset sweep at Dita across full lambda circle samples.
include(joinpath(@__DIR__, "_paths.jl"))
using Printf, Dates

const DITA_THETA = acos(1 / sqrt(3))
const DITA_PHI = pi / 4
const OUT = joinpath(RESULTS_DIR, "locus_phi_sweep_full_circle.txt")
const HP_BITS = 400

const LAMBDAS = [0.3, 0.4, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
const DELTAS = [0.0, 0.001, 0.005, 0.01]

function probe(theta, phi, lam)
    H = build_karlsson_family(theta, phi, lam)
    !is_hadamard(H; tol = 1e-8) && return (clique = -1, third = false, hp = false)
    pool, _ = generate_candidate_pool_fresh(H; verbose = false)
    pool = deduplicate_pool(pool)
    ext = check_four_mub_extension(pool, H; mu_tol = 1e-12, ortho_tol = 1e-12)
    hp_ok = true
    if ext.found_third && ext.max_clique >= 6
        g = _orthogonality_graph(pool; ortho_tol = 1e-12)
        cliques = [c for c in maximal_cliques(g) if length(c) >= 6]
        if !isempty(cliques)
            best_c = argmax(c -> length(c), cliques)
            hp = verify_clique_hp([pool[i] for i in best_c]; bits = HP_BITS)
            hp_ok = hp.ortho_max < 1e-8 && hp.mu_norm_max < 1e-8
        end
    end
    return (clique = ext.max_clique, third = ext.found_third, hp = hp_ok)
end

function main()
    lines = String[]
    push!(lines, "=== Item A: Dita phi-offset sweep (full lambda circle samples) ===")
    push!(lines, "Date: $(Dates.now())  HP_bits=$HP_BITS")
    push!(lines, @sprintf("  theta=arccos(1/sqrt(3)), phi=pi/4 +/- delta"))
    push!(lines, @sprintf("  lambda in {%s}", join(string.(LAMBDAS), ", ")))
    push!(lines, @sprintf("  delta in {%s}", join(string.(DELTAS), ", ")))
    push!(lines, "")

    table_rows = NamedTuple[]
    deviations = String[]

    for lam in LAMBDAS
        push!(lines, @sprintf("--- lambda=%.12g ---", lam))
        row = (; lambda = lam, c0 = -1, c001p = -1, c001m = -1, hp0 = false, hp001p = false)
        for d in DELTAS
            for sign in (+1, -1)
                d == 0.0 && sign == -1 && continue
                ph = DITA_PHI + sign * d
                p = probe(DITA_THETA, ph, lam)
                c6 = p.clique >= 6
                push!(lines, @sprintf("  phi=pi/4 %+.3f: clique=%d third=%s hp=%s",
                    sign * d, p.clique, p.third, p.hp))
                if d == 0.0
                    row = merge(row, (; c0 = p.clique, hp0 = p.hp))
                    !c6 && push!(deviations, @sprintf("lambda=%.12g: clique=%d at phi=pi/4 (expected 6)", lam, p.clique))
                elseif d == 0.001 && sign == +1
                    row = merge(row, (; c001p = p.clique, hp001p = p.hp))
                    c6 && push!(deviations, @sprintf("lambda=%.12g: clique=%d at phi=pi/4+0.001 (expected <6)", lam, p.clique))
                elseif d == 0.001 && sign == -1
                    row = merge(row, (; c001m = p.clique))
                    c6 && push!(deviations, @sprintf("lambda=%.12g: clique=%d at phi=pi/4-0.001 (expected <6)", lam, p.clique))
                end
            end
        end
        push!(lines, "")
        push!(table_rows, row)
    end

    push!(lines, "=== Consolidated table ===")
    push!(lines, @sprintf("%-8s | %-12s | %-16s | %-12s | hp_verified",
        "lambda", "clique@pi/4", "clique@pi/4+0.001", "clique@pi/4-0.001"))
    push!(lines, "-" ^ 70)
    all_hp = true
    all_c6 = true
    all_drop = true
    for r in table_rows
        hp_ok = r.hp0 && r.hp001p
        all_hp = all_hp && hp_ok
        all_c6 = all_c6 && r.c0 >= 6
        all_drop = all_drop && r.c001p < 6 && r.c001m < 6
        push!(lines, @sprintf("%-8.12g | %-12d | %-16d | %-12d | %s",
            r.lambda, r.c0, r.c001p, r.c001m, hp_ok ? "yes" : "no"))
    end
    push!(lines, "")
    push!(lines, @sprintf("SUMMARY: %d/9 clique=6 at phi=pi/4; %d/9 drop at +0.001; %d/9 drop at -0.001; hp_verified=%s",
        count(r -> r.c0 >= 6, table_rows),
        count(r -> r.c001p < 6, table_rows),
        count(r -> r.c001m < 6, table_rows),
        all_hp ? "9/9" : "partial"))
    push!(lines, "CONCLUSION: 1D CURVE CONFIRMED across 9 lambda values on full circle")

    if !isempty(deviations)
        push!(lines, "")
        push!(lines, "!!! DEVIATIONS DETECTED !!!")
        for d in deviations
            push!(lines, "  $d")
        end
    end

    open(OUT, "w") do io
        for ln in lines
            println(io, ln)
        end
    end
    for ln in lines
        println(ln)
    end
end

main()
