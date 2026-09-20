# E4: exact-enumeration probe at D0 (BW style), not Groebner.
#
# (A) φ_D flat-vector search against textbook D(0) and against D_bc.
# (B) Seeded exact third basis F_D for {I, D_bc} (Bengtsson eqs. 78–80).
# (C) Pairwise MU among any distinct 6-cliques found.
#
# Groebner is only the fallback if this enumeration fails to decide.
# This script does NOT run Groebner / Macaulay2.
#
# Usage: julia --project=. --compiled-modules=no scripts/julia/enumerate_d0_cliques.jl
# Output: results/enumerate_d0_cliques.txt

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(ROOT, "src", "brierley_weigert_notes.jl"))

using LinearAlgebra
using Printf
using Dates

const RESULTS_DIR = joinpath(ROOT, "results")
const OUT = joinpath(RESULTS_DIR, "enumerate_d0_cliques.txt")
const MU_TOL = 1e-10
const ORTHO_TOL = 1e-10

function mu_to_H_ok(z, H; tol=MU_TOL)
    for k in 1:6
        s = zero(ComplexF64)
        @inbounds for j in 1:6
            s += conj(H[j, k]) * z[j]
        end
        abs(abs2(s) - 6) < tol || return false
    end
    return true
end

ortho_ok(u, v; tol=ORTHO_TOL) = abs(dot(u, v)) < tol
mu_vec_ok(u, v; tol=MU_TOL) = abs(abs2(dot(u, v)) - 1 / 6) < tol

function list_cliques6(adj)
    n = size(adj, 1)
    cliques = Vector{Vector{Int}}()
    function rec(R, P)
        length(R) == 6 && (push!(cliques, copy(R)); return)
        length(R) + length(P) < 6 && return
        for (idx, v) in enumerate(P)
            newP = Int[u for u in P[(idx + 1):end] if adj[v, u]]
            push!(R, v)
            rec(R, newP)
            pop!(R)
        end
    end
    rec(Int[], collect(1:n))
    return cliques
end

function bases_pairwise_mu(cliques, vecs)
    nB = length(cliques)
    pairs = mu_pairs = 0
    for a in 1:nB, b in (a + 1):nB
        pairs += 1
        ok = true
        for i in cliques[a], j in cliques[b]
            if !mu_vec_ok(vecs[i], vecs[j])
                ok = false
                break
            end
        end
        ok && (mu_pairs += 1)
    end
    return pairs, mu_pairs
end

function vec_key(v)
    return (round.(real.(v); digits=10), round.(imag.(v); digits=10))
end

"""Dephased unique φ_D pool MU to H (all 6 components in φ_D)."""
function phi_D_pool(H, phases)
    nph = length(phases)
    z = ones(ComplexF64, 6)
    pool = Vector{Vector{ComplexF64}}()
    n_tested = 0
    @inbounds for a in 1:nph, b in 1:nph, c in 1:nph, d in 1:nph, e in 1:nph, f in 1:nph
        n_tested += 1
        z[1] = phases[a]; z[2] = phases[b]; z[3] = phases[c]
        z[4] = phases[d]; z[5] = phases[e]; z[6] = phases[f]
        mu_to_H_ok(z, H) || continue
        push!(pool, copy(z .* conj(z[1])))
    end
    seen = Set{Any}()
    uniq = Vector{Vector{ComplexF64}}()
    for v in pool
        k = vec_key(v)
        if k ∉ seen
            push!(seen, k)
            push!(uniq, v)
        end
    end
    return n_tested, uniq
end

function clique_report(pool_z, lines; label="")
    vecs = [v ./ sqrt(6) for v in pool_z]
    n = length(vecs)
    adj = falses(n, n)
    for i in 1:n, j in (i + 1):n
        ortho_ok(vecs[i], vecs[j]) && (adj[i, j] = adj[j, i] = true)
    end
    n_edges = count(adj) ÷ 2
    push!(lines, "  orthogonality edges = $n_edges")
    if n_edges > 8000
        push!(lines, "  ERROR: graph too dense; skip clique listing.")
        return (n_pool=n, n_cliques=0, pairs=0, mu_pairs=0, cliques=Vector{Vector{Int}}())
    end
    raw = list_cliques6(adj)
    sets = Set(Set(c) for c in raw)
    unique_cliques = [sort(collect(s)) for s in sets]
    pairs, mu_pairs = bases_pairwise_mu(unique_cliques, vecs)
    push!(lines, "  6-cliques (distinct sets) = $(length(unique_cliques))")
    push!(lines, "  pairs of distinct third bases = $pairs")
    push!(lines, "  pairs that are pairwise MU = $mu_pairs")
    return (n_pool=n, n_cliques=length(unique_cliques), pairs=pairs, mu_pairs=mu_pairs,
            cliques=unique_cliques, vecs=vecs)
end

"""Verify an unnormalised CHM as third MUB for {I, H}."""
function verify_third_basis(B, H, lines; name="B3")
    dB = chm_defect(B)
    dMU = mu_defect_chm(B, H)
    # ONB of columns (up to 1/√6 scale): B'B = 6 I for a CHM
    cols = [B[:, j] / sqrt(6) for j in 1:6]
    ortho = maximum(abs(dot(cols[i], cols[j])) for i in 1:6 for j in (i + 1):6)
    muI = maximum(abs(abs2(cols[j][k]) - 1 / 6) for j in 1:6 for k in 1:6)
    ok = dB.ok && dMU < MU_TOL && ortho < ORTHO_TOL && muI < MU_TOL
    push!(lines, @sprintf(
        "  %-12s CHM=%s  MU-to-H=%.3e  col-ortho=%.3e  MU-to-I=%.3e  ok=%s",
        name, dB.ok, dMU, ortho, muI, ok))
    return ok, cols
end

function main()
    t0 = time()
    D0 = dita_D0()
    Dbc = dita_D_bc()
    FD = third_mub_F_D_at_D0()
    phases = phase_set_phi_D()
    nph = length(phases)
    lines = String[]
    push!(lines, "=== E4: exact-enumeration probe at D0 (BW style, not Groebner) ===")
    push!(lines, "Date: $(Dates.format(now(), "yyyy-mm-ddTHH:MM:SS"))")
    push!(lines, "Phases: |φ_D| = $nph  (24th roots + ±α, tan α = 2).")
    push!(lines, "Certificate style: enumerate pool, form 6-cliques, test pairwise MU.")
    push!(lines, "Groebner is NOT run.")
    push!(lines, "")

    # --- (B) seeded exact third basis for D_bc ---
    push!(lines, "--- Seeded exact third basis F_D for {I, D_bc} (Bengtsson eq. 78) ---")
    ok_FD, FD_cols = verify_third_basis(FD, Dbc, lines; name="F_D")
    ok_Dbc = chm_defect(Dbc).ok
    push!(lines, "  D_bc CHM ok=$ok_Dbc")
    push!(lines, "")

    # --- (A1) φ_D × textbook D0 ---
    push!(lines, "--- (A1) φ_D pool vs textbook D0 = D(0) ---")
    n_tested, pool_D0 = phi_D_pool(D0, phases)
    push!(lines, "  tested = $n_tested")
    push!(lines, "  MU pool size = $(length(pool_D0))  (BW reports 120)")
    r0 = clique_report(pool_D0, lines; label="D0")
    push!(lines, "")

    # --- (A2) φ_D × D_bc, optionally seeding F_D columns ---
    push!(lines, "--- (A2) φ_D pool vs D_bc (block-circulant ≈ D0) ---")
    _, pool_Dbc = phi_D_pool(Dbc, phases)
    # Seed with dephased F_D columns (exact third basis)
    for j in 1:6
        z = FD[:, j]                    # unimodular CHM column
        zd = z .* conj(z[1])
        k = vec_key(zd)
        if !any(vec_key(v) == k for v in pool_Dbc)
            push!(pool_Dbc, zd)
        end
    end
    push!(lines, "  MU pool size after seeding F_D columns = $(length(pool_Dbc))")
    rbc = clique_report(pool_Dbc, lines; label="D_bc")
    # Check that F_D columns appear as a clique in the seeded pool
    fd_idx = Int[]
    for j in 1:6
        zd = FD[:, j] .* conj(FD[1, j])
        hit = findfirst(v -> vec_key(v) == vec_key(zd), pool_Dbc)
        hit === nothing || push!(fd_idx, hit)
    end
    push!(lines, "  F_D column indices in pool = $fd_idx")
    if length(fd_idx) == 6
        fd_set = Set(fd_idx)
        found = any(Set(c) == fd_set for c in rbc.cliques)
        push!(lines, "  F_D recovered as a 6-clique among pool = $found")
    end
    push!(lines, "")

    # --- Verdict ---
    push!(lines, "=== VERDICT ===")
    if ok_FD
        push!(lines, "Exact third basis exists: F_D is MU to {I, D_bc} (pipeline check PASS).")
    else
        push!(lines, "FAIL: F_D is not a verified third basis for D_bc.")
    end
    if rbc.n_cliques >= 1 && rbc.mu_pairs == 0
        push!(lines, "Among $(rbc.n_cliques) third bases from the D_bc φ_D+seed pool, none are pairwise MU.")
        if rbc.n_cliques == 10 && rbc.n_pool == 120
            push!(lines, "Counts match BW 2009 (120/10) — corollary, not a new theorem.")
        else
            push!(lines, "Counts ($(rbc.n_pool) vectors, $(rbc.n_cliques) cliques) differ from BW 120/10;")
            push!(lines, "not claimed complete. Groebner was NOT run.")
        end
    elseif rbc.mu_pairs > 0
        push!(lines, "UNEXPECTED: $(rbc.mu_pairs) pairwise-MU third-basis pairs at D_bc — check tolerances.")
    else
        push!(lines, "D_bc φ_D+seed pool: $(rbc.n_pool) vectors, $(rbc.n_cliques) six-cliques.")
        push!(lines, "Textbook D0 φ_D probe: $(r0.n_pool) vectors, $(r0.n_cliques) six-cliques.")
        push!(lines, "Enumeration incomplete vs BW 120/10; seeded F_D still certifies one third basis.")
        push!(lines, "Groebner was NOT run (preferred only if enumeration fails to decide).")
    end
    push!(lines, @sprintf("elapsed = %.1f s", time() - t0))
    push!(lines, "No claim on Zauner, all of K6^(3), or dim I=-1 over CC.")

    mkpath(RESULTS_DIR)
    open(OUT, "w") do io
        println.(Ref(io), lines)
    end
    foreach(println, lines)
    println("Wrote $OUT")
    ok_FD || error("E4 seeded F_D verification failed")
end

main()
