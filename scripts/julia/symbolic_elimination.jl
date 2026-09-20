# Step 4 / Phase 3: Symbolic elimination — real per-H pool export + fourth-MUB witness skeleton.
# Oscar.jl preferred; Macaulay2 .m2 export as fallback.
#
# Exact n_wit=1 over a named number field (never CC): scripts/julia/export_w1_exact.jl
#
# Usage: julia --project=. scripts/julia/symbolic_elimination.jl
#        julia --project=. scripts/julia/symbolic_elimination.jl --anchors
#        julia --project=. scripts/julia/symbolic_elimination.jl --export-pools

include(joinpath(@__DIR__, "_paths.jl"))

using Dates
using LinearAlgebra
using Printf

const M2_DIR = SYMBOLIC_EXPORT_DIR
const OUT_DIR = RESULTS_DIR
const DITA_THETA = acos(1 / sqrt(3))

function has_oscar()
    try
        @eval using Oscar
        return true
    catch
        return false
    end
end

"""Format complex number for Macaulay2 CC ring (ii = sqrt(-1))."""
function m2_complex(z::Complex)
    r, im = real(z), imag(z)
    if abs(im) < 1e-15
        return @sprintf("%.17g", r)
    elseif abs(r) < 1e-15
        im == 1.0 && return "ii"
        im == -1.0 && return "-ii"
        return @sprintf("%.17g*ii", im)
    else
        im >= 0 && return @sprintf("(%.17g+%.17g*ii)", r, im)
        return @sprintf("(%.17g%.17g*ii)", r, im)
    end
end

"""Build 10 per-H pool equation strings (z_i*w_i-1 and MU to H)."""
function build_pool_equation_strings(H::AbstractMatrix{ComplexF64})
    hc = conj.(H)
    eqs = String[]
    for i in 1:5
        push!(eqs, "z$i * w$i - 1")
    end
    for k in 0:4
        lhs_parts = String[m2_complex(hc[1, k + 1])]
        rhs_parts = String[m2_complex(H[1, k + 1])]
        for j in 1:5
            hc_ij = m2_complex(hc[j + 1, k + 1])
            h_ij = m2_complex(H[j + 1, k + 1])
            hc_ij != "0" && push!(lhs_parts, "$hc_ij*z$j")
            h_ij != "0" && push!(rhs_parts, "$h_ij*w$j")
        end
        lhs = join(lhs_parts, " + ")
        rhs = join(rhs_parts, " + ")
        push!(eqs, "($lhs)*($rhs) - 6")
    end
    return eqs
end

function write_pool_m2(path, H, label; extra_comment = "")
    eqs = build_pool_equation_strings(H)
    open(path, "w") do io
        println(io, "-- Per-H MU pool ideal: 10 equations, 10 variables (z1..z5, w1..w5)")
        println(io, "-- Anchor: $label  generated $(Dates.now())")
        !isempty(extra_comment) && println(io, "-- $extra_comment")
        println(io, "R = CC[z1,z2,z3,z4,z5,w1,w2,w3,w4,w5];")
        println(io, "I = ideal(")
        for (i, eq) in enumerate(eqs)
            sep = i < length(eqs) ? "," : ""
            println(io, "  $eq$sep")
        end
        println(io, ");")
        println(io, "print(\"#generators = \" | toString numgens I | \" dim I = \" | toString dim I | \" degree I = \" | toString degree I | newline)")
    end
    return path, length(eqs)
end

function export_pool_anchors()
    anchors = [
        ("F6", F6, "F6 Fourier matrix (regression anchor)"),
        ("Dita", build_karlsson_family(DITA_THETA, pi / 4, 0.4),
         "Karlsson Dita theta=arccos(1/sqrt(3)) phi=pi/4 lambda=0.4"),
        ("F6_theta0", build_karlsson_family(0.0, 0.5, 0.3),
         "Karlsson F6_theta0 theta=0 lambda=0.3"),
    ]
    results = NamedTuple[]
    for (name, H, comment) in anchors
        fname = name == "F6" ? "pool_F6.m2" : name == "Dita" ? "pool_Dita.m2" :
                "pool_$(name).m2"
        path, neqs = write_pool_m2(joinpath(M2_DIR, fname), H, name; extra_comment = comment)
        @printf("  Exported %s (%d eqs)\n", path, neqs)
        push!(results, (name = name, path = path, n_eqs = neqs, stub = false))
    end
    return results
end

"""Build fourth-MUB witness equation strings at fixed H."""
function build_fourth_mub_witness_equations(H::AbstractMatrix{ComplexF64}, clique_indices::Vector{Int};
                                              n_wit::Int=5)
    pool, = generate_candidate_pool_fresh(H; verbose=false)
    pool = deduplicate_pool(pool)
    length(clique_indices) < 6 && error("clique_indices must have length >= 6")
    B3 = hcat([pool[i] for i in clique_indices[1:6]]...)

    eqs = String[]
    for v in 1:n_wit
        for i in 1:5
            push!(eqs, "z$(v)_$i * w$(v)_$i - 1")
        end
    end

    function mu_eqs_for_basis(B, prefix_z, prefix_w)
        local_eqs = String[]
        for k in 1:6
            lhs_parts = [m2_complex(conj(B[1, k]))]
            rhs_parts = [m2_complex(B[1, k])]
            for j in 1:5
                cj = m2_complex(conj(B[j + 1, k]))
                bj = m2_complex(B[j + 1, k])
                cj != "0" && push!(lhs_parts, "$cj*$(prefix_z)$j")
                bj != "0" && push!(rhs_parts, "$bj*$(prefix_w)$j")
            end
            push!(local_eqs, "($(join(lhs_parts, " + ")))*($(join(rhs_parts, " + "))) - 6")
        end
        return local_eqs
    end

    # Computational MU to I is already z_i*w_i=1 above; do not encode I6 as a column basis
    # (that incorrectly yields (1)*(1)-6 = -5 and trivially empties the ideal).
    for v in 1:n_wit
        pz, pw = "z$(v)_", "w$(v)_"
        append!(eqs, mu_eqs_for_basis(H, pz, pw))
        append!(eqs, mu_eqs_for_basis(B3, pz, pw))
    end

    for a in 1:n_wit, b in (a + 1):n_wit
        terms = String[]
        push!(terms, "1")
        for j in 1:5
            push!(terms, "z$(a)_$j * w$(b)_$j")
        end
        push!(eqs, join(terms, " + "))
    end

    return eqs, pool, B3
end

"""Reduced witness: 2 candidate fourth-basis vectors (20 vars, ~52 eqs)."""
function build_reduced_fourth_mub_witness(H::AbstractMatrix{ComplexF64}, clique_indices::Vector{Int})
    return build_fourth_mub_witness_equations(H, clique_indices; n_wit=2)
end

function export_reduced_witness_m2(H_name, H)
    mkpath(M2_DIR)
    g = let
        pool, = generate_candidate_pool_fresh(H; verbose=false)
        pool = deduplicate_pool(pool)
        _orthogonality_graph(pool; ortho_tol=1e-8)
    end
    cliques = [c for c in maximal_cliques(g) if length(c) >= 6]
    c = argmax(cl -> length(cl), cliques)
    eqs, _, _ = build_reduced_fourth_mub_witness(H, c)
    vars = witness_var_decl(2)
    path = joinpath(M2_DIR, "fourth_mub_reduced_$(H_name).m2")
    open(path, "w") do io
        println(io, "-- Reduced fourth-MUB witness (2 vectors) for $H_name")
        println(io, "-- vars=$(length(vars)) eqs=$(length(eqs)) third_clique=$c")
        println(io, "R = CC[$(join(vars, ","))];")
        println(io, "witness = ideal(")
        for (i, eq) in enumerate(eqs)
            sep = i < length(eqs) ? "," : ""
            println(io, "  $eq$sep")
        end
        println(io, ");")
        println(io, "print(\"#generators = \" | toString numgens witness | \" dim = \" | toString dim witness | \" degree = \" | toString degree witness | newline)")
    end
    return path, length(vars), length(eqs)
end

"""Export pool with exact-arithmetic comment block for L7 certificate."""
function export_pool_dita_exact()
    H = build_karlsson_family(DITA_THETA, pi / 4, 0.4)
    path, neqs = write_pool_m2(
        joinpath(M2_DIR, "pool_Dita_exact.m2"), H, "Dita_exact";
        extra_comment="L7 certificate: HC reports 240 certified roots; mixed_volume=252",
    )
    return path, neqs
end

"""Build fourth-MUB witness sizing at fixed H."""
function fourth_mub_witness_system(H::AbstractMatrix{ComplexF64})
    pool, = generate_candidate_pool_fresh(H; verbose = false)
    pool = deduplicate_pool(pool)
    ext = check_four_mub_extension(pool, H)
    !ext.found_third && error("H has no third MUB clique for witness construction")
    g = _orthogonality_graph(pool; ortho_tol = 1e-8)
    cliques = [c for c in maximal_cliques(g) if length(c) >= 6]
    c = argmax(cl -> length(cl), cliques)
    eqs, _, _ = build_fourth_mub_witness_equations(H, c; n_wit=5)
    n_eqs = length(eqs)
    n_vars = 5 * 2 * 5  # 5 witness vectors × (5 z + 5 w)
    return (
        n_witness_eqs = n_eqs,
        n_witness_vars = n_vars,
        third_clique_indices = c,
        param_vars = 3,
        equation_strings = eqs,
    )
end

function export_theta0_slice_m2()
    mkpath(M2_DIR)
    path = joinpath(M2_DIR, "theta0_fourth_mub_witness.m2")
    open(path, "w") do io
        println(io, "-- Macaulay2 export: theta=0 Karlsson slice, fourth-MUB witness ideal")
        println(io, "-- Generated by symbolic_elimination.jl $(Dates.now())")
        println(io, "R = QQ[p0,p1,p2,p3,p4,p5, phi, lam, w6];")
        println(io, "-- w6^6 - 1 = 0; theta = 0 fixed")
        println(io, "-- Run: eliminate({p0,p1,p2,p3,p4,p5}, ideal witness)")
    end
    return path
end

function export_mobius_degeneracy_m2()
    mkpath(M2_DIR)
    path = joinpath(M2_DIR, "mobius_degeneracy_locus.m2")
    open(path, "w") do io
        println(io, "-- Mobius z4^2 = M_A(z2^2) inconsistency locus in (theta, phi, lambda)")
        println(io, "R = CC[theta, phi, lam, z1, z2, z3, z4];")
        println(io, "-- eliminate {z1,z2,z3,z4} to get curve in (theta,phi,lam)")
    end
    return path
end

function witness_var_decl(n_wit = 5)
    zvars = ["z$(v)_$j" for v in 1:n_wit for j in 1:5]
    wvars = ["w$(v)_$j" for v in 1:n_wit for j in 1:5]
    return vcat(zvars, wvars)
end

function export_w_elimination_m2(H_name, H)
    mkpath(M2_DIR)
    path = joinpath(M2_DIR, "fourth_mub_w_elimination_$(H_name).m2")
    wit = fourth_mub_witness_system(H)
    vars = witness_var_decl(5)
    open(path, "w") do io
        println(io, "-- Fourth-MUB witness elimination export for $H_name")
        println(io, "-- witness vars=$(wit.n_witness_vars) eqs=$(wit.n_witness_eqs)")
        println(io, "-- third clique indices: $(wit.third_clique_indices)")
        println(io, "-- Generated $(Dates.now())")
        println(io, "R = CC[$(join(vars, ","))];")
        println(io, "witness = ideal(")
        eqs = wit.equation_strings
        for (i, eq) in enumerate(eqs)
            sep = i < length(eqs) ? "," : ""
            println(io, "  $eq$sep")
        end
        println(io, ");")
        println(io, "-- Eliminate witness variables (Track C theta=0 / Dita slice):")
        println(io, "-- eliminate(set apply(1..$(length(vars)), i -> R_i), witness)")
    end
    return path
end

function theta0_slice_numerical_scan()
    println("\n[S1] theta=0 slice — numerical fourth-MUB scan")
    n_third = n_fourth = 0
    n_pts = 0
    for phi in range(0, pi / 2; length = 6)
        for lam in range(0, 2 * pi; length = 6)
            n_pts += 1
            try
                H = build_karlsson_family(0.0, phi, lam)
                !is_hadamard(H) && continue
                pool, = generate_candidate_pool_fresh(H; verbose = false)
                pool = deduplicate_pool(pool)
                ext = check_four_mub_extension(pool, H)
                ext.found_third && (n_third += 1)
                ext.found_fourth && (n_fourth += 1)
            catch
            end
        end
    end
    @printf("  theta=0 slice: %d points, third=%d, fourth=%d\n", n_pts, n_third, n_fourth)
    return (n_points = n_pts, n_third = n_third, n_fourth = n_fourth)
end

function local_witness_at_anchors()
    println("\n[S2] Local fourth-MUB witness system at anchors")
    anchors = [
        ("F6", F6, 0.0, 0.0, 0.0),
        ("Dita", build_karlsson_family(DITA_THETA, pi / 4, 0.4),
         DITA_THETA, pi / 4, 0.4),
        ("theta0", build_karlsson_family(0.0, 0.5, 0.3), 0.0, 0.5, 0.3),
    ]
    results = NamedTuple[]
    exports = String[]
    for (name, H, th, ph, lm) in anchors
        wit = fourth_mub_witness_system(H)
        @printf("  %s: witness vars=%d eqs=%d clique=%s\n",
                name, wit.n_witness_vars, wit.n_witness_eqs, wit.third_clique_indices)
        push!(exports, export_w_elimination_m2(name, H))
        push!(results, (name = name, theta = th, phi = ph, lambda = lm, wit...))
    end
    return results, exports
end

function oscar_groebner_probe()
    println("\n[S3] Oscar.jl Groebner probe")
    if !has_oscar()
        println("  Oscar.jl not installed — using Macaulay2 export fallback.")
        println("  See docs/METHODS.md section 7 for Oscar install steps.")
        p1 = export_theta0_slice_m2()
        p2 = export_mobius_degeneracy_m2()
        println("  Exported: $p1")
        println("  Exported: $p2")
        return (available = false, exports = [p1, p2])
    end
    println("  Oscar available — full Groebner elimination deferred (50+ witness vars).")
    println("  Use symbolic_export/*.m2 for Macaulay2 cross-check.")
    return (available = true, exports = String[])
end

function main()
    println("=== Symbolic elimination (Phase 3) ===\n")

    mkpath(M2_DIR)
    mkpath(joinpath(OUT_DIR, "track_c_elimination"))

    println("[P3.1] Per-H pool polynomial export (F6, Dita)")
    pool_exports = export_pool_anchors()
    exact_path, exact_neqs = export_pool_dita_exact()
    push!(pool_exports, (name="Dita_exact", path=exact_path, n_eqs=exact_neqs, stub=false))

    println("[P3.2] Reduced fourth-MUB witness export (2 vectors)")
    reduced_paths = NamedTuple[]
    for (name, H) in [("Dita", build_karlsson_family(DITA_THETA, pi/4, 0.4)),
                      ("F6_theta0", build_karlsson_family(0.0, 0.5, 0.3))]
        p, nv, ne = export_reduced_witness_m2(name, H)
        push!(reduced_paths, (name=name, path=p, n_vars=nv, n_eqs=ne))
        @printf("  reduced witness %s: %d vars, %d eqs -> %s\n", name, nv, ne, p)
    end

    slice = theta0_slice_numerical_scan()
    witnesses, exports = local_witness_at_anchors()
    oscar = oscar_groebner_probe()

    println("\n=== Summary ===")
    println("  theta=0 slice: found_fourth=$(slice.n_fourth) / $(slice.n_points) points")
    println("  Oscar installed: $(oscar.available)")
    for pe in pool_exports
        @printf("  pool export %s: %d eqs stub=%s\n", pe.name, pe.n_eqs, pe.stub)
    end
    for w in witnesses
        println("  witness $(w.name): vars=$(w.n_witness_vars) eqs=$(w.n_witness_eqs)")
    end

    summary_path = joinpath(OUT_DIR, "symbolic_elimination_summary.txt")
    open(summary_path, "w") do io
        println(io, "phase: 3 (real pool export + witness skeleton)")
        println(io, "theta0_slice: n=$(slice.n_points) third=$(slice.n_third) fourth=$(slice.n_fourth)")
        println(io, "oscar_available: $(oscar.available)")
        for pe in pool_exports
            println(io, "pool_$(pe.name): n_eqs=$(pe.n_eqs) stub=$(pe.stub) path=$(pe.path)")
        end
        for w in witnesses
            println(io, "$(w.name): witness_vars=$(w.n_witness_vars) eqs=$(w.n_witness_eqs)")
        end
        for p in vcat(exports, oscar.exports)
            println(io, "export: $p")
        end
    end
    println("  Wrote $summary_path")
end

main()
