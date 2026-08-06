# Phase 3 validation for the per-H pipeline (Step 3 — must pass before sweep).
#
# Usage:  julia --project=. validate_perH_pipeline.jl

include(joinpath(@__DIR__, "_paths.jl"))

function check(label, ok, detail="")
    status = ok ? "PASS" : "FAIL"
    println("  [$status] $label" * (detail == "" ? "" : " — $detail"))
    return ok
end

function main()
    println("=== Per-H pipeline validation ===")
    println("Question tested: vectors MU to {I, H}; then clique search for 3rd/4th MUB.\n")

    all_ok = true

    println("[V1] F6 regression anchor (156 raw / 48 pool)")
    pool_f6, stats_f6 = generate_candidate_pool_fresh(F6; verbose=true)
    all_ok &= check("raw = 156", stats_f6.n_raw == 156, "got $(stats_f6.n_raw)")
    all_ok &= check("pool = 48", length(deduplicate_pool(pool_f6)) == 48,
                    "got $(length(pool_f6))")

    println("\n[V2] Homotopy + fresh cross-check (drop detection)")
    tracker = init_parametric_pool_tracker(verbose=true)
    cases = [
        ("F6", F6, 48),
        ("theta=0", build_karlsson_family(0.0, 0.5, 0.3), 48),
        ("generic", build_karlsson_family(0.3, 0.5, 0.2), 48),
        ("Dita", build_karlsson_family(acos(1 / sqrt(3)), pi / 4, 0.4), nothing),
    ]
    for (name, H, expected) in cases
        pool, _, meta = generate_candidate_pool(H; tracker=tracker, cross_check=true, verbose=true)
        if expected !== nothing
            all_ok &= check("$name pool size", length(pool) == expected, "got $(length(pool))")
        else
            all_ok &= check("$name pool in {48,72,120}", length(pool) in (48, 72, 120),
                            "got $(length(pool))")
        end
        cmp = meta.cross_check
        if meta.drop_suspected && cmp.only_in_b > 0
            all_ok &= check("$name fresh fallback restored pool", length(pool) >= cmp.matched,
                            "missed $(cmp.only_in_b)")
        else
            check("$name tracking", true, "no drop or fresh agreed")
        end
    end

    println("\n[V3] Extension semantics at generic Karlsson point")
    H = build_karlsson_family(0.3, 0.5, 0.2)
    pool, _ = generate_candidate_pool_fresh(H)
    ext = check_four_mub_extension(pool, H)
    all_ok &= check("pool > 0 (MU vectors exist)", ext.n_pool > 0, "n=$(ext.n_pool)")
    all_ok &= check("no third MUB (max_clique < 6)", !ext.found_third,
                    "max_clique=$(ext.max_clique)")
    all_ok &= check("no fourth MUB", !ext.found_fourth)
    println("  (Literature: no closed-form count for generic H; Conjecture 8.1 evidence")
    println("   suggests Karlsson triples are rare/absent — clique failure is expected.)")

    println("\n=== Verdict ===")
    if all_ok
        println("  ALL CHECKS PASSED")
    else
        println("  FAILURES — fix before sweeping")
        exit(1)
    end
end

main()
