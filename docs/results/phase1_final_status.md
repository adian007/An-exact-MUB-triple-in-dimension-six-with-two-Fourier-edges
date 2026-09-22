# Phase 1 Definitions Suite — Final Status

## Run Summary

| Metric | Value |
|---|---|
| **Exit code** | 0 (success) |
| **Passed** | 45 |
| **Failed** | 0 |
| **Total checks** | 45 |
| **Runtime** | ~6.9–11.2 s (varies by machine load; SymPy 1.14.0, Python 3.13.14, Windows 11) |
| **Suite file** | `research/definitions/phase1_definitions_suite.py` |
| **Report file** | `research/definitions/phase1_definitions_suite.out.txt` |

## Verification Performed (finalization only — no test logic altered)

1. **Dead code after `return` in `main()` — NONE found.** Line 819 is
   `return 0 if n_fail == 0 else 1`. The only following lines are 822–823:
   `if __name__ == "__main__":` / `raise SystemExit(main())`, which is the
   standard module entry point (unreachable *inside* `main()` but is the sole
   invocation path for the script — not dead code). No statements exist between
   line 819 and line 822. The file was already clean; no deletion was needed.

2. **T2.2 / T2.3 `closure_ok` bool conversion — verified correct, no change:**
   - Line 416: `closure_ok = bool(chm_unit(DM))` — `chm_unit(H, n=6)` (line 230)
     returns `is_zero_exact(H * dag(H) - n * sp.eye(n))`, which is annotated
     `-> bool` and returns a genuine Python `bool`. The `bool()` wrapper is a
     safe no-op that guarantees a strict `bool` is passed to `record()`.
   - Line 420–426 (T2.3 record call): uses `closure_ok` as the pass/fail boolean
     consistently with the `record(tag, ok, note)` signature. The note string
     honestly reports that modulus multiset preservation is *observed* (not
     assumed) — it reports `False` for the concrete row-phase choice on D_bc,
     exactly as the research plan requires.

## T6 Notes — Honesty Confirmation

All three T6 notes report **non-equivalence** findings honestly:

- **T6.1**: `dephased H_Dita(i) != dephased D_bc` — records that H_Dita(i) is
  NOT monomial-equivalent to D_bc (different root-of-unity families: 4th roots
  vs 3rd roots). Marked "FINDING (expected)".
- **T6.2**: `dephased B3(i) != dephased F_D` — records that B3(i) is NOT
  monomial-equivalent to F_D (row of all 1s vs sqrt(10)*(1-3i)/10 entries).
  Marked "FINDING (expected)".
- **T6.3**: Consequence record — H_Dita(i) is NOT CHM-equivalent to D_bc and
  B3(i) is NOT CHM-equivalent to F_D; they are independent candidate anchors.
  Marked "CONSEQUENCE (recorded, not assumed)".

No test logic or pass/fail criteria were changed. The non-equivalence findings
are the intended results and are reported as such in the notes.

## Remaining Issues

None. The suite is clean:
- No dead code.
- No logic changes.
- All 45 checks pass with exit code 0.
- T6 notes are honest about non-equivalence.
- The run-directory copy was created automatically at
  `runs/2026-09-17T170232Z_phase1_definitions/`.