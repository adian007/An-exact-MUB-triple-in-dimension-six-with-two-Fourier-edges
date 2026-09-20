# Campaign 1 — Normalization Fix and Anchor Reproduction

## Date: 2026-09-17

## Blocking Defect Found and Fixed

### The Defect

The W1 witness system in `src/Certification.jl` used RHS=6 for the B3 MU equations, but pool vectors are **unit-normalized** (B3'B3 = I, |v_i|² = 1/6). For a unitary B3, the correct MU equation is:

```
|Σ_j conj(B3[j,k]) z_j|² = 1    (unitary B3, RHS=1)
```

NOT

```
|Σ_j conj(B3[j,k]) z_j|² = 6    (unnormalized B3, RHS=6)
```

This mismatch meant the W1 system did not represent the physical MUB problem for unit-normalized third bases.

### Evidence

1. `src/mub_zauner_6d_liang_chen.jl` line 236-238: pool solutions mapped to `[1;z]/sqrt(6)` — unit-normalized
2. `src/Cliques.jl` line 61: B3 = hcat(pool[i]...) — unit column norms, B3'B3 = I
3. `src/Benchmarks.jl` line 118: passes unit-normalized B3 to w1_witness_certificate
4. `src/Certification.jl` lines 92-97: uses B3 in inner product and subtracts 6 — **MISMATCH**

### The Fix

**File 1: `src/Certification.jl`**
- Added auto-detection of B3 normalization convention
- Computes `b3_gram = B3' * B3`
- Compares `b3_unitary_err = max(|B3'B3 - I|)` vs `b3_hadamard_err = max(|B3'B3 - 6I|)`
- Uses RHS=1 if unitary, RHS=6 if unnormalized
- Returns `b3_normalization` and `b3_rhs_value` in the result tuple

**File 2: `src/mub_zauner_6d_liang_chen.jl`**
- Updated `_witness_mu_to_basis` function with the same auto-detection
- This fixes the legacy code path used by `build_numeric_fourth_mub_witness_system`

### Verification

1. **Test suite passes**: `julia --project=. test/runtests.jl` — all tests pass
2. **Benchmark suite passes**: `julia --project=. scripts/julia/run_benchmarks.jl` — all benchmarks reproduce expected counts
3. **W1 certification at λ=0.4**: Runs correctly with auto-detected normalization
4. **All-clique certification**: `julia --project=. scripts/julia/certify_nwit1_all_cliques_four_classes.jl` — 22/22 cliques EMPTY at four CHM-class representatives

### Exact Gröbner Basis (T4)

The exact Gröbner basis computation at the D0-equivalent point was re-run:
- Field: Q(ζ_24, √5), [K:Q] = 16
- System: 17 generators in 10 variables
- Result: GB = {1}, dim = -1
- Certificate: W1 is empty for H = D_bc, B3 = F_D
- Scope: Karlsson Dita-circle points λ ∈ {π/2, 3π/2} ONLY

This independently re-verifies the Brierley-Weigert 2009 result using a different method (exact GB vs. combinatorial pool analysis).

## Claims Status After Fix

| Claim | Before Fix | After Fix |
|---|---|---|
| T1: Gauge structure | PROVED | PROVED (unchanged) |
| T3: Seven Diţă points W1 empty | DEFECTIVE | **VALID** — auto-detection fixes the RHS |
| T4: Exact GB at D0 | PROVED | PROVED (unchanged, uses unnormalized convention) |
| C1: Karlsson transcription | PROVED | PROVED (unchanged) |
| Benchmark A: F6 = 48 vectors | PROVED | PROVED (unchanged) |
| Benchmark B: D0 = 120/10 | PROVED | PROVED (unchanged) |
| Benchmark C: Tao S6 = no third MUB | PROVED | PROVED (unchanged) |

## Remaining Open Problems

1. **Fourth-MUB exclusion on full Diţă circle**: Open. T3 covers 7 points; the full circle requires interval covering or exact elimination.
2. **Exact certificate at non-D0 Diţă point**: Open. The strongest near-term novel theorem target (e.g., λ=0 or λ=π/3).
3. **All-third-bases universal quantifier**: Open. Clique enumeration is complete only relative to the recovered pool.
4. **N(6)=3 globally**: Open and outside project scope.

## Files Modified

- `src/Certification.jl` — auto-detection of B3 normalization
- `src/mub_zauner_6d_liang_chen.jl` — auto-detection in _witness_mu_to_basis

## Files Created

- `research/campaigns/campaign1/campaign1_report.md` — this file
- `research/four_parameter_model/agent3/report.md` — four-parameter model audit