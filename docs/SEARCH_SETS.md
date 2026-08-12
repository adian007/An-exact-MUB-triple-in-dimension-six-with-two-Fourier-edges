# Search sets for Karlsson \(K_6^{(3)}\) (Tracks C and D)

Formal definitions used by `search_special_loci.jl`, `csv_reconciliation.jl`, and the Track C impossibility claim on a **covering region** \(\mathcal{R}\).

---

## Notation

| Symbol | Meaning |
|--------|---------|
| \(K_6^{(3)}\) | Karlsson three-parameter CHM family \(H(\theta,\phi,\lambda)\) |
| Primary queue | Anchors + degeneracy candidates before adaptive refinement |
| `pool_complete_flag` | `n_certified == n_tracked` and `n_verified > 0` (corrected gate) |
| HP clique-6 | `max_clique ≥ 6` with 256-bit (or 400-bit) clique re-verification |

---

## \(\mathcal{S}_{588}\) — certified search with degen cap 80

**Definition.** The **588-row target** from the original `search_special_loci.jl` design:

\[
\mathcal{S}_{588} = \mathcal{P}_{\mathrm{primary}}^{(80)} \;\cup\; \mathcal{R}_{\mathrm{ref}}^{(1)}
\]

where:

- **Anchors** (20 points): `F6_cyclic`, `F6_theta0`, `Dita`, `Mobius_audit_fail`, five `theta_*` slices at \(\phi=\pi/4\), \(\lambda=\pi/6\), and twelve `lambda_kpi6_k=*` (\(k=0..11\)).
- **Degeneracy load:** unique \((\theta,\phi,\lambda)\) from `results/degeneracy_candidates.csv`, **capped at 80** (historical default), deduplicated at 10-digit coordinate keys.
- **One-level refinement:** each primary point with `max_clique ≥ 6`, incomplete pool, or cert/tracked mismatch spawns a \(5^3\) (clique-6) or \(3^3\) local grid; second-level `ref_ref_*` is blocked by `can_spawn_refinement`.

**Expected row count:** \(20 + 80 + \text{refinements} \approx 588\) (verify via `--verify-only`).

**Purpose:** Reproducible baseline for “no fourth MUB on certified set” before scaling degen load.

---

## \(\mathcal{S}^\*\) — full degen primary (Track D completeness)

**Definition.** Same anchor set, but degeneracy candidates loaded with:

1. **Severity sort** by descending `z4_dev` (Mobius inconsistency first).
2. **Full `Float64` parse** of CSV coordinates (no extra truncation on load).
3. Configurable cap via `--degen-cap N` (default 80 for \(\mathcal{S}_{588}\); phased rollout uses 200, 500, or 1845).

\[
\mathcal{S}^\* = \mathcal{P}_{\mathrm{primary}}^{(N)} \;\cup\; \mathcal{R}_{\mathrm{ref}}^{(1)}
\]

**Batch strategy (self-directed):**

| Batch | Content | Expected clique-6 yield |
|-------|---------|-------------------------|
| 1 | Top 200 by `z4_dev` + all `mobius_z4` | Dita + high-dev loci |
| 2 | \(\theta=0\) Fourier slice rows | F6 arc neighbors |
| 3 | Remaining `circulant_match` | Mostly clique=2 (negative evidence) |

Use `--no-refine` for bulk primary-only passes (skip \(5^3\) spawn).

---

## \(\mathcal{R}\) — covering region for Track C

**Initial definition** (refined as Track D discovers components):

\[
\mathcal{R} = \bigl\{ (\theta,\phi,\lambda) \in K_6^{(3)} :
  \text{queued by `search_special_loci` with full degen load } (\mathcal{S}^\*)
  \;\text{OR}\; \text{on Dita } \lambda\text{-circle at } \phi=\pi/4
  \;\text{OR}\; \text{F6 } \theta=0 \text{ bounded } \lambda\text{-arc}
  \;\text{OR}\; \text{anchor set} \bigr\}
\]

**Operational membership:**

| Subregion | Script / artifact | Status |
|-----------|-------------------|--------|
| Anchors + capped degen + ref | `special_loci_search.csv` | Partial (\(\mathcal{S}_{588}\) incomplete) |
| Dita \(\lambda\)-circle | `lambda_periodicity_dita.jl`, `dita_lambda_fourth_dense.jl` | 126/126 + dense probe |
| F6 \(\theta=0\) \(\lambda\)-arc | `f6_boundary_reverify.jl` | Width \(\approx 0.118\) rad |
| Full degen \(\mathcal{S}^\*\) | `search_special_loci.jl --degen-cap 1845 --no-refine` | **1863/1865** primary complete; see coverage table |

**Coverage (2026-08-06, no-refine primary queue):**

| Artifact | `degen_cap` requested | `degen_cap` actual | Rows | Distinct \((\theta,\phi,\lambda)\) | Pool-complete | Clique-6 | Fourth MUB |
|----------|----------------------|-------------------|------|-------------------------------------|---------------|----------|------------|
| `special_loci_degen200.csv` | 200 | 200 | 219 | 219 (subset of degen500) | 217 | 4 | 0 |
| `special_loci_degen500.csv` | 500 | 500 | 519 | 519 (subset of degen1845) | 517 | 4 | 0 |
| `special_loci_degen500_run2.csv` | 1845 | **500** (mislabeled) | 519 | identical to degen500 | 517 | 4 | 0 |
| `special_loci_degen1845.meta.txt` | 1845 | **1845** | 1863 | **1863** (authoritative; log-backed) | 1852 | 21 loci | 0 |
| `special_loci_degen1845_CORRUPTED592.csv` | — | — | 592 | partial accidental rerun | — | — | do not cite |

Caps are **cumulative by severity sort** (top-\(N\) degen + 20 anchors), not disjoint partitions. Do not add row counts across batches. The design target is \(20 + 1845 = 1865\); two coordinate keys deduplicate at queue build (1863 distinct).

**Clique-6 / third-MUB count (verified 2026-08-06):** summary and CSV field `found_third` both mean `max_clique ≥ 6` at that candidate \((\theta,\phi,\lambda)\) via `check_four_mub_extension` (one hit per queue point, not per distinct third-basis solution). The degen1845 run reports **21** hits: all **21** distinct coordinate keys; **4** overlap the degen500 set (identical `max_clique`/`found_third` on all 519 shared keys); **17** are new at \(\theta=0\) `degen_circulant_match` points. Both runs used `--no-refine`. Independent re-evaluation of 4 spot coordinates confirmed `max_clique=6`.

**Track C theorem skeleton (honest scope):**

> For all \((\theta,\phi,\lambda) \in \mathcal{R}\) with `pool_complete_flag=true`, no fourth MUB exists.

Evidence: (i) exhaustive certified search with `found_fourth=false`; (ii) symbolic empty witness ideal on \(\theta=0\) slice (Phase 3); (iii) uniform clique bound on Dita circle.

This is **not** global \(N(6)=3\) for all CHMs in \(\mathbb{C}^6\).

---

## Certified subset \(\mathcal{S}\)

\[
\mathcal{S} = \{ p \in \mathcal{R} : \texttt{pool\_complete\_flag=true},\; n_{\mathrm{certified}} = n_{\mathrm{tracked}} \}
\]

Negative fourth-MUB claims apply to \(\mathcal{S}\) and \(\mathcal{R}\) only, not all of \(K_6^{(3)}\).

---

## Cross-reference

| File | Role |
|------|------|
| `scripts/julia/search_special_loci.jl` | Queue builder + certified search |
| `scripts/julia/csv_reconciliation.jl` | Row-count gap analysis |
| `scripts/julia/locus_classification.jl` | Track D component catalog |
| `results/locus_classification.json` | Classification deliverable |
| `REPRODUCE.md` | One-command regeneration |
| `results/final_honest_status.md` | Claim ledger |
