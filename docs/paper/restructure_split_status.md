# Option A split — titles locked

**Titles chosen 2026-08-16:**
- **Paper 1 (P1-B):** Gauge Structure and Certified Non-Extendability for Selected Karlsson Hadamards in Dimension Six
- **Paper 2 (P2-B):** Computational Identification and Audit of Third-MUB Loci Within Karlsson's Three-Parameter Family

Companion bibitems updated to match. Overleaf: upload bundled `main_theorems.tex` / `methods_audit.tex` as `main.tex`.

---

Frozen pre-split monolith (for diffs only): `mub6_karlsson_ieee.tex`, `main.tex`, `main_standalone.tex`, old `proofs/{gauge_and_locus,dita_third_mub,fourth_mub_obstruction}.tex`.

---

## (a) Title candidates

### Paper 1 — theorems (`main_theorems.tex`)

| ID | Title | One-line justification |
|----|-------|------------------------|
| **P1-A** | Certified Fourth-MUB Obstruction at Seven Dita Points in Karlsson's Family | Lead with T3; affirmative loci stay in Paper 2. |
| **P1-B** | Gauge Structure and Certified Non-Extendability for Selected Karlsson Hadamards in Dimension Six | Names T1+T3; “selected” hedges pointwise scope. |
| **P1-C** | Exact and Certified Witness Emptiness for Fourth MUBs in Karlsson's Three-Parameter Family | Emphasizes T3/T4 epistemic contrast (certified vs exact); T4 still demoted in abstract. |

### Paper 2 — methods/audit (`methods_audit.tex`)

| ID | Title | One-line justification |
|----|-------|------------------------|
| **P2-A** | A Reproducible Per-$H$ Pipeline for Third-MUB Loci in Karlsson's Family | Methods/reproducibility register; pipeline as product. |
| **P2-B** | Computational Identification and Audit of Third-MUB Loci Within Karlsson's Three-Parameter Family | Matches “identification not classification”; audit in the title. |
| **P2-C** | Search-Set Completeness, Grid Misses, and HP Re-Verification for MUBs in Karlsson's $K_6^{(3)}$ | Foregrounds S*/grid-miss honesty as the citable contribution. |

Reply with e.g. `P1-A + P2-B` (or a free-form title) before any PDF build.

---

## (b) File-by-file: what moved where

### Added (active sources)

| File | Role |
|------|------|
| `paper/preamble_common.tex` | Shared packages, theorems, `\repo` |
| `docs/paper/main_theorems.tex` | **Paper 1** body |
| `paper/methods_audit.tex` | **Paper 2** body |
| `docs/paper/proofs/gauge_structure.tex` | T1 lemmas + theorem (from old `gauge_and_locus` algebraic half) |
| `docs/paper/proofs/fourth_mub_theorems.tex` | T3 + T4 + remarks (from old `fourth_mub_obstruction`; L7 removed) |
| `docs/paper/proofs/locus_geometry.tex` | `lem:L5-semicont` + Findings L5/L6 + **FLAG comment** |
| `docs/paper/proofs/dita_third_mub_methods.tex` | Finding T2 + Finding L7; companion cites for gauge facts |
| `paper/restructure_split_status.md` | This checklist |

### Content map

| Material | Paper |
|----------|-------|
| T1 gauge algebra | 1 |
| T3 certified seven-point obstruction (**lead**) | 1 |
| T4 exact D0 GB (**pipeline validation billing**) | 1 |
| C1 A-matrix audit + `tab:audit` | 1 |
| New Scope and Logical Status | both (document-specific) |
| Minimal methods for W1 / certify / GB | 1 |
| Full per-H pipeline, V1–V4, implementation | 2 |
| Track D loci + rename “Computational Identification…” | 2 |
| Findings Track D, T2, L5, L6, L7 | 2 |
| `lem:L5-semicont` (flag only; content unchanged) | 2 |
| S*, R, Open Problem `op:trackc` | 2 |
| 512-grid miss, `tab:sweep`, `tab:reaudit` | 2 |
| Pool completeness `tab:poolcomplete`, probe tables, reproduce app | 2 |
| Pipeline figure | 2 |

### Fixes 1–7 status

1. **Abstract verbs (P1):** T3 = numerical-interval / `certify()`; T4 = exact ideal-membership / GB `{1}`; dense hedges → Scope §.  
2. **T4 billing (P1):** subordinate in abstract, contributions, results heading, conclusion, Remark `rem:T4-scope`.  
3. **Classification → identification (P2):** section title + prose; CHM “classification” in literature sense left alone where it referred to McNulty–Weigert surveys (none in P2 intro beyond family context). Script name `locus_classification.jl` / JSON path unchanged (repo artifacts).  
4. **Lemma flag:** `%% FLAG FOR AUTHOR REVIEW:` immediately before `lem:L5-semicont` in `proofs/locus_geometry.tex`. Lemma body untouched.  
5. Titles: candidates above; placeholders in `.tex`.  
6. **Booktabs:** all result tables in both new docs use `\toprule/\midrule/\bottomrule`.  
7. Cross-refs: see (c).

### Intentionally not updated yet

- `bundle_paper_tex.py` / Overleaf `main.tex` still target the monolith.  
- Old `proofs/*.tex` unchanged (diff baseline).  
- No PDF compile this pass.

---

## (c) Cross-reference pass — resolved vs ambiguous

### Resolved (internal `\ref` → companion `\cite`)

| Former internal target | Now |
|------------------------|-----|
| P1 → Finding Track D, T2, `op:trackc`, `app:reproduce`, probe tables | `\cite{JaupiMethods2026}` |
| P1 T3 remark on 628/628 dense probe | cite methods companion |
| P2 → T1 / L1–L4 / φ-gauge / λ-periodicity / Dita inequivalence | `\cite{JaupiTheorems2026}` |
| P2 → T3 certificate | `\cite{JaupiTheorems2026}` |
| Finding L7 | moved entirely into Paper 2 |

Companion bibkeys (titles pending): `JaupiMethods2026`, `JaupiTheorems2026`.

### Ambiguous — need your call (listed, not guessed)

1. **Open Problem on $\mathcal{R}$.** Placed formally in Paper 2 (`op:trackc`). Paper 1 only mentions the obstacle in Scope/Methods. Alternative: keep a short restated conjecture in Paper 1 for readers who never open Paper 2.

2. **Symbolic M2 / inexact-field material (`rem:m2-inexact`, `sec:symbolic`).** Kept in Paper 1 (needed for T3/T4 epistemology). Paper 2 still lists the M2 log path in the reproduce appendix. Alternative: duplicate a short “central obstacle” subsection in Paper 2.

3. **Gauge-equivalence testing prose** that sat in monolith Results (analytic φ-gauge + numerical CHM search). Algebraic half is Theorem T1 in Paper 1; the long numerical CHM-search narrative was folded into Paper 2 locus discussion via companion cite rather than copied verbatim. Say if you want that subsection restored in full in Paper 1 or 2.

4. **Deprecated fixed-pool approach** paragraph: currently Paper 2 only. Harmless either way; confirm OK.

5. **Phase C / Liang–Chen future-work bullet** from monolith Conclusion: **omitted** from both new conclusions (not in the Option A keep-lists). Restore to Paper 1, Paper 2, or drop?

6. **Script/artifact name `locus_classification.*`.** Left as-is (filesystem truth). Rename in prose only already done; renaming files would be a separate repo change.

7. **Numbered “Lemma 7” vs label `lem:L5-semicont`.** Flag is on upper-semicontinuity of $\kappa$ (the continuity/threshold tension). Separate Finding L7 (pool zero-dimensionality) is unrelated and sits in Paper 2 without a flag. Confirm that matches your intent.

8. **Inter-paper citation until arXiv IDs exist.** Bibitems are “companion manuscript, 2026 (title pending).” After titles: update both bibitems; optionally add `note={to appear}` / arXiv ids.

---

## Compile (after title approval)

```powershell
cd paper
pdflatex main_theorems.tex
pdflatex main_theorems.tex
pdflatex methods_audit.tex
pdflatex methods_audit.tex
```

---

**STOP:** reply with title picks (and optional answers to ambiguities 1–7) before final pass / PDFs.
