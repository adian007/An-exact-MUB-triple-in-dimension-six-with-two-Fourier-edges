# Step 0 — Restructuring Proposal (no LaTeX changes yet)

**Source of truth for edits (once approved):** `docs/paper/main_theorems.tex` + `docs/paper/proofs/*.tex`.
`main.tex` / `main_standalone.tex` are bundled regenerations — update the multifile sources, then rebundle.

**Current architecture (compressed):**
| Block | Contents |
|-------|----------|
| Front | Abstract (T1/T3/T4 + Finding T2 + C1/C2 + S*) |
| Intro | Scope/prior; contributions C1–C4 |
| Background → Methods → Implementation | Pipeline |
| Results | C1 audit; coarse-grid miss; re-audit; S_588/S*; **“Third-MUB Locus Classification”**; Track C/A + T3; pool tables; reproducibility |
| Discussion / Conclusion | Epistemic summary + caveats scattered |
| Appendices | Proofs (T1, L5–L6 findings, T2, T3, T4); probe tables; reproduce logs |

---

## Option (A) — Split into two documents

### Paper 1 — Main (“theorem-first”)
**Target:** short, sharp; one clear answer to “what is this paper about?”

**Keep:**
- **T1** — gauge structure (Lemmas L1–L4 + Theorem T1; `proofs/gauge_and_locus.tex` algebraic half)
- **T3** — seven-point certified fourth-MUB obstruction (`proofs/fourth_mub_obstruction.tex` through T3)
- **T4** — D0 exact re-verification, **billed as pipeline validation / independent BW check**, not co-equal headline
- **C1** — Karlsson vs McNulty–Weigert A-matrix correction (Table `tab:audit`)
- Minimal Methods: problem formulation, `certify()` / W1, enough for T3/T4 to be readable
- Dedicated **“Scope and Logical Status of Results”** (new)
- Slim Discussion/Conclusion focused on T1/T3 (+ T4 as validation) and C1
- Proof appendix for T1/T3/T4 only

**Move out to Paper 2 (cite Paper 2; one-paragraph pointer in Paper 1):**
- Finding Track D / Finding T2 / prop L5–L6 (third-MUB loci)
- Entire Track D Results subsection; locus probe tables (`tab:locusdim`, `tab:phisweep`)
- S* sweep, 512-point grid-miss story, coarse/re-audit tables (`tab:sweep`, `tab:reaudit`)
- Pool-completeness tables beyond what T3 needs (`tab:poolcomplete` → methods paper unless one row needed for T3)
- Reproducibility appendices, CSV/JSON inventory, pipeline-overview figure as primary artifact
- Open Problem on all of R (can stay as one sentence in Paper 1 pointing to Paper 2)

**Gained:** Referee-proof focus; T3 as the negative certificate headline; C1 as clean side result; T4 demoted by design.
**Lost:** Single citable object; Paper 1 cannot “own” the affirmative third-MUB loci discovery without cross-cite; some narrative glue (why Dita matters for T3) must be summarized, not proved computationally in-place. Risk: Paper 2 looks like an appendix that escaped.

**Suggested Paper 1 section map:**
1. Intro (+ contributions: T1, T3, C1; T4 as validation note)
2. Scope and Logical Status
3. Background (short)
4. Methods (minimal for T3/T4)
5. Results: C1; T1 summary; T3; T4-as-validation
6. Discussion / Conclusion
7. Appendix: proofs T1/T3/T4

---

### Paper 2 — Methods / audit (“computational methodology report”)
**Stand-alone, citable.**

**Keep as primary content:**
- Full per-H pipeline, completeness gate, planted-clique validation
- **Computational identification of third-MUB loci** (rename; Finding Track D + L5/L6 + T2)
- S* / R covering region; 512-point grid-miss as methodology lesson
- Pool completeness; reproducibility; symbolic-elimination failure log (central obstacle context)
- Cross-reference Paper 1 for T3 certificate and T1 gauge facts used in locus geometry

**Gained:** Honest home for audit trail; “classification → identification” language fits without competing with theorems.
**Lost:** Methods paper without T3 may feel incomplete unless Paper 1 is cited early; duplicate Background/Methods boilerplate.

---

## Option (B) — Single paper, restructured internally

**Spine (three contributions only):**
1. **Main theorem:** T3 (seven-point certified W1 emptiness)
2. **Main computational discovery:** Dita/F6 third-MUB loci (Finding Track D + supporting HP findings; **not** “classification”)
3. **Methodological contribution:** per-H pipeline (incl. completeness gate, planted validation)

**Subordinate (keep, demote billing):**
- T1 — algebraic prerequisite / gauge facts (needed for loci + inequivalence)
- T4 — pipeline validation / independent BW re-verification (not headline)
- C1 — literature correction (short Results subsection or dedicated short section)

**New / relocated sections:**
| New home | What moves there |
|----------|------------------|
| **§ Scope and Logical Status of Results** | Dense abstract qualifiers; scattered “not a theorem on R / N(6)” caveats; Proved vs certify() vs HP vs conjectured ledger |
| **Main Results** | T3 first (or after brief T1); then loci identification; then pipeline value |
| **Back-matter / Audit appendix** | 512-grid miss narrative depth; S* CSV detail; pool-completeness tables; reproduce logs; superseded coarse sweep |
| Front Results | Keep C1 + high-level loci + T3 summary; push audit depth back |

**Current sections → B layout (sketch):**
- Intro contributions list → rewrite to three-spine + demoted T4/C1
- Results § “Third-MUB Locus Classification” → rename + keep as discovery pillar
- Results § Track C/A + long S* prose → shorten in main text; detail → Audit appendix
- Discussion caveats list → absorb into Scope section; Discussion focuses on interpretation
- Appendices proofs unchanged in mathematical content; Lemma `lem:L5-semicont` flagged only (item 4)

**Gained:** One arXiv object; loci discovery stays visible; epistemic hedges centralized; audit trail stop competing with T3.
**Lost:** Still longer than a pure theorem note; referee may still want a split later; discipline required so Audit appendix does not leak into the abstract again.

---

## Fixes that apply under either choice (deferred until you pick A or B)

1. **Abstract verbs:** T3 = numerical-interval / `certify()` certificate at seven points; T4 = exact ideal-membership (GB = {1}) at one D0-equivalent pair. Trim qualifier chain → Scope section.
2. **T4 billing:** subordinate pipeline-validation throughout abstract / intro / conclusion; proof untouched.
3. **“Classification” → “Computational Identification of Third-MUB Loci”** (headers, finding titles, cross-refs; keep CHM-classification literature sense untouched).
4. **Lemma `lem:L5-semicont` (Upper semicontinuity of κ):** FLAG ONLY — proof sketch claims orthogonality-edge stability under small perturbations, but hard-threshold clique membership is discontinuous; paper’s own F6 `\|Δλ\|≈0.055/0.063` and Dita `\|Δϕ\|≥10^{-3}` drops (Table `tab:phisweep`) document sharp drops. No lemma text change.
5. **Title options** (for the doc that holds T1/T3/T4) — show, don’t pick:
   - *Certified Fourth-MUB Obstructions at Dita Points in Karlsson’s Family* (negative/certificate-led)
   - *Third-MUB Loci in Karlsson’s Three-Parameter Family: Identification and Certified Non-Extendability* (affirmative then negative)
   - *Gauge Structure and Certified Non-Extendability for Selected Karlsson Hadamards in Dimension Six* (T1+T3; loci left to subtitle/Paper 2)
6. **Tables → booktabs** (`\toprule/\midrule/\bottomrule`); content/numbers unchanged. (`booktabs` already loaded; several tables already use it — audit remaining.)
7. **Cross-ref integrity** after moves: all `\label`/`\ref`/`\cite`; regenerate bundle if using `main.tex`.

---

## Constraint reminder (for the execution pass)
- No change to numerical results, theorem/lemma *statements*, or proof bodies except: abstract trim, T4 billing, classification wording, Lemma flag comment, table rules, structural moves.
- After structural edits: file change summary for your approval before any “final” cleanup pass; keep diffs reviewable section-by-section.

---

## Recommendation (non-binding)
- Choose **(A)** if the primary goal is a sharp theorem note for a theory-leaning venue and a separate citable methods artifact.
- Choose **(B)** if you want one arXiv preprint that still foregrounds the Dita/F6 loci discovery alongside T3.

**STOP — awaiting your choice: A or B** (and optionally a preferred title from §5, or “defer title”).
