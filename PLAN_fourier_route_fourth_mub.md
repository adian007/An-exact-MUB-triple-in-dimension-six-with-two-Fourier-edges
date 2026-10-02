# Shortest Fourier-Family Route Toward a Fourth-MUB Result in Dimension 6

## Objective and scope

Use the exact Fourier-family nonextension theorem of Jaming, Matolcsi, Móra,
Szöllősi, and Weiner to prove the strongest result available from the current
\(\lambda=\pi/3\) candidate with the least additional computation.

The first target is a theorem about one **fixed exact triple**
\(\{I,H,B_3\}\): prove that it cannot be extended to a fourth MUB. This is
a meaningful local theorem, but it does not settle whether *any* four MUBs
exist in \(\mathbb C^6\). Keep those claims separate throughout.

## Logical route

The published theorem says that a pair \(\{I,F(a,b)\}\), where \(F(a,b)\)
is a member of the order-six Fourier family, cannot be extended to a MUB
quartet. The same conclusion applies after a common unitary change of
coordinates, and the project has a checked transpose/orientation lemma.

Therefore, for the exact triple \(\{I,H,B_3\}\), it is enough to prove that
**one pair inside the triple** has a transition matrix exactly equivalent to
a Fourier-family member. Check all three pairs:

- \(\{I,B_3\}\), represented by \(B_3\);
- \(\{I,H\}\), represented by \(H\);
- \(\{H,B_3\}\), represented by \(H^\dagger B_3\) or its adjoint,
  depending on the matrix convention.

If any one pair is exactly Fourier-family, the published theorem rules out a
quartet containing that pair, and hence rules out a fourth basis extending
the fixed triple. Do not infer this from a small numerical fitting residual.

## Starting evidence in this repository

- The exact \(\pi/3\) candidate basis and phase-root isolating intervals are
  recorded in [the exact B3 ansatz](symbolic_export/pi3_exact_B3_ansatz.md).
- Its supplied ansatz has exact checks for all 15 pairwise orthogonality
  equations and all 18 Hadamard-unbiasedness equations:
  [the exact remainder report](docs/research/reports/phase2_pi3_exact_remainder_verification_2026-10-01.md).
- The candidate was reconstructed from one clique in a 72-vector pool. That
  graph has four six-cliques, but only within the recovered pool; the pool is
  not known to be complete.
- Existing Fourier-family fits at \(\lambda=\pi/3\) have residuals near
  machine precision, but are numerical fits:
  [the Fourier-family audit](docs/results/k3_family_membership_2026-09-29.md).
- The exact \(W_1\) elimination attempt timed out without a conclusion, and
  its joint coefficient-field component still needs validation:
  [the \(W_1\) Gröbner report](docs/research/reports/phase2_w1_groebner_run_2026-10-01.md).

## Work plan

### Gate 0 — Freeze and validate the exact input

1. Select the saved \(\lambda=\pi/3\) clique used by the exact ansatz.
2. Confirm that the root intervals for \(a,b,c\) isolate the intended complex
   roots, not merely roots of the same phase polynomial.
3. Verify that the six exact columns, with the stated \(1/\sqrt 6\)
   normalization, match the chosen numerical clique after permitted column
   phases and permutations.
4. Verify exactly that \(H\) is the intended Diţă matrix at
   \(\lambda=\pi/3\), and that \(B_3\) is unbiased to both \(I\) and \(H\).

**Pass condition:** one named algebraic embedding represents a valid exact
triple and is tied to the stored candidate. If the exact root selection or
normalization fails, repair this before attempting a family fit.

### Gate 1 — Find the easiest Fourier chart

Build the exact transition matrices for the three pairs listed above from the
selected algebraic embedding. For each, test:

- direct and transposed Fourier-family conventions;
- the row/column phase and permutation equivalences allowed by the cited
  theorem and repository conventions;
- whether the candidate is Fourier-family rather than merely another
  Hadamard family.

Use the existing numerical fitter only to identify a likely chart, parameter
values, and row/column relabeling. Rank the candidate charts by how few
independent algebraic quantities must be reconstructed. Start with the
simplest chart; keep the other pair matrices as independent cross-checks.

**Pass condition:** a precise proposed identity of the form
\(M=D_1P_1F(x,y)P_2D_2\) or the documented transpose/conjugate version,
where \(M\) is one exact pair transition matrix, \(D_i\) are diagonal unitary
matrices, and \(P_i\) are permutation matrices.

### Gate 2 — Replace the fit by an exact membership certificate

1. Recover exact algebraic values for the two Fourier parameters \(x,y\)
   from entries or entry ratios of the exact matrix. Do not infer them from
   rounded angles alone.
2. Express all quantities in an explicit number field or a rigorously
   specified quotient with a selected prime component.
3. Use the isolating intervals/complex rectangles for the selected roots to
   identify the correct embedding.
4. Verify the proposed matrix identity exactly entry by entry after applying
   the stated phases, permutations, and orientation.
5. Independently verify the exact transition matrix is Hadamard and that the
   identity uses the same indexing and Fourier convention as the theorem.

A field defined by equations is not enough unless the component and embedding
used by the numerical candidate are specified. Separate resultant
compatibility checks do not by themselves prove that the three phase choices
lie in one common field component.

**Pass condition:** a reproducible exact certificate (for example, zero
remainders in a named number field plus root-isolation data) establishing
Fourier-family membership. A floating-point residual, however small, does not
pass this gate.

### Gate 3 — Apply the published theorem with the correct orientation

1. Write down the exact pair of bases whose transition matrix is the Fourier
   member.
2. State the allowed row/column transformations and show they preserve the
   extension question. Row operations must be justified as a common ambient
   unitary on both bases; column operations relabel or rephase basis vectors.
3. If using the transpose or conjugate convention, invoke and verify the
   repository's common-unitary/conjugation equivalence argument.
4. Apply Theorem 1.4 of Jaming et al. to that pair.

**Pass condition:** a formal contradiction from the assumption that the fixed
triple extends to a quartet.

**Permitted conclusion:** “This specified exact \(\{I,H,B_3\}\) does not
extend to a fourth MUB.”

Do not conclude from this one clique that no fourth MUB exists in dimension
six, or even that the selected \(H\) has no other possible third basis.

### Gate 4 — Expand the result only as far as completeness is proved

To exclude every fourth basis extending the selected pair \(\{I,H\}\):

1. Prove that every possible third basis for this fixed \(H\) is represented
   by the recovered pool, or replace pool enumeration with a complete exact
   classification.
2. Prove exact symmetry/orbit coverage of the four recovered six-cliques.
3. Prove Fourier-family membership for one representative of every exact
   orbit, and show the transformations preserve the pair-extension property.

The current four-clique graph is exhaustive only inside its 72-vector input
pool. Until pool completeness is established, the conclusion remains
candidate-specific.

### Gate 5 — State the boundary to the global problem

A proof for the \(\pi/3\) triple, or even for all triples at one fixed \(H\),
does not settle global nonexistence. The global claim requires ruling out every
possible set of four bases, including triples whose pair transition matrices
are not in the Fourier family and Hadamard matrices outside the currently
studied Karlsson/Diţă region.

If the aim changes from a local theorem to \(N(6)=3\), first state and prove
the missing classification/coverage theorem. Do not extrapolate from samples
over the Diţă circle or from a classification that is only a preprint claim.

## Tooling and execution order

- **Python/SymPy:** generate and reduce exact equations, derive entry-ratio
  identities, manage root-isolation data, and emit reproducible inputs.
  Avoid using floating-point fitting as the final proof.
- **Macaulay2 or Singular:** exact ideal membership, elimination, prime
  components, and field/quotient checks. First solve the smaller
  Fourier-membership problem. Do not rerun the ten-variable \(W_1\)
  Gröbner job until the coefficient component is explicit.
- **Julia/HomotopyContinuation.jl:** use for candidate discovery, high-
  precision refinement, and certification of isolated nonsingular roots.
  It is supporting numerical algebraic geometry, not the final global
  argument.
- **Reproducibility:** every pass gate should record the input hash, exact
  command, software versions, field presentation, root isolations, output
  hash, and claim scope.

Benchmark Macaulay2 and Singular on the *same reduced exact membership
problem* if runtime is a concern. A different CAS may speed up the computation,
but it cannot repair an ambiguous field component or turn a numerical fit
into an exact identity.

## Stop / redirect criteria

- If no pair transition has a plausible exact Fourier chart, stop this route
  and inspect whether the numerical family fit was a false positive or the
  matrix is in a different Hadamard family.
- If the candidate is exactly in another family, do not apply the Fourier
  theorem; identify a theorem that covers that family or return to exact
  \(W_1\) elimination for the fixed triple.
- If exact Fourier membership passes but the recovered pool is incomplete,
  report the fixed-triple theorem and continue pool-completeness work
  separately.
- If membership is proved only numerically, retain the result as a conjectural
  bridge and do not claim nonextension.

## Key references

- Jaming, Matolcsi, Móra, Szöllősi, and Weiner, “A generalized Pauli problem
  and an infinite family of MUB-triplets in dimension 6,” Theorem 1.4,
  [arXiv:0902.0882](https://arxiv.org/abs/0902.0882).
- Project Fourier fit and theorem-convention audit:
  [the Fourier-family report](docs/results/k3_family_membership_2026-09-29.md).
- Exact \(\pi/3\) basis candidate:
  [the exact B3 ansatz](symbolic_export/pi3_exact_B3_ansatz.md).
- Current unresolved proof status:
  [the \(W_1\) run report](docs/research/reports/phase2_w1_groebner_run_2026-10-01.md).
