# K3-to-Fourier structural theorem target

## Definitions
Let H_D(λ) denote the exact Laurent-polynomial representative of the Dita circle inside Karlsson's K_6^(3) family.

Let V_3 be the incidence variety

  V_3 = {(λ,B_3) : B_3 is a third MUB for {I,H_D(λ)}}.

For a six-vector basis B_3, define the transition matrix

  T = H_D(λ)^† B_3.

Then T is complex Hadamard whenever B_3 is MU to H_D(λ).

## Target statement
The first target is

  Every relevant irreducible component of V_3 has a representative whose transition satisfies Theorem 1's three-distinct-columns condition.

The theorem concludes transposed Fourier or 2-circulant membership. The
Fourier-family target below is stronger and requires separately resolving
that alternative.

or, more ambitiously,

  the three transition matrices of every relevant MUB triplet are, up to permutational unitary equivalence, one member of X(α), one member of F(a,b), and one member of F^T(a,b).

## Why this is enough
Jaming–Matolcsi–Móra–Szöllősi–Weiner proved that no pair {I,F(a,b)} extends to a MUB quartet. This rules out a quartet only if the Fourier-family branch is independently established; Theorem 1's 2-circulant alternative does not by itself supply that implication.

## Computational subgoals
1. Use the exact I3 ideal and A4 quotient to enumerate the third-MUB incidence components.
2. Introduce Theorem 1 chart equations for T = H_D^† B_3, then analyze the
   transposed-Fourier and 2-circulant alternatives separately. Fourier-family
   incidence equations address only the Fourier-specific subgoal.
3. Eliminate B_3 variables to compare the I3 component with the Fourier incidence locus.
4. At the fold λ*, use a deflated/augmented system rather than ordinary nonsingular certification.
5. Prove component containment by exact elimination or exact polynomial identity, not by sample fitting.
6. Handle both Fourier and transposed-Fourier orientations and all relevant permutations/rephasings.

## Stop condition
Do not claim the theorem from a finite list of numerical samples. The required proof must establish algebraic containment, certified continuation across every relevant component, or an equivalent exact argument.
