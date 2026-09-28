# K3-to-Fourier structural theorem target

## Definitions
Let H_D(λ) denote the exact Laurent-polynomial representative of the Dita circle inside Karlsson's K_6^(3) family.

Let V_3 be the incidence variety

  V_3 = {(λ,B_3) : B_3 is a third MUB for {I,H_D(λ)}}.

For a six-vector basis B_3, define the transition matrix

  T = H_D(λ)^† B_3.

Then T is complex Hadamard whenever B_3 is MU to H_D(λ).

## Target statement
The preferred target is

  Every irreducible component of V_3 containing the recovered Dita-circle branches has a representative for which T is in the Fourier family F(a,b),

or, more ambitiously,

  the three transition matrices of every relevant MUB triplet are, up to permutational unitary equivalence, one member of X(α), one member of F(a,b), and one member of F^T(a,b).

## Why this is enough
Jaming–Matolcsi–Móra–Szöllősi–Weiner proved that no pair {I,F(a,b)} extends to a MUB quartet. Therefore a proof that a hypothetical quartet containing one of our K3 triplets necessarily contains a Fourier transition matrix would rule out the quartet.

## Computational subgoals
1. Use the exact I3 ideal and A4 quotient to enumerate the third-MUB incidence components.
2. Introduce Fourier-family incidence equations for T = H_D^† B_3.
3. Eliminate B_3 variables to compare the I3 component with the Fourier incidence locus.
4. At the fold λ*, use a deflated/augmented system rather than ordinary nonsingular certification.
5. Prove component containment by exact elimination or exact polynomial identity, not by sample fitting.
6. Handle both Fourier and transposed-Fourier orientations and all relevant permutations/rephasings.

## Stop condition
Do not claim the theorem from a finite list of numerical samples. The required proof must establish algebraic containment, certified continuation across every relevant component, or an equivalent exact argument.
