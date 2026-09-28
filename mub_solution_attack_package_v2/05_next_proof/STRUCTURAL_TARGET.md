# Structural theorem target

Let V3(lambda) be the complete third-MUB incidence variety for H_D(lambda),
and let C be a complete B3 clique. The target is to establish, modulo the
verified A4 action and natural MUB equivalences, that a transition associated
with every relevant component satisfies Theorem 1's three-distinct-columns
condition. Then determine whether its conclusion is transposed Fourier or
2-circulant.

Only a separately established Fourier-family branch lets the
Jaming--Matolcsi--Mora--Szollosi--Weiner theorem rule out such a quartet; the
2-circulant alternative needs separate treatment.

The difficult steps are to prove the theorem condition on the entire relevant
component, resolve its family alternative, and show that no additional B3
components escape the A4/known-component classification.

A particularly promising implementation is:

1. Represent the B3 transition matrix by dephased polynomial variables.
2. For the Fourier branch only, add Fourier-family incidence equations (two
   unit-modulus parameters plus a finite row/column permutation chart).
   Theorem 1's condition is not itself Fourier-family incidence, and the
   2-circulant branch must also be handled.
3. Compute the elimination ideal of the B3 incidence variety intersected with the complement of the Fourier subvariety.
4. Show the resulting component is empty, ideally by Groebner elimination or exact numerical algebraic geometry.
5. Repeat for the transpose-Fourier chart if necessary.
6. Handle the fold lambda* by working on the algebraic incidence variety rather than continuation through the singular parameter.

This is a considerably stronger and more focused target than searching directly for a fourth vector.
