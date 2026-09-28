# Structural theorem target

Let V3(lambda) be the complete third-MUB incidence variety for H_D(lambda). Let C be a complete B3 clique. We want to prove that, modulo the verified A4 action and the natural MUB equivalences, at least one transition matrix among

1. H_D(lambda),
2. B3(lambda),
3. H_D(lambda)^* B3(lambda)

belongs to the Fourier family F(x,y).

A successful proof would imply that any fourth basis extending the triple would extend a pair of bases whose transition is Fourier-family. The Jaming--Matolcsi--Mora--Szollosi--Weiner theorem rules out such a quartet for every Fourier-family parameter value.

The difficult step is to prove the Fourier-family membership algebraically for the entire relevant component and to prove that no additional B3 components escape the A4/known-component classification.

A particularly promising implementation is:

1. Represent the B3 transition matrix by dephased polynomial variables.
2. Add the Fourier-family incidence equations (two unit-modulus parameters plus a finite row/column permutation chart).
3. Compute the elimination ideal of the B3 incidence variety intersected with the complement of the Fourier subvariety.
4. Show the resulting component is empty, ideally by Groebner elimination or exact numerical algebraic geometry.
5. Repeat for the transpose-Fourier chart if necessary.
6. Handle the fold lambda* by working on the algebraic incidence variety rather than continuation through the singular parameter.

This is a considerably stronger and more focused target than searching directly for a fourth vector.
