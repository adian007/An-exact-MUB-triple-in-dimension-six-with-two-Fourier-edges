# Symbolic I3 Third-MUB Branch Ideal — Diţă Circle

This bundle is the first concrete algebraic object in the fourth-MUB elimination program.

## Exact Hadamard family

We use

`H_D(z) = 1/sqrt(6) * [[1,1,1,1,1,1], [1,-1,z,-z,i,-i], [1,-i,i,i,-i,-1], [1,i,-z,z,-1,-i], [1,z^-1,-i,-1,-z^-1,i], [1,-z^-1,-1,-i,z^-1,i]]`

with `z*t-1=0` and `t=z^{-1}`.

## Single-vector ideal I3(z)

A flat vector is gauged as

`v=(1,x1,...,x5)/sqrt(6)`

and inverse variables `y_j` satisfy `x_j*y_j-1=0`.

The six exact MU equations are

`F_k=(v^* h_k)(h_k^* v)-6=0`, `k=1,...,6`,

where conjugation is represented algebraically by `x_j -> y_j`, `z -> t`, `i -> -i`.

Thus

`I3_single = < z*t-1, x_j*y_j-1, F_1,...,F_6 >`.

The equations are quadratic in the x/y variables; after retaining z,t as independent generators, F2-F5 have total degree 4 because of products containing z*t. Modulo `z*t-1`, they have the expected Laurent structure.

## Six-vector third-MUB branch ideal

For six gauged vectors `v_a`, take six copies of the single-vector MU ideal and add, for every `a<b`,

`<v_a,v_b>=0` and its conjugate.

This gives a direct algebraic model of the third-MUB incidence fiber. A branch is a six-vector fiber over z with all pairwise orthogonality equations satisfied.

## Exact symmetries found

### 1. Complex conjugation

`z -> t`, `x_j -> y_j`, `y_j -> x_j`, `i -> -i`.

Each F_k maps exactly to itself. This is an involution and sends a solution at z to a solution at z^{-1}.

### 2. Parameter shift lambda -> lambda + pi

`z -> -z`, `t -> -t`, `x4 <-> x5`, `y4 <-> y5`.

The exact equations transform as

`F1 -> F1, F2 -> F2, F3 <-> F4, F5 -> F5, F6 -> F6`.

At the matrix level this comes from the exact identity that `H_D(-z)` is obtained from `H_D(z)` by swapping rows 5 and 6 and columns 3 and 4. The row operation acts on MU vectors; the column operation only permutes Hadamard basis vectors.

### 3. Generic fixed-z automorphism group

A numerical discovery followed by exact symbolic verification finds exactly 12 monomial row automorphisms of `H_D(z)` at generic z. Their row-permutation parts form a closed group with cycle types

- identity: 1
- double transpositions: 3
- products of two 3-cycles: 8

so the permutation group is isomorphic to `A4`.

Two generators are:

`G2 = (0,3,2,1,5,4)`

and

`G3 = (1,4,3,5,0,2)`

in zero-based row-index permutation notation. The corresponding monomial matrices U satisfy `U H_D(z) Q = H_D(z)` for exact monomial Q, symbolically over Q(i,z,z^{-1}).

The induced action on MU vectors is `v -> U v`, followed by the usual global-phase gauge fixing of the first coordinate.

Therefore the generic fixed-z symmetry group is order 12, not 24. The previously observed 24 singular roots should therefore be investigated as possible two-orbit structure under A4, not as evidence for a 24-element symmetry group.

## Files

- `build_i3.py`: constructs the exact single-vector equations and tests the parameter symmetries.
- `I3_single_vector.txt`: expanded exact generators F1,...,F6.
- `I3_equations.py`: machine-readable SymPy equations.
- `i3_branch_definition.md`: definition of the six-vector branch ideal.
- `symmetry_test.txt`: exact equation-level tests of conjugation and z -> -z.
- `fixed_symmetry_exact.txt`: exact monomial automorphisms and generators.
- `fixed_symmetry.txt`: generic numerical discovery of the 12 automorphisms.
- `search_zminus_sym.py`: numerical discovery of the lambda -> lambda+pi matrix equivalence.
- `extract_fixed_sym.py`: exact verification of the fixed-z automorphisms.

## Status

This is an exact formulation and exact symmetry verification. It does NOT yet prove that every third-MUB branch has been classified, nor does it prove fourth-MUB non-extendability. The next object is the A4 orbit decomposition of the six-vector branch and then the reduced fourth-vector incidence ideal on one representative orbit.
