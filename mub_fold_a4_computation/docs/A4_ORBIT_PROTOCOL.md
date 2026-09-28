# A4 orbit protocol

The fixed-z monomial automorphism search previously found 12 compatible row monomial automorphisms at a generic numerical sample, each subsequently checked symbolically. Their row permutations form a group of order 12 with cycle types identity, three double transpositions, and eight products of two 3-cycles, consistent with `A4`.

The two stored generators are:

`G2 = (0,3,2,1,5,4)`

`G3 = (1,4,3,5,0,2)`

with the corresponding monomial matrices in the source-context files.

For an isolated MU vector `v`, apply every group element, regauge so the first component is one, and match the result against the isolated solution list. Record:

- orbit size;
- stabilizer size;
- all matched indices;
- maximum matching residual;
- whether the generated orbit closes inside the numerical pool.

Do not infer a 24 = 12 + 12 decomposition until all 24 singular roots have been isolated independently of the symmetry action.
