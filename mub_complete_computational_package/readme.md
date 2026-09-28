# Dimension-6 MUB Computational Research Package

This package collects the computational artifacts developed for the dimension-6 mutually unbiased bases (MUB) project, centered on Karlsson's three-parameter family, its Diţă subfamily/circle, third-MUB pool geometry, exact symbolic ideals, fixed-parameter symmetries, and the fourth-vector extension problem.

## Scope

The package is a research workspace, not a claim that all open problems have been solved. It separates:

- **exact algebraic constructions**: symbolic ideals, Laurent-polynomial Diţă representation, exact symmetry identities;
- **certified numerical computations**: where HomotopyContinuation interval certification was actually used;
- **ordinary numerical evidence**: continuation, singular-value measurements, optimization probes, and finite searches;
- **open hypotheses**: especially the global fourth-MUB problem over the whole Karlsson family and the interpretation of the 24 singular roots.

## Main mathematical object

For the Diţă circle, with z = exp(i lambda), the exact unnormalized Hadamard representative is

H_D(z) =
[[1,1,1,1,1,1],
 [1,-1,z,-z,i,-i],
 [1,-i,i,i,-i,-1],
 [1,i,-z,z,-1,-i],
 [1,z^-1,-i,-1,-z^-1,i],
 [1,-z^-1,-1,-i,z^-1,i]].

The normalized matrix is H_D(z)/sqrt(6).

A third-MUB vector is represented in gauge v=(1,x1,...,x5)/sqrt(6), with inverse variables y_j satisfying x_j y_j=1. The exact single-vector ideal is

I3^(1) = < z t - 1, x_j y_j - 1, F_1,...,F_6 >,

where F_k=(v^* h_k)(h_k^* v)-6 for the unnormalized H_D(z).

The six-vector third-MUB branch is built from six copies of this vector ideal plus pairwise orthogonality equations and one shared parameter relation zt-1.

## Current known computational status

1. The symbolic single-vector ideal was generated exactly.
2. The six-vector I3 branch has 97 generators before redundancy reduction and 62 algebraic variables in the explicit inverse-variable formulation.
3. Complex conjugation was checked exactly as z<->t and x_j<->y_j with i->-i.
4. A fixed-z exhaustive numerical search over 720 row permutations found 12 compatible monomial automorphisms at a generic numerical sample; each found automorphism was then verified symbolically. Their row permutations form A4. This establishes the tested 12-element symmetry set at that sample, not a global completeness theorem for symbolic generic z.
5. Earlier computations found a numerical MU-pool transition 120 -> 72 near lambda*=0.11148022437796655429..., with 24 singular MU-vector roots recovered in an augmented singularity computation. The interpretation of those 24 roots as two A4 orbits of size 12 is a hypothesis until explicitly checked.
6. Seven Dita-circle parameter values have certified fourth-MUB obstruction in the theorem work. That result is pointwise over the recovered third-MUB cliques at those parameters, not a theorem over the full Karlsson family.

## Recommended next computation

Run an A4 orbit decomposition on the actual third-MUB cliques at lambda=0.4, pi/3, and 2pi/3, then test the 24 singular roots at lambda*. If the roots split into two 12-element A4 orbits, record this as a verified orbit statement. Only after this should the reduced fourth-vector incidence ideal I4 be constructed for inequivalent B3 representatives.

## Reproducibility

The canonical exact I3 implementation is integrated at
`scripts/python/dita_i3.py`. From the main repository root:

    python scripts/python/dita_i3.py
    python -m unittest discover -s test -p "test_dita_i3.py" -v

Generated equations, symmetry checks, term counts, and a run manifest are
written to `symbolic_export/I3/`. The original
`01_symbolic_I3/build_i3.py` path remains a compatibility wrapper. Python
requires SymPy. For numerical algebraic geometry, the project uses Julia +
HomotopyContinuation.jl. See the official documentation:
https://www.juliahomotopycontinuation.org/HomotopyContinuation.jl/stable/

The package is designed for Windows/Anaconda + Julia workflows, but the symbolic scripts are platform-independent.

## Important research hygiene

Do not promote numerical observations into theorem statements. In particular:

- a finite parameter sweep does not prove absence on a continuum;
- numerical continuation does not prove branch completeness at singularities;
- a small optimization residual is not a non-existence certificate;
- an orbit pattern must be checked explicitly before being used in a proof;
- a Gröbner basis equal to {1} is exact only for the exact ideal actually supplied to the CAS.
