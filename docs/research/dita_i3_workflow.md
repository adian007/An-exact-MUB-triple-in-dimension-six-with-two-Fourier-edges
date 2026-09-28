# Exact Diţă-circle I₃ workflow

## Purpose and scope

[`scripts/python/dita_i3.py`](../../scripts/python/dita_i3.py) is the canonical
SymPy implementation of the exact Diţă representative
\(H_D(z)\), its single-vector MU equations, the six-vector third-basis
incidence ideal, and the uploaded fixed-\(z\) symmetry certificates. It keeps
the algebraic parameter variables \(z,t\) and records the Laurent relation
\(zt-1=0\), rather than replacing the parameter with a floating-point sample.

The single-vector gauge is
\[
v=(1,x_1,\ldots,x_5)/\sqrt 6,\qquad x_jy_j=1,
\]
and the six MU equations are
\[
F_k=(v^*h_k)(h_k^*v)-6=0.
\]
The conjugate variables are represented algebraically by \(x_j\leftrightarrow
y_j\), \(z\leftrightarrow t\), and \(i\mapsto-i\).

The six-vector branch system shares one \(zt-1\) relation, adds five inverse
relations and six MU equations for each of six vectors, and adds both
orthogonality equations for every pair. Before redundancy reduction this is
62 variables and 97 generators (30 of them pairwise orthogonality equations).

## Reproduce

From the repository root, with Python and SymPy available:

```powershell
python scripts/python/dita_i3.py
python -m unittest discover -s test -p "test_dita_i3.py" -v
```

Generated exact equations, symmetry checks, term counts, and a run manifest
are written to `symbolic_export/I3/`. Choose a separate output directory with
`--out-dir PATH` to compare or preserve an independent run. The generated
`I3_equations.py` file is importable; it defines `z`, `t`, `x`, `y`,
`relations`, and `F`.

The uploaded package's original entry point remains as a compatibility
wrapper:

```powershell
python mub_complete_computational_package/01_symbolic_I3/build_i3.py
```

`a4_vector_orbit` and `a4_clique_orbits` in the same module apply the two exact
monomial generators numerically and identify projective vector/basis orbits.
The clique routine fails explicitly when a generator image is absent from the
provided set; it does not turn incomplete input into a complete orbit claim.

## Interface boundary with the existing pool solver

The project pool/clique pipeline builds a Karlsson representative; the
uploaded exact certificates are stated for the explicit \(H_D(z)\) matrix.
The exact algebra and orbit routines therefore do not silently consume the
current pool's vector coordinates. Before applying the action to a pool
result, establish and test the row/column permutation, dephasing, and vector
coordinate map between those representatives. Until then, pass only vectors
expressed in the documented \(H_D(z)\) row gauge to the orbit functions.

For follow-up data, use `enumerate_all_third_mub_bases` from
[`src/Cliques.jl`](../../src/Cliques.jl) to enumerate verified pool-relative
cliques, then convert them with an explicit equivalence certificate before
calling `a4_clique_orbits`. A result remains relative to the supplied pool;
the exact ideal itself does not prove that the pool is complete.

## Evidence boundaries

- The \(H_D(z)\) formulas, equations, \(zt=1\) reductions, and two generator
  identities are exact symbolic checks.
- The two fixed-\(z\) generators produce a 12-element row-permutation group,
  and the numerical orbit helper describes their action on supplied vectors.
- This does not prove generic completeness of all possible monomial
  automorphisms, classify all third-MUB branches, or prove fourth-MUB
  nonexistence over the full Karlsson family.
- Singular-root counts and their proposed orbit decomposition remain
  numerical until the actual roots and complete corresponding clique set are
  checked under the exact representative map.

The longer-term path is: recover cliques with the existing solver, establish
the representative/gauge conversion, audit A₄ orbits at the selected
parameters and near the fold, then pass inequivalent validated representatives
to the existing fourth-vector witness and certification workflow.

## Connected structural campaign

The fold/A₄ results now feed the repository-native
[`i3_singular_locus`](../../research/campaigns/i3_singular_locus/protocol.md)
campaign. The Fourier-family work is tracked separately in the
[`k3_fourier_structure` campaign](../../research/campaigns/k3_fourier_structure/protocol.md).
It is a parallel analysis of the same validated B₃ cliques:

```text
complete pool -> all B3 cliques
                       |-> W1 fourth-vector certification
                       `-> Fourier/X/F^T transition audit
```

The second branch is currently diagnostic. It becomes a proof route only after
the gauge map, all relevant permutation/orientation charts, and exact
component containment have been established.

## Critical-parameter singular-locus campaign

[`scripts/julia/i3_singular_locus.jl`](../../scripts/julia/i3_singular_locus.jl)
uses the explicit \(H_D(z)\) rows and solves against its **columns**, matching
the symbolic equations above. The pool stage measures the regular physical
roots at \(\lambda_\ast\pm2\cdot10^{-6}\) and at \(\lambda_\ast\). At the
critical value, ordinary path tracking can lose singular endpoints: its
physical pool count is the regular-root count, not the expected total after
the singular roots are restored.

The augmented stage uses the project's square chart of ten variables
\((x_1,\ldots,x_5,y_1,\ldots,y_5)\): five inverse equations and five
independent MU-column equations (the sixth follows from Parseval on the
physical locus). It introduces \(z,t\), imposes \(zt=1\), and adds a null
vector \(u\) for the \(10\times10\) Jacobian in the vector variables, with
the algebraic normalization \(\sum_j u_j^2=1\). The resulting polynomial
system has 22 equations in 22 variables. Because \(z\) is solved as a
variable, the solve can refine the discriminant parameter near the supplied
decimal \(\lambda_\ast\); roots are filtered to the physical unit-circle and
conjugate locus and then deduplicated by vector and parameter.

From the repository root (WSL with the project Julia environment):

```bash
julia --project=. scripts/julia/i3_singular_locus.jl counts
julia --project=. scripts/julia/i3_singular_locus.jl singular
julia --project=. scripts/julia/i3_singular_locus.jl singular --certify
julia --project=. scripts/julia/i3_singular_locus.jl cliques
python scripts/python/analyze_i3_singular.py prepare-certification
julia --project=. scripts/julia/i3_singular_locus.jl certify-saved
python scripts/python/analyze_i3_singular.py singular
python scripts/python/analyze_i3_singular.py cliques
```

Outputs are under `results/campaigns/i3_singular_locus/`. The `singular`
command computes the mixed volume before path tracking and records every
finite augmented-system root, physical-locus residuals, and roots near the
target parameter. `--certify` requests HomotopyContinuation certification;
uncertified numerical roots must not be described as isolated/certified.
Alternatively, `prepare-certification` exports the saved approximations and
`certify-saved` certifies them without repeating the full augmented solve.
The Python analyzer tests the supplied roots for closure under the two
discovered A₄ generators and reports orbit sizes and stabilizer orders without
assuming a \(12+12\) split. A missing image fails the analysis.

The `cliques` command independently recovers and 128-bit verifies the
pool-relative third bases at \(\lambda=0.4,\pi/3,2\pi/3\); the analyzer reports
their \(A_4\) orbits only when the recovered clique set is closed. These are
claims about the supplied recovered pools, not a proof of pool completeness.
The expected 120/96/72 transition is validated only when the regular counts
and augmented singular-root count are combined; a raw 72 at the critical
parameter alone does not confirm the expected total of 96.

## Recorded computation

The seeded \(H_D\)-column pool solves recovered **120** regular physical
vectors at \(\lambda_\ast-2\cdot10^{-6}\), **72** at the supplied critical
parameter, and **72** at \(\lambda_\ast+2\cdot10^{-6}\). The ordinary critical
pool solve reports the regular-root count; it does not include the
singular-root representatives.

The Jacobian-null augmented system recovered 24 distinct physical vectors
near the target parameter. Their mean refined value was
\(\lambda=0.11148022437796644\) at Float64 output precision, and the maximum
augmented residual among the 24 is below \(4\cdot10^{-15}\). Interval
certification reported all 1920 recovered nonsingular augmented-system roots
certified, including all 24 selected physical roots. Their computed A₄
decomposition is **12+12**, with stabilizer order 1 for each orbit. The
regular/singular reconciliation is \(72+24=96\), giving the requested
120/96/72 counts.

This does **not** yet prove that the 24-root list is complete. The augmented
solve tracked 13,096 paths and reported 2,425 singular endpoints; the
explicit mixed-volume estimate was 12,224. The singular endpoints were not
all resolved and certified as distinct physical roots. Certification proves
the recovered roots are isolated solutions of the augmented system, but the
unresolved endpoints/path-count discrepancy prevents a completeness claim
for the physical singular locus.

At each of \(\lambda=0.4,\pi/3,2\pi/3\), fresh solves recovered 72 pool vectors
and four third-MUB cliques. Each four-clique set forms one A₄ orbit of size 4
with stabilizer order 3. The reported cliques passed independent basis checks
and 128-bit orthogonality, normalization, and MU-to-\(H_D\) checks. This
classifies the verified cliques in each recovered pool, not all possible
third-MUB bases.

Raw outputs, certification records, and orbit analyses are preserved in
`results/campaigns/i3_singular_locus/`.
