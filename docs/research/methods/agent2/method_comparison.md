# Agent 2: numerical-algebraic-geometry design

Status: **INCOMPLETE / design only**.  No large solve, parameter sweep, Gröbner computation, or witness decomposition was run by this agent.  The only repository inspection was read-only; the prototype below was intentionally not executed.

## Exact target and the indispensable real-locus distinction

For an unnormalised Hadamard matrix `H` and third basis `B`, fix a fourth-vector phase by `z_0=1`, write `v=z/sqrt(6)`, and introduce independent complex variables `w_1,...,w_5`.  The repository's fixed-triple system is

```
z_i*w_i - 1 = 0                                      i=1,...,5
(sum_j conj(H[j,k])*z_j)(sum_j H[j,k]*w_j) - 6 = 0   k=1,...,6
(sum_j conj(B[j,k])*z_j)(sum_j B[j,k]*w_j) - 6 = 0   k=1,...,6,
```

with `z_0=w_0=1`; it has 17 equations in 10 complex variables.  Physical solutions are additionally on the anti-holomorphic locus `w_i=conj(z_i)`.  Treating `w` as independent is a **complexification**: it can introduce nonphysical roots, and no complex solution count equals a physical count without the subsequent conjugacy test.  Conversely, complex emptiness implies physical emptiness.  If `H` or `B` is given only by binary floating point, this is a polynomial system over approximate coefficients, not an exact system for the intended MUB triple.

For the pair-pool system, the local code uses the first five `H` equations plus the five torus equations.  Its omission of the sixth `H` equation is justified only on the physical locus using unitarity; this reduction is not automatically valid on the complexification.  For a complex-algebraic completeness statement retain all six equations, or prove the omitted equation is in the relevant saturated ideal.  This is a required audit item before any theorem claim.

Logical scopes:

* `W1(H,B)=empty` proves that this exact fixed triple has no fourth MUB, since any fourth basis provides a W1 vector.
* A W1 vector proves neither a fourth basis nor a fourth MUB.
* Fixed-parameter/selected-`B` evidence says nothing about other third bases or other CHM parameters.

## Review of present implementation

`src/mub_zauner_6d_liang_chen.jl` builds the pool as a 10-by-10 complexified square system. `src/Certification.jl:w1_witness_certificate` uses a fixed-seed random 10-by-17 matrix to form a square subsystem `A F=0`, mixed-volume path tracking, `HomotopyContinuation.certify`, and an ordinary floating-point full-residual threshold.

The reduction `F=0 => AF=0` is sound.  But a floating full-residual threshold at an approximate square-system root is not, by itself, a certified exclusion of `F=0`. A defensible numerical certificate needs (i) complete accounting of the square-system paths, (ii) an existence/uniqueness enclosure for every finite square root, and (iii) interval evaluation of at least one original equation on each enclosure that excludes zero (or a validated interval/Krawczyk exclusion of `F` on the enclosure). Singular roots and paths at infinity must be separately resolved.  Therefore existing `EMPTY_CERTIFIED` wording should be treated as **not yet established by this agent's standard** until its residual gate is upgraded. This is a methods finding, not a claim that the underlying result is false.

## Methods and selection criteria

The companion CSV is the compact comparison. Sources: HomotopyContinuation's current documentation describes polynomial solving, parameter homotopy, endgame tracking, certification, witness sets, numerical irreducible decomposition and trace tests [HC docs](https://www.juliahomotopycontinuation.org/HomotopyContinuation.jl/v2.14/) and [witness-set documentation](https://github.com/JuliaHomotopyContinuation/HomotopyContinuation.jl/blob/main/docs/src/witness_sets.md).  The numerical-algebraic geometry reference is Sommese--Wampler, *The Numerical Solution of Systems of Polynomials Arising in Engineering and Science* (World Scientific, 2005, ISBN 9789814480888).  The interval certification basis for HC is Breiding, Rose and Telen, “Certifying zeros of polynomial systems using interval arithmetic,” arXiv:2011.05000 (2021).  PHCpack documents polyhedral endgames and extended precision in [release 2.4.93 documentation](https://homepages.math.uic.edu/~jan/PHCpack.pdf); Bertini documents regeneration, deflation and numerical irreducible decomposition in its [manual/examples](https://bertini.nd.edu/BertiniExamples/).  The independently maintained [CertifiedHomotopyTracking.jl project](https://github.com/klee669/CertifiedHomotopyTracking.jl) describes Krawczyk-based certified paths; availability has not been tested here.

`Hom4PS` was not found in the repository or PATH and has no assigned role yet.  It is a possible cross-check only after its exact version, licence, and input conventions are recorded.  Likewise Bertini/Bertini2, PHCpack and alphaCertified were not found/used by this agent.

## Recommended staged pipeline

1. **Exact input boundary.** Store `H,B` over an explicit number field or rational real-coordinate representation. Independently verify `HH^†=BB^†=6I` and `|(H^†B)jk|^2=1` before building W1. For a parameter, introduce real variables with circle constraints or use a separately proved rational parametrisation; do not substitute machine `exp` values into a claimed exact computation.
2. **Geometry diagnosis.** Homogenize/projectivize as appropriate, use regeneration/witness sets, monodromy and trace tests to determine dimensions/components of the complexified variety. Do not use a generic fibre result at a discriminant fibre.
3. **Finite fixed instance.** Form `AF` with recorded seed and exact/certified `A`. Compute BKK bound, track every path with adaptive precision/endgames; record finite, infinite, failed and singular endpoints. Deflate/endgame every unresolved endpoint. A path failure means `INCOMPLETE`.
4. **Certified full-system gate.** For each isolated `AF` root, use an interval root enclosure. Interval-evaluate all 17 equations over the enclosure; to rule out a W1 root, one equation's interval must exclude zero. If an enclosure is inconclusive, subdivide/raise precision or mark incomplete. Separately test the physical conjugation condition only when classifying nonempty roots.
5. **Cross-system reproducibility.** Repeat on an independently encoded system (e.g. PHCpack/Bertini or exact real coordinates), different gamma/random-square seeds, and independently computed BKK support. Disagreement is a stop condition.
6. **Parameter study.** First compute discriminant/critical loci or validated continuation charts. Parameter homotopy is an exploration accelerator, not a special-fibre or universal statement. A continuous theorem requires exact elimination/saturation or certified parameter-box covering, including boundaries.

## Prototype benchmark (not executed)

`prototype/mub_h2_w1.jl` is a 2D physical anchor: for `I,H2`, after `z0=1`, the complexified system is `z*w-1=0`, `(1+z)(1+w)-2=0`, whose two roots are `z=+i,-i` and lie on `w=conj(z)`. It checks the known count and invokes `certify`. It does not test completeness of a dimension-six system; it is an API/infrastructure smoke test only.

Run command, deliberately **NOT_RUN**:

```
julia --project=. research/methods/agent2/prototype/mub_h2_w1.jl
```

The expected mathematical answer is two roots, but no output, solver version, path count, precision, tolerance, or certificate has been asserted because it was not run. `prototype/NOT_RUN.log` is the raw-log placeholder and `prototype/command.txt` records the command.

## Reproducibility record and limitations

Inspected local package declaration: `HomotopyContinuation = "2.22"`, Julia compatibility `1.12` in `Project.toml`; installed binary versions were not established because no executable version command was run. No random seed was used by this agent; a future random square must log its seed, coefficient generation, start-path count, endpoint classification, precision escalation, tolerances, Julia/solver versions, Git commit and input hashes.

This report does not establish root counts, W1 emptiness, a valid four-parameter model, completeness of any third-basis list, or any dimension-six nonexistence result.
