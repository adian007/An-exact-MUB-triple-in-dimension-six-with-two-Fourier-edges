# Applicability assessment: methods for a reproducible MUB program

## Scope and status

Status: LITERATURE_RESULT / planning assessment only. No solver or interval computation was run by Agent 9. The open question is not assumed to have a negative answer.

The physical MUB constraints are conjugate-polynomial in complex entries. For numerical algebraic geometry (NAG), a proposed formulation must first be made honest: either split every complex variable into real and imaginary parts with real polynomial equations, or introduce barred variables only together with equations enforcing the physical conjugation relation in the real-coordinate interpretation. Treating barred symbols as independent complex variables produces a complexification, useful for upper bounds or candidate discovery, but not automatically equivalent to the physical problem.

## Method comparison

| Method | Appropriate input and output | Completeness condition | Main failure modes for MUB work | Recommended role |
|---|---|---|---|---|
| Total-degree/polyhedral homotopy (HC.jl/PHCpack) | Square polynomial system at one fixed exact parameter; approximations to isolated complex roots. | A justified start count (Bézout/mixed volume), every path disposition, and a policy for roots at infinity/toric boundary. | Huge mixed volume; paths diverge; singular endpoints; physical real slice is only a subset. | Fixed-pair MU-vector enumeration and small fixed-triple witness systems. |
| Parameter homotopy | A family with a solved generic base point and fixed support. | Certified generic fiber and tracked paths to every target; exceptional/discriminant targets separately handled. | Root paths merge/diverge at discriminant; a generic result does not cover special fibers. | Transport candidate pools along a rigorously defined selected orbit; never as a blanket universal proof. |
| Witness sets + monodromy + trace test | Positive-dimensional complex algebraic set. | Numerical irreducible decomposition with trace/membership evidence, then separate physical-real analysis. | Monodromy is heuristic unless its stopping/trace evidence is recorded; singular/reducible components. | Map third-basis or MU-vector loci and identify branches for later exact/interval treatment. |
| Regeneration / random squaring | Overdetermined polynomial systems or incremental equations. | Same component/root accounting after randomization; generic linear combinations must be logged. | Random squaring can alter interpretation; omitted components; nonphysical complex roots. | Prototype only, paired with original-equation residual checks. |
| Endgames and deflation | Singular or multiple endpoints detected by tracking. | Deflation system and multiplicity/rank evidence are recorded and independently checked. | Numerical rank thresholds; deflation may change geometry; accuracy loss. | Mandatory contingency, not evidence of absence after failure. |
| Certified path tracking (Beltrán--Leykin; interval a-posteriori tools) | A known regular start solution and a specified homotopy path. | Certificate succeeds for the full path under stated hypotheses. | Cannot create missing start roots; singular paths lie outside regular-path assumptions. | Upgrade selected numerical traces to rigorous local continuation evidence. |
| Alpha theory / alphaCertified | Approximate isolated root of square polynomial system. | Alpha criterion proves Newton convergence to a root; distinctness/realness are separately checked as appropriate. | Local only; does not count all roots or prove emptiness; overdetermined support is limited/heuristic. | Certify a discovered candidate or endpoints returned by a complete fixed-fiber solve. |
| Interval Newton / Krawczyk B&B | Compact real box after exact real-coordinate formulation. | Outward-rounded interval arithmetic; every box excluded, contracted to a unique root, or subdivided; boundaries included. | Dependency blowup, conditioning, noncompact phase variables, near-singular gauges. | Rigorous exclusion/existence on small gauge-fixed parameter boxes and boundary strata. |
| Real-algebraic QE/CAD | Exact finite real polynomial equalities/inequalities with quantified variables on compact or explicitly bounded domains. | Correct first-order formula over a real closed field; exact arithmetic and output certificate/replay. | Doubly/exponentially difficult in dimensions/degree; algebraic coefficients and unit-circle parametrization enlarge formulas. | Small fixed triple or low-dimensional restricted orbit; theoretical gold standard, not default. |
| SDP/SOS / moment relaxations | Polynomial optimization/feasibility encoded in moments with compactness/Archimedean assumptions. | A rigorously verified dual Positivstellensatz/SOS certificate with exact/rational reconstruction. | Floating SDP infeasibility alone is not proof; hierarchy may not terminate; rank extraction may fail. | Discovery bounds and potential exact certificate search after independent symbolic verification. |
| Riemannian/unitary or phase-torus optimization | Smooth residual objective over unitary matrices/phases. | None for a global conclusion unless linked to a separate global certificate. | Local minima, gauge duplicates, parametrization bias, finite restarts. | Positive-candidate discovery and benchmark checks only. |

The NAG items are supported by [Breiding--Timme](https://arxiv.org/abs/1711.10911), [PHCpack](https://homepages.math.uic.edu/~jan/PHCpack/phcpack.html), [Bertini](https://bertini.nd.edu/BertiniUsersManual.pdf), and [Sommese--Verschelde--Wampler](https://doi.org/10.1007/978-3-662-05148-1_6). Local certification limits follow [Beltrán--Leykin](https://doi.org/10.1080/10586458.2011.606184) and [Hauenstein--Sottile](https://arxiv.org/abs/1011.1091). Interval and exact-real methods follow [Krawczyk--Neumaier](https://doi.org/10.1016/0022-247X(86)90303-3), [Neumaier](https://doi.org/10.1137/1.9780898717716.ch8), and [Basu--Pollack--Roy](https://doi.org/10.1007/978-3-662-05355-3).

## Suggested staged pipeline (not yet executed)

1. **Specification gate.** Fix one pair or triple of *exact* unnormalized matrices. Independently check \(MM^\dagger=6I\), \(NN^\dagger=6I\), and \(|(M^\dagger N)_{jk}|^2=1\). Record coefficient field, gauges, and physical real-variable equations.
2. **Anchor gate.** On dimensions 2--4 and a known fixed dimension-six anchor, compare two independent implementations. A discrepancy stops the campaign.
3. **Geometry gate.** For an apparently positive-dimensional system, compute witness data before pretending an isolated-root solver is complete. Use monodromy/trace diagnostics and retain the raw path data.
4. **Fixed-fiber gate.** For a selected fixed triple, use a justified sparse/total-degree bound or exact elimination. Path failures, divergent paths, and singular endpoints are findings requiring analysis, not zero roots.
5. **Certification gate.** Certify discovered isolated roots locally with alpha theory/interval methods. To prove emptiness on a compact real domain, use an exhaustive validated cover, exact elimination, or an independently checkable algebraic certificate.
6. **Parameterized-orbit gate.** Prove the orbit itself is a valid MUB triple, compute/test discriminant and denominator strata, and certify continuation only on regular cells. Directly treat every exceptional fiber. The resulting theorem is about that orbit, not all third bases.
7. **Universal-quantifier gate.** Do not promote an orbit theorem to a slice/family theorem until the representation is proved exhaustive or the all-third-bases formula is actually decided.

## Promising directions with precise cautions

- **Sparse/Laurent formulation.** Unit-modulus phase variables naturally suggest a torus description. Polyhedral homotopy can exploit sparsity, but the physical condition \(\bar z=z^{-1}\) must be reflected correctly and toric-boundary solutions cannot be silently discarded.
- **Numerical irreducible decomposition of third-basis loci.** This can reveal whether candidate third bases lie on curves/surfaces rather than in finite pools. It is a discovery map, not a proof that all physical third bases were found, unless numerical decomposition and real-locus analysis are completed.
- **Validated continuation on gauge-fixed compact cells.** This is a plausible bridge from point results to a selected arc/box. It needs an invertible branch Jacobian on each regular cell, explicit overlap/coverage data, and separate singular/boundary treatment.
- **Exact/fixed-field certificates.** If entries and a chosen triple lie in a modest cyclotomic field, Gröbner saturation/elimination may give a replayable fixed-triple emptiness proof. A generic calculation must not be reused at special parameter values without specialization or valid saturation.
- **SDP/SOS.** The MUB-SDP literature makes this reasonable for lower bounds and candidate restrictions. A floating-point infeasibility report is only exploratory; pursue it toward an exact dual certificate only if the chosen relaxation has an explicit theorem linking infeasibility to the target scope.
- **Frame/Gram/operator formulations.** Gram matrices and rank/PSD constraints expose symmetry and may yield lower-dimensional relaxations, but any relaxation must retain enough rank and phase information; otherwise an infeasibility result may concern a stronger surrogate and a feasibility result may concern a weaker one.

## Literature-based no-go statements

- No source above verifies that Karlsson's three-parameter family has a justified four-parameter extension.
- No source above proves either existence or nonexistence of four MUBs in \(\mathbb C^6\).
- The numerical optimization approaches found are useful for conjecture generation and regression tests, not standalone fourth-MUB nonexistence proofs.
- Neither numerical genericity nor a generic fiber certificate covers a discriminant/exceptional fiber without additional argument.
