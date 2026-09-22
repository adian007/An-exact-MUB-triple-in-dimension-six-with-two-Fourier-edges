# Certificate design and claim tiers

## Valid certificate types

### Exact algebraic

A number-field polynomial system has an exact unit ideal, or an equivalent exact contradiction. This is the strongest fixed-instance certificate and should include the field and saturation conventions.

### Certified numerical

A complete homotopy computation plus interval enclosures and full-system residual filtering can certify the stated numerical object. In the current witness method, the logical chain is:

`W1(H,B3) = empty` implies no fourth MUB extending the fixed triple.

The square subsystem must be generic, all start paths must be tracked, certified endpoints must be accounted for, and the full overdetermined equations must be checked.

### Numerical support

High precision, reseeding, tolerance sweeps, and independent implementations provide strong evidence but do not replace a proof or an interval cover.

## Required metadata

Every certificate record should state:

- exact parameter formula or decimal precision;
- dephasing and equivalence convention;
- pool-completeness evidence;
- all-clique count;
- witness system dimensions and mixed volume;
- random seed and software versions;
- residual thresholds and precision;
- whether the third basis is exact or numerically reconstructed.

## Important logical caveat

A numerical B3 reconstructed from approximate pool vectors introduces an epsilon-transfer issue. W1 emptiness for the computed matrix is not automatically W1 emptiness for the exact intended B3 unless the perturbation is bounded and the exclusion has a margin. The current referee notes correctly identify this as a point that should be explicit in theorem statements.

## Recommended wording

Use `certified numerical obstruction for the recorded reconstructed triple` unless the exact-to-numerical transfer has been proved. Reserve `proved` for analytic or exact algebraic arguments.
