# Certified fold singular roots — 2026-10-01

The saved augmented-system certification completed in WSL with Julia 1.12.6
and HomotopyContinuation 2.22.1. All 1,920 supplied isolated-root
approximations were certified at up to 512-bit precision. A coordinate-level
cross-match maps all 192 stored physical candidates uniquely to that certified
root list with zero difference in the serialized coordinates.

Exactly 24 physical candidates have parameter angle within `1e-8` of the
recorded fold value
`lambda_* = 0.1114802243779665542913031975274172717685818097497045`;
their maximum angular difference is `5.42e-16`. The independent refined
singular-root dataset also contains 24 candidates, and all 24 match distinct
certified roots within maximum phase error `3.45e-10`. The two recorded A4
orbits each map to 12 distinct certified roots; their source matching errors
are below `8.51e-10`.

Evidence is in `results/campaigns/i3_singular_locus/root_certification.json`
and `results/phase2_fold_certified_singular_roots.json`. This certifies
isolated roots of the augmented Jacobian-null system and supports the 24-root
count near the numerically specified fold parameter. It does not prove exact
equality with the rounded decimal `lambda_*`, completeness of the physical
pool, or maximality of the A4 symmetry group.
