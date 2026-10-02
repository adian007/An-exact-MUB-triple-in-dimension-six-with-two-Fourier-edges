# Canonical pi/3 pool rebuild — 2026-10-01

The repository's Julia/HomotopyContinuation path was rerun for `lambda=pi/3`
inside `scripts/julia/i3_singular_locus.jl cliques`. All 252 mixed-volume
paths were tracked. The solver returned 240 finite solutions, retained 72
verified physical MU vectors, and verified four third-basis cliques with
128-bit checks. The maximum recorded Hadamard MU defect was `7.45e-16`.

The 24 vectors across the four rebuilt cliques match the corresponding 24
vectors in the archived pi/3 pool within maximum component error
`5.82e-15`. Machine-readable cross-check:
`results/phase2_pi3_julia_rebuild_audit.json`.

This independently corroborates the restored-pool graph result and the four
known cliques. Homotopy tracking is numerical; this does not certify that the
72-vector pool is complete or exclude other third bases.
