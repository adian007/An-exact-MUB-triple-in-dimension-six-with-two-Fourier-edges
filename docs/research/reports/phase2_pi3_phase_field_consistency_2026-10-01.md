# pi/3 phase-field consistency check — 2026-10-01

Exact univariate elimination was applied to the three pairs `(P(a),L0(q,a))`,
`(P(b),L1(q,b))`, and `(P(c),L2(q,c))`. Each resultant in `q` has degree 24.
Their successive gcd with `Phi_36(q) = q^12-q^6+1` has degree 12, so the
phase ideal is consistent over the algebraic closure for each primitive
36th-root embedding.

The machine-readable result is
`results/phase2_pi3_phase_field_ideal.json`. This establishes consistency of
the phase relations used by the corrected export; it does not analyze the W1
witness variables, prove ideal dimension, or establish W1 emptiness.
