# pi/3 exact remainder verification — 2026-10-01

The corrected distinct-phase ansatz was checked directly with exact SymPy
polynomial arithmetic. For each of the six signed columns, all six cleared
Hadamard-unbiasedness equations reduce to zero modulo:

1. the cyclotomic relation `q^12 - q^6 + 1`,
2. the common degree-24 phase polynomial `P(u)`, and
3. the corresponding phase-specific linear-in-`q` relation `L0`, `L1`, or
   `L2`.

The result is recorded in
`results/phase2_pi3_exact_remainder_verification.json`: 36 of 36 equations
reduce exactly to zero.

The independent pairwise check also reduces all 15 inner products
`O01,O02,...,O45` to zero modulo `q^12-q^6+1`, with no specialization of
`a`, `b`, or `c`.

This verifies the exact MUB identities for the supplied algebraic candidate.
It is not a classification theorem and does not establish that the candidate
is the only solution. The external Macaulay2 Gröbner run remains useful for
ideal dimension, elimination, and component analysis when a CAS runtime is
available.
