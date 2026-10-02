# W1 exact Gröbner run — 2026-10-01

## Outcome

The saved exact \(W_1\) witness ideal was run with Macaulay2 1.26.06 in WSL, using a 30-minute wall-clock limit. The process remained CPU-active and reached the timeout without producing the script's first status print or any Gröbner-basis result. The durable log `results/track_c_elimination/m2_w1_pi3_exact_abc_wsl.log` contains only `interrupted, stopping`. This is an inconclusive timeout, not a unit-ideal result and not evidence that a witness exists.

The script's equations were unchanged; its final output was narrowed to the basis element count, the remainder of 1, and ideal dimension rather than dumping the full basis. The input is `symbolic_export/w1_pi3_exact_abc.m2`.

## Scope and prerequisite

The intended calculation concerns one exact, selected triple at \(\lambda=\pi/3\): whether its witness ideal is the unit ideal. Even a valid unit-ideal result would exclude a fourth flat vector only for that fixed triple, not prove a global no-fourth-MUB theorem.

At the time of the timed-out run, the coefficient ring constructed with `toField(S / ideal(...))` had not been verified. A later exact Singular audit established that the joint phase quotient is a degree-96 field over \(\mathbb Q\), using the three irreducible quadratic gcds and independence of their square classes. See `phase2_w1_coefficient_field_2026-10-01.md` for the proof and logs.

## Next step

The coefficient-field gate has now passed, and `symbolic_export/w1_pi3_exact_abc.m2` has been reformulated over that field. The remaining step is to run its Macaulay2 smoke test and then retry the ten-variable witness Gröbner computation. The smoke test could not start because WSL reported insufficient host paging-file resources. The Gröbner computation is not checkpointable; a timeout restarts it from the beginning.

## Artifacts

- Exact input: `symbolic_export/w1_pi3_exact_abc.m2`
- WSL log: `results/track_c_elimination/m2_w1_pi3_exact_abc_wsl.log`
- Setup status: `results/phase2_w1_setup.json`
- Related separate-phase consistency check: `results/phase2_pi3_phase_field_ideal.json`
