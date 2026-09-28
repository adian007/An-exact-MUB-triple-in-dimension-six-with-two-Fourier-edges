# MUB Solution Attack Package: Karlsson/Diţă -> MUB Triplet Structure -> Fourier Transition

## Purpose
This package records the strongest valid computational route reached so far toward the dimension-6 MUB problem. It is intentionally conservative: it does **not** claim a proof of N(6)=3.

## Main current result
At lambda = 0.4, pi/3, and 2*pi/3 on the exact Diţă slice of Karlsson K6^(3), the recovered third-MUB pool contains 72 MU vectors and 4 complete B3 cliques. An independent graph reconstruction gives exactly four size-6 cliques at each tested parameter. The previously verified order-12 monomial symmetry group acts transitively on the four cliques at these tested points.

For lambda = 0.4 and lambda = 2*pi/3, exhaustive finite permutation matching of a selected B3 transition matrix against the two-parameter Fourier family produced numerical equivalence residuals of approximately 2.5e-15 and 5.4e-15, respectively. At lambda = pi/3, the tested representative did not pass that particular Fourier-family orientation test; another transition must be checked.

This is evidence consistent with the 2025/2026 Hadamard-cube program: a MUB triplet should have transition matrices related to a Szollosi X-family matrix, a Fourier-family matrix, and a transposed Fourier-family matrix. That statement remains conjectural in the cited 2026 work.

## Critical literature status
A 2021 theorem claiming an exact 2x2-Hadamard-count restriction for H2-reducible matrices and using it to exclude Dita-family cases was challenged in a 2025 Comment. The Reply did not restore the disputed theorem and explicitly states that Theorem 2 and Lemma 1 could not be saved at present. Therefore this package does not use that theorem as a proof bridge.

The 2026 review still treats the dimension-6 four-MUB problem as open. The 2026 classification of order-six complex Hadamard matrices supplies an exhaustive classification of CHM classes, but does not by itself solve MUB extension.

## What is now the cleanest proof target?
Prove, for the relevant K6^(3) third-MUB component, that the complete B3 incidence variety is contained in the known Fourier/X/FT triplet structure. Then invoke the rigorous theorem that a Fourier-family transition matrix cannot occur in a MUB quartet. The missing step is the K6^(3)-specific structural proof, not another random numerical fourth-vector search.

## Exact vs numerical
EXACT: the Diţă reduction of Karlsson, the polynomial I3 ideal, the verified order-12 A4 automorphisms, and the algebraic definitions of the incidence systems.
NUMERICAL: pool counts, clique enumeration, fold roots, singular-root orbits, and Fourier-family membership tests of sampled matrices.
CERTIFIED: prior pointwise fourth-vector obstructions at the seven previously certified lambda points from the earlier research package.
NOT PROVED: continuum non-extendability of the full Diţă circle; global non-extendability of all K6^(3); N(6)=3.

## Provenance
See the ZIP files in 04_provenance/ for complete earlier source trees and result logs.
