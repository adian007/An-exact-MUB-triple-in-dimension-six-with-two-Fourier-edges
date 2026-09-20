# Known Results Table

| Result | Source | Exact statement | Method | Scope | Status |
|---|---|---|---|---|---|
| Universal bounds | Bengtsson et al., Sec. 1; Grassl, Sec. 4 | `3 <= N(6) <= 7` | General construction and dimension bound | Global | Rigorous |
| Zauner conjecture | Grassl, Sec. 4; Zauner thesis is UNVERIFIED here | `N(6)=3` | Conjecture | Global | Open |
| Heisenberg pair maximality | Grassl, Theorem 2 and Corollary 1 | 48 vectors, 16 third bases, no vector unbiased to any resulting triple | Exact MAGMA algebra | Heisenberg/Jacobi-equivalent pairs | Rigorous local theorem |
| Fourier pair maximality | Brierley--Weigert, Sec. 3.3; Grassl | 48 vectors and 16 third bases; no two third bases MU | Gröbner/algebraic computation | `{I,F6}` and equivalent cases | Rigorous for that pair |
| D0 nonextendability | Brierley--Weigert, Sec. 4.1, Table 1, Eq. (20) | 120 vectors, 10 third bases, no pair MU; no four MUBs containing D0 | Exact closed phase set and algebra | `{I,D0}` | Rigorous for D0 |
| Diţă family sampled behavior | Brierley--Weigert, Sec. 4.2, Table 2 | Sampled Diţă matrices have 48/72/120 vectors and four third bases; no four in the sampled study | Exact/approximate depending on point | Sampled Diţă parameter set | Not a full-family theorem |
| Karlsson family constraint | Maxwell--Brierley, Theorem 4.1 | `g_j(rho)=0` for specified permutations `rho` | Analytic calculation | Karlsson family | Rigorous structural result |
| Karlsson-only LP attempt | Maxwell--Brierley, Secs. 5-6 | LP did not yield contradiction | Numerical LP | Complete sets built from Karlsson matrices | Inconclusive |
| K6^(3) coverage | Karlsson, Theorem 11, Sec. 7 | Every H2-reducible order-six CHM is equivalent to a member of K6^(3) | Parametrization theorem | H2-reducible subclass | Rigorous, restricted |
| Full CHM classification | Wuttig--Tindall v2, Theorem 6 and Corollary 25 | Claims an exhaustive finite-corner classification and `H6 = G6^(4) union K6^(3) union T6` | Exact/Lean-audited preprint claim | All order-six CHMs; independently unverified in this repository | Current 2026 classification claim |
| Repository T3 | `paper/proofs/fourth_mub_obstruction.tex`; result logs | Seven listed Dita points, all recovered-pool 6-cliques: `W1` certified empty, 40/40 | Certified numerical homotopy/interval workflow | Finite local set in K6^(3) | Local certified-numerical |
| Repository T4 | `results/track_c_elimination/w1_D0_groebner.log` | `W1(D_bc,F_D)` has exact unit ideal over `Q(zeta24,sqrt5)` | Exact Macaulay2 Gröbner basis | One D0-equivalent representative and one third basis | Exact local re-verification, not new D0 mathematics |
