# Verified method bibliography

This is a source-verified starting bibliography for campaigns on MUBs in \(\mathbb C^6\). It contains no new mathematical result. A row's `certification_scope` is intentionally narrower than its use case: numerical algebraic geometry and optimization do not turn an unproved universal quantifier into a theorem.

The machine-readable version is [method_bibliography.csv](method_bibliography.csv). Search provenance is [raw_search_notes.md](raw_search_notes.md).

| ID | Directly verified source | What it supports for this program | What it does **not** establish |
|---|---|---|---|
| MUB-01 | [Brierley--Weigert (2009)](https://arxiv.org/abs/0901.4051) | Algebraic-computational construction of matrices MU to a prescribed Hadamard. | A global exclusion in dimension six. |
| MUB-02 | [Bengtsson et al. (2006)](https://arxiv.org/abs/quant-ph/0610161) | Numerical study of MUB triplets/CHM landscape. | An established four-parameter CHM classification. The paper describes evidence for such a family. |
| MUB-03 | [Tadej--Życzkowski (2006)](https://doi.org/10.1007/s11080-006-8220-2) | CHM definitions, equivalence, and catalogue context. | Completeness of the order-six catalogue. |
| MUB-04 | [Brierley--Weigert (2010)](https://arxiv.org/abs/1006.0093) | SDP, discretization, and Gröbner perspectives. | That an SDP/finite discretization solves the general \(d=6\) question. |
| MUB-05 | [Pál--Vértesi--Navascués (2022)](https://arxiv.org/abs/2203.09429) | Optimization-based discovery objectives and useful low-dimensional checks. | Nonexistence inferred from numerical optima. |
| NAG-01 | [Breiding--Timme (2018)](https://arxiv.org/abs/1711.10911) and [current HC.jl docs](https://www.juliahomotopycontinuation.org/HomotopyContinuation.jl/stable/) | Parameter homotopies, witness sets, monodromy, endgames, and package interfaces. | Completeness without path-count, start-system, and singularity evidence. |
| NAG-03 | [Beltrán--Leykin (2012)](https://doi.org/10.1080/10586458.2011.606184) | Certified tracking of a regular homotopy path. | Enumeration of all paths/components or treatment of singular endpoints by itself. |
| NAG-04 | [Hauenstein--Sottile (2012)](https://arxiv.org/abs/1011.1091) | Alpha-theory certification of approximate roots for square polynomial systems. | Global emptiness; and exact certification of arbitrary overdetermined formulations. |
| NAG-05 | [Sommese--Verschelde--Wampler (2003)](https://doi.org/10.1007/978-3-662-05148-1_6) | Numerical irreducible decomposition/witness-set route for positive-dimensional loci. | Symbolic proof of all real physical points or coverage of a parameter family. |
| NAG-06 | [Bertini manual](https://bertini.nd.edu/BertiniUsersManual.pdf) | Regeneration, witness sets, monodromy/trace test, and deflation facilities. | That a particular user execution is valid without preserved logs/data. |
| NAG-07 | [Verschelde / PHCpack](https://homepages.math.uic.edu/~jan/PHCpack/phcpack.html) | Mixed-volume/polyhedral homotopy and Laurent polynomial workflows. | Automatic handling of every toric-boundary or nongeneric solution. |
| NAG-08 | [CertifiedHomotopyTracking.jl docs](https://klee669.github.io/CertifiedHomotopyTracking.jl/dev/) | Interval/Krawczyk certified tracking and a-posteriori checking. | A substitute for independently benchmarked completeness. |
| INT-01 | [Krawczyk--Neumaier (1986)](https://doi.org/10.1016/0022-247X(86)90303-3) | Interval root-existence/uniqueness/exclusion operators. | A global infeasibility theorem without a finite validated cover. |
| INT-02 | [Neumaier interval Newton chapter](https://doi.org/10.1137/1.9780898717716.ch8) | Formal conditions and limitations for interval Newton methods. | Feasible scaling in high-dimensional near-singular systems. |
| REAL-01 | [Basu--Pollack--Roy](https://doi.org/10.1007/978-3-662-05355-3) | Exact real-algebraic quantifier elimination and semialgebraic decision theory. | Practical feasibility for a full \(6\)-vector/all-third-bases formula. |

## Citation discipline for downstream agents

- Cite MUB-02 only as a numerical study and historical conjectural context, never as verification of a four-parameter order-six CHM family.
- Cite MUB-01/MUB-05 only with their explicit prescribed-matrix and algorithmic scopes.
- A numerically complete zero-dimensional solve needs: the exact polynomial system, homogenization/toric boundary policy, start/root count, tracked path disposition, precision/tolerances, singular-path policy, deduplication, and an independent check.
- A certified tracked path or alpha certificate is local. It is not an emptiness certificate.
- An interval exclusion becomes global only after a rigorously specified compact domain is covered, including boundary strata, by certified exclusion/existence boxes.
