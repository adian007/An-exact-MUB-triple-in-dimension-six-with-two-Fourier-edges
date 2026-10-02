# pi/3 restored-pool graph audit — 2026-10-01

The restored pool contains 72 unit-modulus, norm-\u221a6 vectors. Using the
orthogonality relation `|<v_i,v_j>| < 1e-9`:

- the orthogonality graph has 226 edges;
- its maximum clique size is exactly 6;
- exhaustive maximal-clique enumeration finds exactly the four stored
  six-cliques;
- the largest orthogonality-edge residual is `9.76e-10`.

The pool graph contains no pairs classified as mutually unbiased under the
`|<v_i,v_j>| = sqrt(6)` threshold. The pool is a candidate set for constructing
orthogonal B3 bases, rather than a union of multiple MUBs.

The result is recorded in `results/phase2_pi3_pool_graph_audit.json`. It applies
only to the restored 72-vector pool and is not a completeness theorem for all
possible MU vectors.
