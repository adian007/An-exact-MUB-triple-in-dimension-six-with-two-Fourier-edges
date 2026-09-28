# Integration into the Main MUB Repository

The exact Diţă-circle I₃ workflow is integrated using the repository's
existing layout:

- canonical symbolic implementation: `scripts/python/dita_i3.py`;
- regression tests: `test/test_dita_i3.py`;
- generated equations, checks, and run manifest: `symbolic_export/I3/`;
- workflow and evidence boundaries: `docs/research/dita_i3_workflow.md`.

From the repository root:

    python scripts/python/dita_i3.py
    python -m unittest discover -s test -p "test_dita_i3.py" -v

The original package entry point
`mub_complete_computational_package/01_symbolic_I3/build_i3.py` delegates to
the canonical builder, so there is only one implementation of the equations.

## Connection to numerical pool/clique results

The existing Julia pool solver constructs a Karlsson representative, whereas
the uploaded exact symmetry identities refer to the displayed \(H_D(z)\)
representative. The new orbit routines intentionally accept vectors/cliques
in the latter gauge and fail if a supplied clique set is not closed under the
generators. Do not apply them directly to pool vectors until the exact
row/column/dephasing equivalence and its induced vector-coordinate map have
been verified.

After that map is established, enumerate pool-relative third bases with
`enumerate_all_third_mub_bases` in `src/Cliques.jl`, convert each clique to the
exact representative, and pass the complete converted collection to
`a4_clique_orbits`. Results remain relative to the completeness evidence for
the input pool.

The 24-root/two-orbit interpretation, generic completeness of the discovered
monomial automorphisms, and global fourth-MUB nonexistence remain unproved.
