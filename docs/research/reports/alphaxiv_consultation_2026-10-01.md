# AlphaXiv consultation — 2026-10-01

## Execution-relevant findings

- Fatah, McLoughlin, and Ghafoor, [A Reproducible Software Workflow for
  Unanchored Approximate MUB Optimization: A Case Study in Dimension Six](https://www.alphaxiv.org/abs/2607.10615), supports the current provenance-first
  workflow: retain seeds, precision, pairwise defects, and configuration
  metadata. It also cautions that fixing too many bases can exclude parts of
  the search space.
- Matolcsi, Matszangosz, Varga, and Weiner, [Triplets of Mutually Unbiased
  Bases](https://www.alphaxiv.org/abs/2503.14752), introduces Hadamard cubes and
  an algebraic identity conjectured to imply the dimension-six upper bound of
  three MUBs. This is a promising later route for classifying the exact
  triple, but it does not replace the current candidate verification.
- Sarkar, [Degree-Four Vector-Coordinate SoS Cannot Detect the MUB Upper
  Bound](https://www.alphaxiv.org/abs/2606.13903), shows that degree-four
  vector-coordinate SoS is too weak for this problem, while a centered
  projector-coordinate Gram formulation recovers the elementary bound.

## Consequence for this repository

The immediate execution order remains:

1. verify the corrected distinct-phase exact candidate;
2. run the exported ideal in a CAS for component/dimension information;
3. only then investigate Hadamard-cube and projector-coordinate formulations
   as classification or non-extension tools.

The first item is now complete: the exact remainder audit reports 18/18
Hadamard equations equal to zero. The second item is still pending because no
local CAS runtime is installed and the Docker daemon is unavailable.
