# Referee risk register

## Scope overstatement

Risk: a result about Karlsson's family is read as a proof of `N(6)=3`.

Response: put the family restriction and the fixed-triple quantifier in every theorem statement.

## Approximate B3 issue

Risk: W1 emptiness is certified for a reconstructed numerical B3 but stated for an exact B3.

Response: provide an epsilon-transfer bound or weaken the wording to the reconstructed triple.

## Incomplete pool issue

Risk: missing homotopy paths are silently treated as no physical vectors.

Response: report path failures, reseed agreement, conjugate filtering, and dedup diagnostics.

## Clique-selection issue

Risk: one arbitrary third basis stands in for all third bases.

Response: enumerate all size-six cliques and certify each, or explicitly state the weaker scope.

## Equivalence issue

Risk: several lambda points are counted as independent despite CHM equivalence.

Response: record the exact row/column/permutation/conjugation map and distinguish parameter classes from CHM classes.

## Literature transcription issue

Risk: a printed A-block is used without verifying `A A^* = 2I`.

Response: retain the local symbolic and numerical transcription audit as a prerequisite.
