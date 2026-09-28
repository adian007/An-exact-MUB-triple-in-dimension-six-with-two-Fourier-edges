# Research Context and Computational Narrative

## Research question

The project investigates how third mutually unbiased bases occur inside Karlsson's three-parameter family K_6^(3) of complex Hadamard matrices of order six, and what can be established about extension of those third MUBs to a fourth MUB.

The computational narrative is

K_6^(3) -> exact MU-vector ideal -> finite MU pool -> six-vector B3 cliques -> third-MUB loci/branches -> singular pool geometry -> fourth-vector incidence problem -> certified pointwise obstructions.

## Diţă circle

At theta_D=arccos(1/sqrt(3)) and phi=pi/4, the Karlsson construction simplifies to an exact Laurent-polynomial family H_D(z), |z|=1, with z=exp(i lambda).

The exact relation to the canonical Diţă parameter is x=lambda/(2pi). Standard Diţă equivalences induce lambda ~ lambda+pi and lambda ~ -lambda+pi/2. A convenient fundamental interval is lambda in [-pi/4,pi/4].

## Certified seven-point obstruction

The theorem work used the parameter set

Lambda_cert={0,0.4,pi/3,2pi/3,pi,4pi/3,5pi/3}.

At each certified parameter, every size-6 third-MUB clique recovered in the relevant pool was tested using a fourth-vector witness ideal. The reported bridge was N_expected=N_tracked=N_certified=252 and N_full-system=0 for the generic square subsystem used for certification. This is a pointwise result over the recovered-pool cliques and does not establish fourth-MUB absence over all K_6^(3).

## Exact D0 cross-check

At a D0-equivalent point, an exact fourth-vector witness ideal over Q(zeta_24,sqrt(5)) was reported to have Groebner basis {1}. This is an exact machine-checkable cross-check, but it is not a stronger theorem than the known complete-pool D0 result.

## Karlsson transcription audit

The corrected Karlsson A-matrix satisfies (AA^dagger-2I)_11=sqrt(3) sin(phi) sin(theta). The printed McNulty-Weigert transcription used in the audit was found not to be unitary on a Zariski-open subset. This audit belongs to the reproducibility/correction layer and should be kept separate from the main mathematical non-extendability claim.

## Singular pool geometry

A numerical singularity was located at

lambda*=0.1114802243779665542913031975274172717685818097497045...

with x*=lambda*/(2pi)=0.01774262876674699005494250338493207496415495655...

The numerical pool count changed from 120 to 72 across a small neighborhood of lambda*. An augmented system recovered 24 simultaneous singular MU-vector roots. A representative root had singular values approximately (11.73,9.70,8.68,4.65,0), and representative fold coefficients w^T F_lambda approximately -3.01808639127 and w^T F_aa[v,v] approximately -1.19270960464.

These are numerical local diagnostics, not an exact proof of a fold for every singular root and not a proof that the 24 roots form two A4 orbits.

## Exact-vs-numerical boundary

Exact:
- H_D(z) Laurent representation;
- Diţă parameter map x=lambda/(2pi);
- symbolic I3 ideal;
- tested exact symmetry identities;
- exact algebraic calculations explicitly recorded as such.

Numerical/certified:
- finite pool extraction and clique recovery unless an explicit certification record says otherwise;
- continuation around lambda;
- singular-value/fold diagnostics;
- finite optimization probes;
- numerical orbit patterns until exact orbit membership has been checked.

Open:
- global classification of third-MUB loci in all K_6^(3);
- fourth-MUB non-extendability over an entire continuum/component without handling all singularities;
- N(6)=3.
