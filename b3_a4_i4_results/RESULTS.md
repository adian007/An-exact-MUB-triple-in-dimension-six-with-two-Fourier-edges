# B3 A4 quotient and reduced I4 computation

## Tested parameters
- lambda = 0.4
- lambda = pi/3
- lambda = 2*pi/3

## Third-MUB pool
At each tested parameter, the numerical fixed-H MU-vector recovery returned 72 distinct vectors. The six-vector orthogonality graph contained exactly four size-6 cliques under the numerical tolerance used by the script.

## A4 quotient
The verified fixed-z monomial generators were applied to every vector in each clique. At each tested parameter, both generators permuted the four recovered B3 cliques, and the generated action is transitive on those four cliques. Thus the four recovered B3 bases form a single numerical A4 orbit at each tested parameter.

Consequently, for the fourth-vector search, one representative B3 per parameter is sufficient within this recovered numerical pool, because the remaining three representatives are related by the tested A4 symmetry.

This is a symmetry reduction of the recovered numerical pool, not a global theorem about all third-MUB bases in K_6^(3).

## Reduced fourth-vector witness
For a representative B3 = {v_1,...,v_6}, the fourth vector is w=(1,q_1,...,q_5)/sqrt(6), |q_j|=1. The numerical witness system used here is

  |w^* h_k(lambda)|^2 = 6,  k=1,...,6,
  |w^* v_j|^2 = 6,       j=1,...,6.

The search used nonlinear least squares over the five phase variables of w from multiple random starts. No zero-residual fourth vector was found in these probes. Best infinity-norm residuals for one representative per parameter were:

- lambda = 0.4:     0.4888464091
- lambda = pi/3:    0.5190991753
- lambda = 2pi/3:   0.6352006073

These residuals are numerical optimization evidence only. They do not prove non-extendability.

## Algebraic I4 formulation
The exact symbolic formulation should introduce

  z t - 1 = 0,
  q_j r_j - 1 = 0,

and the MU equations to the columns of H_D(z) and to each of the six B3 vectors. For exact B3 coefficients this gives the incidence ideal I4. The present lambda=0.4, pi/3, 2pi/3 B3 vectors are stored numerically, so the generated witness systems are numerical instances rather than exact coefficient ideals.

## What this accomplishes
The main computational gain is a symmetry quotient: at the three tested Dita-circle parameters, four recovered third-MUB bases collapse to one A4 orbit. Therefore the fourth-MUB extension problem needs only one representative B3 per parameter within the recovered pool.

The next rigorous step is to replace the numerical B3 coefficients by exact/algebraic descriptions (where available) or certified high-precision representations, then perform parameter homotopy / elimination on the reduced I4 incidence system.
