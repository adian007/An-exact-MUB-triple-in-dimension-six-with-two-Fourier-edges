# Interval arithmetic and parameter covering

## Fixed point versus region

Interval certification at isolated parameter values does not cover an interval of parameters. To cover a region R, one needs either a symbolic positivity/inequality argument or a finite interval cover whose boxes collectively exclude W1.

## Practical covering scheme

For a parameter box X and vector variables Y:

1. validate that the Karlsson construction is defined throughout X;
2. bound CHM and B3 coefficient variation over X;
3. subdivide X near poles, root collisions, and clique changes;
4. use interval Newton/Krawczyk exclusion or a certified homotopy continuation from a box interior;
5. retain a numerical margin separating every candidate endpoint from the full witness equations;
6. recurse until every box is excluded or promoted to an exact/algebraic exceptional stratum.

## Why the current W1 approach helps

The witness is overdetermined: 17 equations in 10 variables. A generic square projection gives a manageable mixed volume of 252 in the local implementation. This is suitable for fixed-point certification and may be suitable for small parameter boxes, but the projection must remain generic over the whole box.

## Main hazards

- interval dependency can become too wide near Dita degeneracies;
- a moving B3 is not naturally an independent fixed coefficient matrix;
- clique identities can change when pool vectors collide;
- parameter denominators create false boxes unless excluded explicitly;
- a zero residual at finite precision is not an interval exclusion.

## Research recommendation

Do not try to cover all `K_6^(3)` immediately. First prove an interval exclusion theorem on a small compact subarc of the Dita circle away from the four special constellation-change windows. This would test whether the local certified point method has a viable uniform upgrade.
