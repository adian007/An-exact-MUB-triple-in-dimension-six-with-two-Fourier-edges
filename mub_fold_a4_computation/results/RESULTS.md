# Preliminary computational results

A 1200-start phase-coordinate least-squares search at the critical parameter found **96 distinct numerical MU-vector solutions**. The smallest singular values of the 6x5 MU Jacobian split sharply: the first 24 candidates have sigma_min below 1.4e-6, while the next candidates are around 1.12 or larger. This is strong numerical evidence for **24 singular/double-root candidates plus 72 regular roots** at the critical parameter.

This is still not a singular-root certificate. The 24 candidates must be refined with a square augmented singular system/deflation and then subjected to an appropriate validation procedure.

The earlier 400-start runs produced incomplete pools (114/93/70), so those files are retained only as exploratory diagnostics. The 1200-start critical run is the current primary numerical result.
