# Brierley-Weigert 2009 reference notes for Dita third-MUB construction.
#
# Paper: S. Brierley and S. Weigert, Phys. Rev. A 79, 062316 (2009)
# "Maximal sets of mutually unbiased quantum states in dimension six"
#
# Key results relevant to T2:
# - Third MUBs exist at special Karlsson/Dita-type CHM parameters
# - Construction uses algebraic (Groebner) methods to find MU vectors to {I, H}
# - Our implementation: extract_third_mub_basis in src/dita_third_mub_construction.jl
#   provides a numerical constructive certificate at each lambda on the Dita slice
#
# Parametric lift strategy:
# 1. H(lambda) depends only on z1=exp(i*lambda) at fixed (theta_D, pi/4) — Theorem T1 L3
# 2. Third clique B3(lambda) extracted from MU-pool at each lambda
# 3. Full-circle HP audit (126/126, 628/628) establishes T2 numerically
# 4. Closed-form B3(lambda) from BW2009 Gröbner basis — open for future symbolic work

const BW2009_REFERENCE = """
Brierley, Weigert, Phys. Rev. A 79, 062316 (2009)
Dita point in Karlsson family: theta=arccos(1/sqrt(3)), phi=pi/4
"""

export BW2009_REFERENCE
