# Closed-form Dita / Brierley–Weigert objects for algebraic T3 work.
#
# Primary source (Groebner / closed-form pool at D0):
#   S. Brierley and S. Weigert, Phys. Rev. A 79, 052316 (2009)
#   arXiv:0901.4051  "Constructing mutually unbiased bases in dimension six"
#
# Numerical-constellation paper (different article; do not conflate):
#   S. Brierley and S. Weigert, Phys. Rev. A 78, 042312 (2008)
#   "Maximal sets of mutually unbiased quantum states in dimension six"
#
# Block-circulant transcription:
#   Bengtsson et al., arXiv:quant-ph/0610161, §7 (eqs. 11, 61, 68, 74, 79–80).
#   BW note that the last two list entries in that §7 must be swapped.
#
# Convention: UNnormalised CHMs (|H_ij|=1, H H† = 6 I), matching the rest of this repo.
# Bengtsson sometimes writes 1/√6; we drop that factor.
#
# This file does not depend on HomotopyContinuation.

using LinearAlgebra

const BW2009_REFERENCE = """
Brierley, Weigert, Phys. Rev. A 79, 052316 (2009); arXiv:0901.4051
Constructing mutually unbiased bases in dimension six
Dita D0: 120 MU vectors, 10 third bases, none pairwise MU.
Phases φ_D = {0, π, ±kπ/12, ±α} with tan α = 2.
"""

const BW2008_CONSTELLATION_REFERENCE = """
Brierley, Weigert, Phys. Rev. A 78, 042312 (2008)
Maximal sets of mutually unbiased quantum states in dimension six
(numerical constellations; not the Groebner/closed-form paper)
"""

"""Unnormalised 3×3 Fourier matrix (entries sixth/third roots)."""
function fourier3_unnorm()
    ω = cis(2π / 3)
    return [ω^(j * k) for j in 0:2, k in 0:2]
end

"""
Diţă affine family D(x), Bengtsson quant-ph/0610161 eq. (11), unnormalised.
z = exp(2π i x). Textbook D0 is D(0) (fourth roots of unity).
"""
function dita_D(x::Real)
    z = cis(2π * x)
    zb = conj(z)
    i = im
    return ComplexF64[
        1  1   1      1   1     1
        1 -1   i     -i  -i     i
        1  i  -1    i*z -i*z   -i
        1 -i  i*zb   -1   i   -i*zb
        1 -i -i*zb    i  -1    i*zb
        1  i  -i   -i*z  i*z   -1
    ]
end

dita_D0() = dita_D(0)

"""α with tan α = 2, taken in (0, π/2). Then e^{iα} = (1+2i)/√5."""
cis_alpha_tan2() = (1 + 2im) / sqrt(5)

"""e^{iθ} with tan θ = -2 in (−π/2, π/2): (1−2i)/√5 (Bengtsson b₂)."""
cis_beta_tan_m2() = (1 - 2im) / sqrt(5)

"""
BW 0901.4051 phase set φ_D: 24th-root phases together with ±α, tan α = 2.
Returned as unimodular ComplexF64 values (the phases themselves).
"""
function phase_set_phi_D()
    roots24 = [cis(k * π / 12) for k in 0:23]
    α = cis_alpha_tan2()
    return vcat(roots24, [α, conj(α)])
end

"""Fourier affine family F(x1,x2), Bengtsson eq. (4), already unimodular."""
function fourier_F(x1::Real, x2::Real)
    q = cis(2π / 6)
    z1 = cis(2π * x1)
    z2 = cis(2π * x2)
    return ComplexF64[
        1  1         1         1   1         1
        1  q*z1      q^2*z2    q^3 q^4*z1    q^5*z2
        1  q^2       q^4       1   q^2       q^4
        1  q^3*z1    z2        q^3 z1        q^3*z2
        1  q^4       q^2       1   q^4       q^2
        1  q^5*z1    q^4*z2    q^3 q^2*z1    q*z2
    ]
end

"""
Twisted-product form F_D, Bengtsson eq. (61), unnormalised CHM:
    [ F3 ,  F3  ]
    [ F3*D, -F3*D ]
D is a 3×3 diagonal unitary.
"""
function fourier_F_D(Ddiag::AbstractVector{<:Number})
    F3 = fourier3_unnorm()
    D = Diagonal(ComplexF64.(Ddiag))
    top = hcat(F3, F3)
    bot = hcat(F3 * D, -(F3 * D))
    return vcat(top, bot)
end

"""
Third MUB for the block-circulant form of D0 (Bengtsson eqs. 78–80):
F_D with D = diag(ω^9 b₂, 1, 1), ω = exp(2π i/24), tan(2π c₂) = −2.
The triplet is {I, F_D, D_bc} with D_bc ≈ D(0).
"""
function third_mub_F_D_at_D0()
    ω = cis(2π / 24)
    b2 = cis_beta_tan_m2()
    return fourier_F_D([ω^9 * b2, 1, 1])
end

"""
Block-circulant matrix equivalent to D(0), Bengtsson eqs. (79)–(80).
ω = exp(2π i/24). Unnormalised.
"""
function dita_D_bc()
    ω = cis(2π / 24)
    C3 = ComplexF64[
        1     ω^6  ω^6
        ω^6   1    ω^6
        ω^6   ω^6  1
    ]
    C4 = ComplexF64[
        ω^15  ω^3  ω^3
        ω^3   ω^15 ω^3
        ω^3   ω^3  ω^15
    ]
    top = hcat(C3, C4)
    bot = hcat(C4, -im * C3')
    return vcat(top, bot)
end

"""
3×3 circulant blocks in Bengtsson eq. (74) at D(1/8).
ar5iv HTML renders C1[3,3] as 0; the matrix is circulant, so that entry is 1.
"""
function bengtsson_C1_C2()
    ω = cis(2π / 24)
    C1 = ComplexF64[
        1      ω^11  ω
        ω      1     ω^11
        ω^11   ω     1
    ]
    C2 = ComplexF64[
        1     ω^5  ω^7
        ω^7   1    ω^5
        ω^5   ω^7  1
    ]
    return C1, C2
end

"""
Four unnormalised CHMs of eq. (68) style: third MUBs for {I, F(1/6,1/12)},
each equivalent to D(1/8). Not claimed to be B3 at D0.
"""
function bengtsson_eq68_third_mubs()
    C1, C2 = bengtsson_C1_C2()
    H1 = vcat(hcat(C1, C2), hcat(C2', -C1'))
    H2 = vcat(hcat(C2, C1), hcat(C1', -C2'))
    H3 = vcat(hcat(C1', C2'), hcat(C2, -C1))
    H4 = vcat(hcat(C2', C1'), hcat(C1, -C2))
    return (H1, H2, H3, H4)
end

function chm_defect(H; tol=1e-10)
    n = size(H, 1)
    unit = norm(H * H' - n * I(n))
    modu = maximum(abs.(abs.(H) .- 1))
    return (unitary=unit, modulus=modu, ok=unit < tol && modu < tol)
end

"""MU defect of two unnormalised CHMs (column bases). Zero iff MU."""
function mu_defect_chm(H1, H2)
    n = size(H1, 1)
    G = (H1' * H2) / n
    return maximum(abs.(abs2.(G) .- 1 / n))
end

export BW2009_REFERENCE, BW2008_CONSTELLATION_REFERENCE,
       dita_D, dita_D0, dita_D_bc,
       fourier_F, fourier_F_D, fourier3_unnorm,
       third_mub_F_D_at_D0, bengtsson_C1_C2, bengtsson_eq68_third_mubs,
       phase_set_phi_D, cis_alpha_tan2, cis_beta_tan_m2,
       chm_defect, mu_defect_chm
