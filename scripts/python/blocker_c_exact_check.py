"""
blocker_c_exact_check.py  --  EXACT-ARITHMETIC verification (no floating point).

Verifies the hand derivation in research/reports/audit_followup_2026-09-26.md.

  (1) |A11|^2, |A12|^2, |B11|^2, |B12|^2 as exact functions of (theta, phi)
  (2) D_A = |alpha_A|^2-|beta_A|^2 =  2*r3*sin(th)*(sin(ph) - r3*cos(th)*cos(ph))
      D_B = |alpha_B|^2-|beta_B|^2 = -2*r3*sin(th)*(sin(ph) + r3*cos(th)*cos(ph))
      i.e. exactly proportional to Karlsson's two published Mobius degeneracy
      factors (arXiv:1003.4177 Sec. 5).
  (3) Dita slice cos(th)=1/sqrt(3), phi=pi/4: A, B, D_A, D_B, A non-Hermitian,
      M_A == -1, M_B(z) == 1/z, z2=z3=i, z4=1/z1, z2-recovery denominator != 0.
  (4) theta=0 (F6 seam): M_A == M_B == 1, D_A = D_B = 0, A non-Hermitian.
  (5) Hermiticity locus of A; the two degeneracy factors at both anchors.

Exact identity check (unit-test level). NOT a new numerical experiment and
NOT an interval-arithmetic certification. No floating point anywhere.
"""
import sympy as sp

I = sp.I
th, ph = sp.symbols('theta phi', real=True)
c, s, cp, sp_ = sp.cos(th), sp.sin(th), sp.cos(ph), sp.sin(ph)
r3 = sp.sqrt(3)

# A11 = (-1/2 + r3/2*s*sp_) + i*(r3/2)(c + s*cp)
# A12 = (-1/2 - r3/2*s*sp_) + i*(r3/2)(-c + s*cp)
# B11 = -1 - A11 ,  B12 = -1 - A12
reA11, imA11 = -sp.Rational(1, 2) + r3 / 2 * s * sp_, r3 / 2 * (c + s * cp)
reA12, imA12 = -sp.Rational(1, 2) - r3 / 2 * s * sp_, r3 / 2 * (-c + s * cp)
T = sp.trigsimp
m11 = T(sp.expand(reA11**2 + imA11**2))
m12 = T(sp.expand(reA12**2 + imA12**2))
mb11 = T(sp.expand((-1 - reA11)**2 + imA11**2))
mb12 = T(sp.expand((-1 - reA12)**2 + imA12**2))
print("== (1) exact entry moduli ==")
for nm, v in (("|A11|^2", m11), ("|A12|^2", m12), ("|B11|^2", mb11), ("|B12|^2", mb12)):
    print("  ", nm, "=", v)
print("   sumA-2 =", T(sp.simplify(m11 + m12 - 2)), " sumB-2 =", T(sp.simplify(mb11 + mb12 - 2)))

# D = |alpha|^2-|beta|^2 = |entry|^4 difference = (x-y)(x+y); here x+y = 2.
D_A = T(sp.simplify(2 * (m12 - m11)))
D_B = T(sp.simplify(2 * (mb12 - mb11)))

# ---- (3) Dita slice ----------------------------------------------------
cD, sD, h = 1 / r3, sp.sqrt(2) / r3, 1 / sp.sqrt(2)
subD = {c: cD, s: sD, cp: h, sp_: h}
A11D = sp.simplify(reA11.subs(subD) + I * imA11.subs(subD))
A12D = sp.simplify(reA12.subs(subD) + I * imA12.subs(subD))
print("== (3) Dita slice (cos th=1/sqrt3, phi=pi/4) ==")
print("   A11 =", A11D, " A12 =", A12D)
print("   A = [[i,-1],[-1,i]] ?", A11D == I and A12D == -1)
print("   B11 =", -1 - A11D, " B12 =", -1 - A12D)
print("   B = [[-1-i,0],[0,1-i]] ?", (-1 - A11D) == -1 - I and (-1 - A12D) == 0)
print("   A Hermitian? A11-conj(A11) =", sp.simplify(A11D - sp.conjugate(A11D)),
      "->", A11D == sp.conjugate(A11D))
alA, beA = sp.simplify(A12D**2), sp.simplify(A11D**2)
alB, beB = sp.simplify((-1 - A12D)**2), sp.simplify((-1 - A11D)**2)
print("   alpha_A =", alA, " beta_A =", beA,
      " D_A =", sp.simplify(alA*sp.conjugate(alA) - beA*sp.conjugate(beA)))
print("   alpha_B =", alB, " beta_B =", beB,
      " D_B =", sp.simplify(alB*sp.conjugate(alB) - beB*sp.conjugate(beB)))
z = sp.symbols('z')
mA = sp.cancel(((alA*z - beA)/(sp.conjugate(beA)*z - sp.conjugate(alA))))
mB = sp.cancel(((alB*z - beB)/(sp.conjugate(beB)*z - sp.conjugate(alB))))
print("   M_A(z) =", mA, " (constant -1) ->", sp.simplify(mA + 1) == 0)
print("   M_B(z)*z =", sp.simplify(sp.expand(mB*z)), " (so M_B(z) = 1/z)")
print("   z3sq = M_A(z1^2) = -1 ; z4sq = M_B(z1^2) = 1/z1^2")
# z2 recovery, karlsson_k6_3.py:107-109, with z3sq = -1
num = sp.simplify(beB + sp.conjugate(alB))
den = sp.simplify(alB + sp.conjugate(beB))
print("   z2-recovery num =", num, " den =", den, " z2sq =", sp.simplify(num/den))
print("   den == 0 ?", den == 0, "(False => determinant does NOT vanish at Dita)")

# ---- (4) theta = 0 (F6 seam) -------------------------------------------
sub0 = {c: sp.Integer(1), s: sp.Integer(0), cp: sp.Symbol('q', real=True), sp_: sp.Symbol('q', real=True)}
bA0 = sp.simplify(reA11.subs(sub0) + I * imA11.subs(sub0))
aA0 = sp.simplify(reA12.subs(sub0) + I * imA12.subs(sub0))
alA0, beA0 = sp.simplify(aA0**2), sp.simplify(bA0**2)
alB0, beB0 = sp.simplify((-1 - aA0)**2), sp.simplify((-1 - bA0)**2)
print("== (4) F6 seam (theta=0) ==")
print("   A11 =", bA0, " A12 =", aA0, " alpha_A =", alA0, " beta_A =", beA0)
print("   alpha_B =", alB0, " beta_B =", beB0)
print("   D_A =", sp.simplify(alA0*sp.conjugate(alA0) - beA0*sp.conjugate(beA0)),
      " D_B =", sp.simplify(alB0*sp.conjugate(alB0) - beB0*sp.conjugate(beB0)))
mA0 = sp.cancel(((alA0*z - beA0)/(sp.conjugate(beA0)*z - sp.conjugate(alA0))))
mB0 = sp.cancel(((alB0*z - beB0)/(sp.conjugate(beB0)*z - sp.conjugate(alB0))))
print("   M_A(z) =", mA0, " M_B(z) =", mB0, " (both constant 1 ->",
      sp.simplify(mA0 - 1) == 0 and sp.simplify(mB0 - 1) == 0, ")")
print("   A Hermitian? A11 =", bA0, "conjugate =", sp.conjugate(bA0), "->", bA0 == sp.conjugate(bA0))

# ---- (5) Hermiticity locus and degeneracy factors ---------------------
print("== (5) strata ==")
print("   A Hermitian <=> A11 real <=> Im(A11) = r3/2*(cos th + sin th*cos ph) = 0")
print("   Im(A11) at Dita =", sp.simplify(cD + sD*h), "(nonzero)   at theta=0:", 1, "(nonzero)")
print("   M_A factor at Dita = sin(th)(sin ph - r3*cos th*cos ph) =",
      sp.simplify(sD*(h - r3*cD*h)), "= 0   <-- Dita ON the M_A stratum")
print("   M_B factor at Dita = sin(th)(sin ph + r3*cos th*cos ph) =",
      sp.simplify(sD*(h + r3*cD*h)), "  <-- Dita NOT on the M_B stratum")
print("   at theta=0: sin(th)=0 so BOTH factors vanish  <-- F6 is a stratum for BOTH")
print("   tan(ph) = r3*cos(th):  at th_D, r3*cos(th_D) =", sp.simplify(r3*cD), "= tan(pi/4)")
print("   D_B on the M_A-stratum curve = -2*r3*s*(2*r3*c*cp) =",
      sp.simplify(-2*r3*sD*(h + r3*cD*h)))
print("   alpha_B = 0 <=> A12 = -1 <=> sin(th)sin(ph) = 1/sqrt(3) and sin(th)cos(ph) = cos(th)")
print("   on the M_A curve this forces 3*cos(th)^2 = 1 -> cos(th) = +/-1/sqrt(3)")

tgt_A = 2 * r3 * s * (sp_ - r3 * c * cp)
tgt_B = -2 * r3 * s * (sp_ + r3 * c * cp)
print("== (2) D-invariants vs Karlsson degeneracy factors ==")
print("   D_A - 2*r3*s*(sin ph - r3*cos th*cos ph) =", T(sp.simplify(D_A - tgt_A)))
print("   D_B + 2*r3*s*(sin ph + r3*cos th*cos ph) =", T(sp.simplify(D_B - tgt_B)))
