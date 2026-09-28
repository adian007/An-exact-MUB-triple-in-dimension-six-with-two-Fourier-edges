"""
d_invariant_check.py

PURPOSE: resolve the D-invariant arc/circle question (task: "Close or demote the
D-invariant arc/circle claim").

REPRODUCTION / re-derivation using the repo's OWN family definition
(`scripts/python/karlsson_k6_3.py`, following McNulty & Weigert Eqs. (7.2)-(7.6),
citing Karlsson arXiv:1003.4177). It runs NO new homotopy solve and produces NO
new clique search. It only evaluates the family and its Mobius data.

Answers:
  Q1  Is D = |alpha|^2 - |beta|^2 independent of lambda?
  Q2  What is the exact z-structure at each locus (z2,z3,z4 as functions of z1)?
  Q3  Which entries of H(lambda) actually depend on lambda at each locus?

Output: results/d_invariant_check_out.txt
"""
import io
import numpy as np
import sympy as sp

OUT = io.open("results/d_invariant_check_out.txt", "w", encoding="utf-8")


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    OUT.write(s + "\n")
    OUT.flush()


th, ph = sp.symbols("theta phi", real=True)
c, s = sp.cos(th), sp.sin(th)
p, q = sp.cos(ph), sp.sin(ph)
I = sp.I

A11 = sp.Rational(-1, 2) + I * sp.sqrt(3) / 2 * (c + sp.exp(-I * ph) * s)
A12 = sp.Rational(-1, 2) + I * sp.sqrt(3) / 2 * (-c + sp.exp(I * ph) * s)
B11 = -1 - A11
B12 = -1 - A12

aA, bA = A12 ** 2, A11 ** 2
aB, bB = B12 ** 2, B11 ** 2

P("=" * 78)
P("Q1  Is D = |alpha|^2 - |beta|^2 independent of lambda?")
P("=" * 78)
P("Structural: A = build_A(theta,phi) and B = -F2 - A contain no lambda.")
P("alpha_x = X12^2, beta_x = X11^2 are therefore functions of (theta,phi) only.")
P("lambda enters ONLY as z1 = exp(i*lambda), used afterwards in M_A, M_B.")
P("Hence dD/dlambda == 0 identically.")
P("")
P("alpha_A free symbols:", sorted(str(x) for x in aA.free_symbols))
P("beta_A  free symbols:", sorted(str(x) for x in bA.free_symbols))
P("alpha_B free symbols:", sorted(str(x) for x in aB.free_symbols))
P("beta_B  free symbols:", sorted(str(x) for x in bB.free_symbols))
P("  -> 'lambda' appears in NONE of them.   [exact symbolic]")
P("")

m11 = sp.simplify(sp.expand(sp.re(A11) ** 2 + sp.im(A11) ** 2))
m12 = sp.simplify(sp.expand(sp.re(A12) ** 2 + sp.im(A12) ** 2))
mB11 = sp.simplify(sp.expand(sp.re(B11) ** 2 + sp.im(B11) ** 2))
mB12 = sp.simplify(sp.expand(sp.re(B12) ** 2 + sp.im(B12) ** 2))
P("Closed forms of the moduli:")
P("  |A11|^2 =", m11)
P("  |A12|^2 =", m12)
P("  |B11|^2 =", mB11)
P("  |B12|^2 =", mB12)

fA = sp.simplify(m12 - m11)
fB = sp.simplify(mB12 - mB11)
P("")
P("Degeneracy discriminants (vanishing <=> the Mobius map collapses):")
P("  |A12|^2 - |A11|^2 =", sp.factor(fA))
P("  |B12|^2 - |B11|^2 =", sp.factor(fB))
P("")
targetA = sp.sqrt(3) * s * (q - sp.sqrt(3) * c * p)
targetB = sp.sqrt(3) * s * (q + sp.sqrt(3) * c * p)
P("Karlsson Sec.5 factors:  sin(th)(sin(ph) - sqrt3 cos(th)cos(ph))")
P("                         sin(th)(sin(ph) + sqrt3 cos(th)cos(ph))")
P("  residual vs |A12|^2-|A11|^2 :", sp.simplify(fA - targetA))
P("  residual vs |B12|^2-|B11|^2 :", sp.simplify(fB + targetB))

thD = sp.acos(1 / sp.sqrt(3))
P("")
P("At the anchors:")
P("  F6   (th=0)          :", sp.simplify(fA.subs(th, 0)), ",",
  sp.simplify(fB.subs(th, 0)))
P("  Dita (thD, phi=pi/4) :",
  sp.simplify(fA.subs({th: thD, ph: sp.pi / 4})), ",",
  sp.simplify(fB.subs({th: thD, ph: sp.pi / 4})))
P("")
P("TEST (i): zero/pole/sign-change at the arc boundaries 0.055044803035 or")
P("          0.0625984193874 (lambda offsets)?")
P("  Answer: NO, and it cannot. D has no lambda dependence at all, so it has")
P("  no structure at any lambda, boundary or interior.   [exact symbolic]")
P("")

# ---------------------------------------------------------------- z-structure
P("=" * 78)
P("Q2  Exact z-structure at each locus")
P("=" * 78)
z1 = sp.symbols("z1", positive=True)


def mobius(z, al, be):
    return (al * z - be) / (sp.conjugate(be) * z - sp.conjugate(al))


def zstruct(thv, phv, label):
    Av11 = sp.simplify(A11.subs({th: thv, ph: phv}))
    Av12 = sp.simplify(A12.subs({th: thv, ph: phv}))
    Bv11 = sp.simplify(B11.subs({th: thv, ph: phv}))
    Bv12 = sp.simplify(B12.subs({th: thv, ph: phv}))
    aAv, bAv = sp.simplify(Av12 ** 2), sp.simplify(Av11 ** 2)
    aBv, bBv = sp.simplify(Bv12 ** 2), sp.simplify(Bv11 ** 2)
    P("")
    P("--- %s ---" % label)
    P("  A11=%s  A12=%s" % (Av11, Av12))
    P("  B11=%s  B12=%s" % (Bv11, Bv12))
    P("  alpha_A=%s beta_A=%s   (|a|=%s |b|=%s)" %
      (aAv, bAv, sp.Abs(aAv), sp.Abs(bAv)))
    P("  alpha_B=%s beta_B=%s   (|a|=%s |b|=%s)" %
      (aBv, bBv, sp.Abs(aBv), sp.Abs(bBv)))
    z3sq = sp.simplify(mobius(z1 ** 2, aAv, bAv))
    z4sq = sp.simplify(mobius(z1 ** 2, aBv, bBv))
    P("  z3^2 = M_A(z1^2) = %s" % z3sq)
    P("  z4^2 = M_B(z1^2) = %s" % z4sq)
    num = sp.simplify(bBv - z3sq * sp.conjugate(aBv))
    den = sp.simplify(aBv - z3sq * sp.conjugate(bBv))
    P("  z2^2 = num/den ;  num=%s  den=%s" % (num, den))
    if den == 0:
        P("     ** DENOMINATOR IDENTICALLY ZERO -> 0/0 form; z2 needs a LIMIT. **")
    else:
        P("     den != 0  ->  z2^2 = %s" % sp.simplify(num / den))
    return z3sq, z4sq


zstruct(0, sp.Rational(1, 2), "F6_theta0  (theta=0, phi=0.5)")
zstruct(thD, sp.pi / 4, "Dita       (theta=arccos(1/sqrt3), phi=pi/4)")

P("")
P("-" * 78)
P("F6 z2: the inversion is a 0/0 form at th=0. Resolve it NUMERICALLY by")
P("approaching th -> 0 along phi = 0.5 and watching z2^2.")
for e in (1e-1, 1e-2, 1e-3, 1e-4, 1e-6):
    Av11 = complex(A11.subs({th: e, ph: sp.Rational(1, 2)}).evalf(30))
    Av12 = complex(A12.subs({th: e, ph: sp.Rational(1, 2)}).evalf(30))
    Bv11, Bv12 = -1 - Av11, -1 - Av12
    aAv, bAv = Av12 ** 2, Av11 ** 2
    aBv, bBv = Bv12 ** 2, Bv11 ** 2
    z3sq_n = sp.N(mobius(z1 ** 2, sp.nsimplify(aAv, rational=False),
                          sp.nsimplify(bAv, rational=False)), 20)
    # numeric: substitute a sample z1 to see the map's value (should be ~const)
    z3v = complex(sp.N(mobius(1.0, sp.nsimplify(aAv, rational=False),
                              sp.nsimplify(bAv, rational=False)), 20))
    z3sq_ev = complex(z3v)
    num_n = bBv - z3sq_ev * (aBv.conjugate())
    den_n = aBv - z3sq_ev * (bBv.conjugate())
    z2sq_n = num_n / den_n if abs(den_n) > 0 else float("nan")
    P("  eps=%-8.0e  |A11|=%.12f |B11|=%.12f  M_A(z1^2)=%s  z2^2=%.10f%+.3fi"
      % (e, abs(Av11), abs(Bv11), ("%.10f%+.3fi" % (z3sq_ev.real, z3sq_ev.imag)),
         z2sq_n.real, z2sq_n.imag))
P("  (M_A, M_B collapse to a constant as th->0; the z2 inversion degenerates.)")
P("")

P("=" * 78)
P("Q3  Which entries of H(lambda) depend on lambda?")
P("=" * 78)
P("From the assembly (karlsson_k6_3.py:128-137) H is built from")
P("  Z1=Zleft(z1), Z2=Zleft(z2), Z3=Zright(z3), Z4=Zright(z4),")
P("with  Zleft(z)=[[1,1],[z,-z]]  and  Zright(z)=[[1,z],[1,-z]].")
P("  F6  : z3=z4=1 constant -> lambda enters ONLY through Z1 (affine in z1).")
P("  Dita: z2=z3=i constant, z4=1/z1 -> lambda enters through Z1 (z1) AND")
P("        Z4 (1/z1): H is Laurent in z1 with terms paired as z1 <-> 1/z1.")
P("")

# ------------------------------------------------------- numeric confirmation
import importlib.util
spec = importlib.util.spec_from_file_location("k6", "scripts/python/karlsson_k6_3.py")
k6 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(k6)

P("-" * 78)
P("NUMERIC ENTRY-LEVEL CHECK (uses the repo builder directly)")
P("  NOTE: comparing |H| is useless (all CHM entries have modulus 1), so we")
P("  compare the COMPLEX entries against the lam=0 reference.")
P("-" * 78)
lams = [0.0, 0.3, 1.0471975512, 2.0, 3.14159265358979, 4.0, 5.235987756, 6.2]
for label, thv, phv in (("F6_theta0", 0.0, 0.5),
                        ("F6_th=1e-9", 1e-9, 0.5),
                        ("Dita", np.arccos(1 / np.sqrt(3)), np.pi / 4)):
    Hs = []
    ok = True
    for lam in lams:
        try:
            H = k6.build_k6(thv, phv, lam)
            if not np.all(np.isfinite(H)):
                ok = False
            Hs.append(H)
        except Exception as e:
            P("  %s lam=%.4f EXCEPTION %s" % (label, lam, e))
            ok = False
    if ok and Hs:
        ref = Hs[0]
        varying = int(np.sum(np.array([np.max(np.abs(H - ref)) for H in Hs]) > 0))
        # count entries that change for at least one lambda
        per_entry = np.zeros((6, 6), dtype=int)
        for H in Hs:
            per_entry += (np.abs(H - ref) > 1e-12).astype(int)
        nchg = int(np.sum(per_entry > 0))
        P("  %-11s built OK at all %d lambdas; ENTRIES changing with lambda: %d / 36"
          % (label, len(lams), nchg))
        P("               Hadamardness max|H H^dag - 6I| = %.3e"
          % np.max(np.abs(ref @ ref.conj().T - 6 * np.eye(6))))
    else:
        P("  %-11s NOT finite -> repo builder DEGENERATE at this (theta,phi)."
          % label)
        P("               (0/0 in the z2 inversion; see the stderr warnings above)")


P("")
P("=" * 78)
P("SUMMARY")
P("=" * 78)
P(" Q1  D is a function of (theta,phi) ONLY. dD/dlambda == 0 identically.")
P("     [exact symbolic]")
P(" Q2  Both loci are lambda-loops at FIXED (theta,phi):")
P("       F6   at (theta,phi) = (0, 0.5)")
P("       Dita at (theta,phi) = (arccos(1/sqrt3), pi/4)")
P("     Therefore D is CONSTANT along each locus.")
P("     => A point INSIDE the F6 arc and a point OUTSIDE it (same locus) have")
P("        IDENTICAL D. D has zero discriminating power for survival.")
P(" Q3  lambda enters only via the z's; the z-structures differ between loci.")
P("     That difference is a genuine, exact structural asymmetry (hypothesis C),")
P("     and it is NOT a variation of D.")
P("")
P(" Evidence tags: Q1/Q2 exact symbolic (sympy). Q3 entry-level counts are")
P(" sampled/numerical only (float64 evaluation of the repo builder).")
OUT.close()


