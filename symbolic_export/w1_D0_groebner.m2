-- Step 1 certificate run for w1_D0_exact.m2 (E3 follow-through, 2026-08-15).
-- Full Groebner basis of the n_wit=1 witness ideal over Q(zt_24, sqrt5)
-- (zt = zeta_24; [K:Q] = 16, verified in verify_field_degree_w1.jl: Phi_24
-- irreducible over Q and 5 not a square in Q(zeta_24), so toField's quotient
-- is genuinely a field and GB arithmetic is valid).
-- Expected from the bounded cost-estimate run: witness = (1), i.e. W1 empty
-- at H = D_bc with B3 = F_D, scoped to Karlsson lambda in {pi/2, 3pi/2} only.

load "w1_D0_exact.m2"
print("numgens witness = " | toString numgens witness);
print("numvars R       = " | toString numgens R);
print("field: Q(zeta_24, sqrt5), [K:Q] = 16");
elapsedTime G = gens gb witness;
print("Groebner basis = " | toString G);
print("GB size        = " | toString numColumns G);
r = 1_R % witness;
print("1 % witness    = " | toString r);
elapsedTime d = dim witness;
print("dim  = " | toString d);
if r == 0 then (
    print("CERTIFICATE: witness = (1) over Q(zeta_24, sqrt5).");
    print("W1 is empty for H = D_bc, B3 = F_D: no sixth flat vector unbiased to");
    print("both. Scope: Karlsson Dita-circle points lambda in {pi/2, 3pi/2} ONLY.");
) else (
    print("NO certificate: 1 not in witness ideal; dim/degree above apply.");
);
print("certificate run done");
