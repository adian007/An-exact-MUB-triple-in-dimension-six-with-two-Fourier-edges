-- Bounded cost-estimate driver for w1_D0_exact.m2 (E3 follow-through, 2026-08-15).
-- Stage A: specialize zeta_24 -> a, sqrt5 -> b in ZZ/241 (both exist since
--          241 = 1 mod 24 and 5 is a QR mod 241) and compute dim/degree fast.
--          This is a signal/cost probe only, NOT a characteristic-0 certificate.
-- Stage B: dim/degree over the exact field Q(zeta_24, sqrt5) (degree 16 over Q).
--          If the witness ideal is the unit ideal here, W1 is empty at D_bc,
--          i.e. lambda in {pi/2, 3pi/2} on the Karlsson Dita circle.
-- Run bounded (external timeout); do not silently escalate.

load "w1_D0_exact.m2"
print("numgens witness = " | toString numgens witness);
print("numvars R       = " | toString numgens R);
print("coefficient field: Q(zeta_24, sqrt5), [K:Q] = 16 (verified separately)");

-- ---------- Stage A: finite-field specialization ----------
p = 241;
kp = ZZ/p;
a = null; b = null;
for i from 2 to p-1 do (c := sub(i, kp); if c^24 == 1 and c^12 != 1 and c^8 != 1 then (a = c; break));
for i from 2 to p-1 do (c := sub(i, kp); if c^2 == 5 then (b = c; break));
assert(a^8 - a^4 + 1 == 0);
assert(b^2 == 5);
print("Stage A: p = 241, zeta -> " | toString a | ", s5 -> " | toString b);
Rp = kp[z1p, z2p, z3p, z4p, z5p, w1p, w2p, w3p, w4p, w5p];
phi = map(Rp, R, (gens Rp) | {a, b});  -- images for z1..w5, then zt, s5
Ip = phi witness;
print("Stage A: computing dim over ZZ/241 ...");
elapsedTime dp = dim Ip;
print("Stage A: dim (mod 241) = " | toString dp);
if dp >= 0 then (
    elapsedTime degp = degree Ip;
    print("Stage A: degree (mod 241) = " | toString degp);
) else (
    print("Stage A: unit ideal mod 241 (signal only, not a char-0 certificate)");
);

-- ---------- Stage B: exact number field ----------
print("Stage B: computing dim over Q(zeta_24, sqrt5) ...");
gbTrace = 1;
elapsedTime d = dim witness;
print("Stage B: dim (exact) = " | toString d);
if d >= 0 then (
    elapsedTime dg = degree witness;
    print("Stage B: degree (exact) = " | toString dg);
) else (
    print("Stage B: witness is the UNIT IDEAL over Q(zeta_24, sqrt5).");
    print("=> W1 empty at D_bc, lambda in {pi/2, 3pi/2} (exact-coefficient certificate).");
);
print("driver done");
