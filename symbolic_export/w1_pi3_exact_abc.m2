-- Exact pi/3 W1 export using the verified multiquadratic phase field.
-- The three degree-2 gcd factors and their independent square classes are
-- certified in results/phase2_w1_coefficient_field_preflight.log and
-- results/phase2_w1_quadratic_square_classes.log.
A = QQ[q,a,b,c] / ideal(
  q^12-q^6+1,
  323*a^2 + 16*q^11 - 35*q^10 + 170*q^9 - 18*q^8 + 120*q^7 + 32*q^6 - 166*q^5 - 157*q^4 + 8*q^3 - 46*q^2 - 48*q - 66,
  323*b^2 - 166*q^11 + 192*q^10 + 170*q^9 - 46*q^8 - 72*q^7 + 32*q^6 + 150*q^5 - 35*q^4 + 8*q^3 + 64*q^2 + 120*q - 66,
  323*c^2 + 150*q^11 - 157*q^10 + 170*q^9 + 64*q^8 - 48*q^7 + 32*q^6 + 16*q^5 + 192*q^4 + 8*q^3 - 18*q^2 - 72*q - 66
);
kk = toField A;
fieldMap = map(kk,A,{A_0,A_1,A_2,A_3});
q = fieldMap(A_0);
a = fieldMap(A_1);
b = fieldMap(A_2);
c = fieldMap(A_3);
ai = 1/a;
bi = 1/b;
ci = 1/c;
R = kk[z1,z2,z3,z4,z5,w1,w2,w3,w4,w5];
witness = ideal(
  w1*z1 - 1,
  w2*z2 - 1,
  w3*z3 - 1,
  w4*z4 - 1,
  w5*z5 - 1,
  w1*z1 + w1*z2 + w1*z3 + w1*z4 + w1*z5 + w1 + w2*z1 + w2*z2 + w2*z3 + w2*z4 + w2*z5 + w2 + w3*z1 + w3*z2 + w3*z3 + w3*z4 + w3*z5 + w3 + w4*z1 + w4*z2 + w4*z3 + w4*z4 + w4*z5 + w4 + w5*z1 + w5*z2 + w5*z3 + w5*z4 + w5*z5 + w5 + z1 + z2 + z3 + z4 + z5 - 5,
  q^1080*w4*z4 - q^1080*w4*z5 - q^1080*w5*z4 + q^1080*w5*z5 - q^1059*w2*z4 + q^1059*w2*z5 + q^1059*w3*z4 - q^1059*w3*z5 - q^1050*w1*z4 + q^1050*w1*z5 + q^1050*z4 - q^1050*z5 - q^345*w4*z2 + q^345*w4*z3 + q^345*w5*z2 - q^345*w5*z3 + q^324*w2*z2 - q^324*w2*z3 - q^324*w3*z2 + q^324*w3*z3 + q^315*w1*z2 - q^315*w1*z3 - q^315*z2 + q^315*z3 - q^30*w4*z1 + q^30*w4 + q^30*w5*z1 - q^30*w5 + q^9*w2*z1 - q^9*w2 - q^9*w3*z1 + q^9*w3 + w1*z1 - w1 - z1 - 5,
  q^324*w2*z2 - q^324*w2*z4 - q^324*w4*z2 + q^324*w4*z4 + q^321*w1*z2 - q^321*w1*z4 - q^321*w3*z2 + q^321*w3*z4 - q^315*w5*z2 + q^315*w5*z4 + q^315*z2 - q^315*z4 + q^219*w2*z1 - q^219*w2*z3 - q^219*w4*z1 + q^219*w4*z3 + q^216*w1*z1 - q^216*w1*z3 - q^216*w3*z1 + q^216*w3*z3 - q^210*w5*z1 + q^210*w5*z3 + q^210*z1 - q^210*z3 - q^9*w2*z5 + q^9*w2 + q^9*w4*z5 - q^9*w4 - q^6*w1*z5 + q^6*w1 + q^6*w3*z5 - q^6*w3 + w5*z5 - w5 - z5 - 5,
  q^324*w2*z2 - q^324*w2*z5 - q^324*w5*z2 + q^324*w5*z5 - q^321*w1*z2 + q^321*w1*z5 + q^321*w3*z2 - q^321*w3*z5 - q^315*w4*z2 + q^315*w4*z5 + q^315*z2 - q^315*z5 - q^219*w2*z1 + q^219*w2*z3 + q^219*w5*z1 - q^219*w5*z3 + q^216*w1*z1 - q^216*w1*z3 - q^216*w3*z1 + q^216*w3*z3 + q^210*w4*z1 - q^210*w4*z3 - q^210*z1 + q^210*z3 - q^9*w2*z4 + q^9*w2 + q^9*w5*z4 - q^9*w5 + q^6*w1*z4 - q^6*w1 - q^6*w3*z4 + q^6*w3 + w4*z4 - w4 - z4 - 5,
  q^1080*w4*z4 - q^1080*w4*z5 - q^1080*w5*z4 + q^1080*w5*z5 - q^1059*w1*z4 + q^1059*w1*z5 + q^1059*w2*z4 - q^1059*w2*z5 + q^1050*w3*z4 - q^1050*w3*z5 - q^1050*z4 + q^1050*z5 - q^345*w4*z1 + q^345*w4*z2 + q^345*w5*z1 - q^345*w5*z2 + q^324*w1*z1 - q^324*w1*z2 - q^324*w2*z1 + q^324*w2*z2 - q^315*w3*z1 + q^315*w3*z2 + q^315*z1 - q^315*z2 + q^30*w4*z3 - q^30*w4 - q^30*w5*z3 + q^30*w5 - q^9*w1*z3 + q^9*w1 + q^9*w2*z3 - q^9*w2 + w3*z3 - w3 - z3 - 5,
  q^324*w1*z1 + q^324*w1*z3 - q^324*w1*z4 - q^324*w1*z5 + q^324*w3*z1 + q^324*w3*z3 - q^324*w3*z4 - q^324*w3*z5 - q^324*w4*z1 - q^324*w4*z3 + q^324*w4*z4 + q^324*w4*z5 - q^324*w5*z1 - q^324*w5*z3 + q^324*w5*z4 + q^324*w5*z5 + q^315*w2*z1 + q^315*w2*z3 - q^315*w2*z4 - q^315*w2*z5 - q^315*z1 - q^315*z3 + q^315*z4 + q^315*z5 + q^9*w1*z2 - q^9*w1 + q^9*w3*z2 - q^9*w3 - q^9*w4*z2 + q^9*w4 - q^9*w5*z2 + q^9*w5 + w2*z2 - w2 - z2 - 5,
  a*ai*q^720*w4*z4 - a*ai*q^701*w2*z4 + a*ai*q^700*w1*z4 - a*ai*q^55*w4*z2 + a*ai*q^36*w2*z2 - a*ai*q^35*w1*z2 + a*ai*q^20*w4*z1 - a*ai*q*w2*z1 + a*ai*w1*z1 + a*q^1000*w4*z5 - a*q^981*w2*z5 + a*q^980*w1*z5 + a*q^930*w4*z3 - a*q^911*w2*z3 + a*q^910*w1*z3 + a*q^20*w4 - a*q*w2 + a*w1 + ai*q^728*w5*z4 + ai*q^726*w3*z4 + ai*q^700*z4 - ai*q^63*w5*z2 - ai*q^61*w3*z2 - ai*q^35*z2 + ai*q^28*w5*z1 + ai*q^26*w3*z1 + ai*z1 + q^1008*w5*z5 + q^1006*w3*z5 + q^980*z5 + q^938*w5*z3 + q^936*w3*z3 + q^910*z3 + q^28*w5 + q^26*w3 - 5,
  a*ai*q^72*w4*z4 + a*ai*q^71*w2*z4 - a*ai*q^70*w1*z4 + a*ai*q^37*w4*z2 + a*ai*q^36*w2*z2 - a*ai*q^35*w1*z2 - a*ai*q^2*w4*z1 - a*ai*q*w2*z1 + a*ai*w1*z1 + a*q^982*w4*z5 + a*q^981*w2*z5 - a*q^980*w1*z5 + a*q^912*w4*z3 + a*q^911*w2*z3 - a*q^910*w1*z3 + a*q^2*w4 + a*q*w2 - a*w1 + ai*q^98*w5*z4 + ai*q^96*w3*z4 + ai*q^70*z4 + ai*q^63*w5*z2 + ai*q^61*w3*z2 + ai*q^35*z2 - ai*q^28*w5*z1 - ai*q^26*w3*z1 - ai*z1 + q^1008*w5*z5 + q^1006*w3*z5 + q^980*z5 + q^938*w5*z3 + q^936*w3*z3 + q^910*z3 + q^28*w5 + q^26*w3 - 5,
  b*bi*q^1152*w4*z4 + b*bi*q^1127*w2*z4 + b*bi*q^1120*w1*z4 + b*bi*q^277*w4*z2 + b*bi*q^252*w2*z2 + b*bi*q^245*w1*z2 + b*bi*q^32*w4*z1 + b*bi*q^7*w2*z1 + b*bi*w1*z1 + b*q^592*w4*z5 + b*q^567*w2*z5 + b*q^560*w1*z5 + b*q^102*w4*z3 + b*q^77*w2*z3 + b*q^70*w1*z3 + b*q^32*w4 + b*q^7*w2 + b*w1 + bi*q^1136*w5*z4 + bi*q^1122*w3*z4 + bi*q^1120*z4 + bi*q^261*w5*z2 + bi*q^247*w3*z2 + bi*q^245*z2 + bi*q^16*w5*z1 + bi*q^2*w3*z1 + bi*z1 + q^576*w5*z5 + q^562*w3*z5 + q^560*z5 + q^86*w5*z3 + q^72*w3*z3 + q^70*z3 + q^16*w5 + q^2*w3 - 5,
  b*bi*q^1152*w4*z4 + b*bi*q^1127*w2*z4 + b*bi*q^1120*w1*z4 + b*bi*q^277*w4*z2 + b*bi*q^252*w2*z2 + b*bi*q^245*w1*z2 + b*bi*q^32*w4*z1 + b*bi*q^7*w2*z1 + b*bi*w1*z1 - b*q^592*w4*z5 - b*q^567*w2*z5 - b*q^560*w1*z5 - b*q^102*w4*z3 - b*q^77*w2*z3 - b*q^70*w1*z3 - b*q^32*w4 - b*q^7*w2 - b*w1 - bi*q^1136*w5*z4 - bi*q^1122*w3*z4 - bi*q^1120*z4 - bi*q^261*w5*z2 - bi*q^247*w3*z2 - bi*q^245*z2 - bi*q^16*w5*z1 - bi*q^2*w3*z1 - bi*z1 + q^576*w5*z5 + q^562*w3*z5 + q^560*z5 + q^86*w5*z3 + q^72*w3*z3 + q^70*z3 + q^16*w5 + q^2*w3 - 5,
  c*ci*q^1116*w2*z2 + c*ci*q^1093*w4*z2 + c*ci*q^1085*w1*z2 + c*ci*q^311*w2*z4 + c*ci*q^288*w4*z4 + c*ci*q^280*w1*z4 + c*ci*q^31*w2*z1 + c*ci*q^8*w4*z1 + c*ci*w1*z1 + c*q^521*w2*z3 + c*q^498*w4*z3 + c*q^490*w1*z3 + c*q^171*w2*z5 + c*q^148*w4*z5 + c*q^140*w1*z5 + c*q^31*w2 + c*q^8*w4 + c*w1 + ci*q^1099*w3*z2 + ci*q^1089*w5*z2 + ci*q^1085*z2 + ci*q^294*w3*z4 + ci*q^284*w5*z4 + ci*q^280*z4 + ci*q^14*w3*z1 + ci*q^4*w5*z1 + ci*z1 + q^504*w3*z3 + q^494*w5*z3 + q^490*z3 + q^154*w3*z5 + q^144*w5*z5 + q^140*z5 + q^14*w3 + q^4*w5 - 5,
  c*ci*q^1116*w2*z2 + c*ci*q^1093*w4*z2 + c*ci*q^1085*w1*z2 + c*ci*q^311*w2*z4 + c*ci*q^288*w4*z4 + c*ci*q^280*w1*z4 + c*ci*q^31*w2*z1 + c*ci*q^8*w4*z1 + c*ci*w1*z1 - c*q^521*w2*z3 - c*q^498*w4*z3 - c*q^490*w1*z3 - c*q^171*w2*z5 - c*q^148*w4*z5 - c*q^140*w1*z5 - c*q^31*w2 - c*q^8*w4 - c*w1 - ci*q^1099*w3*z2 - ci*q^1089*w5*z2 - ci*q^1085*z2 - ci*q^294*w3*z4 - ci*q^284*w5*z4 - ci*q^280*z4 - ci*q^14*w3*z1 - ci*q^4*w5*z1 - ci*z1 + q^504*w3*z3 + q^494*w5*z3 + q^490*z3 + q^154*w3*z5 + q^144*w5*z5 + q^140*z5 + q^14*w3 + q^4*w5 - 5
);
print("numgens witness = " | toString numgens witness);
elapsedTime G = gens gb witness;
print("Groebner basis element count = " | toString numgens source G);
print("1 % witness = " | toString (1_R % witness));
print("dim = " | toString dim witness);
