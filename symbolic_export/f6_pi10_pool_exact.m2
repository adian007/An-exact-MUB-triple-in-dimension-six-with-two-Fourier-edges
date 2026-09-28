-- Exact F6(theta=0, lambda=pi/10) per-vector MU ideal
-- Coefficient field Q(zeta_60), no floating-point coefficients
-- One MU equation is omitted by the unitary sum identity
kk = toField(QQ[zz]/ideal(zz^16 + zz^14 - zz^10 - zz^8 - zz^6 + zz^2 + 1));
R = kk[x1,x2,x3,x4,x5,y1,y2,y3,y4,y5];
I = ideal(
  x1*y1-1,
  x2*y2-1,
  x3*y3-1,
  x4*y4-1,
  x5*y5-1,
  ((1)+(1)*x1+(1)*x2+(1)*x3+(1)*x4+(1)*x5)*((1)+(1)*y1+(1)*y2+(1)*y3+(1)*y4+(1)*y5)-6,
  ((1)+(-1)*x1+(1)*x2+(-1)*x3+(1)*x4+(-1)*x5)*((1)+(-1)*y1+(1)*y2+(-1)*y3+(1)*y4+(-1)*y5)-6,
  ((1)+((1)*zz^1 + (1)*zz^3 + (-1)*zz^9 + (-1)*zz^11 + (1)*zz^15)*x1+((-1)*zz^10)*x2+((-1)*zz^1 + (-1)*zz^3 + (1)*zz^7 + (1)*zz^9 + (1)*zz^11 + (-1)*zz^15)*x3+(-1 + (1)*zz^10)*x4+((-1)*zz^7)*x5)*((1)+((1)*zz^3)*y1+(-1 + (1)*zz^10)*y2+((-1)*zz^13)*y3+((-1)*zz^10)*y4+((-1)*zz^3 + (1)*zz^13)*y5)-6,
  ((1)+((-1)*zz^1 + (-1)*zz^3 + (1)*zz^9 + (1)*zz^11 + (-1)*zz^15)*x1+((-1)*zz^10)*x2+((1)*zz^1 + (1)*zz^3 + (-1)*zz^7 + (-1)*zz^9 + (-1)*zz^11 + (1)*zz^15)*x3+(-1 + (1)*zz^10)*x4+((1)*zz^7)*x5)*((1)+((-1)*zz^3)*y1+(-1 + (1)*zz^10)*y2+((1)*zz^13)*y3+((-1)*zz^10)*y4+((1)*zz^3 + (-1)*zz^13)*y5)-6,
  ((1)+((-1)*zz^10)*x1+(-1 + (1)*zz^10)*x2+(-1 + (1)*zz^10)*x3+((-1)*zz^10)*x4+(1)*x5)*((1)+(-1 + (1)*zz^10)*y1+((-1)*zz^10)*y2+((-1)*zz^10)*y3+(-1 + (1)*zz^10)*y4+(1)*y5)-6
);
print("Starting exact Groebner calculation" | newline);
G = gb I;
print("Exact MU-pool ideal computed; dimension=" | toString dim I | newline);
