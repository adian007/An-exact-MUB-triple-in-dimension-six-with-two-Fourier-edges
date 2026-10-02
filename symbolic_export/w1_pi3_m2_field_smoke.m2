A = QQ[q,a,b,c] / ideal(
  q^12-q^6+1,
  323*a^2 + 16*q^11 - 35*q^10 + 170*q^9 - 18*q^8 + 120*q^7 + 32*q^6 - 166*q^5 - 157*q^4 + 8*q^3 - 46*q^2 - 48*q - 66,
  323*b^2 - 166*q^11 + 192*q^10 + 170*q^9 - 46*q^8 - 72*q^7 + 32*q^6 + 150*q^5 - 35*q^4 + 8*q^3 + 64*q^2 + 120*q - 66,
  323*c^2 + 150*q^11 - 157*q^10 + 170*q^9 + 64*q^8 - 48*q^7 + 32*q^6 + 16*q^5 + 192*q^4 + 8*q^3 - 18*q^2 - 72*q - 66
);
kk = toField A;
fieldMap = map(kk,A,{A_0,A_1,A_2,A_3});
qK = fieldMap(A_0);
aK = fieldMap(A_1);
bK = fieldMap(A_2);
cK = fieldMap(A_3);
print("dimension of coefficient algebra = " | toString dim A);
print("degree of coefficient algebra = " | toString degree A);
print("q relation residual = " | toString (qK^12-qK^6+1));
print("phase inverse check = " | toString ((1/aK)*aK-1) | ", " | toString ((1/bK)*bK-1) | ", " | toString ((1/cK)*cK-1));
