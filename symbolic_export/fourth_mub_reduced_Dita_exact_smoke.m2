-- Exact-arithmetic probe: cyclotomic-6 coefficient field (expected to fail coefficient match for Dita)
needsPackage "Cyclotomic"
K = cyclotomicField 6
R = K[z1_1,z1_2,z1_3,z1_4,z1_5,w1_1,w1_2,w1_3,w1_4,w1_5]
z = K_0  -- primitive 6th root of unity
I = ideal(z1_1*w1_1 - 1, z1_2*w1_2 - 1, z1_3 - z)
print("cyclotomic6 smoke: dim I = " | toString dim I | newline)
