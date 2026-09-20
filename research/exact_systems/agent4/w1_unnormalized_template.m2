-- NON-EXECUTED TEMPLATE. Status NOT_RUN.
-- Replace K, H, Hbar, T, Tbar by exact data; H and T MUST be unnormalized.
K = QQ; -- placeholder only
R = K[z1,z2,z3,z4,z5,w1,w2,w3,w4,w5];
z = {1,z1,z2,z3,z4,z5}; w = {1,w1,w2,w3,w4,w5};
-- q(Cbar,C,k) := (sum(0..5,j -> Cbar_(j,k)*z#j))*(sum(0..5,j -> C_(j,k)*w#j))-6;
-- witness = ideal(apply(1..5,i -> z#i*w#i-1) |
--   apply(0..5,k -> q(Hbar,H,k)) | apply(0..5,k -> q(Tbar,T,k)));

