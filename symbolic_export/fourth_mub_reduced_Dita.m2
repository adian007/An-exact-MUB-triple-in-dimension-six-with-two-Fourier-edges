-- Reduced fourth-MUB witness (2 vectors) for Dita
-- Regenerate via: julia --project=. scripts/julia/symbolic_elimination.jl
-- Expected: ~20 vars, ~52 eqs (2 witness vectors x MU to I, H, B3 + orthogonality)
-- Status: export skeleton; run symbolic_elimination.jl for full equations

R = CC[z1_1,z1_2,z1_3,z1_4,z1_5,w1_1,w1_2,w1_3,w1_4,w1_5,z2_1,z2_2,z2_3,z2_4,z2_5,w2_1,w2_2,w2_3,w2_4,w2_5];
witness = ideal(0);
print "#generators = "; print numgens witness; print " dim = "; print dim witness; println ""
