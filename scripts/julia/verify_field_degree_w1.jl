# Step 0.1 of the w1_D0_exact.m2 elimination attempt:
# verify [Q(zeta_24, sqrt(5)) : Q] = 16, i.e. that s5^2 - 5 stays irreducible
# over the cyclotomic field Q(zeta_24) (Phi_24 = x^8 - x^4 + 1, phi(24) = 8).
# Output: results/track_c_elimination/w1_D0_field_degree.txt
using Nemo

const OUT = joinpath(normpath(joinpath(@__DIR__, "..", "..")),
                     "results", "track_c_elimination", "w1_D0_field_degree.txt")
mkpath(dirname(OUT))
const LOG = String[]
function logln(args...)
    s = join(string.(args))
    push!(LOG, s)
    println(s)
end

Qx, x = polynomial_ring(QQ, "x")
phi24 = x^8 - x^4 + 1
logln("Phi_24 = ", phi24)
logln("Phi_24 irreducible over Q: ", is_irreducible(phi24))

K, zeta = number_field(phi24, "zeta")
logln("[Q(zeta_24):Q] = ", degree(K))

# t^2 - 5 is irreducible over K iff 5 is not a square in K.
# Residue-field certificate: take p = 73 (p ≡ 1 mod 24, so Phi_24 splits
# completely mod p and every residue field of K above p is F_p itself).
# If 5 were a square in K, its image would be a square in F_73.
p = 73
@assert p % 24 == 1
F = Nemo.Native.GF(p)
Fy, y = polynomial_ring(F, "y")
phi24_p = y^8 - y^4 + 1
fac = factor(phi24_p)
degs = sort([degree(g) for (g, e) in fac])
logln("Phi_24 factor degrees mod $p: ", degs, "  (all 1 => p splits completely)")
@assert all(==(1), degs)
sq5_modp = is_square(F(5))
logln("5 is a square mod $p: ", sq5_modp)
@assert !sq5_modp
logln("=> 5 is NOT a square in Q(zeta_24); t^2 - 5 is irreducible over it.")
logln("[Q(zeta_24, sqrt5):Q] = 16")

open(OUT, "w") do io
    for s in LOG
        println(io, s)
    end
end
println("Wrote $OUT")
