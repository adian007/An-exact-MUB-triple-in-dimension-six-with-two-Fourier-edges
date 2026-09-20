import sympy as sp

OUT = open("results/phase1_simplify_probe.txt", "w", encoding="utf-8")


def p(*a):
    print(*a)
    print(*a, file=OUT)


I = sp.I

e = sp.exp(-2 * sp.pi * I / 6)
p("raw exp  e^3 :", sp.simplify(e ** 3))
S = sum(e ** (1 * j) for j in range(6))
p("geom sum (k=1):", sp.simplify(S), "| expand:", sp.simplify(sp.expand(S)))

p("expand_complex exp(i pi/12):", sp.expand_complex(sp.exp(I * sp.pi / 12)))
p("cos(pi/12) =", sp.cos(sp.pi / 12), " sin(pi/12) =", sp.sin(sp.pi / 12))

w = sp.expand_complex(sp.exp(-2 * sp.pi * I / 6))
p("w6 radical:", w)
S2 = sum(w ** (1 * j) for j in range(6))
p("geom sum with radical (k=1):", sp.simplify(sp.expand(S2)))

z24 = sp.expand_complex(sp.exp(2 * sp.pi * I / 24))
p("z24 radical:", z24)

F6 = sp.Matrix(6, 6, lambda k, j: sp.expand(sp.expand_complex(sp.exp(-2 * sp.pi * I / 6)) ** (k * j)))
D = sp.simplify(sp.expand(F6 * F6.H - 6 * sp.eye(6)))
p("F6 F6^dag - 6I is zero matrix:", D.is_zero_matrix)
p("entries of D:", list(D))

# Try the alternative: exp-based with powsimp / trigsimp
F6e = sp.Matrix(6, 6, lambda k, j: sp.exp(-2 * sp.pi * I * k * j / 6))
De = F6e * F6e.H
p("exp-based F6*F6^dag row0:", [sp.simplify(sp.expand(De[0, c])) for c in range(6)])
p("exp-based trigsimp row0:", [sp.trigsimp(sp.expand_complex(De[0, c])) for c in range(6)])
OUT.close()