import sympy as sp
from pathlib import Path

I=sp.I
z,t=sp.symbols('z t')
x=sp.symbols('x1:6')
y=sp.symbols('y1:6')

H=sp.Matrix([
[1,1,1,1,1,1],
[1,-1,z,-z,I,-I],
[1,-I,I,I,-I,-1],
[1,I,-z,z,-1,-I],
[1,t,-I,-1,-t,I],
[1,-t,-1,-I,t,I],
])

v=sp.Matrix([1,*x])
vbar=sp.Matrix([1,*y])

# Unbiasedness: |v^* h_k|^2 = 6 in unnormalized coordinates.
# y_j are formal inverses of x_j; t=z^{-1}.
raw=[]
for k in range(6):
    a=(vbar.T*H[:,k])[0]
    Hbar=H.conjugate().xreplace({sp.conjugate(z):t, sp.conjugate(t):z})
    b=(v.T*Hbar[:,k])[0]
    b=sp.expand(b)
    # x variables are independent, so conjugate(x) should not occur
    raw.append(sp.expand(a*b-6))

# Remove duplicate/dependent equations by polynomial linear-algebra rank over Q(i,z,t) heuristically.
# For the object, retain all six equations and identify the first two parameter-independent ones.
rels=[z*t-1]
rels += [x[j]*y[j]-1 for j in range(5)]

out=Path('/mnt/data/mub_i3')
out.joinpath('I3_single_vector.txt').write_text('')
with open(out/'I3_single_vector.txt','w') as f:
    f.write('Exact single-vector MU ideal I3(z) for the Diţă circle\n')
    f.write('Variables: z,t,x1..x5,y1..y5; relations z*t-1 and xj*yj-1.\n')
    f.write('H_D(z) is the exact Laurent representative from the theorem work.\n\n')
    for j,r in enumerate(raw,1):
        f.write(f'F{j} = {sp.sstr(r)}\n\n')

# Save expanded equations as a SymPy-loadable module-ish repr.
with open(out/'I3_equations.py','w') as f:
    f.write('import sympy as sp\n')
    f.write('I=sp.I\n')
    f.write('z,t=sp.symbols("z t")\n')
    f.write('x=sp.symbols("x1:6"); y=sp.symbols("y1:6")\n')
    f.write('F=[\n')
    for r in raw: f.write('    '+repr(r)+',\n')
    f.write(']\n')

# Basic structural statistics
with open(out/'I3_stats.txt','w') as f:
    f.write('6 MU equations generated.\n')
    for j,r in enumerate(raw,1):
        f.write(f'F{j}: terms={len(sp.Poly(r, *([z,t]+list(x)+list(y))).terms())}, degree={sp.Poly(r, *([z,t]+list(x)+list(y))).total_degree()}\n')
    f.write('Parameter dependence:\n')
    for j,r in enumerate(raw,1):
        f.write(f'F{j}: z={r.has(z)}, t={r.has(t)}\n')

# Exact conjugation involution on variables: z<->t, x_j<->y_j, i->-i.
def conj_expr(e):
    ee=e.xreplace({sp.conjugate(z):t,sp.conjugate(t):z})
    return sp.conjugate(ee).xreplace({sp.conjugate(z):t,sp.conjugate(t):z, **{sp.conjugate(x[j]):y[j] for j in range(5)}, **{sp.conjugate(y[j]):x[j] for j in range(5)}})

# Instead construct involution by simultaneous substitution through dummy symbols.
def invol(e):
    return sp.expand(e.xreplace({z:t,t:z, **{x[j]:y[j] for j in range(5)}, **{y[j]:x[j] for j in range(5)}}).subs({I:-I}))
with open(out/'symmetry_test.txt','w') as f:
    f.write('Conjugation candidate: z<->t, x_j<->y_j, i->-i.\n')
    for j,r in enumerate(raw,1):
        q=invol(r)
        matches=[]
        for k,s2 in enumerate(raw,1):
            if sp.simplify(q-s2)==0: matches.append(k)
            elif sp.simplify(q+s2)==0: matches.append(-k)
        f.write(f'F{j} -> {matches}\n')


# Exact parameter involution z -> -z with x4 <-> x5, y4 <-> y5.
def shift_pi(e):
    sub={z:-z,t:-t,x[3]:x[4],x[4]:x[3],y[3]:y[4],y[4]:y[3]}
    return sp.expand(e.xreplace(sub))
with open(out/'symmetry_test.txt','a') as f:
    f.write('\nParameter symmetry z->-z, x4<->x5, y4<->y5:\n')
    for j,r in enumerate(raw,1):
        q=shift_pi(r); matches=[]
        for k,s2 in enumerate(raw,1):
            if sp.simplify(q-s2)==0: matches.append(k)
            elif sp.simplify(q+s2)==0: matches.append(-k)
        f.write(f'F{j} -> {matches}\n')
# Construct the exact projective third-basis fiber ideal template.
with open(out/'I3_branch_definition.md','w') as f:
    f.write('# Symbolic I3 branch ideal\n\n')
    f.write('For six vectors v_a=(1,x_{a1},...,x_{a5})/sqrt(6), introduce inverse variables y_{aj}.\n')
    f.write('The branch ideal consists of six copies of the single-vector MU ideal plus pairwise orthogonality equations.\n\n')
    f.write('For a<b, impose <v_a,v_b>=0 and its conjugate. In algebraic form:\n\n')
    f.write('`sum_j y[a,j]*x[b,j] + 1 = 0` and `sum_j x[a,j]*y[b,j] + 1 = 0`, with x*y=1.\n\n')
    f.write('A third-MUB branch is a six-vector fiber of this ideal whose projection to z is nonempty.\n')

print('generated', out)
