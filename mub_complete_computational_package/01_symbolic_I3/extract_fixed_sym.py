import numpy as np,itertools, sympy as sp
I=sp.I
z=sp.symbols('z', nonzero=True); t=1/z
H=sp.Matrix([[1,1,1,1,1,1],[1,-1,z,-z,I,-I],[1,-I,I,I,-I,-1],[1,I,-z,z,-1,-I],[1,t,-I,-1,-t,I],[1,-t,-1,-I,t,I]])
# numerical finder gives 12 pairs, but rederive with exact candidate using numeric matching
zn=np.exp(0.371j)
Hn=np.array(H.subs(z,zn)).astype(complex)
found=[]; tol=2e-8
for r in itertools.permutations(range(6)):
 Ar=Hn[list(r),:]
 for j0 in range(6):
  phases=1/Ar[:,j0]; C=phases[:,None]*Ar
  src=[C[:,j]/C[0,j] for j in range(6)]; tgt=[Hn[:,j]/Hn[0,j] for j in range(6)]
  used=set(); perm=[None]*6; ok=True
  for tj,tcol in enumerate(tgt):
   matches=[sj for sj,scol in enumerate(src) if sj not in used and np.max(np.abs(scol-tcol))<tol]
   if not matches: ok=False; break
   sj=matches[0]; used.add(sj); perm[tj]=sj
  if ok:
   # exact check: U=diag(phases) P_r; V diagonal phases after source->target perm
   U=sp.zeros(6)
   for ii,ri in enumerate(r): U[ii,ri]=1/H[ri,j0]
   # C=U H. target col tj = source col sj * v_tj, where v from first row
   V=sp.zeros(6)
   for tj,sj in enumerate(perm):
    v=H[0,tj]/(U*H[:,sj])[0]
    V[sj,sj]=v # source column gets v? We need U H V, target ordering via permutation.
   # Instead build Q with Q[sj,tj]=v_tj such that U H Q = H.
   Q=sp.zeros(6)
   for tj,sj in enumerate(perm):
    v=H[0,tj]/(U*H[:,sj])[0]
    Q[sj,tj]=v
   diff=sp.simplify(U*H*Q-H)
   if diff==sp.zeros(6):
    found.append((r,j0,tuple(perm),U,Q));
   break
print('exact fixed-z automorphisms:',len(found))
out='/mnt/data/mub_i3/fixed_symmetry_exact.txt'
with open(out,'w') as f:
 for n,(r,j0,perm,U,Q) in enumerate(found,1):
  f.write(f'G{n}: rows={r}, j0={j0}, target_to_source={perm}\n')
  f.write('U=\n'+sp.sstr(U)+'\n')
  f.write('Q=\n'+sp.sstr(Q)+'\n')
  # gauge action: v' = U v / (U v)_0
  f.write('gauge first component multiplier = '+sp.sstr(U[0,:])+'\n\n')
