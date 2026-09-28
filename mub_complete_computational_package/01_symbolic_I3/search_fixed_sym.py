import numpy as np,itertools

def H(z):
 return np.array([[1,1,1,1,1,1],[1,-1,z,-z,1j,-1j],[1,-1j,1j,1j,-1j,-1],[1,1j,-z,z,-1,-1j],[1,1/z,-1j,-1,-1/z,1j],[1,-1/z,-1,-1j,1/z,1j]],complex)
z=np.exp(0.371j); A=H(z); B=A
tol=2e-8
found=[]
for r in itertools.permutations(range(6)):
 Ar=A[list(r),:]
 for j0 in range(6):
  phases=1/Ar[:,j0]
  C=phases[:,None]*Ar
  src=[C[:,j]/C[0,j] for j in range(6)]
  tgt=[B[:,j]/B[0,j] for j in range(6)]
  used=set(); perm=[None]*6; ok=True
  for tj,tcol in enumerate(tgt):
   matches=[sj for sj,scol in enumerate(src) if sj not in used and np.max(np.abs(scol-tcol))<tol]
   if not matches: ok=False; break
   sj=matches[0]; used.add(sj); perm[tj]=sj
  if ok:
   found.append((r,tuple(perm)))
   break
print('row automorphisms',len(found))
for a in found: print(a)
