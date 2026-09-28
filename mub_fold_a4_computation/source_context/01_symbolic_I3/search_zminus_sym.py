import numpy as np,itertools
I=1j

def H(z):
 return np.array([
[1,1,1,1,1,1],
[1,-1,z,-z,I,-I],
[1,-I,I,I,-I,-1],
[1,I,-z,z,-1,-I],
[1,1/z,-I,-1,-1/z,I],
[1,-1/z,-1,-I,1/z,I]],complex)

z=np.exp(0.371j); A=H(z); B=H(-z)
# Search monomial equivalence A -> B. For a row permutation r, set row phases from first column to keep first col 1; then normalized row vectors.
# Column phases/permutation are determined by matching column 0 first: choose a column j of A mapped to B column 0, but first-column constraint suggests j=0 under possible column permutation after row phases. Try all source columns j0.
tol=1e-8
for r in itertools.permutations(range(6)):
 Ar=A[list(r),:]
 # choose row phases so source column j0 becomes 1
 for j0 in range(6):
  phases=1/Ar[:,j0]
  C=phases[:,None]*Ar
  # target B has first column 1. C[:,j0]=1. Now map source j0 to target col0.
  # normalized target columns after row permutation must match up to a phase; since target col0=1, require source columns normalized by first row?
  src=[]
  for j in range(6):
   c=C[:,j]
   # normalize by first entry
   src.append(tuple(np.round(c/c[0],8)))
  tgt=[tuple(np.round(B[:,j]/B[0,j],8)) for j in range(6)]
  used=set(); perm=[None]*6; ok=True
  # target normalized columns need equal source normalized columns
  for tj,tcol in enumerate(tgt):
   matches=[sj for sj,scol in enumerate(src) if sj not in used and np.max(np.abs(np.array(scol)-np.array(tcol)))<tol]
   if not matches: ok=False; break
   sj=matches[0]; used.add(sj); perm[tj]=sj
  if ok:
   print('FOUND rows',r,'j0',j0,'perm target->source',perm)
   # derive U phases and V phases such that B = diag(target?) P A diag
   # Our C=diag(phases) P A. target col tj is C[:,sj]*v_tj; v= B[0,tj]/C[0,sj]
   v=[]
   for tj,sj in enumerate(perm): v.append(B[0,tj]/C[0,sj])
   print('row phases',phases)
   print('col phases target order',v)
   print('res',np.max(np.abs(C[:,perm]*np.array(v)[None,:]-B)))
   raise SystemExit
print('none')
