import numpy as np,itertools,math,os,time
base='/mnt/data/mub_b3_work/b3_a4_i4_results'
perms=list(itertools.permutations(range(6)))
F6=np.array([[np.exp(2j*np.pi*r*c/6) for c in range(6)] for r in range(6)],complex)
var={(r,c):('a' if c in (1,4) else 'b') for r in (1,3,5) for c in (1,2,4,5)}
fixed=np.ones((6,6),bool)
for r,c in var: fixed[r,c]=False

def test(M,tol_fixed=1e-7):
 best=(1e9,None)
 M=np.asarray(M)
 for p in perms:
  # do rows as array
  Mr=M[list(p),:]
  for q in perms:
   X=Mr[:,list(q)]
   d = X*X[0,0]/(X[:,[0]]*X[[0],:])
   err=np.max(np.abs(d[fixed]-F6[fixed]))
   if err<best[0]:
    # extract a/b robustly from slots
    aa=[];bb=[]
    for r in (1,3,5):
     for c in (1,4): aa.append(d[r,c]/F6[r,c])
     for c in (2,5): bb.append(d[r,c]/F6[r,c])
    A=np.mean(aa);B=np.mean(bb)
    fit=F6.copy()
    for r,c in var: fit[r,c]*=A if var[(r,c)]=='a' else B
    total=np.max(np.abs(d-fit))
    best=(err,(p,q,A,B,total))
   if best[0]<tol_fixed:
    # still continue? first candidate may not best total
    pass
 return best

def getcl(V):
 n=len(V); G=np.abs(V@V.conj().T);adj=[set(np.where((G[i]<1e-7)&(np.arange(n)!=i))[0]) for i in range(n)];cl=set()
 def bk(R,P):
  if len(R)==6:cl.add(tuple(sorted(R)));return
  while P:
   v=P.pop();bk(R+[v],[u for u in P if u in adj[v]])
 bk([],list(range(n)));return sorted(cl)

for fname in ['pool_0.400000000000.npz','pool_1.047197551197.npz','pool_2.094395102393.npz']:
 V=np.load(os.path.join(base,fname))['V']
 C=getcl(V)[0]
 B=V[list(C)].T
 print('\n',fname,'testing')
 t=time.time(); best=test(B); print('time',time.time()-t,'best',best)
