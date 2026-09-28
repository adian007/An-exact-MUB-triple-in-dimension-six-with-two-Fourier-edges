import numpy as np,itertools,math,os,time
base='/mnt/data/mub_b3_work/b3_a4_i4_results'
perms=np.array(list(itertools.permutations(range(6))),dtype=int)
F6=np.array([[np.exp(2j*np.pi*r*c/6) for c in range(6)] for r in range(6)],complex)
varA=[(r,c) for r in (1,3,5) for c in (1,4)]
varB=[(r,c) for r in (1,3,5) for c in (2,5)]
fixed=np.ones((6,6),bool)
for rc in varA+varB: fixed[rc]=False

def test(M):
 best=(1e9,None)
 # normalize M already unit modulus; dephase each candidate
 for p in perms:
  Mr=M[p,:]
  for q in perms:
   X=Mr[:,q]
   d=X*X[0,0]/(X[:,None,0]*X[0,None,:])
   # fixed pattern quick
   e=np.max(np.abs(d[fixed]-F6[fixed]))
   if e>2e-8: continue
   a=np.mean([d[r,c]/F6[r,c] for r,c in varA])
   b=np.mean([d[r,c]/F6[r,c] for r,c in varB])
   fit=F6.copy()
   for r,c in varA: fit[r,c]*=a
   for r,c in varB: fit[r,c]*=b
   total=np.max(np.abs(d-fit))
   if total<best[0]: best=(total,(tuple(p),tuple(q),a,b,e))
 return best

def getcl(V):
 n=len(V);G=np.abs(V@V.conj().T);adj=[set(np.where((G[i]<1e-7)&(np.arange(n)!=i))[0]) for i in range(n)];cl=set()
 def bk(R,P):
  if len(R)==6:cl.add(tuple(sorted(R)));return
  while P:
   v=P.pop();bk(R+[v],[u for u in P if u in adj[v]])
 bk([],list(range(n)));return sorted(cl)
for fname in ['pool_0.400000000000.npz','pool_1.047197551197.npz','pool_2.094395102393.npz']:
 V=np.load(os.path.join(base,fname))['V']; print('\n',fname)
 for idx,C in enumerate(getcl(V)):
  B=V[list(C)].T
  t=time.time();best=test(B)
  print('clique',idx,'time',round(time.time()-t,2),'best',best)
