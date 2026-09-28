import numpy as np,itertools,math,os,time,json
base='/mnt/data/mub_b3_work/b3_a4_i4_results'
perms=list(itertools.permutations(range(6)))
F6=np.array([[np.exp(2j*np.pi*r*c/6) for c in range(6)] for r in range(6)],complex)
varA=[(r,c) for r in (1,3,5) for c in (1,4)]
varB=[(r,c) for r in (1,3,5) for c in (2,5)]
fixed=np.ones((6,6),bool)
for rc in varA+varB: fixed[rc]=False

def fourier_fit(M):
 best=(1e9,None)
 for p in perms:
  Mr=M[list(p),:]
  for q in perms:
   X=Mr[:,list(q)]
   d=X*X[0,0]/(X[:,None,0]*X[0,None,:])
   if np.max(np.abs(d[fixed]-F6[fixed]))>2e-8: continue
   aa=np.array([d[r,c]/F6[r,c] for r,c in varA]); bb=np.array([d[r,c]/F6[r,c] for r,c in varB])
   a=aa.mean(); b=bb.mean()
   consistency=max(np.max(np.abs(aa-a)),np.max(np.abs(bb-b)))
   fit=F6.copy()
   for r,c in varA: fit[r,c]*=a
   for r,c in varB: fit[r,c]*=b
   total=np.max(np.abs(d-fit))
   score=max(total,consistency,abs(abs(a)-1),abs(abs(b)-1))
   if score<best[0]: best=(score,(tuple(p),tuple(q),complex(a),complex(b),float(total),float(consistency)))
 return best

def getcl(V):
 n=len(V);G=np.abs(V@V.conj().T);adj=[set(np.where((G[i]<1e-7)&(np.arange(n)!=i))[0]) for i in range(n)];cl=set()
 def bk(R,P):
  if len(R)==6:cl.add(tuple(sorted(R)));return
  while P:
   v=P.pop(); bk(R+[v],[u for u in P if u in adj[v]])
 bk([],list(range(n)));return sorted(cl)

def HD(lam):
 z=np.exp(1j*lam)
 return np.array([[1,1,1,1,1,1],[1,-1,z,-z,1j,-1j],[1,-1j,1j,1j,-1j,-1],[1,1j,-z,z,-1,-1j],[1,z.conjugate(),-1j,-1,-z.conjugate(),1j],[1,-z.conjugate(),-1,-1j,z.conjugate(),1j]],complex)/np.sqrt(6)

outs={}
for lam,fname in [(0.4,'pool_0.400000000000.npz'),(math.pi/3,'pool_1.047197551197.npz'),(2*math.pi/3,'pool_2.094395102393.npz')]:
 V=np.load(os.path.join(base,fname))['V']; C=getcl(V)[0]; B=V[list(C)].T/np.sqrt(6); H=HD(lam)
 for label,M in [('B',B),('HdagB',H.conj().T@B)]:
  t=time.time(); res=fourier_fit(M); print(lam,label,'time',time.time()-t,'score',res)
  outs[f'{lam:.15f}_{label}']=res
print(json.dumps({k:{'score':float(v[0]),'fit':None if v[1] is None else {'p':v[1][0],'q':v[1][1],'a':str(v[1][2]),'b':str(v[1][3]),'total':v[1][4],'consistency':v[1][5]}} for k,v in outs.items()},indent=2))
