import json, numpy as np
R=json.load(open('/mnt/data/mub_fold_a4_computation/results/refined_singular_roots.json'))
R=[x for x in R if x['F_norm']<1e-7 and x['J_u_norm']<1e-7 and x['norm_err']<1e-8]
LAM=0.1114802243779665542913031975274172717685818097497045

def H(lam):
 z=np.exp(1j*lam); zi=np.conj(z)
 return np.array([[1,1,1,1,1,1],[1,-1,z,-z,1j,-1j],[1,-1j,1j,1j,-1j,-1],[1,1j,-z,z,-1,-1j],[1,zi,-1j,-1,-zi,1j],[1,-zi,-1,-1j,zi,1j]],complex)
def F(a,lam):
 v=np.r_[1.,np.exp(1j*a)]; inn=np.conj(H(lam)).T@v
 return (np.abs(inn)**2-6).real

def J(a,lam):
 h=2e-5; out=np.empty((6,5))
 for j in range(5):
  ap=a.copy();am=a.copy();ap[j]+=h;am[j]-=h;out[:,j]=(F(ap,lam)-F(am,lam))/(2*h)
 return out

def diag(a,u):
 J0=J(a,LAM); J5=J0[:5,:]
 U,S,Vh=np.linalg.svd(J0,full_matrices=False); u2=Vh[-1]; u2*=np.sign(np.dot(u2,u))
 U5,S5,Vh5=np.linalg.svd(J5,full_matrices=True); w=U5[:,-1]
 hp=2e-5
 Fl=(F(a,LAM+hp)-F(a,LAM-hp))/(2*hp)
 q=(F(a+hp*u2,LAM)-2*F(a,LAM)+F(a-hp*u2,LAM))/hp**2
 return {'sigma':S.tolist(),'sigma5':S5.tolist(),'wFl':float(np.dot(w,Fl[:5])),'wQ':float(np.dot(w,q[:5])),'abs_wFl':float(abs(np.dot(w,Fl[:5]))),'abs_wQ':float(abs(np.dot(w,q[:5])))}
D=[diag(np.array(x['a']),np.array(x['u'])) for x in R]
print('wFl abs min/max',min(x['abs_wFl'] for x in D),max(x['abs_wFl'] for x in D))
print('wQ abs min/max',min(x['abs_wQ'] for x in D),max(x['abs_wQ'] for x in D))
for i,x in enumerate(D): print(i,x['wFl'],x['wQ'],x['sigma'][-1])
json.dump(D,open('/mnt/data/mub_fold_a4_computation/results/fold_diagnostics_24.json','w'),indent=2)
