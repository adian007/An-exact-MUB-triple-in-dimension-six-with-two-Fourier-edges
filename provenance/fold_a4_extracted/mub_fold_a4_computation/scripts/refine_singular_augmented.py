import numpy as np, json, os
from scipy.optimize import least_squares
LAM=0.1114802243779665542913031975274172717685818097497045

def H(lam):
 z=np.exp(1j*lam); zi=np.conj(z)
 return np.array([[1,1,1,1,1,1],[1,-1,z,-z,1j,-1j],[1,-1j,1j,1j,-1j,-1],[1,1j,-z,z,-1,-1j],[1,zi,-1j,-1,-zi,1j],[1,-zi,-1,-1j,zi,1j]],complex)

def F(a,lam=LAM):
 v=np.r_[1.,np.exp(1j*a)]; M=H(lam); inn=np.conj(M).T@v
 return (np.abs(inn)**2-6).real

def J(a,lam=LAM):
 # complex-step impossible through real exp; central diff
 h=1e-6; out=np.empty((6,5))
 for j in range(5):
  ap=a.copy(); am=a.copy(); ap[j]+=h; am[j]-=h
  out[:,j]=(F(ap,lam)-F(am,lam))/(2*h)
 return out

def obj(y):
 a=y[:5]; u=y[5:]
 ff=F(a)
 jj=J(a)
 return np.r_[ff, jj@u, np.dot(u,u)-1]

cand=np.loadtxt('/mnt/data/mub_fold_a4_computation/results/critical_singular_candidates.csv',delimiter=',')
rows=[]
for k,a0 in enumerate(cand):
 u=np.linalg.svd(J(a0),full_matrices=False)[2][-1]
 for s in [1,-1]:
  y0=np.r_[a0,s*u]
  r=least_squares(obj,y0,max_nfev=4000,xtol=1e-14,ftol=1e-14,gtol=1e-14,verbose=0)
  a=(r.x[:5]+np.pi)%(2*np.pi)-np.pi
  ff=np.linalg.norm(F(a)); j=J(a); sv=np.linalg.svd(j,compute_uv=False); null=np.linalg.norm(j@r.x[5:]); normerr=abs(np.dot(r.x[5:],r.x[5:])-1)
  rows.append(dict(seed=k,sign=s,cost=float(r.cost),optimality=float(r.optimality),success=bool(r.success),F_norm=float(ff),J_u_norm=float(null),norm_err=float(normerr),sigma_min=float(sv[-1]),a=a.tolist(),u=r.x[5:].tolist()))
# dedupe a
uniq=[]
for row in rows:
 a=np.array(row['a'])
 if not any(np.linalg.norm(np.angle(np.exp(1j*(a-np.array(q['a'])))))<1e-6 for q in uniq): uniq.append(row)
print('raw',len(rows),'unique',len(uniq))
print('sigma min range',min(r['sigma_min'] for r in uniq),max(r['sigma_min'] for r in uniq))
print('F max',max(r['F_norm'] for r in uniq),'Ju max',max(r['J_u_norm'] for r in uniq),'cost max',max(r['cost'] for r in uniq))
os.makedirs('/mnt/data/mub_fold_a4_computation/results',exist_ok=True)
with open('/mnt/data/mub_fold_a4_computation/results/refined_singular_roots.json','w') as f: json.dump(uniq,f,indent=2)
np.savetxt('/mnt/data/mub_fold_a4_computation/results/refined_singular_roots.csv',np.array([r['a'] for r in uniq]),delimiter=',')
