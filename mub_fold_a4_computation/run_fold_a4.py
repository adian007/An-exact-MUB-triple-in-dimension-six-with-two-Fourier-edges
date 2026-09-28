"""Numerical fold/A4 analysis for the Dita circle MU-vector pool.
Uses column-wise MU equations and phase coordinates x_j=exp(i a_j), z=exp(i lambda).
This is an exploratory numerical driver; certification of singular roots is separate.
"""
import numpy as np, json, csv, os
from scipy.optimize import least_squares

LAM=0.1114802243779665542913031975274172717685818097497045

def H(lam):
    z=np.exp(1j*lam); zi=np.conj(z)
    return np.array([
        [1,1,1,1,1,1],
        [1,-1,z,-z,1j,-1j],
        [1,-1j,1j,1j,-1j,-1],
        [1,1j,-z,z,-1,-1j],
        [1,zi,-1j,-1,-zi,1j],
        [1,-zi,-1,-1j,zi,1j]
    ],complex)

def residual(a,lam):
    v=np.r_[1.0,np.exp(1j*np.asarray(a))]
    M=H(lam)
    # columns, matching the symbolic I3 construction
    inn=np.conj(M).T @ v
    return (np.abs(inn)**2-6.0).real

def jac_fd(a,lam):
    h=2e-6
    J=np.empty((6,5))
    for j in range(5):
        ap=np.array(a,float); am=np.array(a,float); ap[j]+=h; am[j]-=h
        J[:,j]=(residual(ap,lam)-residual(am,lam))/(2*h)
    return J

def solve_many(lam,n=400,seed=20260928):
    rng=np.random.default_rng(seed)
    sols=[]
    for k in range(n):
        a0=rng.uniform(-np.pi,np.pi,5)
        r=least_squares(lambda a: residual(a,lam)[:5], a0, max_nfev=1200, xtol=1e-12, ftol=1e-12, gtol=1e-12)
        rr=np.linalg.norm(residual(r.x,lam))
        if rr<1e-7:
            a=(r.x+np.pi)%(2*np.pi)-np.pi
            if not any(np.linalg.norm(np.angle(np.exp(1j*(a-s))))<1e-5 for s in sols):
                sols.append(a)
    return np.array(sols)

def A4_maps(a,lam):
    """Generators from the verified monomial automorphisms, applied to v then regauged."""
    z=np.exp(1j*lam)
    G2=np.array([[1,0,0,0,0,0],[0,0,0,1,0,0],[0,0,1,0,0,0],[0,1,0,0,0,0],[0,0,0,0,0,1],[0,0,0,0,1,0]],complex)
    G3=np.array([[0,-1,0,0,0,0],[0,0,0,0,z,0],[0,0,0,-1j,0,0],[0,0,0,0,0,-z],[1,0,0,0,0,0],[0,0,1j,0,0,0]],complex)
    v=np.r_[1,np.exp(1j*a)]
    outs=[]
    for U in (G2,G3):
        w=U@v; w=w/w[0]
        outs.append(np.angle(w[1:]))
    return outs

def close(a,b,tol=2e-5):
    return np.linalg.norm(np.angle(np.exp(1j*(a-b))))<tol

def orbit(a,sols,lam):
    # Generate group by BFS using the two generators numerically.
    seen=[]; queue=[a]
    while queue and len(seen)<12:
        x=queue.pop(0)
        if any(close(x,y) for y in seen): continue
        seen.append(x)
        for y in A4_maps(x,lam): queue.append(y)
    inds=[]
    for x in seen:
        ds=[np.linalg.norm(np.angle(np.exp(1j*(x-s)))) for s in sols]
        inds.append(int(np.argmin(ds)) if min(ds)<2e-4 else -1)
    return seen,inds

def main():
    os.makedirs('results',exist_ok=True)
    for label,lam in [('below',LAM-2e-6),('critical',LAM),('above',LAM+2e-6)]:
        sols=solve_many(lam)
        np.savetxt(f'results/{label}_solutions.csv',sols,delimiter=',')
        with open(f'results/{label}_summary.json','w') as f: json.dump({'lambda':lam,'n_distinct_found':len(sols)},f,indent=2)
        print(label,lam,len(sols))
    # Use critical solutions as the pool; inspect A4 closure.
    sols=np.loadtxt('results/critical_solutions.csv',delimiter=',')
    if sols.ndim==1: sols=sols[None,:]
    rows=[]; assigned=set(); orbit_id=0
    for i,a in enumerate(sols):
        if i in assigned: continue
        pts,inds=orbit(a,sols,LAM)
        valid=[j for j in inds if j>=0]
        orbit_id+=1
        for j in valid: assigned.add(j)
        rows.append({'orbit':orbit_id,'seed':i,'size':len(set(valid)),'matched_indices':sorted(set(valid)),'closure_ok':len(valid)==len(set(valid)) and len(valid)>0})
    with open('results/a4_orbits.json','w') as f: json.dump(rows,f,indent=2)
    print('orbits',[(r['orbit'],r['size'],r['closure_ok']) for r in rows])

if __name__=='__main__': main()
