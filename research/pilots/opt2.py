import numpy as np, gzl
from scipy.optimize import minimize
def A_tau(tau):
    A=np.array([[1.0, tau.real],[0.0, tau.imag]]); return A/np.sqrt(tau.imag)
def f(p,nu):
    x,y=p
    if x<0 or x>0.5 or x*x+y*y<1 or y>3: return 1e9
    return gzl.zeta_circle([nu]*3, A_tau(x+1j*y))
for nu in [2.5,3.0,3.25,3.5,3.75,4.0,5.0]:
    res=[]
    for st in [(0.5,0.87),(0.0,1.0),(0.25,1.0),(0.0,1.3),(0.45,0.9),(0.1,1.05)]:
        r=minimize(f,st,args=(nu,),method='Nelder-Mead',options={'xatol':1e-5,'fatol':1e-12,'maxiter':400})
        res.append((r.fun,r.x[0],r.x[1]))
    res.sort()
    h=f((0.5,np.sqrt(3)/2),nu); s=f((0,1),nu)
    print('nu=%.2f best %.10f at tau=%.4f+%.4fi | hex %.10f square %.10f'%(nu,res[0][0],res[0][1],res[0][2],h,s))
# brute force sanity check: square lattice nu=6 via direct summation
x=np.arange(-40,41); X,Y=np.meshgrid(x,x); P=np.stack([X.ravel(),Y.ravel()],1).astype(float)
r=np.linalg.norm(P,axis=1); m=r>0; P=P[m]; r=r[m]
sel=r<12; Ps=P[sel]; rs=r[sel]
D=np.linalg.norm(Ps[:,None,:]-Ps[None,:,:],axis=2); np.fill_diagonal(D,np.inf)
print('direct triangle sum square nu=6 (|x|,|y|<12):', np.sum(np.outer(rs**-6,rs**-6)*D**-6), ' GZL:', gzl.zeta_circle([6.0]*3,np.eye(2)))
