import numpy as np, gzl, time
def A_tau(tau):
    A=np.array([[1.0, tau.real],[0.0, tau.imag]]); return A/np.sqrt(tau.imag)
# 2D: scan fundamental domain for triangle (3-body) Riesz energy with nu_e=nu on all edges, unit covolume
for nu in [2.5, 4.0, 6.0]:
    best=None
    for x in np.linspace(0,0.5,11):
        for y in np.linspace(np.sqrt(1-x*x), 2.0, 16):
            v=gzl.zeta_circle([nu]*3, A_tau(x+1j*y))
            if best is None or v<best[0]: best=(v,x,y)
    hexv=gzl.zeta_circle([nu]*3, A_tau(0.5+1j*np.sqrt(3)/2)); sq=gzl.zeta_circle([nu]*3, A_tau(1j))
    print('2D nu=%.1f  grid-min %.12g at tau=%.3f+%.3fi | hex %.12g | square %.12g'%(nu,best[0],best[1],best[2],hexv,sq))
# 3D: FCC vs BCC vs SC at unit covolume, three-body triangle zeta with equal exponents; also 2-body Epstein for reference
fcc=np.array([[0,1,1],[1,0,1],[1,1,0]],float); bcc=np.array([[-1,1,1],[1,-1,1],[1,1,-1]],float); sc=np.eye(3)
def norm(A): return A/abs(np.linalg.det(A))**(1/3)
for nu in [3.5,4.5,6.0,9.0]:
    t=time.time(); r={n:gzl.zeta_circle([nu]*3,norm(A)) for n,A in [('FCC',fcc),('BCC',bcc),('SC',sc)]}
    print('3D 3-body nu=%.1f'%nu, {k:round(v,10) for k,v in r.items()}, ' %.1fs'%(time.time()-t))
