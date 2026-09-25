import numpy as np, gzl
from scipy.optimize import brentq
def A_tau(tau):
    A=np.array([[1.0, tau.real],[0.0, tau.imag]]); return A/np.sqrt(tau.imag)
H=0.5+1j*np.sqrt(3)/2
g=lambda nu: gzl.zeta_circle([nu]*3,A_tau(H))-gzl.zeta_circle([nu]*3,A_tau(1j))
nus=brentq(g,3.75,4.0,xtol=1e-6); print('nu* (hex=square) =',nus)
# energy along boundary arc |tau|=1 from square (theta=90deg) to hex (60deg) and along x=0 line at nu*
for nu in [nus-0.05,nus,nus+0.05]:
    arc=[(th,gzl.zeta_circle([nu]*3,A_tau(np.exp(1j*np.radians(th))))) for th in np.linspace(60,90,13)]
    line=[(y,gzl.zeta_circle([nu]*3,A_tau(1j*y))) for y in [1.0,1.05,1.1,1.2]]
    print('nu=%.4f arc:'%nu,' '.join('%.0f:%.6f'%a for a in arc)); print('         rect:',' '.join('%.2f:%.6f'%a for a in line))
