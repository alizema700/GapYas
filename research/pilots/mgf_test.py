import numpy as np, gzl, epsteinlib, mpmath as mp, time
def A_tau(tau): return np.array([[1.0, tau.real],[0.0, tau.imag]])
def E(s,tau):
    A=A_tau(tau); z=epsteinlib.epstein_zeta(2*s,A,np.zeros(2),np.zeros(2)).real
    return (tau.imag/np.pi)**s*z
def C(avec,tau):
    A=A_tau(tau); w=sum(avec)
    return (tau.imag/np.pi)**w*gzl.zeta_circle([2*a for a in avec],A)
zeta5=float(mp.zeta(5)); zeta3=float(mp.zeta(3))
for tau in [1j, 0.5+np.sqrt(3)/2*1j, 0.21+1.37j, -0.4+0.95j]:
    rhs=0.4*E(5,tau)+zeta5/30
    print('tau',tau,'  RHS 2/5E5+z5/30 =',repr(rhs))
    for eps in [1e-1,1e-2,1e-3,1e-4,1e-6]:
        t=time.time(); lhs=C([2,2,1+eps/2],tau); dt=time.time()-t
        print('   eps',eps,'C_{2,2,1+eps/2}=',repr(lhs),' rel.diff',abs(lhs-rhs)/rhs,' %.2fs'%dt)
    # E_3 + zeta3 test for C_{1,1,1} via limit
    lhs=C([1+1e-6,1+1e-6,1+1e-6],tau); print('   C_{1,1,1}(eps=1e-6)',lhs,' E3+z3=',E(3,tau)+zeta3, 'rel',abs(lhs-E(3,tau)-zeta3)/(E(3,tau)+zeta3))
