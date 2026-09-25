import numpy as np, gzl, time, mpmath as mp
exec(open('mgf_test.py').read().split('zeta5=')[0])
zeta5=float(mp.zeta(5)); zeta3=float(mp.zeta(3))
def rich(f, h=1e-3):
    # quadratic Richardson in eps using h,2h,4h
    f1,f2,f4=f(h),f(2*h),f(4*h)
    return (8*f1-6*f2+f4)/3
for tau in [1j, 0.21+1.37j, -0.4+0.95j]:
    r=0.4*E(5,tau)+zeta5/30; l=rich(lambda e: C([2,2,1+e],tau))
    r3=E(3,tau)+zeta3; l3=rich(lambda e: C([1+e,1+e,1+e],tau))
    print(tau,' C221 extrap',repr(l),' rel err vs 2/5E5+z5/30:',abs(l-r)/r,' | C111 extrap rel err vs E3+z3:',abs(l3-r3)/r3)
# tetrahedral K4 MGF with all a=2 (nu=4) on tau=i and generic tau
edges=np.array([[0,1],[0,2],[0,3],[1,2],[1,3],[2,3]])
for tau in [1j,0.21+1.37j]:
    A=A_tau(tau)
    for n in [32,48,64]:
        t=time.time()
        try:
            v=gzl.evaluate_graph(edges, np.full(6,4.0), A, n_points=n)
            print('K4 nu=4 tau',tau,'n',n,'zeta=',repr(float(v)),' MGF value',repr(float(v)*(tau.imag/np.pi)**12),' %.1fs'%(time.time()-t))
        except Exception as ex: print('ERR',type(ex).__name__,ex); break
