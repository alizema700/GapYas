"""k-body planar steep-decay limit for k = 2..6: which lattice maximises the smallest product of all
k(k-1)/2 mutual distances of k distinct lattice points?  (k = 2: hexagonal; 3: rectangle y_inf; 4: hexagonal)"""
import itertools, sys, numpy as np
from scipy.optimize import minimize
k = int(sys.argv[1])
LAB = np.array([(m, n) for m in range(-4, 5) for n in range(-4, 5) if (m, n) != (0, 0)], float)
def Qk(x, y, K=16):
    B = np.array([[1, x], [0, y]]) / np.sqrt(y); V = LAB @ B.T
    V = V[np.argsort(np.linalg.norm(V, axis=1))[:K]]
    P = np.vstack([np.zeros(2), V])
    D = np.linalg.norm(P[:, None] - P[None, :], axis=2); np.fill_diagonal(D, 1.0)
    L = np.log(D); best = np.inf
    for c in itertools.combinations(range(1, K + 1), k - 1):
        idx = (0,) + c; s = sum(L[a, b] for a, b in itertools.combinations(idx, 2))
        best = min(best, s)
    return np.exp(best)
grid = [(Qk(x, y), x, y) for x in np.linspace(0, 0.5, 21) for y in np.linspace(max(np.sqrt(max(1 - x * x, 0)), 0.8661), 1.7, 25)]
grid.sort(reverse=True)
cands = []
for _, x0, y0 in grid[:5]:
    rr = minimize(lambda p: -Qk(p[0], p[1]), [x0, y0], method="Nelder-Mead", options={"xatol": 1e-9, "fatol": 1e-12})
    cands.append((-rr.fun, rr.x[0], rr.x[1]))
best = max(cands)
r = minimize(lambda p: -Qk(p[0], p[1]), [best[1], best[2]], method="Nelder-Mead", options={"xatol": 1e-10, "fatol": 1e-13})
print("k=%d: max Q = %.8f at tau = %.6f + %.6f i ;  hex %.8f  square %.8f" % (k, -r.fun, r.x[0], r.x[1], Qk(0.5, np.sqrt(3) / 2), Qk(0, 1)))
