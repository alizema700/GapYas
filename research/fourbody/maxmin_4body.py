"""Planar k-body steep-decay limit, k = 4: maximise over unit-covolume lattices
Q(L) = min over 4-point sets {0,x,y,z} of distinct lattice points of the product of the 6 mutual distances."""
import itertools, json, numpy as np
ax = range(-3, 4)
LAB = np.array([(m, n) for m in ax for n in ax if (m, n) != (0, 0)], float)
def Q(x, y, K=18):
    B = np.array([[1, x], [0, y]]) / np.sqrt(y)
    V = LAB @ B.T; r = np.linalg.norm(V, axis=1); idx = np.argsort(r)[:K]; V = V[idx]; r = r[idx]
    best = np.inf; arg = None
    for i, j, k in itertools.combinations(range(K), 3):
        p = r[i] * r[j] * r[k] * np.linalg.norm(V[i] - V[j]) * np.linalg.norm(V[i] - V[k]) * np.linalg.norm(V[j] - V[k])
        if p < best: best, arg = p, (i, j, k)
    return best, [LAB[idx[a]].tolist() for a in arg]
res = []
for x in np.linspace(0, 0.5, 51):
    for y in np.linspace(np.sqrt(max(1 - x * x, 0)), 1.8, 60):
        q, a = Q(x, y); res.append((q, x, y))
res.sort(reverse=True)
print("top grid points (Q, x, y):", [tuple(round(v, 5) for v in t) for t in res[:6]])
from scipy.optimize import minimize
best = None
for q, x0, y0 in res[:10]:
    r = minimize(lambda p: -Q(p[0], p[1])[0], [x0, y0], method="Nelder-Mead", options={"xatol": 1e-10, "fatol": 1e-13, "maxiter": 4000})
    if best is None or -r.fun > best[0]: best = (-r.fun, *r.x)
q, conf = Q(best[1], best[2])
named = {"hex": Q(0.5, np.sqrt(3) / 2)[0], "square": Q(0.0, 1.0)[0]}
print("named:", named); print("maximum Q = %.8f at tau = %.8f + %.8f i; minimising set %s" % (best[0], best[1], best[2], conf))
json.dump(dict(best=best, named=named), open("maxmin_4body.json", "w"), indent=1)
