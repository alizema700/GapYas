"""3D: maximise P(L) = min triangle product over unit-covolume 3D lattices (nu -> infinity limit)."""
import numpy as np, sys
from scipy.optimize import minimize
sys.path.insert(0, "../verification")
from lattice_sums import FCC, BCC, SC, unit_covolume
from v3_3d_lattices import lll, A_from_params, params_of

def P_of(A, R=3, rmax=3.2):
    A = lll(unit_covolume(A))
    ax = np.arange(-R, R + 1)
    m = np.stack(np.meshgrid(ax, ax, ax, indexing="ij"), -1).reshape(-1, 3)
    m = m[np.any(m != 0, axis=1)]
    X = m @ A.T; r = np.linalg.norm(X, axis=1)
    k = r < rmax; X, r = X[k], r[k]
    D = np.linalg.norm(X[:, None] - X[None], axis=-1); np.fill_diagonal(D, np.inf)
    return float(np.min(np.outer(r, r) * D))

for n, A in [("FCC", FCC), ("BCC", BCC), ("SC", SC)]:
    print(n, "P =", P_of(A))
rng = np.random.default_rng(3)
f = lambda p: -P_of(A_from_params(p)) if np.all(np.abs(p[:2]) < 1.2) else 0.0
res = []
starts = [params_of(BCC), params_of(FCC)] + [np.concatenate([rng.normal(0, .25, 2), rng.normal(0, .5, 3)]) for _ in range(60)]
for s in starts:
    r = minimize(f, s, method="Nelder-Mead", options=dict(xatol=1e-9, fatol=1e-12, maxiter=4000))
    r = minimize(f, r.x, method="Nelder-Mead", options=dict(xatol=1e-10, fatol=1e-13, maxiter=4000))
    res.append((-r.fun, r.x))
res.sort(key=lambda z: -z[0])
for P, x in res[:8]:
    A = lll(unit_covolume(A_from_params(x))); G = A.T @ A
    print("P = %.10f   Gram(reduced) diag %s  offdiag %s" % (P, np.round(np.diag(G), 5), np.round([G[0,1], G[0,2], G[1,2]], 5)))
np.save("maxmin_3d_best.npy", A_from_params(res[0][1]))
