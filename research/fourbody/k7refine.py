"""Refine the k-body max-min lattice (nonsmooth): maximise t s.t. log-product of each near-active set >= t, SLP in (x, y)."""
import itertools, sys, json, numpy as np
from scipy.optimize import linprog
k = int(sys.argv[1]); x, y = float(sys.argv[2]), float(sys.argv[3])
LAB = np.array([(m, n) for m in range(-5, 6) for n in range(-5, 6) if (m, n) != (0, 0)], float)
def sets_and_vals(x, y, K=18, slack=1e-3):
    B = np.array([[1, x], [0, y]]) / np.sqrt(y); V = LAB @ B.T; idx = np.argsort(np.linalg.norm(V, axis=1))[:K]
    lab = LAB[idx]; P = np.vstack([np.zeros(2), V[idx]]); D = np.linalg.norm(P[:, None] - P[None, :], axis=2); np.fill_diagonal(D, 1); L = np.log(D)
    out = []
    for c in itertools.combinations(range(1, K + 1), k - 1):
        idxs = (0,) + c; s = sum(L[a, b] for a, b in itertools.combinations(idxs, 2)); out.append((s, [tuple(lab[i - 1]) for i in c]))
    m = min(v for v, _ in out); return m, [(v, S) for v, S in out if v <= m + slack]
def logprod(S, x, y):
    pts = [(0.0, 0.0)] + list(S); s = 0.0
    for a, b in itertools.combinations(range(k), 2):
        m, n = pts[a][0] - pts[b][0], pts[a][1] - pts[b][1]; s += 0.5 * np.log(((m + n * x) ** 2) / y + n * n * y)
    return s
tr = 0.01
for it in range(200):
    m, act = sets_and_vals(x, y)
    A, b = [], []
    for v, S in act:
        f0 = logprod(S, x, y); h = 1e-7
        gx = (logprod(S, x + h, y) - logprod(S, x - h, y)) / (2 * h); gy = (logprod(S, x, y + h) - logprod(S, x, y - h)) / (2 * h)
        A.append([-gx, -gy, 1.0]); b.append(f0)
    res = linprog([0, 0, -1], A_ub=A, b_ub=b, bounds=[(-tr, tr), (-tr, tr), (None, None)], method="highs")
    xn, yn = x + res.x[0], y + res.x[1]
    mn, _ = sets_and_vals(xn, yn)
    if mn > m + 1e-15: x, y = xn, yn; tr = min(tr * 1.5, 0.02)
    else: tr *= 0.5
    if tr < 1e-13: break
m, act = sets_and_vals(x, y, slack=1e-9)
print("k=%d: optimum tau = %.12f + %.12f i, Q = %.10f, active sets %d" % (k, x, y, np.exp(m), len(act)))
for v, S in act[:12]: print("   ", S)
json.dump(dict(k=k, x=x, y=y, Q=float(np.exp(m)), active=[[list(p) for p in S] for v, S in act]), open("k%d_opt.json" % k, "w"), indent=1)
