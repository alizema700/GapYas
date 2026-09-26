"""Which lattice in R^d (unit covolume) maximises P(L) = min |x||y||x-y| over distinct nonzero x, y in L?
d = 2: rectangle (Theorem 1), d = 3: BCC (Theorem 3).  Here: d = 4 (and d = 2, 3 as a check).
Method: sequential linear programming (SLP) on the Gram matrix with a trust region, many random starts.
P is scale-normalised as P / sqrt(det G) (covolume 1)."""
import itertools, sys, json, numpy as np
from scipy.optimize import linprog
d = int(sys.argv[1]) if len(sys.argv) > 1 else 4
nstart = int(sys.argv[2]) if len(sys.argv) > 2 else 60
rng = np.random.default_rng(int(sys.argv[3]) if len(sys.argv) > 3 else 1)
C = np.array([v for v in itertools.product(range(-2, 3), repeat=d) if any(v)], float)   # coefficient vectors

def lll_gram(G):
    """cheap reduction: pairwise size reduction + sorting (Minkowski-ish) to keep short vectors in C"""
    B = np.linalg.cholesky(G).T; n = d
    for _ in range(50):
        changed = False
        for i in range(n):
            for j in range(n):
                if i == j: continue
                mu = np.round(B[:, i] @ B[:, j] / (B[:, j] @ B[:, j]))
                if mu != 0:
                    B[:, i] -= mu * B[:, j]; changed = True
        idx = np.argsort(np.sum(B * B, 0)); B = B[:, idx]
        if not changed: break
    return B.T @ B

def triangles(G, slack=1.0005):
    """all (x, y) coefficient pairs with product <= slack * P (candidates from C), plus the values"""
    L = np.einsum('ij,jk,ik->i', C, G, C)
    lam2 = L.min()
    # sides must satisfy |x| <= P / lam1^2; estimate P from a quick pass
    order = np.argsort(L); keep = order[:min(len(C), 400)]
    Ck = C[keep]; Lk = L[keep]
    D = Ck[:, None, :] - Ck[None, :, :]
    LD = np.einsum('abj,jk,abk->ab', D, G, D)
    prod2 = Lk[:, None] * Lk[None, :] * LD
    np.fill_diagonal(prod2, np.inf); prod2[LD < 1e-12] = np.inf
    P2 = prod2.min()
    ii, jj = np.where(prod2 <= P2 * slack ** 2)
    return P2, [(Ck[a], Ck[b]) for a, b in zip(ii, jj) if a < b]

def Pnorm(G):
    P2, _ = triangles(G, 1.0)
    return np.sqrt(P2) / np.linalg.det(G) ** 0.5

iu = np.triu_indices(d)
def vec2G(g): G = np.zeros((d, d)); G[iu] = g; return G + G.T - np.diag(np.diag(G))
def G2vec(G): return G[iu]

def slp(G, iters=300):
    tr = 0.05
    for it in range(iters):
        G = lll_gram(G); G = G / np.linalg.det(G) ** (1 / d)
        P2, act = triangles(G, 1.02)
        g0 = G2vec(G); n = len(g0)
        # maximise t subject to log(prod_t(G + dG)) >= t (linearised), log det(G + dG) = 0 (linearised), |dg| <= tr
        Ginv = np.linalg.inv(G)
        def dlog_l(v):  # gradient of log(v^T G v) wrt upper-tri entries
            M = np.outer(v, v); M = M + M.T - np.diag(np.diag(M)); return M[iu] / (v @ G @ v)
        A_ub, b_ub = [], []
        for x, y in act:
            grad = dlog_l(x) + dlog_l(y) + dlog_l(x - y)
            val = np.log((x @ G @ x) * (y @ G @ y) * ((x - y) @ G @ (x - y)))
            A_ub.append(np.concatenate([-grad, [1.0]])); b_ub.append(val)
        gdet = (2 * Ginv - np.diag(np.diag(Ginv)))[iu]
        res = linprog(np.concatenate([np.zeros(n), [-1.0]]), A_ub=np.array(A_ub), b_ub=np.array(b_ub),
                      A_eq=np.array([np.concatenate([gdet, [0.0]])]), b_eq=[0.0],
                      bounds=[(-tr, tr)] * n + [(None, None)], method="highs")
        if not res.success: tr *= 0.5; continue
        Gn = vec2G(g0 + res.x[:n])
        if np.min(np.linalg.eigvalsh(Gn)) <= 0: tr *= 0.5; continue
        Gn = Gn / np.linalg.det(Gn) ** (1 / d)
        if Pnorm(Gn) > Pnorm(G) + 1e-13: G = Gn; tr = min(tr * 1.5, 0.2)
        else: tr *= 0.5
        if tr < 1e-10: break
    G = lll_gram(G); return G / np.linalg.det(G) ** (1 / d)

def named():
    out = {}
    out["Z^d"] = np.eye(d)
    if d >= 3:
        Dn = np.array([[1, -1] + [0] * (d - 2)] + [[0] * i + [1, -1] + [0] * (d - 2 - i) for i in range(1, d - 1)] + [[1, 1] + [0] * (d - 2)], float)
        out["D_d"] = Dn @ Dn.T
        Bs = np.eye(d); Bs[-1] = 0.5                       # D_d^* = Z^d + (1/2,...,1/2)
        out["D_d^*"] = Bs @ Bs.T
    A = np.eye(d + 1)[:, :] ; roots = np.array([np.eye(d + 1)[i] - np.eye(d + 1)[i + 1] for i in range(d)])
    out["A_d"] = roots @ roots.T
    out["A_d^*"] = np.linalg.inv(out["A_d"])
    return {k: v / np.linalg.det(v) ** (1 / d) for k, v in out.items()}

res = {"d": d, "named": {k: float(Pnorm(lll_gram(G))) for k, G in named().items()}}
print("named lattices:", {k: round(v, 6) for k, v in res["named"].items()}, flush=True)
best = []
for s in range(nstart):
    M = rng.normal(size=(d, d)); G = M @ M.T + 0.3 * np.eye(d)
    G = slp(G / np.linalg.det(G) ** (1 / d))
    P = Pnorm(G); ev = np.linalg.eigvalsh(G)
    best.append((P, G)); print("start %3d  P = %.8f  Gram eigenvalues %s" % (s, P, np.round(ev, 4)), flush=True)
best.sort(key=lambda t: -t[0])
Gb = best[0][1]; P2, act = triangles(Gb, 1.0 + 1e-7)
res["best_P"] = float(best[0][0]); res["best_G"] = Gb.tolist(); res["top10"] = [float(p) for p, _ in best[:10]]
res["n_active_triangles"] = len(act)
L = np.sort(np.einsum('ij,jk,ik->i', C, Gb, C)); res["shortest_norms2"] = L[:30].tolist()
json.dump(res, open("maxmin_d%d.json" % d, "w"), indent=1)
print("BEST P = %.10f ; active triangles %d ; shortest |v|^2: %s" % (best[0][0], len(act), np.round(L[:24], 5)))
