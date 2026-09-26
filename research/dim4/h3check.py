"""Check (H3) at a max-min lattice: active gradients of log(product) in the tangent space of unit-covolume
Gram matrices must (a) span it and (b) admit a strictly positive combination equal to 0."""
import itertools, sys, numpy as np
from scipy.optimize import linprog
def check(G, k=3, cr=2, tol=1e-7, label=""):
    d = len(G); iu = np.triu_indices(d)
    C = np.array([v for v in itertools.product(range(-cr, cr + 1), repeat=d) if any(v)], float)
    L = np.einsum('ij,jk,ik->i', C, G, C); keep = np.argsort(L)[:300]; C = C[keep]
    sets = []; best = np.inf
    for combo in itertools.combinations(range(len(C)), k - 1):
        pts = [np.zeros(d)] + [C[i] for i in combo]
        p = 1.0
        for a, b in itertools.combinations(range(k), 2):
            v = pts[a] - pts[b]; p *= v @ G @ v
        if p < best * (1 - tol): best = p; sets = [pts]
        elif p <= best * (1 + tol): sets.append(pts)
    Ginv = np.linalg.inv(G)
    def grad(pts):
        g = np.zeros(len(iu[0]))
        for a, b in itertools.combinations(range(k), 2):
            v = pts[a] - pts[b]; M = np.outer(v, v); M = M + M.T - np.diag(np.diag(M)); g += M[iu] / (v @ G @ v)
        return g
    A = np.array([grad(s) for s in sets])
    # project out the scaling direction: tangent space {dG : tr(G^-1 dG) = 0}; use metric-free projection on the vector of entries
    n = (2 * Ginv - np.diag(np.diag(Ginv)))[iu]; n = n / np.linalg.norm(n)
    Ap = A - np.outer(A @ n, n)
    rank = np.linalg.matrix_rank(Ap, tol=1e-8)
    m = len(sets)
    lp = linprog(np.zeros(m), A_eq=Ap.T, b_eq=np.zeros(Ap.shape[1]), bounds=[(1, None)] * m, method="highs")
    U = np.unique(np.round(Ap, 6), axis=0)
    print("%s: min product^(1/2) = %.8f, active sets %d, distinct gradients %d, rank %d / %d, strictly positive combination: %s"
          % (label, np.sqrt(best), m, len(U), rank, len(iu[0]) - 1, lp.status == 0))
if __name__ == "__main__":
    d = 4
    Dn = np.array([[1, -1, 0, 0], [0, 1, -1, 0], [0, 0, 1, -1], [1, 1, 0, 0]], float); G = Dn @ Dn.T; G /= np.linalg.det(G) ** (1 / 4)
    check(G, 3, 2, label="D4, three-body")
    B = np.array([[1, 0.5], [0, np.sqrt(3) / 2]]); G = B.T @ B; G /= np.sqrt(np.linalg.det(G))
    check(G, 4, 3, label="hexagonal, four-body")
    Bs = np.eye(3); Bs[-1] = 0.5; G = Bs @ Bs.T; G /= np.linalg.det(G) ** (1 / 3)
    check(G, 3, 2, label="BCC, three-body (reference)")
