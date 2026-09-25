"""V3: 3D three-body lattice energy T_nu over ALL 3D Bravais lattices (unit covolume).

(a) accurate FCC / BCC / SC values vs nu (reference sums),
(b) global search: random starts + Nelder-Mead over the 5-dim space of lattice shapes,
(c) local stability (5x5 Hessian) of FCC and BCC,
(d) Bain path fcc <-> bcc (body-centred tetragonal family).
"""
import json, sys, time
import numpy as np
from scipy.optimize import minimize
from lattice_sums import triangle_sum, triangle_sum_box, FCC, BCC, SC, unit_covolume

def lll(B, delta=0.75):
    B = np.array(B, float).T.copy()   # rows = basis vectors
    n = len(B)
    def gs(B):
        Bs = np.zeros_like(B); mu = np.zeros((n, n))
        for i in range(n):
            Bs[i] = B[i]
            for j in range(i):
                mu[i, j] = B[i] @ Bs[j] / (Bs[j] @ Bs[j]); Bs[i] -= mu[i, j] * Bs[j]
        return Bs, mu
    k = 1
    Bs, mu = gs(B)
    while k < n:
        for j in range(k - 1, -1, -1):
            q = round(mu[k, j])
            if q:
                B[k] -= q * B[j]; Bs, mu = gs(B)
        if Bs[k] @ Bs[k] >= (delta - mu[k, k - 1] ** 2) * (Bs[k - 1] @ Bs[k - 1]):
            k += 1
        else:
            B[[k, k - 1]] = B[[k - 1, k]]; Bs, mu = gs(B); k = max(k - 1, 1)
    return B.T

def A_from_params(p):
    # upper-triangular basis with log-diagonal, det = 1
    l1, l2, b, c, e = p
    a, d = np.exp(l1), np.exp(l2); f = 1.0 / (a * d)
    return np.array([[a, b, c], [0, d, e], [0, 0, f]])

def T(A, nu, fast=False):
    A = lll(unit_covolume(A))
    if fast:
        return triangle_sum_box(A, (nu,) * 3, 20)
    return triangle_sum(A, (nu,) * 3, Rs=(16, 24, 32, 40, 48))

def shells(A, nmax=3):
    A = lll(unit_covolume(A)); ax = np.arange(-4, 5)
    m = np.stack(np.meshgrid(ax, ax, ax, indexing="ij"), -1).reshape(-1, 3)
    r = np.linalg.norm(m @ A.T, axis=1); r = np.sort(r[r > 0])
    out = []; i = 0
    while len(out) < nmax:
        j = i
        while j < len(r) and abs(r[j] - r[i]) < 1e-6: j += 1
        out.append((round(float(r[i]), 5), j - i)); i = j
    return out

def params_of(A):
    # QR -> upper triangular with positive diagonal
    A = unit_covolume(A); Q, Rm = np.linalg.qr(A)
    S = np.diag(np.sign(np.diag(Rm))); Rm = S @ Rm
    return np.array([np.log(Rm[0, 0]), np.log(Rm[1, 1]), Rm[0, 1], Rm[0, 2], Rm[1, 2]])

def main():
    res = {}
    nus = [float(a) for a in sys.argv[1:]] or [3.2, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 9.0, 12.0]
    rng = np.random.default_rng(7)
    for nu in nus:
        t0 = time.time()
        ref = {n: T(A, nu) for n, A in [("FCC", FCC), ("BCC", BCC), ("SC", SC)]}
        # global search
        f = lambda p: T(A_from_params(p), nu, fast=True) if np.all(np.abs(p[:2]) < 1.5) else 1e9
        starts = [params_of(FCC), params_of(BCC), params_of(SC)] + \
                 [np.concatenate([rng.normal(0, 0.25, 2), rng.normal(0, 0.5, 3)]) for _ in range(13)]
        found = []
        for s in starts:
            r = minimize(f, s, method="Nelder-Mead", options=dict(xatol=1e-5, fatol=1e-11, maxiter=2500, maxfev=2500))
            found.append((r.fun, r.x))
        found.sort(key=lambda z: z[0])
        bestA = A_from_params(found[0][1])
        best_acc = T(bestA, nu)
        distinct = []
        for v, x in found:
            if all(abs(v - w) > 1e-6 * v for w, _ in distinct):
                distinct.append((v, x))
        res[nu] = dict(reference=ref, best_found=best_acc, best_shells=shells(bestA),
                       local_minima=[(float(v), shells(A_from_params(x))) for v, x in distinct[:6]])
        print("nu=%.2f  FCC %.12f  BCC %.12f  SC %.12f | best found %.12f shells %s | (%.0fs)"
              % (nu, ref["FCC"], ref["BCC"], ref["SC"], best_acc, shells(bestA), time.time() - t0), flush=True)
        for v, sh in res[nu]["local_minima"]:
            print("     local min  %.10f  shells %s" % (v, sh), flush=True)
    json.dump(res, open("results/v3_3d_%s.json" % "_".join(str(n) for n in nus), "w"), indent=1, default=float)


if __name__ == '__main__':
    main()
