"""Hessian of T_nu at BCC in the Gram chart q = (g22,g33,g12,g13,g23), g11 = 1 (scale invariant form):
smallest eigenvalue vs nu. Direct summation over labels in [-R,R]^3 (log-sum-exp scaled by (3/2)^nu)."""
import itertools, numpy as np, sys
R = int(sys.argv[1]) if len(sys.argv) > 1 else 5
labs = np.array([v for v in itertools.product(range(-R, R + 1), repeat=3) if any(v)], float)
def gram(q): g22, g33, g12, g13, g23 = q; return np.array([[1, g12, g13], [g12, g22, g23], [g13, g23, g33]])
def T(q, nu):
    G = gram(q); det = np.linalg.det(G)
    l = np.einsum('ij,jk,ik->i', labs, G, labs); ll = np.log(l)
    tot = 0.0
    for s in range(0, len(labs), 400):
        D = labs[s:s + 400, None, :] - labs[None, :, :]
        ld = np.einsum('abj,jk,abk->ab', D, G, D)
        with np.errstate(divide='ignore'): L = 0.5 * (ll[s:s + 400, None] + ll[None, :] + np.log(ld) - np.log(det)) - np.log(1.5)
        L[~np.isfinite(L)] = np.inf
        tot += np.sum(np.exp(-nu * L))
    return tot
q0 = np.array([1, 1, 1/3, 1/3, -1/3])
for nu in [float(v) for v in sys.argv[2:]]:
    h = 1e-3; E = np.eye(5); H = np.zeros((5, 5)); f0 = T(q0, nu)
    g = np.array([(T(q0 + h * E[i], nu) - T(q0 - h * E[i], nu)) / (2 * h) for i in range(5)])
    for i in range(5):
        for j in range(i, 5):
            H[i, j] = H[j, i] = (T(q0 + h*E[i] + h*E[j], nu) - T(q0 + h*E[i] - h*E[j], nu) - T(q0 - h*E[i] + h*E[j], nu) + T(q0 - h*E[i] - h*E[j], nu)) / (4 * h * h)
    ev = np.linalg.eigvalsh(H / f0)
    print("nu=%5.1f  |grad|/T=%.1e  eig(Hess)/T: %s" % (nu, np.linalg.norm(g) / f0, np.array2string(ev, precision=3)), flush=True)
