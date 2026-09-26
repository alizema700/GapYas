import itertools, numpy as np, sys
LAB = np.array([(m, n) for m in range(-6, 7) for n in range(-6, 7) if (m, n) != (0, 0)], float)
def Qk(x, y, k, K):
    B = np.array([[1, x], [0, y]]) / np.sqrt(y); V = LAB @ B.T
    V = V[np.argsort(np.linalg.norm(V, axis=1))[:K]]
    P = np.vstack([np.zeros(2), V]); D = np.linalg.norm(P[:, None] - P[None, :], axis=2); np.fill_diagonal(D, 1.0); L = np.log(D)
    best = np.inf
    for c in itertools.combinations(range(1, K + 1), k - 1):
        idx = (0,) + c; s = 0.0
        for a, b in itertools.combinations(idx, 2):
            s += L[a, b]
            if s >= best: break
        best = min(best, s)
    return np.exp(best)
import sys
KK = [int(v) for v in sys.argv[1:]] or [16, 20, 24]
for K in KK:
    print("K=%d: oblique %.5f   square %.5f   hex %.5f" % (K, Qk(0.194050, 1.248035, 7, K), Qk(0, 1, 7, K), Qk(0.5, np.sqrt(3) / 2, 7, K)), flush=True)
