"""Numerical sanity checks for the nu -> infinity theorem in 2D.

P(L) = min over triples of distinct lattice points of the product of the three mutual distances.
Claims to check:
  (i)   P(L) <= G(a) = min(2a^3, sqrt(a^2 + a^-2)) for every unit-covolume lattice (a = shortest vector)
  (ii)  max_L P(L) = P* = 2((1+sqrt17)/8)^(3/4), attained only by the rectangle with aspect sqrt((sqrt17-1)/2)
  (iii) T_nu(L)^(1/nu) -> 1/P(L) and minimisers of T_nu approach the rectangle
"""
import numpy as np, itertools, sys
sys.path.insert(0, "../verification")
from lattice_sums import A_tau

def P_of(A, R=5):
    ax = np.arange(-R, R + 1)
    m = np.stack(np.meshgrid(ax, ax, indexing="ij"), -1).reshape(-1, 2)
    m = m[np.any(m != 0, axis=1)]
    X = m @ A.T
    r = np.linalg.norm(X, axis=1)
    keep = r < 4.0
    X, r = X[keep], r[keep]
    D = np.linalg.norm(X[:, None] - X[None], axis=-1)
    np.fill_diagonal(D, np.inf)
    return float(np.min(np.outer(r, r) * D))

a_star = ((1 + np.sqrt(17)) / 8) ** 0.25
P_star = 2 * a_star ** 3
y_inf = np.sqrt((np.sqrt(17) - 1) / 2)
print("a* = %.12f  P* = %.12f  y_inf = %.12f  1/a*^2 = %.12f" % (a_star, P_star, y_inf, 1 / a_star**2))
print("P(rectangle y_inf) =", P_of(A_tau(1j * y_inf)))

# (i) random lattices in the fundamental domain
rng = np.random.default_rng(0); worst = -1
for _ in range(20000):
    x = rng.uniform(0, 0.5); y = rng.uniform(np.sqrt(1 - x * x), 3.0)
    A = A_tau(x + 1j * y); a = np.linalg.norm(A[:, 0])
    G = min(2 * a**3, np.sqrt(a**2 + a**-2))
    P = P_of(A)
    worst = max(worst, P - G)
    assert P <= G + 1e-12, (x, y, P, G)
print("(i) max P - G over 20000 random lattices: %.2e (must be <= 0)" % worst)

# (ii) grid maximisation of P over the fundamental domain
best = (0, None)
for x in np.linspace(0, 0.5, 201):
    for y in np.linspace(np.sqrt(1 - x * x), 2.0, 401):
        P = P_of(A_tau(x + 1j * y), R=4)
        if P > best[0]: best = (P, (x, y))
print("(ii) grid max P = %.10f at tau = %.4f + %.4fi  (P* = %.10f, y_inf = %.4f)" % (best[0], *best[1], P_star, y_inf))
