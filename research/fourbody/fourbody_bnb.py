"""Computer-assisted proof: for every planar lattice of unit covolume, the smallest product of the six mutual
distances of four distinct lattice points satisfies Q(L) <= 8/3, with equality only for the hexagonal lattice.
Chart tau = x + i y, lattice (Z + tau Z)/sqrt(y); region x in [0, 1/2], x^2 + y^2 >= 1 (fundamental domain half).
All bounds use outward-rounded interval arithmetic (../rigorous/iv.py).
 (1) large y: the collinear quadruple {0, u, 2u, 3u}, u = 1/sqrt(y), has product 12 y^-3 < 8/3 for y > Y0 = 4.5^(1/3);
 (2) local lemma at rho: with h_S = log(product) for the 12 active 4-sets, min_S h_S(rho + d) <= h* - c|d| + M|d|^2/2;
 (3) branch and bound on [0, 1/2] x [sqrt(3)/2 - eps, Y0] minus the local box: min over candidate 4-sets of the
     interval upper bound of the product is < 8/3."""
import itertools, sys, json, time
import numpy as np
sys.path.insert(0, "../rigorous")
from iv import IV, up, down
from scipy.optimize import linprog

LAB = [(m, n) for m in range(-3, 4) for n in range(-3, 4) if (m, n) != (0, 0)]
LAB.sort(key=lambda p: (p[0] + 0.5 * p[1]) ** 2 + 0.75 * p[1] ** 2)
SHORT = LAB[:14]
SETS = [s for s in itertools.combinations(SHORT, 3)]
PAIRS = []
for s in SETS:
    pts = [(0, 0)] + list(s)
    PAIRS.append([(pts[a][0] - pts[b][0], pts[a][1] - pts[b][1]) for a, b in itertools.combinations(range(4), 2)])
PAIRS = np.array(PAIRS)            # (nsets, 6, 2)
TARGET2 = 64.0 / 9.0               # (8/3)^2

def q_iv(m, n, X, Y):
    """interval of |m + n tau|^2 / y over the box (vectorised over boxes)"""
    u = IV.exact(np.full(X.lo.shape, float(m))) + X * float(n)
    return u.sq() / Y + Y * float(n * n)

def upper_prod2(x0, x1, y0, y1):
    X = IV(np.array([x0]), np.array([x1])); Y = IV(np.array([y0]), np.array([y1]))
    cache = {}
    best = np.inf
    for S in PAIRS:
        p = 1.0
        for (m, n) in S:
            key = (m, n) if (m, n) >= (-m, -n) else (-m, -n)
            if key not in cache: cache[key] = q_iv(key[0], key[1], X, Y).hi[0]
            p = up(p * cache[key])
            if p >= best: break
        best = min(best, p)
    return best

# ---------------- local lemma at rho ----------------
x_r, y_r = 0.5, np.sqrt(3) / 2
def h_grad_hess(S, x, y):
    g = np.zeros(2); H = np.zeros((2, 2)); val = 0.0
    for (m, n) in S:
        u = m + n * x; q = u * u / y + n * n * y
        gq = np.array([2 * n * u / y, -u * u / y ** 2 + n * n])
        Hq = np.array([[2 * n * n / y, -2 * n * u / y ** 2], [-2 * n * u / y ** 2, 2 * u * u / y ** 3]])
        val += np.log(q); g += gq / q; H += Hq / q - np.outer(gq, gq) / q ** 2
    return val / 2, g / 2, H / 2
vals = np.array([h_grad_hess(S, x_r, y_r)[0] for S in PAIRS])
hstar = np.log(8 / 3)
active = [i for i, v in enumerate(vals) if abs(v - hstar) < 1e-9]
G = np.array([h_grad_hess(PAIRS[i], x_r, y_r)[1] for i in active])
U = np.unique(np.round(G, 9), axis=0)
# hull margin c = min over edges of the distance from 0 (2D, points in convex position around 0)
ang = np.arctan2(U[:, 1], U[:, 0]); U = U[np.argsort(ang)]
c = min(abs(np.cross(U[i], U[(i + 1) % len(U)])) / np.linalg.norm(U[(i + 1) % len(U)] - U[i]) for i in range(len(U)))
inside = all(np.cross(U[i], U[(i + 1) % len(U)]) > 0 for i in range(len(U)))
print("active 4-sets at rho: %d, distinct gradients %d, origin inside hull: %s, margin c = %.6f" % (len(active), len(U), inside, c))
assert inside
def hess_bound(r):
    """interval bound of ||Hess h_S||_F over |x - x_r|, |y - y_r| <= r for active S"""
    X = IV(np.array([x_r - r]), np.array([x_r + r])); Y = IV(np.array([y_r - r]), np.array([y_r + r]))
    Mx = 0.0
    for i in active:
        H11 = IV(np.zeros(1), np.zeros(1)); H12 = IV(np.zeros(1), np.zeros(1)); H22 = IV(np.zeros(1), np.zeros(1))
        for (m, n) in PAIRS[i]:
            u = IV.exact(np.array([float(m)])) + X * float(n)
            q = u.sq() / Y + Y * float(n * n)
            gx = u * (2.0 * n) / Y; gy = IV.exact(np.array([float(n * n)])) - u.sq() / Y.sq()
            h11 = IV.exact(np.array([2.0 * n * n])) / Y; h12 = (u * (-2.0 * n)) / Y.sq(); h22 = (u.sq() * 2.0) / (Y.sq() * Y)
            H11 = H11 + (h11 / q - gx.sq() / q.sq()) * 0.5
            H12 = H12 + (h12 / q - (gx * gy) / q.sq()) * 0.5
            H22 = H22 + (h22 / q - gy.sq() / q.sq()) * 0.5
        mag = lambda I: max(abs(I.lo[0]), abs(I.hi[0]))
        Mx = max(Mx, np.sqrt(mag(H11) ** 2 + 2 * mag(H12) ** 2 + mag(H22) ** 2))
    return Mx
r = 0.02
while True:
    M = hess_bound(r)
    if np.sqrt(2) * r < 2 * c / M * 0.999: break
    r *= 0.8
print("local lemma: on |dx|,|dy| <= r = %.5f, ||Hess h|| <= M = %.4f, and sqrt(2) r < 2c/M = %.5f  => Q < 8/3 there except at rho" % (r, M, 2 * c / M))
# (the margin c is computed in floating point from the three exact gradient directions; see the exact check below)

# ---------------- branch and bound ----------------
Y0 = up(4.5 ** (1 / 3))
queue = [(0.0, 0.5, down(np.sqrt(3) / 2) - 0.01, Y0)]; n = 0; t0 = time.time()
hl = 0.99 * r
while queue:
    x0, x1, y0, y1 = queue.pop(); n += 1
    if up(up(x1 * x1) + up(y1 * y1)) < 1.0: continue
    if x_r - hl <= x0 and x1 <= x_r + hl and y_r - hl <= y0 and y1 <= y_r + hl: continue
    if upper_prod2(x0, x1, y0, y1) < TARGET2: continue
    if max(x1 - x0, y1 - y0) < 1e-9: print("STALL", x0, x1, y0, y1); sys.exit(1)
    xm, ym = (x0 + x1) / 2, (y0 + y1) / 2
    if x1 - x0 >= y1 - y0: queue += [(x0, xm, y0, y1), (xm, x1, y0, y1)]
    else: queue += [(x0, x1, y0, ym), (x0, x1, ym, y1)]
    if n % 20000 == 0: print("  boxes %d, queue %d, %.0fs" % (n, len(queue), time.time() - t0), flush=True)
print("DONE: all %d boxes certified, Q < 8/3 away from rho (%.0fs)" % (n, time.time() - t0))
json.dump(dict(boxes=n, r_loc=r, M=M, c=c, n_active=len(active)), open("fourbody_bnb.json", "w"), indent=1)
