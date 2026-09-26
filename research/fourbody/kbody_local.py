"""Rigorous localisation of a non-symmetric k-body max-min lattice (used for k = 7).
Given a floating-point candidate tau0 = x0 + i y0:
 (1) Q(tau0) >= Q_low, by exhaustive pruned enumeration of all k-point sets containing 0, with outward rounding;
 (2) branch and bound: Q(tau) < Q_low for all tau in {0 <= x <= 1/2, |tau| >= 1, y <= Y0} outside the box
     |x - x0|, |y - y0| <= w; for y > Y0 the collinear k-set gives Q < Q_low;
 hence every maximiser of Q lies in that box, which meets none of the symmetry lines x = 0, x = 1/2, |tau| = 1:
 the max-min lattice is oblique.  Usage: kbody_local.py k x0 y0 w"""
import itertools, sys, json, time, math
import numpy as np
sys.path.insert(0, "../rigorous")
from iv import IV, up, down
k = int(sys.argv[1]); x0c, y0c, w = float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])
npair = k * (k - 1) // 2

def q_scalar_iv(m, n, x, y):
    X = IV(np.array([x]), np.array([x])); Y = IV(np.array([y]), np.array([y]))
    u = IV.exact(np.array([float(m)])) + X * float(n); v = u.sq() / Y + Y * float(n * n)
    return v.lo[0], v.hi[0]

# ---------- (1) lower bound of Q(tau0)^2 by exhaustive enumeration ----------
R = 12
QL = {p: q_scalar_iv(p[0], p[1], x0c, y0c) for p in itertools.product(range(-2 * R, 2 * R + 1), repeat=2) if p != (0, 0)}
pts_all = sorted([(m, n) for m in range(-R, R + 1) for n in range(-R, R + 1) if (m, n) != (0, 0)], key=lambda p: QL[p][0])
lam2 = QL[pts_all[0]][0]
def qmin_inf_sphere(x, y):
    """lower bound of min of q(m, n) = (m + n x)^2 / y + n^2 y over max(|m|, |n|) = 1 (the four edges; exact minimisers)"""
    cands = [1 / y + 0.0, x * x / y + y]
    for s in (-1, 1):
        n = -s * x / (x * x + y * y)                       # edge m = s: minimise (s + n x)^2 / y + n^2 y over |n| <= 1
        for nn in (max(min(n, 1), -1), -1.0, 1.0): cands.append(((s + nn * x) ** 2) / y + nn * nn * y)
        for mm in (max(min(-s * x, 1), -1), -1.0, 1.0): cands.append(((mm + s * x) ** 2) / y + y)   # edge n = s
    return down(min(cands) * (1 - 1e-12))
edge = down((R + 1) ** 2 * qmin_inf_sphere(x0c, y0c))   # every lattice point outside the R-box has q >= edge (homogeneity)   # every point outside the R-box has q >= edge
                                                                          # (edge taken on the ring |.|_inf = R+1; q is
                                                                          # convex and the ring separates the box)
def enum_min(j, prev_low):
    """[lo, hi] of min over j-point sets {0, p_1..p_{j-1}} of the squared product.  Completeness: a set containing a
    point p has squared product >= q(p) * lam2^(j-2) * prev_low (drop p: remaining j-1 points form a set whose
    squared product is >= prev_low after translation; the j-2 other distances from p are >= lam2)."""
    best_hi = [math.inf]; best_lo = [math.inf]; cnt = [0]
    def rec(chosen, start, plo, phi):
        if len(chosen) == j - 1:
            cnt[0] += 1; best_hi[0] = min(best_hi[0], phi); best_lo[0] = min(best_lo[0], plo); return
        rem = j - 1 - len(chosen)
        for i in range(start, len(pts_all)):
            p = pts_all[i]; qp = QL[p][0]
            if qp * lam2 ** (j - 2) * prev_low > best_hi[0] * (1 + 1e-9): break      # no set containing p (or later) can win
            flo, fhi = QL[p]
            for c in chosen:
                d = QL[(p[0] - c[0], p[1] - c[1])]; flo = down(flo * d[0]); fhi = up(fhi * d[1])
            if down(plo * flo) <= best_hi[0] * (1 + 1e-9): rec(chosen + [p], i + 1, down(plo * flo), up(phi * fhi))
    rec([], 0, 1.0, 1.0)
    assert edge * lam2 ** (j - 2) * prev_low > best_hi[0] * 1.001, ("enumeration radius too small", j)
    return best_lo[0], best_hi[0], cnt[0]
prev = 1.0                      # j = 1: empty product
for j in range(2, k + 1):
    lo, hi, cnt = enum_min(j, prev)
    print("  j=%d: min squared product in [%.12g, %.12g]  (%d sets)" % (j, lo, hi, cnt), flush=True)
    prev = lo
Q2_low, Q2_hi = lo, hi
best_hi = [Q2_hi]
print("k=%d, tau0 = %.12f + %.12f i: Q(tau0)^2 in [%.12f, %.12f]" % (k, x0c, y0c, Q2_low, Q2_hi), flush=True)

# ---------- (2) branch and bound ----------
LAB = [(m, n) for m in range(-3, 4) for n in range(-3, 4) if (m, n) != (0, 0)]
LAB.sort(key=lambda p: (p[0] + x0c * p[1]) ** 2 / y0c + p[1] ** 2 * y0c)
SHORT = LAB[:13]
PAIRS = []
for s in itertools.combinations(SHORT, k - 1):
    P = [(0, 0)] + list(s)
    PAIRS.append([(P[a][0] - P[b][0], P[a][1] - P[b][1]) for a, b in itertools.combinations(range(k), 2)])
def upper_prod2(x0, x1, y0, y1):
    X = IV(np.array([x0]), np.array([x1])); Y = IV(np.array([y0]), np.array([y1])); cache = {}; best = math.inf
    for S in PAIRS:
        p = 1.0
        for (m, n) in S:
            key = (m, n) if (m, n) >= (-m, -n) else (-m, -n)
            if key not in cache:
                u = IV.exact(np.array([float(key[0])])) + X * float(key[1]); cache[key] = (u.sq() / Y + Y * float(key[1] ** 2)).hi[0]
            p = up(p * cache[key])
            if p >= best: break
        best = min(best, p)
    return best
COLL = math.prod(j - i for i, j in itertools.combinations(range(k), 2))
Y0 = up((COLL ** 2 / Q2_low) ** (2 / npair))       # COLL^2 y^(-npair) < Q2_low for y > Y0
queue = [(0.0, 0.5, down(math.sqrt(3) / 2) - 0.01, Y0)]; n = 0; t0 = time.time()
while queue:
    x0, x1, y0, y1 = queue.pop(); n += 1
    if up(up(x1 * x1) + up(y1 * y1)) < 1.0: continue
    if x0c - w <= x0 and x1 <= x0c + w and y0c - w <= y0 and y1 <= y0c + w: continue
    if upper_prod2(x0, x1, y0, y1) < Q2_low: continue
    if max(x1 - x0, y1 - y0) < 1e-10: print("STALL", x0, x1, y0, y1); sys.exit(1)
    xm, ym = (x0 + x1) / 2, (y0 + y1) / 2
    if x1 - x0 >= y1 - y0: queue += [(x0, xm, y0, y1), (xm, x1, y0, y1)]
    else: queue += [(x0, x1, y0, ym), (x0, x1, ym, y1)]
    if n % 5000 == 0: print("  boxes %d, queue %d, %.0fs" % (n, len(queue), time.time() - t0), flush=True)
sym = (x0c - w > 0) and (x0c + w < 0.5) and ((x0c - w) ** 2 + (y0c - w) ** 2 > 1)
print("DONE k=%d: all %d boxes certified (%.0fs). Every maximiser of Q lies in |x-%.9f|,|y-%.9f| <= %g; box avoids all symmetry lines: %s"
      % (k, n, time.time() - t0, x0c, y0c, w, sym))
json.dump(dict(k=k, x0=x0c, y0=y0c, w=w, Q2=[Q2_low, best_hi[0]], boxes=n, oblique=sym), open("kbody_local_k%d.json" % k, "w"), indent=1)
