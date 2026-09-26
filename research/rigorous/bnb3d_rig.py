"""Rigorous branch and bound for Theorem 3 (BCC maximises P).  Same algorithm as theorem3d/bnb.py, but every
bound is computed with directed (outward) rounding, see iv.py:
  upper(l_v) on a box   : c0 + sum_k c_k x_k with x_k = hi or lo by the sign of the (small integer) c_k,
                          every product and partial sum rounded up;
  lower(det G) on a box : each monomial bounded in the safe direction and rounded accordingly;
  bound Q <= sqrt(up(l_m l_n l_d)) / sqrt(down(det)), rounded up, compared with 3/2 exactly.
Pruning by region constraints is done only when a constraint is violated on the whole box (also with
directed rounding).  The candidate triangles are chosen heuristically; this affects speed, not validity."""
import itertools, json, sys, time
import numpy as np
from iv import up, down

TARGET = 1.5
V1 = [m for m in itertools.product(range(-1, 2), repeat=3) if any(m)]
VECS = np.array(V1 + [tuple(2 * np.array(m)) for m in V1])
n1 = len(V1)
PAIRS = np.array([(i, j) for i in range(n1) for j in range(n1) if i < j] + [(i, n1 + i) for i in range(n1)])
def coeffs(v):
    m1, m2, m3 = v
    return np.array([m1 * m1, m2 * m2, m3 * m3, 2 * m1 * m2, 2 * m1 * m3, 2 * m2 * m3], float)
CV = np.array([coeffs(v) for v in VECS]); CD = np.array([coeffs(VECS[i] - VECS[j]) for i, j in PAIRS])

def lin_up(C, lo, hi):
    """rigorous upper bound of c0 + sum_k c_k x_k over boxes; C (k,6) small integers (exact), lo/hi (b,5)"""
    acc = np.broadcast_to(C[None, :, 0], (lo.shape[0], C.shape[0])).astype(float).copy()
    for k in range(5):
        ck = C[None, :, k + 1]
        x = np.where(ck > 0, hi[:, k:k + 1], lo[:, k:k + 1])
        acc = up(acc + up(ck * x))
    return acc
def lin_up_sel(Csel, lo, hi):
    """same for per-box coefficient rows Csel (b,k,6)"""
    acc = Csel[..., 0].copy()
    for k in range(5):
        ck = Csel[..., k + 1]
        x = np.where(ck > 0, hi[:, k:k + 1], lo[:, k:k + 1])
        acc = up(acc + up(ck * x))
    return acc

def det_down(lo, hi):
    g22l, g33l, g12l, g13l, g23l = lo.T; g22h, g33h, g12h, g13h, g23h = hi.T
    t1 = down(g22l * g33l)                                       # g22, g33 >= 1 > 0
    t2 = up(np.maximum(up(g23l * g23l), up(g23h * g23h)))
    t3 = up(up(g12h * g12h) * g33h)                              # g12, g13 >= 0 on the region
    t4 = up(up(g13h * g13h) * g22h)
    t5 = np.where(g23l < 0, -up(up(up(2 * g12h) * g13h) * (-g23l)), down(down(down(2 * g12l) * g13l) * g23l))
    return down(down(down(down(t1 - t2) - t3) - t4) + t5)

def infeasible(lo, hi):
    bad = hi[:, 1] < lo[:, 0]
    bad |= down(lo[:, 0] * lo[:, 1]) > up(64.0 / 9.0)
    bad |= lo[:, 4] > up(hi[:, 0] / 2)
    bad |= hi[:, 4] < -up(hi[:, 0] / 2)
    for e1 in (1, -1):
        for e2 in (1, -1):
            C = np.array([[1, 1, 0, 2 * e1 * e2, 2 * e1, 2 * e2]], float)
            bad |= lin_up(C, lo, hi)[:, 0] < 0
    return bad

def pick(x, k=24):
    lv = CV[None, :, 0] + x @ CV[:, 1:].T; ld = CD[None, :, 0] + x @ CD[:, 1:].T
    prod = lv[:, PAIRS[:, 0]] * lv[:, PAIRS[:, 1]] * ld
    return np.argpartition(prod, k, axis=1)[:, :k]

def upper_bound(lo, hi, k=24):
    sel = pick((lo + hi) / 2, k)
    vmax = lin_up(CV, lo, hi)
    b = np.arange(len(lo))[:, None]
    dmax = lin_up_sel(CD[sel], lo, hi)
    num = up(up(vmax[b, PAIRS[sel, 0]] * vmax[b, PAIRS[sel, 1]]) * dmax)
    num = np.min(num, axis=1)
    dm = det_down(lo, hi)
    with np.errstate(invalid="ignore", divide="ignore"):
        ub = up(up(np.sqrt(num)) / down(np.sqrt(np.where(dm > 0, dm, np.nan))))
    return np.where(dm > 0, ub, np.inf)

def run(center, hw, batch=50000, log=print):
    lo0 = np.array([[1.0, 1.0, 0.0, 0.0, -up(4 / 3)]]); hi0 = np.array([[up(8 / 3), up(64 / 9), 0.5, 0.5, up(4 / 3)]])
    queue = [(lo0, hi0)]; n = 0; t0 = time.time(); c = np.asarray(center); h = np.asarray(hw)
    elo, ehi = down(c - h), up(c + h)
    while queue:
        lo, hi = queue.pop()
        if len(lo) > batch: queue.append((lo[batch:], hi[batch:])); lo, hi = lo[:batch], hi[:batch]
        n += len(lo)
        inside = np.all((lo >= elo) & (hi <= ehi), axis=1)
        keep = (~infeasible(lo, hi)) & (~inside)
        lo, hi = lo[keep], hi[keep]
        if not len(lo): continue
        ub = upper_bound(lo, hi)
        op = ~(ub < TARGET)
        lo, hi = lo[op], hi[op]
        if not len(lo): continue
        w = hi - lo
        if np.min(np.max(w, axis=1)) < 1e-9: log("STALL", lo[0], hi[0]); return False, n
        ax = np.argmax(w / np.array([1.67, 6.1, 0.5, 0.5, 2.7]), axis=1); r = np.arange(len(lo))
        mid = (lo[r, ax] + hi[r, ax]) / 2
        lo2, hi2 = lo.copy(), hi.copy(); hi[r, ax] = mid; lo2[r, ax] = mid
        queue.append((np.vstack([lo, lo2]), np.vstack([hi, hi2])))
        if n // 500000 != (n - len(lo)) // 500000: log("boxes %d, open %d, %.0fs" % (n, 2 * len(lo), time.time() - t0))
    log("DONE: all boxes certified (rigorous rounding), %d boxes, %.0fs" % (n, time.time() - t0))
    return True, n

if __name__ == "__main__":
    hw = float(sys.argv[1]) if len(sys.argv) > 1 else 0.007
    ok, n = run([1, 1, 1 / 3, 1 / 3, -1 / 3], [hw] * 5)
    json.dump(dict(ok=ok, boxes=n, excluded_halfwidth=hw), open("bnb3d_rig_result.json", "w"))
