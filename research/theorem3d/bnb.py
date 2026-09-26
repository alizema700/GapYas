"""Branch and bound: show Q(G) = P(G)/sqrt(det G) < 3/2 on the reduced region R minus a small box around BCC.

Region R (g11 = 1): 1 <= g22 <= g33, g22*g33 <= 64/9, 0 <= g12, g13 <= 1/2, |g23| <= g22/2, and the four
Minkowski conditions 1 + g22 + 2(e1 e2 g12 + e1 g13 + e2 g23) >= 0.
Upper bound on a box: for triangles t = (0, m, n) (plus collinear ones, included since n = 2m is allowed),
Q <= sqrt(l_m * l_n * l_{m-n}) / sqrt(det), l_v = v^T G v linear in the parameters; take box maxima of the
linear forms and a box minimum of det.  All bounds are inflated by a relative safety factor for rounding.
"""
import itertools, sys, time, json
import numpy as np

SAFE = 1 + 1e-9
TARGET = 1.5
# triangle candidates: (0, m, n) with m, n having coefficients in {-1,0,1}, plus collinear (0, m, 2m)
V1 = [m for m in itertools.product(range(-1, 2), repeat=3) if any(m)]
VECS = np.array(V1 + [tuple(2 * np.array(m)) for m in V1])
n1 = len(V1)
PAIRS = np.array([(i, j) for i in range(n1) for j in range(n1) if i < j] + [(i, n1 + i) for i in range(n1)])

def coeffs(v):   # l_v = c0 + c22 g22 + c33 g33 + c12 g12 + c13 g13 + c23 g23
    m1, m2, m3 = v
    return np.array([m1 * m1, m2 * m2, m3 * m3, 2 * m1 * m2, 2 * m1 * m3, 2 * m2 * m3], float)
CV = np.array([coeffs(v) for v in VECS])                       # (nv, 6)
CD = np.array([coeffs(VECS[i] - VECS[j]) for i, j in PAIRS])  # (np, 6)

def lin_max(C, lo, hi):
    # C: (k,6) coefficient rows; lo, hi: (b,5) -> (b,k) maxima of c0 + c.x over the boxes
    pos = np.clip(C[:, 1:], 0, None); neg = np.clip(C[:, 1:], None, 0)
    return C[None, :, 0] + hi @ pos.T + lo @ neg.T

def lin_eval(C, x):
    return C[None, :, 0] + x @ C[:, 1:].T

def det_min(lo, hi):
    g22l, g33l, g12l, g13l, g23l = lo.T; g22h, g33h, g12h, g13h, g23h = hi.T
    t1 = g22l * g33l
    t2 = np.maximum(g23l ** 2, g23h ** 2)
    t3 = g12h ** 2 * g33h
    t4 = g13h ** 2 * g22h
    t5 = np.where(g23l < 0, 2 * g12h * g13h * g23l, 2 * g12l * g13l * g23l)
    return (t1 - t2 - t3 - t4 + t5) / SAFE - 1e-12   # absolute margin for rounding (all terms are O(10))

def det_at(x):
    g22, g33, g12, g13, g23 = x.T
    return g22 * g33 - g23 ** 2 - g12 ** 2 * g33 - g13 ** 2 * g22 + 2 * g12 * g13 * g23

def infeasible(lo, hi):
    bad = hi[:, 1] < lo[:, 0]                                    # g33 >= g22
    bad |= lo[:, 0] * lo[:, 1] > 64 / 9                           # g22 g33 <= 64/9
    bad |= lo[:, 4] > hi[:, 0] / 2                                 # |g23| <= g22/2
    bad |= hi[:, 4] < -hi[:, 0] / 2
    for e1 in (1, -1):
        for e2 in (1, -1):
            C = np.array([[1, 1, 0, 2 * e1 * e2, 2 * e1, 2 * e2]], float)
            bad |= lin_max(C, lo, hi)[:, 0] < 0
    return bad

def pick_triangles(x, k=24):
    """indices of the k smallest-Q triangles at points x (b,5) -> (b,k) pair indices"""
    lv = lin_eval(CV, x)                                          # (b, nv)
    ld = lin_eval(CD, x)                                          # (b, np)
    prod = lv[:, PAIRS[:, 0]] * lv[:, PAIRS[:, 1]] * ld
    return np.argpartition(prod, k, axis=1)[:, :k]

def upper_bound(lo, hi, k=24):
    ctr = (lo + hi) / 2
    sel = pick_triangles(ctr, k)                                  # (b,k)
    vmax = lin_max(CV, lo, hi)                                    # (b,nv)
    b = np.arange(len(lo))[:, None]
    i, j = PAIRS[sel, 0], PAIRS[sel, 1]
    dmax = np.empty(sel.shape)
    # linear forms of the differences for the selected pairs
    Cd = CD[sel]                                                  # (b,k,6)
    pos = np.clip(Cd[..., 1:], 0, None); neg = np.clip(Cd[..., 1:], None, 0)
    dmax = Cd[..., 0] + np.einsum("bkj,bj->bk", pos, hi) + np.einsum("bkj,bj->bk", neg, lo)
    num = vmax[b, i] * vmax[b, j] * dmax * SAFE
    dm = det_min(lo, hi)
    ub = np.sqrt(np.min(num, axis=1)) / np.sqrt(np.where(dm > 0, dm, np.nan))
    return np.where(dm > 0, ub, np.inf)

def run(excl_center, excl_halfwidth, max_boxes=4 * 10**7, batch=50000, log=print):
    lo0 = np.array([[1.0, 1.0, 0.0, 0.0, -4 / 3]]); hi0 = np.array([[8 / 3, 64 / 9, 0.5, 0.5, 4 / 3]])
    queue = [(lo0, hi0)]; processed = 0; t0 = time.time(); maxub_left = 0
    ec = np.asarray(excl_center); eh = np.asarray(excl_halfwidth)
    while queue:
        lo, hi = queue.pop()
        if len(lo) > batch:
            queue.append((lo[batch:], hi[batch:])); lo, hi = lo[:batch], hi[:batch]
        processed += len(lo)
        keep = ~infeasible(lo, hi)
        inside = np.all((lo >= ec - eh) & (hi <= ec + eh), axis=1)   # fully inside the excluded box
        keep &= ~inside
        lo, hi = lo[keep], hi[keep]
        if len(lo) == 0: continue
        ub = upper_bound(lo, hi)
        open_ = ub >= TARGET
        lo, hi = lo[open_], hi[open_]
        if len(lo) == 0: continue
        w = hi - lo
        if np.min(np.max(w, axis=1)) < 1e-7:
            log("STALL: box too small", lo[np.argmin(np.max(w, axis=1))]); return False, processed
        ax = np.argmax(w / np.array([1.67, 6.1, 0.5, 0.5, 2.7]), axis=1)
        mid = (lo[np.arange(len(lo)), ax] + hi[np.arange(len(lo)), ax]) / 2
        lo2, hi2 = lo.copy(), hi.copy()
        hi[np.arange(len(lo)), ax] = mid; lo2[np.arange(len(lo)), ax] = mid
        queue.append((np.vstack([lo, lo2]), np.vstack([hi, hi2])))
        if processed // 200000 != (processed - len(lo)) // 200000:
            log("processed %d boxes, open %d, queue chunks %d, %.0fs" % (processed, 2 * len(lo), len(queue), time.time() - t0))
        if processed > max_boxes:
            log("ABORT: box budget exceeded"); return False, processed
    log("DONE: all boxes certified, processed %d boxes in %.0fs" % (processed, time.time() - t0))
    return True, processed

if __name__ == "__main__":
    hw = float(sys.argv[1]) if len(sys.argv) > 1 else 0.02
    ok, n = run([1, 1, 1 / 3, 1 / 3, -1 / 3], [hw] * 5)
    json.dump(dict(ok=ok, boxes=n, excluded_halfwidth=hw), open("bnb_result_hw%g.json" % hw, "w"))
