"""Rigorous computer-assisted proof (Theorem 4): the square (tau0 = i) or hexagonal (tau0 = rho) lattice is the
unique global minimiser of T_nu on the fundamental domain.  All quantities are rigorous enclosures (jets_rig.py);
all scalar comparisons are done in Arb.

(1) large y: T(tau) >= y^(3s) K_s with K_s = sum_{j != k, j,k != 0} (|j||k||j-k|)^(-2s) (finite partial sum);
(2) local convexity lemma on the box |A|,|B| <= r around the expansion point (x0, y0) (a float point within
    1e-16 of tau0):  the Hessian of T in (A, B) is bounded using the Taylor coefficients up to order 3 and the
    majorant remainder  rho2 = T * sum_{k>=4} g_k k^2 r^(k-2),  g_k = sum_{i+j=k} Ghat_ij, computed in closed form
    from g(t) = (1 - h(t))^(-3s), h(t) = (2t + 2t^2)/(1 - t):  sum_k g_k k^2 t^k = t g'(t) + t^2 g''(t).
    Positive definiteness on the box + tau0 critical (by symmetry) => tau0 is the unique minimiser in the box;
(3) branch and bound on [0, 1/2] x [0.8, Y0] minus the local box: second-order Taylor at the centre minus the
    majorant remainder of order >= 3 must exceed the upper bound of T(tau0).
"""
import json, sys, time
import numpy as np
from flint import arb
from iv import up, down
import jets_rig as jr

nu = float(sys.argv[1]); POINT = sys.argv[2]
s = arb(nu) / 2; E3 = 3 * s
X0, Y0p = (0.0, 1.0) if POINT == "square" else (0.5, float(np.sqrt(3) / 2))
R = 40
logf = open("prove_rig_%s_nu%g.log" % (POINT, nu), "w")
def say(*a):
    m = " ".join(str(x) for x in a); print(m, flush=True); logf.write(m + "\n"); logf.flush()

# ---------------- expansion at tau0 (order 3) ----------------
c, BR, ctr = jr.T_jet(X0, Y0p, nu, R=R, K=3, return_trunc=True)
T0tr_up = ctr[(0, 0)].upper()
T0 = c[(0, 0)]; T0_up = T0.upper(); T0_hi = arb(T0_up)
say("nu=%g %s: T(tau0) in [%s, %s]   tail bound %s" % (nu, POINT, T0.lower(), T0.upper(), BR.upper()))
for k in [(1, 0), (0, 1), (1, 1), (3, 0), (1, 2)]:
    assert c[k].lower() <= 0 <= c[k].upper(), ("symmetry-zero coefficient does not contain 0", k, c[k])
say("  coefficients: c20=%s c02=%s c21=%s c03=%s (symmetric ones contain 0)" % (c[(2, 0)], c[(0, 2)], c[(2, 1)], c[(0, 3)]))

def mag(a): return max(abs(a.lower()), abs(a.upper()))
def hess_ok(r):
    r = arb(r); t = r
    h = (2 * t + 2 * t * t) / (1 - t)
    if h.upper() >= 1: return False
    hp = (2 + 4 * t - 2 * t * t) / (1 - t) ** 2; hpp = 8 / (1 - t) ** 3
    g = (1 - h) ** (-E3); gp = E3 * (1 - h) ** (-E3 - 1) * hp
    gpp = E3 * (E3 + 1) * (1 - h) ** (-E3 - 2) * hp * hp + E3 * (1 - h) ** (-E3 - 1) * hpp
    G = jr.majorant_coeffs(E3, 3)
    low = sum((sum((G[i][k - i] for i in range(k + 1)), arb(0)) * k * k * t ** k for k in range(4)), arb(0))
    rho2 = T0_hi * (t * gp + t * t * gpp - low) / (t * t)
    H11 = 2 * c[(2, 0)] - 6 * mag(c[(3, 0)]) * r - 2 * mag(c[(2, 1)]) * r - rho2
    H22 = 2 * c[(0, 2)] - 2 * mag(c[(1, 2)]) * r - 6 * mag(c[(0, 3)]) * r - rho2
    H12 = mag(c[(1, 1)]) + 2 * mag(c[(2, 1)]) * r + 2 * mag(c[(1, 2)]) * r + rho2
    return H11.lower() > 0 and (H11 * H22 - H12 * H12).lower() > 0
def direct_ok(r):
    """square point only (x0 = 0, y0 = 1 exact): c10 = c01 = c11 = c30 = c12 = 0 exactly by symmetry, and
    T(i+d) - T(i) >= A^2 (c20 - |c21| r - rho) + B^2 (c02 - |c03| r - rho), rho = T sum_{i+j>=4} Ghat_ij r^(i+j-2)"""
    r = arb(r)
    try:
        G3 = jr.majorant_coeffs(E3, 3)
        low = sum((G3[i][j] * r ** (i + j) for i in range(4) for j in range(4 - i)), arb(0))
        rho = T0_hi * (jr.majorant_value(E3, r, r) - low) / (r * r)
    except AssertionError:
        return False
    return (c[(2, 0)] - arb(mag(c[(2, 1)])) * r - rho).lower() > 0 and (c[(0, 2)] - arb(mag(c[(0, 3)])) * r - rho).lower() > 0
test = direct_ok if POINT == "square" else hess_ok
rloc = 0.0
for r in np.geomspace(1e-4, 0.1, 200):
    if test(r): rloc = float(r)
    else: break
say("  local lemma (%s) holds for max(|A|,|B|) <= r_loc = %.5f" % ("direct, exact symmetry" if POINT == "square" else "convexity", rloc))
assert rloc > 0

# ---------------- large y ----------------
K = arb(0)
for j in range(-40, 41):
    for k in range(-40, 41):
        if j and k and j != k: K += (arb(abs(j) * abs(k) * abs(j - k))) ** (-2 * s)
Ycut = (T0_hi / K) ** (1 / E3)
Yc = up(float(Ycut.upper()))
say("  large-y cut: T >= y^%s K_s > T(tau0) for y > Y0 = %.6f" % (E3, Yc))

# ---------------- branch and bound ----------------
def box_certified(x0, x1, y0, y1):
    xc, yc = (x0 + x1) / 2, (y0 + y1) / 2
    hx, hy = up(max(xc - x0, x1 - xc)), up(max(yc - y0, y1 - yc))
    cj, _, cjt = jr.T_jet(xc, yc, nu, R=R, K=2, return_trunc=True)
    Ah, Bh = arb(hx) / arb(yc), arb(hy) / arb(yc)
    Ah, Bh = arb(Ah.upper()), arb(Bh.upper())
    P = cj[(0, 0)] - arb(mag(cj[(1, 0)])) * Ah - arb(mag(cj[(0, 1)])) * Bh \
        - arb(max(0.0, -cj[(2, 0)].lower())) * Ah * Ah - arb(mag(cj[(1, 1)])) * Ah * Bh - arb(max(0.0, -cj[(0, 2)].lower())) * Bh * Bh
    try:
        rem = arb(cj[(0, 0)].upper()) * jr.majorant_tail(E3, Ah, Bh, 2)
    except AssertionError:
        return False
    if (P - rem).lower() > T0_up: return True
    # near tau0: compare cube-truncated sums; the tails differ by at most B_R(tau0) (G(A0,B0) - 1), where
    # (A0, B0) bound the scaled distance of the box to the expansion point (term-wise majorant, t/t0 >= 2 - G)
    Pt = cjt[(0, 0)] - arb(mag(cjt[(1, 0)])) * Ah - arb(mag(cjt[(0, 1)])) * Bh \
        - arb(max(0.0, -cjt[(2, 0)].lower())) * Ah * Ah - arb(mag(cjt[(1, 1)])) * Ah * Bh - arb(max(0.0, -cjt[(0, 2)].lower())) * Bh * Bh
    remt = arb(cjt[(0, 0)].upper()) * jr.majorant_tail(E3, Ah, Bh, 2)
    A0 = arb(up(max(abs(x0 - X0), abs(x1 - X0)))) / arb(Y0p); B0 = arb(up(max(abs(y0 - Y0p), abs(y1 - Y0p)))) / arb(Y0p)
    try:
        dtail = BR * (jr.majorant_value(E3, A0, B0) - 1)
    except AssertionError:
        return False
    return (Pt - remt - dtail).lower() > T0tr_up

hl = 0.99 * rloc * Y0p           # excluded box strictly inside the local box
queue = [(0.0, 0.5, 0.8, Yc)]; n = 0; t0 = time.time()
while queue:
    x0, x1, y0, y1 = queue.pop(); n += 1
    if up(up(x1 * x1) + up(y1 * y1)) < 1.0: continue                  # entirely inside the unit disc
    if X0 - hl <= x0 and x1 <= X0 + hl and Y0p - hl <= y0 and y1 <= Y0p + hl: continue
    if box_certified(x0, x1, y0, y1): continue
    if max(x1 - x0, y1 - y0) < 1e-7: say("STALL", (x0, x1, y0, y1)); sys.exit(1)
    xm, ym = (x0 + x1) / 2, (y0 + y1) / 2
    if (x1 - x0) >= (y1 - y0): queue += [(x0, xm, y0, y1), (xm, x1, y0, y1)]
    else: queue += [(x0, x1, y0, ym), (x0, x1, ym, y1)]
    if n % 100 == 0: say("  boxes %d, queue %d, %.0fs" % (n, len(queue), time.time() - t0))
say("DONE nu=%g: all %d boxes certified in %.0fs (rigorous arithmetic) -> %s lattice is the unique global minimiser" % (nu, n, time.time() - t0, POINT))
json.dump(dict(point=POINT, nu=nu, T0=[str(T0.lower()), str(T0.upper())], r_loc=rloc, Y0=Yc, boxes=n), open("prove_rig_%s_nu%g.json" % (POINT, nu), "w"), indent=1)
