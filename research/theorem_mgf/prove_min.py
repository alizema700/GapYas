"""Computer-assisted proof: for nu = 2s in {4, 6, 8} the square lattice (tau = i) is the unique global
minimiser of T_nu on the fundamental domain, equivalently C_{s,s,s}(tau) >= C_{s,s,s}(i) with equality only at i.

Ingredients (see jets.py / tjet.py for the rigorous enclosures):
 (1) large y:  T(tau) >= y^(3s) K_s  (collinear triples along the shortest vector 1), K_s from a finite sum.
 (2) local lemma at i: by the symmetries x -> -x and tau -> 1/conj(tau) the Taylor coefficients
     c10 = c01 = c11 = c30 = c12 = 0; with rigorous c20, c02, c21, c03 and the majorant remainder
     T(i + d) - T(i) >= A^2 (c20 - |c21| r - rho) + B^2 (c02 - |c03| r - rho) > 0 for 0 < max(|A|,|B|) <= r,
     rho = T(i) * sum_{i+j>=4} Ghat_ij r^(i+j-2).
 (3) branch and bound on [0, 1/2] x [y_lo, Y0] minus the local box: second-order Taylor at the box centre with
     rigorous coefficient enclosures and majorant remainder; a box is certified if its lower bound exceeds
     the upper bound of T(i).
"""
import json, sys, time
import numpy as np
from tjet import T_jet
from jets import majorant_power, majorant_eps, majorant_tail

nu = float(sys.argv[1]) if len(sys.argv) > 1 else 4.0
POINT = sys.argv[2] if len(sys.argv) > 2 else "square"
PX, PY = (0.0, 1.0) if POINT == "square" else (0.5, np.sqrt(3) / 2)
s = nu / 2; R = 40
log = open("prove_%s_nu%g.log" % (POINT, nu), "w")
def say(*a):
    msg = " ".join(str(x) for x in a); print(msg, flush=True); log.write(msg + "\n"); log.flush()

# ---- T(i) and local lemma -------------------------------------------------------------------
c3, e3, BR = T_jet(PX, PY, nu, R=R, K=3)
Ti_up = c3[0, 0] + e3[0, 0]; Ti_lo = c3[0, 0] - e3[0, 0]
say("nu=%g  point %s: T in [%.12f, %.12f]" % (nu, POINT, Ti_lo, Ti_up))
for (i, j) in [(1, 0), (0, 1), (1, 1), (3, 0), (1, 2)]:
    say("  symmetric-zero coefficient c%d%d enclosure: %.3e +- %.3e (set to 0 by symmetry)" % (i, j, c3[i, j], e3[i, j]))
c20 = c3[2, 0] - e3[2, 0]; c02 = c3[0, 2] - e3[0, 2]
c21 = abs(c3[2, 1]) + e3[2, 1]; c03 = abs(c3[0, 3]) + e3[0, 3]
say("  c20 >= %.6f  c02 >= %.6f  |c21| <= %.4f  |c03| <= %.4f" % (c20, c02, c21, c03))
G4 = majorant_power(majorant_eps(12), 3 * s, 12)
def rho(r):  # T(i) * sum_{i+j>=4} Ghat_ij r^(i+j-2), via closed form minus the low orders
    low = sum(G4[i, j] * r ** (i + j) for i in range(4) for j in range(4 - i))
    from jets import majorant_value
    return Ti_up * (majorant_value(3 * s, r, r) - low) / r ** 2
rloc = 0.0
for r in np.geomspace(1e-5, 0.2, 400):
    try:
        if c20 - c21 * r - rho(r) > 0 and c02 - c03 * r - rho(r) > 0: rloc = r
        else: break
    except AssertionError:
        break
say("  local lemma holds for max(|A|,|B|) <= r_loc = %.5f" % rloc)
if rloc == 0: sys.exit("local lemma failed")

# ---- large y ------------------------------------------------------------------------------
K = 0.0
for j in range(-60, 61):
    for k in range(-60, 61):
        if j and k and j != k: K += (abs(j) * abs(k) * abs(j - k)) ** (-nu)
Y0 = (Ti_up / K) ** (1 / (3 * s)) * (1 + 1e-9)
say("  large-y cut: T(tau) >= y^%g * %.6f > T(i) for y > Y0 = %.5f" % (3 * s, K, Y0))

# ---- branch and bound -----------------------------------------------------------------------
Gc = majorant_power(majorant_eps(12), 3 * s, 12)
def box_lower(xc, yc, hx, hy):
    c, e, _ = T_jet(xc, yc, nu, R=R, K=2)
    Ah, Bh = hx / yc, hy / yc
    lo = c - np.sign(c) * e                           # move each coefficient towards zero... (conservative below)
    P = (c[0, 0] - e[0, 0]) - (abs(c[1, 0]) + e[1, 0]) * Ah - (abs(c[0, 1]) + e[0, 1]) * Bh \
        - max(0.0, -(c[2, 0] - e[2, 0])) * Ah**2 - (abs(c[1, 1]) + e[1, 1]) * Ah * Bh - max(0.0, -(c[0, 2] - e[0, 2])) * Bh**2
    try:
        rem = (c[0, 0] + e[0, 0]) * majorant_tail(3 * s, Ah, Bh, 2, Gc[:3, :3] if False else None)
    except AssertionError:
        return -np.inf
    return P - rem

x_lo, x_hi, y_lo = 0.0, 0.5, 0.8
queue = [(x_lo, x_hi, y_lo, Y0)]
n_boxes = 0; t0 = time.time(); rl = rloc   # local box: |x| <= rl, |y - 1| <= rl (A = x, B = y - 1 at y0 = 1)
while queue:
    x0, x1, y0, y1 = queue.pop()
    n_boxes += 1
    if x1 ** 2 + y1 ** 2 < 1.0: continue                             # entirely below the unit circle
    if PX - rl * PY <= x0 and x1 <= PX + rl * PY and PY - rl * PY <= y0 and y1 <= PY + rl * PY: continue  # local-lemma box
    xc, yc, hx, hy = (x0 + x1) / 2, (y0 + y1) / 2, (x1 - x0) / 2, (y1 - y0) / 2
    lb = box_lower(xc, yc, hx, hy)
    if lb > Ti_up: continue
    if max(hx, hy) < 1e-7: say("STALL at box", (x0, x1, y0, y1), "lb", lb); sys.exit(1)
    if hx / yc >= hy / yc: queue += [(x0, xc, y0, y1), (xc, x1, y0, y1)]
    else: queue += [(x0, x1, y0, yc), (x0, x1, yc, y1)]
    if n_boxes % 200 == 0: say("  boxes %d, queue %d, %.0fs" % (n_boxes, len(queue), time.time() - t0))
say("DONE nu=%g: all %d boxes certified in %.0fs -> %s lattice is the unique global minimiser" % (nu, n_boxes, time.time() - t0, POINT))
json.dump(dict(point=POINT, nu=nu, T_i=[Ti_lo, Ti_up], r_loc=rloc, Y0=Y0, boxes=n_boxes), open("prove_%s_nu%g.json" % (POINT, nu), "w"), indent=1)
