"""Rigorous proof that the square lattice (tau = i) is the unique minimiser of T_nu for ALL nu in [NU_A, NU_B].

Structure (all quantities are rigorous enclosures from jets_rig.py; nu enters as an Arb ball where needed):
 (L) local region  |x| <= r1, |y - 1| <= r1  (A = x, B = y - 1 at tau0 = i):
     for each nu-slice N (adaptively bisected) the Taylor coefficients c_ij (i + j <= 12) of T at i are enclosed for all
     nu in N at once; c_ij = 0 for i odd (reflection x -> -x) and c_01 = 0 (criticality).  Then for nu in N
        T(i + d) - T(i) >= sum_{2 <= i+j <= 12} c_ij A^i B^j - T(i) * sum_{i+j > 12} Ghat_ij |A|^i |B|^j
     (T(i) cancels exactly).  Inner square max(|A|,|B|) <= r_loc: the direct form of the local lemma (as in
     prove_min_rig.py); ring r_loc <= max <= r1: interval evaluation of the polynomial bound on sub-boxes.
 (O) outer region: branch and bound over boxes (x, y, nu) of the half fundamental domain minus the unit disc and minus
     the local region, with nu in [NU_A, NU_B].  With nu0 = mid(N), eta = rad(N) and D(nu) = T(tau; nu) - T(i; nu):
        D(nu) >= D(nu0) - eta * sup_{N} (|d_nu T(tau)| + |d_nu T(i)|),   |d_nu T| <= (T_{nu-1} + T_{nu+1}) / e,
     because |log p| p^-nu <= p^-(nu-1)/e for p >= 1 and <= p^-(nu+1)/e for p <= 1.  D(nu0) is bounded below by the
     second-order Taylor bound at the box centre (as in prove_min_rig.py) at the single exponent nu0.
 (Y) large y: T(tau) >= y^(3s) K_s > T(i) for y > Y0, with s ranging over N/2.
Usage: prove_nuint_rig.py NU_A NU_B"""
import json, math, sys, time
import numpy as np
from flint import arb
from iv import up, down
import jets_rig as jr

NU_A, NU_B = float(sys.argv[1]), float(sys.argv[2])
DEBUG = False
R = 20; R1 = 0.06; KL = 12
tag = "square_nu%g-%g" % (NU_A, NU_B)
logf = open("prove_nuint_%s.log" % tag, "w")
def say(*a):
    m = " ".join(str(x) for x in a); print(m, flush=True); logf.write(m + "\n"); logf.flush()
def ball(a, b): return arb((a + b) / 2, (b - a) / 2) + arb(0, up((b - a) * 1e-15))      # contains [a, b]
def mag(a): return max(abs(float(a.lower())), abs(float(a.upper())))

# ---------------------------------------------------------------- (L) local region, per nu-slice
def local_ok(a, b):
    nu = ball(a, b)
    c, BR = jr.T_jet(0.0, 1.0, nu, R=R, K=KL)
    # the majorant coefficients Ghat_ij(E) are increasing in E (positive series in E), so the remainder bounds are
    # evaluated at the upper end E3 = 3 nu_max / 2 as a point (a ball would destroy the cancellation in rho)
    E3 = arb(float((3 * nu / 2).upper()))
    T0_hi = arb(c[(0, 0)].upper())
    for (i, j) in list(c):
        if i % 2 == 1: c[(i, j)] = arb(0)                 # exact: reflection symmetry
    c[(0, 1)] = arb(0)                                    # exact: tau = i is critical (tau -> -1/tau symmetry)
    if c[(2, 0)].lower() <= 0 or c[(0, 2)].lower() <= 0: return None
    # inner square: direct form with the order-3 coefficients
    G3 = jr.majorant_coeffs(E3, 3)
    def direct_ok(r):
        r = arb(r)
        try:
            low = sum((G3[i][j] * r ** (i + j) for i in range(4) for j in range(4 - i)), arb(0))
            rho = T0_hi * (jr.majorant_value(E3, r, r) - low) / (r * r)
        except AssertionError:
            return False
        return (c[(2, 0)] - arb(mag(c[(2, 1)])) * r - rho).lower() > 0 and (c[(0, 2)] - arb(mag(c[(0, 3)])) * r - rho).lower() > 0
    rloc = 0.0
    for r in np.geomspace(1e-4, 0.05, 120):
        if direct_ok(r): rloc = float(r)
        else: break
    if rloc == 0:
        if DEBUG: print('inner lemma fails')
        return None
    # ring: polynomial lower bound on sub-boxes of [-R1, R1]^2 outside the inner square
    keys = [(i, j) for i in range(KL + 1) for j in range(KL + 1 - i) if i + j >= 2 and i % 2 == 0]
    def box_ok(a0, a1, b0, b1):
        A = arb((a0 + a1) / 2, (a1 - a0) / 2); B = arb((b0 + b1) / 2, (b1 - b0) / 2)
        Am, Bm = arb(max(abs(a0), abs(a1))), arb(max(abs(b0), abs(b1)))
        try: rem = T0_hi * jr.majorant_tail(E3, Am, Bm, KL)
        except AssertionError: return False
        # split the quadratic part exactly (A^2, B^2 >= known lower bounds on the box)
        A2lo = 0.0 if a0 <= 0 <= a1 else min(a0 * a0, a1 * a1); B2lo = 0.0 if b0 <= 0 <= b1 else min(b0 * b0, b1 * b1)
        q = c[(2, 0)].lower() * down(A2lo) + c[(0, 2)].lower() * down(B2lo)
        pA = [arb(1)]; pB = [arb(1)]
        for _ in range(KL): pA.append(pA[-1] * A); pB.append(pB[-1] * B)          # explicit powers (0-containing ball ** 0 is nan)
        rest = sum((c[k] * pA[k[0]] * pB[k[1]] for k in keys if k not in ((2, 0), (0, 2))), arb(0))
        return (arb(q) + rest - rem).lower() > 0
    hl = 0.99 * rloc; stack = [(-R1, R1, -R1, R1)]; nb = 0
    while stack:
        a0, a1, b0, b1 = stack.pop(); nb += 1
        if -hl <= a0 and a1 <= hl and -hl <= b0 and b1 <= hl: continue
        if box_ok(a0, a1, b0, b1): continue
        if max(a1 - a0, b1 - b0) < 1e-5 or nb > 200000:
            if DEBUG: print('ring fail at', (a0, a1, b0, b1), 'rloc', rloc)
            return None
        am, bm = (a0 + a1) / 2, (b0 + b1) / 2
        if a1 - a0 >= b1 - b0: stack += [(a0, am, b0, b1), (am, a1, b0, b1)]
        else: stack += [(a0, a1, b0, bm), (a0, a1, bm, b1)]
    return rloc, nb

t0 = time.time()
slices = [(NU_A, NU_B)]; done_slices = []
while slices:
    a, b = slices.pop()
    res = local_ok(a, b)
    if res is None:
        if b - a < 1e-4: say("LOCAL FAIL", a, b); sys.exit(1)
        m = (a + b) / 2; slices += [(m, b), (a, m)]; continue
    done_slices.append((a, b, res[0], res[1]))
    say("  local region certified for nu in [%.6f, %.6f]: r_loc = %.5f, ring boxes %d  (%.0fs)" % (a, b, res[0], res[1], time.time() - t0))
say("(L) local region |x|, |y-1| <= %g certified on %d nu-slices" % (R1, len(done_slices)))

# ---------------------------------------------------------------- (Y) large-y cut, uniform in nu
nuB = ball(NU_A, NU_B); sB = nuB / 2
cB, _ = jr.T_jet(0.0, 1.0, nuB, R=R, K=1)
T0_all_hi = arb(cB[(0, 0)].upper())
K = arb(0)
for j in range(-40, 41):
    for k in range(-40, 41):
        if j and k and j != k: K += (arb(abs(j) * abs(k) * abs(j - k))) ** (-2 * sB)
Klow = arb(K.lower())
# y^(3s) K > T0 for all s in the ball: y > (T0/K)^(1/(3s)); take the max over the ball (T0/K > 1 or < 1 handled by ball)
Yc = up(float(((T0_all_hi / Klow) ** (1 / (3 * sB))).upper())); Yc = max(Yc, 1.0 + R1 + 1e-3)
say("(Y) large-y cut: T > T(i) for all nu in the interval when y > Y0 = %.6f" % Yc)

# ---------------------------------------------------------------- (O) outer branch and bound over (x, y, nu)
_c0 = {}
def at_i(a, b):
    """point data at tau = i for the slice [a, b]: T(i; nu0) upper, sup_N (T_{nu-1} + T_{nu+1})(i) / e"""
    key = (a, b)
    if key not in _c0:
        nu0 = (a + b) / 2
        c, _ = jr.T_jet(0.0, 1.0, nu0, R=R, K=1)
        cm, _ = jr.T_jet(0.0, 1.0, ball(a, b) - 1, R=R, K=1)
        cp, _ = jr.T_jet(0.0, 1.0, ball(a, b) + 1, R=R, K=1)
        _c0[key] = (float(c[(0, 0)].upper()), float(((cm[(0, 0)] + cp[(0, 0)]) / math.e).upper()))
    return _c0[key]
def T_lower_box(x0, x1, y0, y1, nu):
    """lower bound of T over the spatial box (second-order Taylor at the centre minus majorant remainder)"""
    xc, yc = (x0 + x1) / 2, (y0 + y1) / 2
    hx, hy = up(max(xc - x0, x1 - xc)), up(max(yc - y0, y1 - yc))
    cj, _ = jr.T_jet(xc, yc, nu, R=R, K=2)
    Ah, Bh = arb(arb(hx) / arb(yc)).upper(), arb(arb(hy) / arb(yc)).upper(); Ah, Bh = arb(Ah), arb(Bh)
    E3 = 3 * arb(nu) / 2
    P = cj[(0, 0)] - arb(mag(cj[(1, 0)])) * Ah - arb(mag(cj[(0, 1)])) * Bh \
        - arb(max(0.0, -float(cj[(2, 0)].lower()))) * Ah * Ah - arb(mag(cj[(1, 1)])) * Ah * Bh - arb(max(0.0, -float(cj[(0, 2)].lower()))) * Bh * Bh
    rem = arb(cj[(0, 0)].upper()) * jr.majorant_tail(E3, Ah, Bh, 2)
    return P - rem, Ah, Bh
def T_upper_box(x0, x1, y0, y1, nuball):
    xc, yc = (x0 + x1) / 2, (y0 + y1) / 2
    hx, hy = up(max(xc - x0, x1 - xc)), up(max(yc - y0, y1 - yc))
    cj, _ = jr.T_jet(xc, yc, nuball, R=R, K=1)
    Ah, Bh = arb(arb(arb(hx) / arb(yc)).upper()), arb(arb(arb(hy) / arb(yc)).upper())
    E3 = arb(float((3 * nuball / 2).upper()))        # majorant increasing in the exponent
    return cj[(0, 0)] + arb(mag(cj[(1, 0)])) * Ah + arb(mag(cj[(0, 1)])) * Bh + arb(cj[(0, 0)].upper()) * jr.majorant_tail(E3, Ah, Bh, 1)
def certified(x0, x1, y0, y1, a, b):
    nu0 = (a + b) / 2; eta = (b - a) / 2
    try:
        Tl, Ah, Bh = T_lower_box(x0, x1, y0, y1, nu0)
    except AssertionError:
        return False, "space"
    T0hi, U0 = at_i(a, b)
    D0 = float(Tl.lower()) - T0hi
    if D0 <= 0: return False, "space"
    if eta == 0: return True, ""
    try:
        Ub = T_upper_box(x0, x1, y0, y1, ball(a, b) - 1) + T_upper_box(x0, x1, y0, y1, ball(a, b) + 1)
    except AssertionError:
        return False, "space"
    slack = D0 - eta * (float((Ub / math.e).upper()) + U0)
    if slack > 0: return True, ""
    return False, ("nu" if eta * (float((Ub / math.e).upper()) + U0) > 0.5 * D0 else "space")

queue = [(0.0, 0.5, 0.8, Yc, NU_A, NU_B)]; n = 0; t1 = time.time()
while queue:
    x0, x1, y0, y1, a, b = queue.pop(); n += 1
    if up(up(x1 * x1) + up(y1 * y1)) < 1.0: continue
    if -R1 <= x0 and x1 <= R1 and 1 - R1 <= y0 and y1 <= 1 + R1: continue          # inside the local region
    ok, why = certified(x0, x1, y0, y1, a, b)
    if ok: continue
    if max(x1 - x0, y1 - y0) < 1e-7 or b - a < 1e-7: say("STALL", (x0, x1, y0, y1, a, b)); sys.exit(1)
    if why == "nu": m = (a + b) / 2; queue += [(x0, x1, y0, y1, a, m), (x0, x1, y0, y1, m, b)]
    else:
        xm, ym = (x0 + x1) / 2, (y0 + y1) / 2
        if x1 - x0 >= y1 - y0: queue += [(x0, xm, y0, y1, a, b), (xm, x1, y0, y1, a, b)]
        else: queue += [(x0, x1, y0, ym, a, b), (x0, x1, ym, y1, a, b)]
    if n % 50 == 0: say("  outer boxes %d, queue %d, %.0fs" % (n, len(queue), time.time() - t1))
say("DONE nu in [%g, %g]: local region on %d slices, %d outer (x,y,nu)-boxes certified in %.0fs -> the square lattice is the "
    "unique global minimiser of T_nu for every nu in the interval" % (NU_A, NU_B, len(done_slices), n, time.time() - t0))
json.dump(dict(nu=[NU_A, NU_B], slices=done_slices, outer_boxes=n, Y0=Yc), open("prove_nuint_%s.json" % tag, "w"), indent=1)
