"""V2: 2D three-body lattice energy T_nu(tau) over the whole fundamental domain.

Independent reference sums only (lattice_sums.py); GZL is used as a cross-check
at the end.  Questions: (a) which lattice is the global minimiser for each nu,
(b) where exactly hex and square exchange (nu*), (c) is the transition first
order (both local minima, barrier in between), (d) spinodals: the nu-range in
which hex / square are local minima.
"""
import json, time
import numpy as np
from scipy.optimize import brentq, minimize
from lattice_sums import triangle_sum, triangle_sum_box, A_tau, HEX, SQUARE
import gzl

def T(tau, nu, fast=False):
    if fast and nu >= 3.0:
        return triangle_sum_box(A_tau(tau), (nu,) * 3, 48)
    Rs = (32, 48, 64) if fast else (48, 64, 96, 128)
    return triangle_sum(A_tau(tau), (nu,) * 3, Rs=Rs)

def in_F(x, y):
    return 0.0 <= x <= 0.5 and x * x + y * y >= 1.0 - 1e-12 and y <= 4.0

out = {}
# (a) global scan over the fundamental domain -------------------------------
nus = [2.2, 2.5, 3.0, 3.5, 3.8, 3.9, 3.95, 4.0, 4.5, 5.0, 6.0, 8.0, 10.0]
scan = {}
for nu in nus:
    t = time.time()
    best = []
    xs = np.linspace(0, 0.5, 26)
    for x in xs:
        y0 = np.sqrt(1 - x * x)
        for y in y0 + np.concatenate([[0], np.geomspace(1e-3, 3.0, 40)]):
            best.append((T(x + 1j * y, nu, fast=True), x, y))
    best.sort()
    # polish the 8 best grid points with Nelder-Mead in (x, y), accurate sums
    pol = []
    for v, x, y in best[:8]:
        f = lambda p: T(p[0] + 1j * p[1], nu) if in_F(*p) else 1e9
        r = minimize(f, (x, y), method="Nelder-Mead",
                     options=dict(xatol=1e-6, fatol=1e-13, maxiter=300))
        pol.append((r.fun, r.x[0], r.x[1]))
    pol.sort()
    h, s = T(HEX, nu), T(SQUARE, nu)
    scan[nu] = dict(min=pol[0][0], x=pol[0][1], y=pol[0][2], hex=h, square=s)
    print("nu=%5.2f  global min %.12f at tau=%.5f+%.5fi | hex %.12f | square %.12f | %.0fs"
          % (nu, pol[0][0], pol[0][1], pol[0][2], h, s, time.time() - t), flush=True)
out["scan"] = scan

# (b) crossing point nu* ------------------------------------------------------
g = lambda nu: T(SQUARE, nu) - T(HEX, nu)
nustar = brentq(g, 3.8, 4.0, xtol=1e-12)
gz = lambda nu: gzl.zeta_circle([nu] * 3, A_tau(SQUARE)) - gzl.zeta_circle([nu] * 3, A_tau(HEX))
nustar_gzl = brentq(gz, 3.8, 4.0, xtol=1e-12)
print("nu* (reference) = %.10f   nu* (GZL) = %.10f   s* = nu*/2 = %.10f" % (nustar, nustar_gzl, nustar / 2))
out["nustar_reference"] = nustar; out["nustar_gzl"] = nustar_gzl

# (c)+(d) local stability of hex and square: Hessian in hyperbolic coordinates
def hessian(tau0, nu, h=2e-3):
    # T is a function on the upper half plane; use local coords (u, v): tau = tau0 + y0*(u + i v)
    y0 = tau0.imag
    f = lambda u, v: T(tau0 + y0 * (u + 1j * v), nu)
    f0 = f(0, 0)
    fuu = (f(h, 0) - 2 * f0 + f(-h, 0)) / h**2
    fvv = (f(0, h) - 2 * f0 + f(0, -h)) / h**2
    fuv = (f(h, h) - f(h, -h) - f(-h, h) + f(-h, -h)) / (4 * h * h)
    return np.linalg.eigvalsh(np.array([[fuu, fuv], [fuv, fvv]]))

stab = {}
for nu in np.round(np.arange(2.2, 8.01, 0.2), 3).tolist() + [3.85, 3.9, 3.95]:
    eh, es = hessian(HEX, nu), hessian(SQUARE, nu)
    stab[nu] = dict(hex_eigs=eh.tolist(), square_eigs=es.tolist())
print("stability (Hessian eigenvalues, >0 = local minimum):")
for nu in sorted(stab):
    print("  nu=%.2f hex %s  square %s" % (nu, np.round(stab[nu]["hex_eigs"], 5), np.round(stab[nu]["square_eigs"], 5)))
out["stability"] = {str(k): v for k, v in stab.items()}

# spinodal of the square lattice: smallest Hessian eigenvalue crosses zero
def sq_min_eig(nu): return hessian(SQUARE, nu)[0]
def hex_min_eig(nu): return hessian(HEX, nu)[0]
nus_grid = sorted(stab)
for name, fun in [("square", sq_min_eig), ("hex", hex_min_eig)]:
    vals = [(nu, stab[nu][name + "_eigs"][0]) for nu in nus_grid]
    roots = []
    for (a, fa), (b, fb) in zip(vals, vals[1:]):
        if fa * fb < 0:
            roots.append(brentq(fun, a, b, xtol=1e-6))
    out[name + "_spinodals"] = roots
    print("%s lattice: local-minimum boundary at nu = %s" % (name, roots))

# barrier along the arc |tau| = 1 at nu*
arc = [(th, T(np.exp(1j * np.radians(th)), nustar)) for th in np.linspace(60, 90, 31)]
out["arc_at_nustar"] = arc
i = int(np.argmax([a[1] for a in arc]))
print("arc at nu*: hex %.10f  max %.10f at theta=%.1f deg  square %.10f"
      % (arc[0][1], arc[i][1], arc[i][0], arc[-1][1]))
json.dump(out, open("results/v2_2d_landscape.json", "w"), indent=1, default=float)
