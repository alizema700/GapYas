"""E3: 2D minimiser of the triangle energy with unequal exponents (nu, nu, mu), i.e. C_{a,a,b} with
a = nu/2, b = mu/2. Global search over the fundamental domain for a grid of (nu, mu)."""
import json, sys, time
import numpy as np
from scipy.optimize import minimize
from common import *

def T(tau, nus, fast):
    if fast: return triangle_sum_box(A_tau(tau), nus, 40)
    return triangle_sum(A_tau(tau), nus, Rs=(48, 64, 96, 128))
def in_F(x, y): return 0 <= x <= 0.5 and x * x + y * y >= 1 - 1e-12 and y <= 4

def classify(x, y):
    if abs(x - 0.5) < 2e-3 and abs(y - np.sqrt(3) / 2) < 2e-3: return "hex"
    if x < 2e-3 and abs(y - 1) < 2e-3: return "square"
    if x < 2e-3: return "rect"
    if abs(x * x + y * y - 1) < 4e-3: return "rhombic"
    if abs(x - 0.5) < 2e-3: return "centred-rect"
    return "oblique"

grid_nu = np.arange(2.5, 12.01, 0.5); grid_mu = np.arange(2.5, 12.01, 0.5)
out = []
for nu in grid_nu:
    for mu in grid_mu:
        if mu < nu: continue          # symmetric in the three exponents only up to relabelling of (nu,nu,mu)? no: keep mu>=nu and mu<nu separately below
for nu in grid_nu:
    for mu in grid_mu:
        nus = (nu, nu, mu)
        if min(nu + nu, nu + mu) <= 2.0: continue
        t = time.time(); pts = []
        for x in np.linspace(0, 0.5, 21):
            y0 = np.sqrt(1 - x * x)
            for y in y0 + np.concatenate([[0], np.geomspace(1e-3, 2.5, 30)]):
                pts.append((T(x + 1j * y, nus, True), x, y))
        pts.sort(); pol = []
        for v, x, y in pts[:4]:
            f = lambda p: T(p[0] + 1j * p[1], nus, False) if in_F(*p) else 1e300
            r = minimize(f, (x, y), method="Nelder-Mead", options=dict(xatol=1e-6, fatol=1e-13, maxiter=250))
            pol.append((r.fun, r.x[0], r.x[1]))
        pol.sort(); v, x, y = pol[0]
        out.append(dict(nu=nu, mu=mu, T=v, x=x, y=y, phase=classify(x, y)))
        print("nu=%.1f mu=%.1f -> %s tau=%.4f+%.4fi T=%.8f (%.0fs)" % (nu, mu, classify(x, y), x, y, v, time.time() - t), flush=True)
        json.dump(out, open("results/e3_unequal_exponents.json", "w"), indent=1)
