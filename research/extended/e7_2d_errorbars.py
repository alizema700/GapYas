"""E7: error bars for the 2D thresholds: square spinodal (~3.806), hex spinodal (~4.278), nu_2 (~8.606).
Each is the zero of a second derivative along a symmetry-adapted direction at a fixed point; we vary the
finite-difference step h and the cut-off set and report mean +- half spread."""
import json
import numpy as np
from scipy.optimize import brentq
from common import triangle_sum, A_tau, HEX, SQUARE

RSETS = {"R48-128": (48, 64, 96, 128), "R64-192": (64, 96, 128, 192)}
def T(tau, nu, Rs): return triangle_sum(A_tau(tau), (nu,) * 3, Rs=Rs)
def d2_square_stretch(nu, h, Rs):   # along tau = i*y at y = 1 (square -> rectangle)
    return (T(1j * np.exp(h), nu, Rs) - 2 * T(1j, nu, Rs) + T(1j * np.exp(-h), nu, Rs)) / h**2
def d2_square_shear(nu, h, Rs):     # along the arc |tau| = 1 (square -> rhombic)
    return (T(np.exp(1j * (np.pi / 2 + h)), nu, Rs) - 2 * T(1j, nu, Rs) + T(np.exp(1j * (np.pi / 2 - h)), nu, Rs)) / h**2
def d2_hex(nu, h, Rs):              # hex: isotropic Hessian; use the arc direction
    th = np.pi / 3
    return (T(np.exp(1j * (th + h)), nu, Rs) - 2 * T(np.exp(1j * th), nu, Rs) + T(np.exp(1j * (th - h)), nu, Rs)) / h**2
jobs = {"square_spinodal": (d2_square_shear, 3.7, 3.9), "hex_spinodal": (d2_hex, 4.2, 4.35), "nu2": (d2_square_stretch, 8.4, 8.8)}
out = {}
for name, (fun, a, b) in jobs.items():
    roots = []
    for rn, Rs in RSETS.items():
        for h in (1e-3, 2e-3, 4e-3):
            r = brentq(lambda nu: fun(nu, h, Rs), a, b, xtol=1e-7); roots.append(r)
            print(name, rn, "h=%g" % h, "root %.6f" % r, flush=True)
    out[name] = dict(roots=roots, mean=float(np.mean(roots)), halfspread=float((max(roots) - min(roots)) / 2))
    print("==> %s = %.5f +- %.5f" % (name, np.mean(roots), (max(roots) - min(roots)) / 2), flush=True)
json.dump(out, open("results/e7_2d_errorbars.json", "w"), indent=1)
