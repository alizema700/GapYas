"""V2b: second transition at large nu (square -> rectangular) and the nu -> infinity limit."""
import json, time
import numpy as np
from scipy.optimize import brentq, minimize, minimize_scalar
from lattice_sums import triangle_sum, triangle_sum_box, A_tau, HEX, SQUARE

def T(tau, nu, fast=False):
    if fast:
        return triangle_sum_box(A_tau(tau), (nu,) * 3, 48)
    return triangle_sum(A_tau(tau), (nu,) * 3, Rs=(48, 64, 96, 128))

def in_F(x, y): return 0.0 <= x <= 0.5 and x * x + y * y >= 1.0 - 1e-12 and y <= 4.0

def global_min(nu):
    pts = []
    for x in np.linspace(0, 0.5, 26):
        y0 = np.sqrt(1 - x * x)
        for y in y0 + np.concatenate([[0], np.geomspace(1e-3, 3.0, 40)]):
            pts.append((T(x + 1j * y, nu, fast=True), x, y))
    pts.sort(); pol = []
    for v, x, y in pts[:8]:
        f = lambda p: T(p[0] + 1j * p[1], nu) if in_F(*p) else 1e300
        r = minimize(f, (x, y), method="Nelder-Mead", options=dict(xatol=1e-7, fatol=1e-15, maxiter=400))
        pol.append((r.fun, r.x[0], r.x[1]))
    pol.sort(); return pol[0]

def hess_square(nu, h=2e-3):
    f = lambda u, v: T(SQUARE + (u + 1j * v), nu)
    f0 = f(0, 0)
    fuu = (f(h, 0) - 2 * f0 + f(-h, 0)) / h**2
    fvv = (f(0, h) - 2 * f0 + f(0, -h)) / h**2
    return fuu, fvv   # (u: shear towards rhombic, v: stretch towards rectangular)

out = {}
# second-order point: curvature of T along the rectangular direction at the square lattice
nus = np.arange(7.6, 9.61, 0.2)
curv = [(nu, *hess_square(nu)) for nu in nus]
for c in curv: print("nu=%.2f  d2T/dx2 (shear) %.6f   d2T/dy2 (stretch) %.6f" % c, flush=True)
g = lambda nu: hess_square(nu)[1]
roots = [brentq(g, a[0], b[0], xtol=1e-6) for a, b in zip(curv, curv[1:]) if a[2] * b[2] < 0]
print("square lattice loses stability towards rectangular at nu =", roots, flush=True)
out["curvature_square"] = curv; out["nu_c2"] = roots

# optimum along the rectangular line and global optimum for larger nu
track = {}
for nu in [8.0, 8.5, 9.0, 9.5, 10.0, 12.0, 15.0, 20.0, 30.0, 50.0, 80.0]:
    t = time.time()
    r = minimize_scalar(lambda y: T(1j * y, nu), bounds=(1.0, 2.0), method="bounded", options=dict(xatol=1e-9))
    v, x, y = global_min(nu)
    track[nu] = dict(rect_y=r.x, rect_T=r.fun, square_T=T(SQUARE, nu), hex_T=T(HEX, nu), glob_T=v, glob_x=x, glob_y=y,
                     rect_side_ratio=r.x)
    print("nu=%5.1f  rect optimum y=%.6f (side ratio b/a=%.6f) T=%.6e | square %.6e | hex %.6e | global %.6e at tau=%.5f+%.5fi | %.0fs"
          % (nu, r.x, r.x, r.fun, track[nu]["square_T"], track[nu]["hex_T"], v, x, y, time.time() - t), flush=True)
out["track"] = track
json.dump(out, open("results/v2b_high_nu.json", "w"), indent=1, default=float)
