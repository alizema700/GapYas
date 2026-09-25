"""V2c: large-nu regime with plain cube sums (R = 48 is converged to machine precision
for nu >= 12; the Richardson fit of V2b becomes ill-conditioned there and is not used)."""
import json, time
import numpy as np
from scipy.optimize import minimize, minimize_scalar
from lattice_sums import triangle_sum_box, A_tau, HEX, SQUARE

T = lambda tau, nu, R=48: triangle_sum_box(A_tau(tau), (nu,) * 3, R)
def in_F(x, y): return 0.0 <= x <= 0.5 and x * x + y * y >= 1.0 - 1e-12 and y <= 4.0
out = {}
for nu in [12.0, 15.0, 20.0, 30.0, 50.0, 80.0, 120.0]:
    t = time.time()
    pts = []
    for x in np.linspace(0, 0.5, 26):
        y0 = np.sqrt(1 - x * x)
        for y in y0 + np.concatenate([[0], np.geomspace(1e-3, 3.0, 40)]):
            pts.append((T(x + 1j * y, nu, 24), x, y))
    pts.sort(); pol = []
    for v, x, y in pts[:8]:
        f = lambda p: T(p[0] + 1j * p[1], nu) if in_F(*p) else 1e300
        r = minimize(f, (x, y), method="Nelder-Mead", options=dict(xatol=1e-8, fatol=0, maxiter=600))
        pol.append((r.fun, r.x[0], r.x[1]))
    pol.sort()
    rect = minimize_scalar(lambda y: T(1j * y, nu), bounds=(1.0, 2.0), method="bounded", options=dict(xatol=1e-10))
    conv = abs(T(1j * rect.x, nu, 64) / rect.fun - 1)
    out[nu] = dict(glob_T=pol[0][0], glob_x=pol[0][1], glob_y=pol[0][2], rect_y=rect.x, rect_T=rect.fun,
                   square_T=T(SQUARE, nu), hex_T=T(HEX, nu), R48_vs_R64=conv)
    print("nu=%6.1f global T=%.10e at tau=%.6f+%.6fi | rect y=%.6f T=%.10e | square %.6e | hex %.6e | conv %.1e | %.0fs"
          % (nu, pol[0][0], pol[0][1], pol[0][2], rect.x, rect.fun, out[nu]["square_T"], out[nu]["hex_T"], conv, time.time() - t), flush=True)
json.dump(out, open("results/v2c_large_nu.json", "w"), indent=1, default=float)
