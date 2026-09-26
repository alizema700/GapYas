"""E2: FCC instability with error bars. Curvature of T_nu along the Bain path at FCC (c/a = sqrt2),
for several finite-difference steps h and cut-off sets; root nu_c of the curvature with spread."""
import json
import numpy as np
from scipy.optimize import brentq
from common import *

def bct(r): return 0.5 * np.array([[-1, 1, 1], [1, -1, 1], [r, r, -r]], float)
RSETS = {"R16-48": (16, 24, 32, 40, 48), "R24-64": (24, 32, 40, 48, 64)}
def T(r, nu, Rs): return triangle_sum(lll(unit_covolume(bct(r))), (nu,) * 3, Rs=Rs)
def curv(nu, h, Rs):
    r0 = np.sqrt(2.0)
    return (T(r0 * np.exp(h), nu, Rs) - 2 * T(r0, nu, Rs) + T(r0 * np.exp(-h), nu, Rs)) / h**2
out = {}
for rname, Rs in RSETS.items():
    for h in (1e-3, 2e-3, 4e-3):
        g = lambda nu: curv(nu, h, Rs)
        grid = [(nu, g(nu)) for nu in (3.6, 3.7, 3.8, 3.9)]
        roots = [brentq(g, a[0], b[0], xtol=1e-6) for a, b in zip(grid, grid[1:]) if a[1] * b[1] < 0]
        out["%s_h%g" % (rname, h)] = dict(grid=grid, roots=roots)
        print(rname, "h=%g" % h, "curvatures", [(round(a, 2), round(b, 5)) for a, b in grid], "root", roots, flush=True)
allroots = [r for v in out.values() for r in v["roots"]]
print("nu_c(FCC unstable along Bain) = %.6f +- %.6f (spread over h and cut-off)" % (np.mean(allroots), (max(allroots) - min(allroots)) / 2))
out["summary"] = dict(mean=float(np.mean(allroots)), halfspread=float((max(allroots) - min(allroots)) / 2))
json.dump(out, open("results/e2_fcc_bain_curvature.json", "w"), indent=1)
