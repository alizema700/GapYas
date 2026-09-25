"""V5b: exponent at which FCC stops being a local minimum (smallest Hessian eigenvalue = 0)."""
import json
import numpy as np
from scipy.optimize import brentq
from lattice_sums import FCC, unit_covolume
from v5_3d_bain_stability import hessian
f = lambda nu: hessian(unit_covolume(FCC), nu)[0]
grid = [(nu, f(nu)) for nu in (3.5, 3.75, 4.0, 4.25, 4.5)]
for g in grid: print("nu=%.2f  min Hessian eigenvalue FCC %.5f" % g, flush=True)
root = [brentq(f, a[0], b[0], xtol=1e-4) for a, b in zip(grid, grid[1:]) if a[1] * b[1] < 0]
print("FCC loses local stability at nu =", root)
json.dump(dict(grid=grid, root=root), open("results/v5b_fcc_threshold.json", "w"), indent=1)
