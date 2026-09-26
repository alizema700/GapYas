"""Data for Fig. (asymptotics of the rectangular minimiser): numerical minimisers of T_nu on the imaginary axis
(direct summation, x = 0 by symmetry) and roots of the balance equation (y sqrt(1+y^2)/2)^nu = 4(y^2-1)/(3(y^2+1))."""
import json, numpy as np, mpmath as mp
from scipy.optimize import minimize_scalar
exec(open("asym2d.py").read().split("a = ((1")[0])
mp.mp.dps = 30
w = (mp.sqrt(17) - 1) / 2; yinf = mp.sqrt(w)
kappa = yinf * (w + 1) / (2 * w + 1) * mp.log(3 * (w + 1) / (4 * (w - 1)))
def bal(nu):
    g = lambda y: nu * mp.log(y * mp.sqrt(1 + y * y) / 2) - mp.log(4 * (y * y - 1) / (3 * (y * y + 1)))
    return float(mp.findroot(g, yinf - kappa / nu))
out = dict(yinf=float(yinf), kappa=float(kappa), num=[], bal=[])
for nu in [9, 10, 11, 12, 14, 16, 18, 20, 25, 30, 40, 50, 70, 100, 150, 200, 300]:
    r = minimize_scalar(lambda y: T(0.0, y, nu, R=6), bounds=(1.0, 1.3), method="bounded", options={"xatol": 1e-12})
    out["num"].append([nu, float(r.x)]); print(nu, r.x, flush=True)
for nu in np.geomspace(9, 400, 120):
    try: out["bal"].append([float(nu), bal(nu)])
    except Exception: pass
json.dump(out, open("asym_data.json", "w"), indent=1)
