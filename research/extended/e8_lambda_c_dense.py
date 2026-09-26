"""E8: lambda_c(nu) on a dense grid for the mixed energy Z_nu + lambda T_nu (unit density).
2D: hex vs square (Z exact via epsteinlib, T via reference sums); 3D: FCC vs BCC.  lambda_c = dZ / (-dT).
Also records the two ingredients separately to explain the shape of the curve."""
import json
import numpy as np
from epsteinlib import epstein_zeta
from common import triangle_sum, A_tau, HEX, SQUARE, FCC, BCC, unit_covolume, lll
def Z(A, nu): A = np.asarray(A, float); return float(epstein_zeta(nu, A, np.zeros(len(A)), np.zeros(len(A))).real)
out = {"2d": [], "3d": []}
for nu in np.arange(4.0, 12.01, 0.25):
    Ah, As = A_tau(HEX), A_tau(SQUARE)
    Th, Ts = triangle_sum(Ah, (nu,) * 3, Rs=(48, 64, 96, 128)), triangle_sum(As, (nu,) * 3, Rs=(48, 64, 96, 128))
    Zh, Zs = Z(Ah, nu), Z(As, nu)
    lc = (Zs - Zh) / (Th - Ts) if Th > Ts else None
    out["2d"].append(dict(nu=nu, Zhex=Zh, Zsq=Zs, Thex=Th, Tsq=Ts, lam_c=lc))
    print("2D nu=%.2f  dZ=%.5f (%.3f%%)  -dT=%.5f (%.3f%%)  lambda_c=%s" % (nu, Zs - Zh, 100 * (Zs / Zh - 1), Th - Ts, 100 * (1 - Ts / Th), lc), flush=True)
for nu in np.arange(3.5, 12.01, 0.5):
    Af, Ab = lll(unit_covolume(FCC)), lll(unit_covolume(BCC))
    Tf, Tb = triangle_sum(Af, (nu,) * 3, Rs=(16, 24, 32, 40, 48)), triangle_sum(Ab, (nu,) * 3, Rs=(16, 24, 32, 40, 48))
    Zf, Zb = Z(Af, nu), Z(Ab, nu)
    out["3d"].append(dict(nu=nu, Zfcc=Zf, Zbcc=Zb, Tfcc=Tf, Tbcc=Tb, lam_c=(Zb - Zf) / (Tf - Tb)))
    print("3D nu=%.2f  lambda_c=%.5f" % (nu, (Zb - Zf) / (Tf - Tb)), flush=True)
json.dump(out, open("results/e8_lambda_c_dense.json", "w"), indent=1)
