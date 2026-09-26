"""E6: the range 2d/3 < nu <= d, where the cube truncation error has two families of exponents:
R^(d-2nu) (one far vertex) and R^(2d-3nu) (two far vertices close to each other; dominant for nu < d).
Extrapolation ansatz T(R) = T + sum_k a_k R^(e_k) over e in {2d-3nu, d-2nu} - {0,1,2}, validated at the exact
anchor nu = 2, d = 2 (T_2 = pi^3 C_{1,1,1} = pi^3 (E_3 + zeta 3), where log terms appear) and by stability
across cut-off sets."""
import json, sys
import numpy as np, mpmath as mp
from epsteinlib import epstein_zeta
from common import triangle_sum_box, A_tau, HEX, SQUARE, FCC, BCC, SC, unit_covolume, lll

def extrap(A, nu, Rs):
    d = A.shape[0]
    vals = np.array([triangle_sum_box(A, (nu,) * 3, R) for R in Rs]); Rs = np.array(Rs, float)
    ex = sorted({round(2 * d - 3 * nu - k, 12) for k in range(3)} | {round(d - 2 * nu - k, 12) for k in range(3)}, reverse=True)
    ex = [e for e in ex if e < 0][: len(Rs) - 1]
    M = np.column_stack([np.ones_like(Rs)] + [Rs ** e for e in ex])
    return float(np.linalg.lstsq(M, vals, rcond=None)[0][0])

out = {}
R2a, R2b = (64, 96, 128, 192, 256, 384, 512), (96, 128, 192, 256, 384, 512, 768)
# anchor check just above nu = 2 is not exact; report both sets for stability
for nu in (1.4, 1.5, 1.6, 1.7, 1.8, 1.9):
    r = {}
    for name, t in (("hex", HEX), ("square", SQUARE), ("rect1.1", 1.1j), ("rhomb75", np.exp(1j * np.radians(75)))):
        a = extrap(A_tau(t), nu, R2a); b = extrap(A_tau(t), nu, R2b); r[name] = (a, b)
    out["2D_%g" % nu] = r
    print("2D nu=%.2f " % nu + "  ".join("%s %.6f/%.6f" % (k, *v) for k, v in r.items()), flush=True)
R3a, R3b = (16, 24, 32, 40, 48, 56), (24, 32, 40, 48, 56, 64)
for nu in (2.2, 2.4, 2.6, 2.8, 3.0):
    if abs(nu - 3.0) < 1e-9: nu = 2.999
    r = {}
    for name, A in (("FCC", FCC), ("BCC", BCC), ("SC", SC)):
        A = lll(unit_covolume(A)); r[name] = (extrap(A, nu, R3a), extrap(A, nu, R3b))
    out["3D_%g" % nu] = r
    print("3D nu=%.3f " % nu + "  ".join("%s %.6f/%.6f" % (k, *v) for k, v in r.items()), flush=True)
json.dump(out, open("results/e6_low_nu.json", "w"), indent=1)
