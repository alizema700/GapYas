"""E4: hcp vs fcc vs bcc, and mixed energies E_lambda = Z_nu + lambda * T_nu at unit density.

Z_nu = pair energy per atom sum_{p != 0} |p|^-nu (Riesz, same decay per bond as the three-body term),
T_nu = three-body energy per atom (triangle sum).  lambda = c3/c2.
3D candidates: fcc, bcc, hcp(c/a), bct(c/a) (Bain path, contains fcc and bcc).  Continuous families are
tabulated on a c/a grid, so that min over the family is available for every lambda at no extra cost.
2D: the whole fundamental domain on a grid (plus the exact hex and square points).
"""
import json, sys
import numpy as np
from multisite import *
from common import triangle_sum, triangle_sum_box, A_tau, unit_covolume, lll
from epsteinlib import epstein_zeta
def epstein_direct(A, nu, Rs=None):
    return float(epstein_zeta(nu, np.asarray(A, float), np.zeros(len(A)), np.zeros(len(A))).real)

def tab3d(nu, Rt=(8, 10, 12, 14), Rz=(16, 24, 32)):
    rows = {}
    for name, (A, S) in [("fcc", fcc_prim()), ("bcc", bcc_prim())]:
        rows[name] = [dict(ca=None, Z=epstein_multisite(A, S, nu, Rz), T=triangle_sum_multisite(A, S, nu, Rt))]
    rows["hcp"] = []
    for ca in np.linspace(1.45, 1.85, 17):
        A, S = hcp(ca); rows["hcp"].append(dict(ca=ca, Z=epstein_multisite(A, S, nu, Rz), T=triangle_sum_multisite(A, S, nu, Rt)))
    rows["bct"] = []
    for r in np.linspace(0.8, 1.6, 33):
        A = unit_covolume(0.5 * np.array([[-1, 1, 1], [1, -1, 1], [r, r, -r]], float)); A = lll(A)
        rows["bct"].append(dict(ca=r, Z=epstein_direct(A, nu, Rs=(16, 24, 32)), T=triangle_sum(A, (nu,) * 3, Rs=(10, 14, 18, 22))))
    return rows

def phase3d(rows, lam):
    best = None
    for name, fam in rows.items():
        for e in fam:
            E = e["Z"] + lam * e["T"]
            if best is None or E < best[0]: best = (E, name, e["ca"])
    return best

def tab2d(nu, nx=26, ny=40):
    pts = []
    for x in np.linspace(0, 0.5, nx):
        y0 = np.sqrt(1 - x * x)
        for y in y0 + np.concatenate([[0], np.geomspace(1e-3, 1.5, ny)]):
            A = A_tau(x + 1j * y)
            pts.append(dict(x=x, y=y, Z=epstein_direct(A, nu, Rs=(32, 48, 64)), T=triangle_sum_box(A, (nu,) * 3, 48)))
    return pts

def classify2d(x, y):
    if abs(x - .5) < 1e-6 and abs(y - np.sqrt(3) / 2) < 1e-6: return "hex"
    if x < 1e-6 and abs(y - 1) < 1e-6: return "square"
    if x < 1e-6: return "rect"
    if abs(x * x + y * y - 1) < 1e-6: return "rhombic"
    return "other"

if __name__ == "__main__":
    which = sys.argv[1]
    nus = [float(a) for a in sys.argv[2:]]
    lams = np.concatenate([[0], np.geomspace(1e-3, 1e3, 121)])
    out = {}
    for nu in nus:
        if which == "3d":
            rows = tab3d(nu)
            ph = [(float(l), *phase3d(rows, l)) for l in lams]
            out[nu] = dict(table=rows, phases=ph)
            seq = []
            for l, E, n, ca in ph:
                if not seq or seq[-1][1] != n: seq.append((l, n, ca))
            print("3D nu=%.1f  pure T: fcc %.8f hcp(best) %.8f bcc %.8f | phase sequence in lambda: %s" % (
                nu, rows["fcc"][0]["T"], min(e["T"] for e in rows["hcp"]), rows["bcc"][0]["T"],
                [(round(l, 4), n, None if ca is None else round(ca, 3)) for l, n, ca in seq]), flush=True)
        else:
            pts = tab2d(nu)
            ph = []
            for l in lams:
                b = min(pts, key=lambda p: p["Z"] + l * p["T"])
                ph.append((float(l), classify2d(b["x"], b["y"]), b["x"], b["y"]))
            out[nu] = dict(points=pts, phases=ph)
            seq = []
            for l, n, x, y in ph:
                if not seq or seq[-1][1] != n: seq.append((round(l, 4), n, round(x, 3), round(y, 3)))
            print("2D nu=%.1f phase sequence in lambda: %s" % (nu, seq), flush=True)
        json.dump(out, open("results/e4_%s_%s.json" % (which, "_".join(sys.argv[2:])), "w"), indent=1, default=float)
