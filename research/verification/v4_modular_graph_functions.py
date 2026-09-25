"""V4: modular graph functions as triangle lattice sums (independent of GZL).

C_{a,b,c}(tau) = sum_{p1+p2+p3=0, p_i != 0} prod_i (tau2/pi)^{a_i} / |p_i|^{2 a_i},  p in Z + tau Z
              = (tau2/pi)^(a+b+c) * T(A_tau_unnormalised; 2a, 2b, 2c)   (dual graph = triangle)
E_s(tau)      = sum_{p != 0} (tau2/pi)^s / |p|^{2s}

Checks the identities C_{1,1,1} = E_3 + zeta(3) (Zagier) and C_{2,2,1} = 2/5 E_5 + zeta(5)/30
with the cube-sum reference (a = 1 edges sit exactly at nu = d = 2, which GZL refuses; the
reference sum converges absolutely and is extrapolated in R), and compares GZL where all nu > d.
"""
import json
import numpy as np, mpmath as mp, gzl
from epsteinlib import epstein_zeta
from lattice_sums import triangle_sum, epstein_direct

def A_raw(tau): return np.array([[1.0, tau.real], [0.0, tau.imag]])
def E(s, tau, direct=False):
    A = A_raw(tau)
    z = epstein_direct(A, 2 * s) if direct else float(epstein_zeta(2 * s, A, np.zeros(2), np.zeros(2)).real)
    return (tau.imag / np.pi) ** s * z
def C(a, tau, Rs):
    return (tau.imag / np.pi) ** sum(a) * triangle_sum(A_raw(tau), [2 * x for x in a], Rs=Rs)

z3, z5 = float(mp.zeta(3)), float(mp.zeta(5))
rows = []
for tau in [1j, 0.5 + 0.5j * np.sqrt(3), 0.21 + 1.37j, -0.4 + 0.95j, 0.1 + 2.0j]:
    e3, e3d = E(3, tau), E(3, tau, direct=True)
    e5 = E(5, tau)
    c111 = C((1, 1, 1), tau, Rs=(128, 192, 256, 384, 512, 768))
    c221 = C((2, 2, 1), tau, Rs=(64, 96, 128, 192, 256, 384))
    c222_ref = C((2, 2, 2), tau, Rs=(48, 64, 96, 128))
    c222_gzl = (tau.imag / np.pi) ** 6 * gzl.zeta_circle([4.0, 4.0, 4.0], A_raw(tau))
    r = dict(tau=str(tau), E3_epsteinlib=e3, E3_direct=e3d,
             C111=c111, E3_plus_z3=e3 + z3, rel111=abs(c111 - e3 - z3) / (e3 + z3),
             C221=c221, rhs221=0.4 * e5 + z5 / 30, rel221=abs(c221 - 0.4 * e5 - z5 / 30) / (0.4 * e5 + z5 / 30),
             C222_ref=c222_ref, C222_gzl=c222_gzl, rel222=abs(c222_ref - c222_gzl) / c222_ref)
    rows.append(r)
    print("tau=%-22s E3 lib/direct rel %.1e | C111 vs E3+z3 rel %.1e | C221 vs 2/5E5+z5/30 rel %.1e | C222 ref vs GZL rel %.1e"
          % (tau, abs(e3 - e3d) / e3, r["rel111"], r["rel221"], r["rel222"]), flush=True)
json.dump(rows, open("results/v4_mgf.json", "w"), indent=1, default=float)
