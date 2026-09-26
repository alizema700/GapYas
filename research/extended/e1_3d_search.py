"""E1: dense global search over 3D Bravais lattices (unit covolume) for T_nu.
64 random starts + FCC/BCC/SC per nu, Nelder-Mead on Richardson-extrapolated small-box sums,
then accurate re-evaluation of the best candidates.  Usage: python e1_3d_search.py nu1 nu2 ..."""
import json, sys, time
import numpy as np
from scipy.optimize import minimize
from common import *

def T(A, nu, fast=True):
    A = lll(unit_covolume(A))
    Rs = (10, 14, 18, 22) if fast else (16, 24, 32, 40, 48)
    return triangle_sum(A, (nu,) * 3, Rs=Rs)

out = {}
for nu in [float(a) for a in sys.argv[1:]]:
    t0 = time.time(); rng = np.random.default_rng(int(nu * 1000))
    f = lambda p: T(A_from_params(p), nu) if np.all(np.abs(p[:2]) < 1.5) else 1e300
    starts = [params_of(FCC), params_of(BCC), params_of(SC)] + \
             [np.concatenate([rng.normal(0, 0.3, 2), rng.normal(0, 0.6, 3)]) for _ in range(64)]
    runs = []
    for s in starts:
        r = minimize(f, s, method="Nelder-Mead", options=dict(xatol=1e-5, fatol=1e-10, maxfev=900))
        runs.append((float(r.fun), r.x.tolist()))
    runs.sort(key=lambda z: z[0])
    # accurate evaluation of distinct best candidates
    cands = []
    for v, x in runs:
        if len(cands) >= 5: break
        sh = shells(A_from_params(np.array(x)))
        if all(sh != c["shells"] for c in cands):
            cands.append(dict(fast=v, accurate=T(A_from_params(np.array(x)), nu, fast=False), shells=sh, params=x))
    ref = {n: T(A, nu, fast=False) for n, A in [("FCC", FCC), ("BCC", BCC), ("SC", SC)]}
    frac_bcc = np.mean([abs(v / ref["BCC"] - 1) < 1e-4 for v, _ in runs])
    out[nu] = dict(reference=ref, candidates=cands, fraction_of_runs_at_BCC=float(frac_bcc), n_runs=len(runs),
                   all_fast_values=[v for v, _ in runs])
    print("nu=%.2f FCC %.10f BCC %.10f SC %.10f | best accurate %.10f shells %s | %.0f%% of runs end at BCC | %.0fs"
          % (nu, ref["FCC"], ref["BCC"], ref["SC"], cands[0]["accurate"], cands[0]["shells"], 100 * frac_bcc, time.time() - t0), flush=True)
    json.dump(out, open("results/e1_3d_search_%s.json" % "_".join(sys.argv[1:]), "w"), indent=1, default=float)
