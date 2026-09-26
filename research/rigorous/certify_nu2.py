"""Rigorous enclosure of nu_2: the curvature c02 of T_nu at tau = i in the stretching direction
(B = Im tau - 1) changes sign; c20 (shear) stays positive.  Coefficients from jets_rig (rigorous)."""
import sys, json
from flint import arb
import jets_rig as jr
out = {}
for nu in [float(a) for a in sys.argv[1:]]:
    c, BR = jr.T_jet(0.0, 1.0, nu, R=40, K=2)
    out[nu] = dict(c02=[str(c[(0, 2)].lower()), str(c[(0, 2)].upper())], c20=[str(c[(2, 0)].lower()), str(c[(2, 0)].upper())])
    print("nu=%.4f  c02 in [%.3e, %.3e]   c20 in [%.4f, %.4f]" % (nu, c[(0, 2)].lower(), c[(0, 2)].upper(), c[(2, 0)].lower(), c[(2, 0)].upper()), flush=True)
json.dump(out, open("certify_nu2.json", "w"), indent=1)
