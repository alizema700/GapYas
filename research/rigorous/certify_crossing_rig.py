"""Rigorous enclosure of the hex/square crossing nu^* (Theorem 4, last statement).

For nu in {3.918364, 3.918366} (the IEEE doubles nearest to these decimals) we enclose T_nu(Z^2) and T_nu(Lambda_rho)
with the rigorous machinery of jets_rig.py (outward-rounded interval arithmetic, Arb powers, Higham bound for the
direct convolution, explicit tail bound B_R) at R = 100.  The hexagonal point is evaluated at the double
(0.5, fl(sqrt(3)/2)); the shift to the exact point rho is bounded with the term-wise majorant,
|T(rho) - T(tau_f)| <= T(tau_f) (G(A, B) - 1), where A = 0, B = |fl(sqrt(3)/2) - sqrt(3)/2| / fl(sqrt(3)/2).
"""
import json
import numpy as np
from flint import arb
import jets_rig as jr

R = 100
out = {}
yf = float(np.sqrt(3) / 2)
dB = abs(arb(yf) - arb(3).sqrt() / 2) / arb(yf)
dB = arb(dB.upper())
for nu in (3.918364, 3.918366):
    E3 = 3 * arb(nu) / 2
    sq, _ = jr.T_jet(0.0, 1.0, nu, R=R, K=1)
    hx, BR = jr.T_jet(0.5, yf, nu, R=R, K=1)
    shift = hx[(0, 0)] * (jr.majorant_value(E3, arb(0), dB) - 1)
    hexT = hx[(0, 0)] + arb(0, 1) * shift.upper()
    D = sq[(0, 0)] - hexT
    out[str(nu)] = dict(square=str(sq[(0, 0)]), hex=str(hexT), D=[str(D.lower()), str(D.upper())], tail=str(BR.upper()))
    print("nu=%r  T(square)=%s  T(hex)=%s  D in [%s, %s]  (tail bound %s)" % (nu, sq[(0, 0)], hexT, D.lower(), D.upper(), BR.upper()), flush=True)
ok = float(arb(out["3.918364"]["D"][0]).lower()) > 0 and float(arb(out["3.918366"]["D"][1]).upper()) < 0
print("CROSSING CERTIFIED in (3.918364, 3.918366)" if ok else "NOT CERTIFIED")
out["ok"] = ok
json.dump(out, open("certify_crossing_rig.json", "w"), indent=1)
