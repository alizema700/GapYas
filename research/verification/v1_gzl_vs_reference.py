"""V1: GZL triangle zeta (zeta_circle, evaluate_graph) vs. independent brute force."""
import time, json
import numpy as np, gzl
from lattice_sums import triangle_sum, A_tau, HEX, SQUARE, FCC, BCC, SC, unit_covolume

rng = np.random.default_rng(1)
cases = []
for tau in [HEX, SQUARE, 0.21 + 1.37j, 0.37 + 0.95j, 0.05 + 2.4j]:
    for nus in [(2.5,) * 3, (3.0,) * 3, (3.918,) * 3, (4.0,) * 3, (6.0,) * 3, (9.0,) * 3, (3.0, 4.5, 6.0)]:
        cases.append(("2D tau=%s" % np.round(tau, 3), A_tau(tau), nus))
Arand = unit_covolume(np.array([[1.0, 0.31, -0.2], [0.0, 1.1, 0.45], [0.0, 0.0, 0.83]]))
for name, A in [("FCC", FCC), ("BCC", BCC), ("SC", SC), ("triclinic", Arand)]:
    for nus in [(3.5,) * 3, (4.5,) * 3, (6.0,) * 3, (9.0,) * 3, (4.0, 5.0, 7.0)]:
        cases.append(("3D " + name, unit_covolume(A), nus))

rows = []
worst = 0.0
for name, A, nus in cases:
    d = A.shape[0]
    Rs = (48, 64, 96, 128, 192) if d == 2 else (16, 24, 32, 40, 48)
    t = time.time(); ref = triangle_sum(A, nus, Rs=Rs); t_ref = time.time() - t
    t = time.time(); gc = gzl.zeta_circle(list(nus), A); t_gzl = time.time() - t
    edges = np.array([[0, 1], [1, 2], [2, 0]])
    ge = float(gzl.evaluate_graph(edges, np.array(nus, float), A))
    rel = abs(gc - ref) / abs(ref); rel2 = abs(ge - gc) / abs(gc)
    worst = max(worst, rel)
    rows.append(dict(case=name, nus=list(nus), reference=ref, gzl_circle=gc, gzl_eval=ge,
                     rel_circle_vs_ref=rel, rel_eval_vs_circle=rel2, t_ref=t_ref, t_gzl=t_gzl))
    print("%-22s nu=%-18s ref=%.15g gzl=%.15g rel=%.1e | eval_graph rel=%.1e | %.2fs/%.2fs"
          % (name, nus, ref, gc, rel, rel2, t_ref, t_gzl))
print("worst relative deviation zeta_circle vs reference: %.2e" % worst)
json.dump(rows, open("results/v1_gzl_vs_reference.json", "w"), indent=1)
