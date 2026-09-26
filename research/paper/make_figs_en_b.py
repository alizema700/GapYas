"""Figures for the extended study: phase map for unequal exponents (E3) and mixed-energy thresholds (E4)."""
import sys, os; sys.path.insert(0, os.getcwd())
import json, glob
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
plt.rcParams.update({"font.size": 9, "figure.dpi": 150, "savefig.bbox": "tight"})
def save(fig, name):
    fig.savefig("../paper/figures/%s.pdf" % name)

# E3: phase map in (nu, mu) for the energy with exponents (nu, nu, mu)
rows = json.load(open("results/e3_unequal_exponents.json"))
phases = ["hex", "square", "rect", "rhombic", "centred-rect", "oblique"]
cols = ["#1f6feb", "#d1242f", "#2da44e", "#bf8700", "#8250df", "#57606a"]
nus = sorted({r["nu"] for r in rows}); mus = sorted({r["mu"] for r in rows})
Z = np.full((len(mus), len(nus)), np.nan)
for r in rows:
    Z[mus.index(r["mu"]), nus.index(r["nu"])] = phases.index(r["phase"])
fig, ax = plt.subplots(figsize=(4.6, 3.8))
ax.pcolormesh(np.array(nus), np.array(mus), Z, cmap=ListedColormap(cols), vmin=-0.5, vmax=len(phases) - 0.5, shading="nearest")
ax.plot([2.5, 12], [2.5, 12], "k--", lw=0.7)
ax.set_xlabel(r"$\nu$ (exponent of the sides $|x|$, $|y|$)"); ax.set_ylabel(r"$\mu$ (exponent of the side $|x-y|$)")
present = sorted({r["phase"] for r in rows}, key=phases.index)
for p in present:
    ax.plot([], [], "s", color=cols[phases.index(p)], label={"hex": "hexagonal", "square": "square", "rect": "rectangular", "rhombic": "rhombic", "centred-rect": "centred rect.", "oblique": "oblique"}[p])
ax.legend(frameon=True, fontsize=7, loc="upper left")
save(fig, "fig_e3_phase_map"); plt.close(fig)

# E4: critical ratio lambda_c for pair + lambda * three-body
d3 = json.load(open("results/e4_3d_3.5_4.5_6.0_9.0.json")); d2 = json.load(open("results/e4_2d_3.0_4.5_6.0_9.0.json"))
fig, ax = plt.subplots(figsize=(4.6, 2.8))
x3 = []; y3 = []
for nu, v in d3.items():
    f, b = v["table"]["fcc"][0], v["table"]["bcc"][0]
    x3.append(float(nu)); y3.append((b["Z"] - f["Z"]) / (f["T"] - b["T"]))
x2 = []; y2 = []
for nu, v in d2.items():
    pts = v["points"]
    hx = min(pts, key=lambda p: (p["x"] - .5) ** 2 + (p["y"] - np.sqrt(3) / 2) ** 2); sq = min(pts, key=lambda p: p["x"] ** 2 + (p["y"] - 1) ** 2)
    if hx["T"] > sq["T"]: x2.append(float(nu)); y2.append((sq["Z"] - hx["Z"]) / (hx["T"] - sq["T"]))
ax.plot(x3, y3, "o-", color="#bf8700", label="3D: FCC → BCC")
ax.plot(x2, y2, "s-", color="#d1242f", label="2D: hexagonal → square")
ax.set_yscale("log"); ax.set_xlabel(r"exponent $\nu$"); ax.set_ylabel(r"critical $\lambda_c=c_3/c_2$")
ax.legend(frameon=False, fontsize=8); save(fig, "fig_e4_lambda_c"); plt.close(fig)
print("ok")
