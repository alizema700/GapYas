"""Figures for the verification appendix (reads results/*.json, computes the heatmaps)."""
import glob, json
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lattice_sums import triangle_sum_box, triangle_sum, A_tau, HEX, SQUARE

plt.rcParams.update({"font.size": 9, "figure.dpi": 150, "savefig.bbox": "tight"})
def save(fig, name):
    fig.savefig("figures/%s.pdf" % name); fig.savefig("figures/%s.png" % name, dpi=200)

C = {"hex": "#1f6feb", "sq": "#d1242f", "rect": "#2da44e", "fcc": "#8250df", "bcc": "#bf8700", "sc": "#57606a"}
v2 = json.load(open("results/v2_2d_landscape.json"))
v2b = json.load(open("results/v2b_high_nu.json"))
nustar = v2["nustar_reference"]; nuc2 = v2b["nu_c2"][0]

# Fig 1: relative energy difference square vs hex, and rectangular branch
nus = np.linspace(2.2, 12, 120)
dh = [(triangle_sum_box(A_tau(SQUARE), (n,) * 3, 64) / triangle_sum_box(A_tau(HEX), (n,) * 3, 64) - 1) for n in nus]
fig, ax = plt.subplots(figsize=(5.2, 3.0))
ax.axhline(0, color="k", lw=0.6)
ax.plot(nus, np.array(dh) * 100, color=C["sq"], label=r"Quadrat: $T_{\mathrm{quad}}/T_{\mathrm{hex}}-1$")
tr = v2b["track"]; kk = sorted(tr, key=float)
ax.plot([float(k) for k in kk if float(k) <= 12], [(min(tr[k]["rect_T"], tr[k]["square_T"]) / tr[k]["hex_T"] - 1) * 100 for k in kk if float(k) <= 12],
        "o", ms=3, color=C["rect"], label=r"bestes Rechteck$/T_{\mathrm{hex}}-1$")
ax.axvline(nustar, color=C["hex"], ls="--", lw=0.8); ax.axvline(nuc2, color=C["rect"], ls="--", lw=0.8)
ax.text(nustar + 0.1, 5, r"$\nu^*=%.4f$" % nustar, fontsize=8); ax.text(nuc2 + 0.1, 5, r"$\nu_2=%.3f$" % nuc2, fontsize=8)
ax.set_xlabel(r"Exponent $\nu$"); ax.set_ylabel("Energiedifferenz in %"); ax.set_ylim(-60, 8)
ax.legend(frameon=False, fontsize=8, loc="lower left")
save(fig, "fig_2d_energy_difference"); plt.close(fig)

# Fig 2: stability (smallest Hessian eigenvalue) of hex and square
st = v2["stability"]; ks = sorted(st, key=float)
fig, ax = plt.subplots(figsize=(5.2, 2.8))
ax.axhline(0, color="k", lw=0.6)
ax.plot([float(k) for k in ks], [st[k]["hex_eigs"][0] for k in ks], "o-", ms=3, color=C["hex"], label="Hexagonal")
ax.plot([float(k) for k in ks], [st[k]["square_eigs"][0] for k in ks], "s-", ms=3, color=C["sq"], label="Quadrat")
cs = v2b["curvature_square"]
ax.plot([c[0] for c in cs], [min(c[1], c[2]) for c in cs], "s-", ms=3, color=C["sq"])
for x in v2["square_spinodals"] + v2["hex_spinodals"] + v2b["nu_c2"]:
    ax.axvline(x, color="0.5", ls=":", lw=0.8)
ax.axvspan(v2["square_spinodals"][0], v2["hex_spinodals"][0], color="0.9", zorder=0)
ax.set_xlabel(r"Exponent $\nu$"); ax.set_ylabel("kleinster Hesse-Eigenwert")
ax.legend(frameon=False, fontsize=8); save(fig, "fig_2d_stability"); plt.close(fig)

# Fig 3: landscape over the fundamental domain for four exponents
fig, axs = plt.subplots(1, 4, figsize=(7.2, 2.6), sharey=True)
xs = np.linspace(0, 0.5, 41); ys = np.linspace(0.85, 1.6, 61)
for ax, nu in zip(axs, [3.0, nustar, 6.0, 10.0]):
    Z = np.full((len(ys), len(xs)), np.nan)
    for i, y in enumerate(ys):
        for j, x in enumerate(xs):
            if x * x + y * y >= 1:
                Z[i, j] = triangle_sum_box(A_tau(x + 1j * y), (nu,) * 3, 40)
    Z = Z / np.nanmin(Z) - 1
    im = ax.pcolormesh(xs, ys, np.log10(Z + 1e-5), shading="auto", cmap="viridis", vmin=-5, vmax=0)
    th = np.linspace(np.pi / 3, np.pi / 2, 50); ax.plot(np.cos(th), np.sin(th), "w-", lw=0.8)
    ax.plot([0.5], [np.sqrt(3) / 2], "o", color="w", ms=3); ax.plot([0], [1], "s", color="w", ms=3)
    ax.set_title(r"$\nu=%.3f$" % nu if nu == nustar else r"$\nu=%g$" % nu, fontsize=9)
    ax.set_xlabel(r"Re $\tau$")
axs[0].set_ylabel(r"Im $\tau$")
cb = fig.colorbar(im, ax=axs, shrink=0.85, pad=0.02); cb.set_label(r"$\log_{10}(T/T_{\min}-1)$")
save(fig, "fig_2d_landscape"); plt.close(fig)

# Fig 4: barrier on the arc |tau|=1 at nu*
arc = np.array(v2["arc_at_nustar"])
fig, ax = plt.subplots(figsize=(4.2, 2.5))
ax.plot(arc[:, 0], arc[:, 1], "-", color="k")
ax.plot([60], [arc[0, 1]], "o", color=C["hex"]); ax.plot([90], [arc[-1, 1]], "s", color=C["sq"])
ax.set_xlabel(r"Winkel $\arg\tau$ in Grad (60° = hex, 90° = Quadrat)"); ax.set_ylabel(r"$T_{\nu^*}(\tau)$")
save(fig, "fig_2d_arc_barrier"); plt.close(fig)

# Fig 5: rectangular aspect ratio of the optimum vs nu
fig, ax = plt.subplots(figsize=(4.2, 2.5))
ax.plot([float(k) for k in kk], [tr[k]["glob_y"] if tr[k]["glob_x"] < 1e-3 else np.nan for k in kk], "o-", ms=3, color=C["rect"])
ax.axvline(nuc2, color=C["rect"], ls="--", lw=0.8)
ax.set_xscale("log"); ax.set_xlabel(r"Exponent $\nu$"); ax.set_ylabel("Seitenverhältnis $b/a$ des Optimums")
save(fig, "fig_2d_rect_ratio"); plt.close(fig)

# Fig 6: 3D FCC/SC relative to BCC (BCC = best lattice found by the global search at every nu)
rows = {}
for f in glob.glob("results/v3_3d_*.json"):
    rows.update(json.load(open(f)))
if rows:
    ks = sorted(rows, key=float)
    fig, ax = plt.subplots(figsize=(5.2, 2.8))
    for name, col in [("FCC", C["fcc"]), ("SC", C["sc"])]:
        ax.plot([float(k) for k in ks], [(rows[k]["reference"][name] / rows[k]["reference"]["BCC"] - 1) * 100 for k in ks], "o-", ms=3, color=col, label=name)
    ax.set_yscale("log"); ax.set_xlabel(r"Exponent $\nu$"); ax.set_ylabel(r"$T/T_{\mathrm{BCC}}-1$ in %")
    ax.legend(frameon=False, fontsize=8); save(fig, "fig_3d_lattices"); plt.close(fig)
print("figures written")
