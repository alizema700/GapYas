"""Fig. 1 of the paper: (a) the extremal rectangular lattice with its two tight triangles,
(b) the 2D phase sequence with the rigorously established points."""
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "serif", "font.size": 8, "savefig.bbox": "tight"})
a = ((1 + np.sqrt(17)) / 8) ** 0.25; b = 1 / a
fig, axs = plt.subplots(1, 2, figsize=(7.0, 2.3), gridspec_kw=dict(width_ratios=[1, 1.6]))
ax = axs[0]
for i in range(-1, 4):
    for j in range(-1, 3):
        ax.plot(i * a, j * b, "o", color="0.25", ms=3)
ax.plot([0, a, 2 * a], [0, 0, 0], "-", color="#d1242f", lw=2)
ax.fill([0, a, 0], [0, 0, b], color="#1f6feb", alpha=0.25); ax.plot([0, a, 0, 0], [0, 0, b, 0], "-", color="#1f6feb", lw=1.5)
ax.annotate("", xy=(3 * a, -0.35), xytext=(2 * a, -0.35), arrowprops=dict(arrowstyle="<->", lw=0.7))
ax.text(2.5 * a, -0.55, r"$a_*$", ha="center")
ax.annotate("", xy=(-0.45 * a, b), xytext=(-0.45 * a, 0), arrowprops=dict(arrowstyle="<->", lw=0.7))
ax.text(-0.55 * a, 0.5 * b, r"$1/a_*$", va="center", ha="right")
ax.set_aspect("equal"); ax.axis("off"); ax.set_title(r"(a) extremal lattice $\Lambda_*$, $2a_*^3=\sqrt{a_*^2+a_*^{-2}}$", fontsize=8)
ax = axs[1]
nustar, nu2 = 3.9183649, 8.6063
ax.axvspan(4 / 3, nustar, color="#1f6feb", alpha=0.25); ax.axvspan(nustar, nu2, color="#d1242f", alpha=0.25); ax.axvspan(nu2, 13, color="#2da44e", alpha=0.25)
ax.text((4 / 3 + nustar) / 2, 0.75, "hexagonal", ha="center"); ax.text((nustar + nu2) / 2, 0.75, "square", ha="center"); ax.text(10.8, 0.75, r"rectangular $\to y_\infty$", ha="center")
ax.plot([2, 3.5], [0.35, 0.35], "v", color="#1f6feb", ms=7, label="proved: hexagonal unique minimiser")
ax.plot([4, 5, 6, 7, 8], [0.35] * 5, "s", color="#d1242f", ms=6, label="proved: square unique minimiser")
ax.axvline(nustar, color="k", lw=0.6, ls="--"); ax.axvline(nu2, color="k", lw=0.6, ls=":")
ax.text(nustar + 0.12, 0.55, r"$\nu^*\approx3.91836$", ha="left", fontsize=7)
ax.text(nu2 + 0.12, 0.55, r"$\nu_2\approx8.606$", ha="left", fontsize=7)
ax.set_xlim(4 / 3, 13); ax.set_ylim(0, 1); ax.set_yticks([]); ax.set_xlabel(r"exponent $\nu$")
ax.legend(loc="lower right", fontsize=6.5, frameon=False)
ax.set_title("(b) phase sequence in $d=2$", fontsize=8)
fig.savefig("figures/fig1_overview.pdf"); fig.savefig("figures/fig1_overview.png", dpi=200)
