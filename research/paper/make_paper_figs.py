"""All figures of the paper in one consistent, print-oriented style.

Style: serif text matching Times (STIX math), 8 pt, only left/bottom spines, thin marks,
recessive light grid, frameless legends or direct labels, validated categorical palette
(slot 1 blue, slot 2 orange, slot 3 aqua; aqua is below 3:1 contrast, so it is always
direct-labelled).  Phases keep one colour everywhere: hexagonal = blue, square = orange,
rectangular = aqua.  Run from research/paper/.
"""
import glob, json, os, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.ticker
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "verification"))
from lattice_sums import triangle_sum_box, A_tau, HEX, SQUARE  # noqa: E402

BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#8a8984", "#e6e5e0"
PHASE = {"hex": BLUE, "square": ORANGE, "rect": AQUA}

plt.rcParams.update({
    "font.family": "serif", "font.serif": ["STIXGeneral", "DejaVu Serif"], "mathtext.fontset": "stix",
    "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7,
    "legend.fontsize": 7, "axes.linewidth": 0.6, "axes.edgecolor": INK2, "axes.labelcolor": INK,
    "xtick.color": INK2, "ytick.color": INK2, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 3, "ytick.major.size": 3, "xtick.direction": "out", "ytick.direction": "out",
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": GRID,
    "grid.linewidth": 0.5, "axes.axisbelow": True, "lines.linewidth": 1.4, "lines.markersize": 4,
    "legend.frameon": False, "savefig.bbox": "tight", "savefig.pad_inches": 0.02, "figure.dpi": 150,
})
COL = 3.35       # one column width in inches (two-column layout)
FULL = 7.0
OUT = os.path.join(HERE, "figures")
V = os.path.join(HERE, "..", "verification", "results")
E = os.path.join(HERE, "..", "extended", "results")

def save(fig, name):
    fig.savefig(os.path.join(OUT, name + ".pdf")); fig.savefig(os.path.join(OUT, name + ".png"), dpi=220); plt.close(fig)

v2 = json.load(open(os.path.join(V, "v2_2d_landscape.json")))
v2b = json.load(open(os.path.join(V, "v2b_high_nu.json")))
NUSTAR, NU2 = v2["nustar_reference"], 8.6063
SQ_SPIN, HEX_SPIN = 3.80632, 4.27848

# ---------------------------------------------------------------- Fig. 1 overview
def fig1():
    a = ((1 + np.sqrt(17)) / 8) ** 0.25; b = 1 / a
    fig = plt.figure(figsize=(FULL, 2.05))
    ax = fig.add_axes([0.0, 0.04, 0.30, 0.84])
    for i in range(-1, 4):
        for j in range(0, 3):
            ax.plot(i * a, j * b, "o", color=MUTED, ms=2.6, zorder=2)
    ax.fill([0, a, 0], [0, 0, b], color=BLUE, alpha=0.18, lw=0, zorder=1)
    ax.plot([0, a, 0, 0], [0, 0, b, 0], color=BLUE, lw=1.3, zorder=3)
    ax.plot([0, 2 * a], [0, 0], color=ORANGE, lw=1.6, zorder=4)
    for p in [(0, 0), (a, 0), (2 * a, 0), (0, b)]:
        ax.plot(*p, "o", color=INK, ms=3.2, zorder=5)
    ax.text(-1.25 * a, -0.62 * b, r"collinear triple:  $\pi=2a_*^3$", color=ORANGE, ha="left", va="center", fontsize=7)
    ax.text(-1.25 * a, -0.95 * b, r"right triangle:  $\pi=\sqrt{a_*^2+a_*^{-2}}$", color=BLUE, ha="left", va="center", fontsize=7)
    ax.annotate("", xy=(3 * a, 2.28 * b), xytext=(2 * a, 2.28 * b), arrowprops=dict(arrowstyle="<->", lw=0.6, color=INK2))
    ax.text(2.5 * a, 2.36 * b, r"$a_*$", ha="center", va="bottom", color=INK2, fontsize=7)
    ax.annotate("", xy=(3.28 * a, 2 * b), xytext=(3.28 * a, b), arrowprops=dict(arrowstyle="<->", lw=0.6, color=INK2))
    ax.text(3.36 * a, 1.5 * b, r"$1/a_*$", va="center", color=INK2, fontsize=7)
    ax.set_xlim(-1.3 * a, 3.9 * a); ax.set_ylim(-1.15 * b, 2.6 * b); ax.set_aspect("equal"); ax.axis("off")
    ax.set_title(r"(a) extremal lattice $\Lambda_*$", loc="left", pad=2)

    ax = fig.add_axes([0.38, 0.22, 0.61, 0.62])
    ax.grid(False); ax.spines["left"].set_visible(False); ax.set_yticks([])
    x0, x1 = 4 / 3, 13.0
    segs = [(x0, NUSTAR, "hex", "hexagonal"), (NUSTAR, NU2, "square", "square"), (NU2, x1, "rect", r"rectangular $\rightarrow y_\infty$")]
    for s0, s1, key, lab in segs:
        ax.add_patch(plt.Rectangle((s0, 0.62), s1 - s0 - 0.06, 0.16, color=PHASE[key], lw=0, alpha=0.9))
        ax.text((s0 + s1) / 2, 0.86, lab, ha="center", va="bottom", color=INK, fontsize=7.5)
    ax.plot([x0, x1], [0.30, 0.30], color=GRID, lw=0.8, zorder=0)
    ax.plot([2, 3.5], [0.30, 0.30], "v", color=BLUE, ms=6, zorder=3)
    ax.plot([4, 5, 6, 7, 8], [0.30] * 5, "s", color=ORANGE, ms=5, zorder=3)
    ax.text(8.85, 0.30, r"$\leftarrow$ proved unique minimiser", color=INK2, fontsize=7, va="center")
    for xv, lab in [(NUSTAR, r"$\nu^*\in(3.918364,\,3.918366)$"), (NU2, r"$\nu_2\approx8.606$")]:
        ax.plot([xv, xv], [0.05, 0.60], color=INK2, lw=0.6, ls=(0, (2, 2)))
        ax.text(xv, 0.08, " " + lab, ha="left", va="bottom", color=INK2, fontsize=7)
    ax.set_xlim(x0, x1); ax.set_ylim(0, 1.05); ax.set_xlabel(r"exponent $\nu$")
    ax.set_title(r"(b) minimiser of $T_\nu$ in $d=2$", loc="left", pad=2)
    save(fig, "fig1_overview")

# ---------------------------------------------------------------- Fig. 2 energy difference
def fig_energy():
    nus = np.linspace(2.2, 12, 160)
    d = np.array([triangle_sum_box(A_tau(SQUARE), (n,) * 3, 64) / triangle_sum_box(A_tau(HEX), (n,) * 3, 64) - 1 for n in nus]) * 100
    tr = v2b["track"]; kk = sorted(tr, key=float)
    xr = [float(k) for k in kk if float(k) > NU2 and float(k) <= 12]
    yr = [(tr[k]["rect_T"] / tr[k]["hex_T"] - 1) * 100 for k in kk if float(k) > NU2 and float(k) <= 12]
    fig, ax = plt.subplots(figsize=(COL, 2.2))
    ax.axhline(0, color=BLUE, lw=1.4)
    ax.plot(nus, d, color=ORANGE)
    ax.plot(xr, yr, "o", color=AQUA, ms=4, mec="white", mew=0.6, zorder=3)
    ax.text(12.1, 1.4, "hexagonal (reference)", color=BLUE, fontsize=7, va="bottom", ha="right")
    ax.text(6.6, -12.5, "square", color=ORANGE, fontsize=7, ha="left")
    ax.text(9.7, -39, "best rectangle", color=INK2, fontsize=7, ha="left")
    for xv, lab in [(NUSTAR, r"$\nu^*$"), (NU2, r"$\nu_2$")]:
        ax.axvline(xv, color=MUTED, lw=0.6, ls=(0, (2, 2))); ax.text(xv + 0.12, -60, lab, ha="left", va="bottom", color=INK2, fontsize=7)
    ax.set_xlim(2, 12.2); ax.set_ylim(-62, 8)
    ax.set_xlabel(r"exponent $\nu$"); ax.set_ylabel(r"$T_\nu/T_\nu(\mathrm{hex})-1$  (%)")
    save(fig, "fig_2d_energy_difference")

# ---------------------------------------------------------------- Fig. 3 stability
def fig_stability():
    st = v2["stability"]; ks = sorted(st, key=float)
    x = np.array([float(k) for k in ks])
    fig, ax = plt.subplots(figsize=(COL, 2.2))
    ax.axvspan(SQ_SPIN, HEX_SPIN, color=GRID, lw=0, zorder=0)
    ax.axhline(0, color=INK2, lw=0.6)
    ax.plot(x, [st[k]["hex_eigs"][0] for k in ks], color=BLUE)
    cs = v2b["curvature_square"]
    xs = list(x) + [c[0] for c in cs if c[0] > x.max()]
    ys = [st[k]["square_eigs"][0] for k in ks] + [min(c[1], c[2]) for c in cs if c[0] > x.max()]
    ax.plot(xs, ys, color=ORANGE)
    ax.text(3.15, 10.4, "hexagonal", color=BLUE, fontsize=7)
    ax.text(6.0, 8.6, "square", color=ORANGE, fontsize=7)
    ax.text((SQ_SPIN + HEX_SPIN) / 2, -9.3, "both\nlocal min.", ha="center", va="bottom", color=INK2, fontsize=6.5)
    ax.axvline(NU2, color=MUTED, lw=0.6, ls=(0, (2, 2))); ax.text(NU2 + 0.1, -9.3, r"$\nu_2$", color=INK2, fontsize=7, va="bottom")
    ax.set_xlim(2.1, 9.7); ax.set_ylim(-10.5, 12.5)
    ax.set_xlabel(r"exponent $\nu$"); ax.set_ylabel("smallest Hessian eigenvalue")
    save(fig, "fig_2d_stability")

# ---------------------------------------------------------------- Fig. 4 landscape (sequential, one hue)
def fig_landscape():
    xs = np.linspace(0, 0.5, 51); ys = np.linspace(0.85, 1.6, 76)
    cmap = plt.get_cmap("Blues_r")
    fig, axs = plt.subplots(1, 4, figsize=(FULL, 2.15), sharey=True, gridspec_kw=dict(wspace=0.12))
    for ax, nu, lab in zip(axs, [3.0, NUSTAR, 6.0, 10.0], [r"$\nu=3$", r"$\nu=\nu^*$", r"$\nu=6$", r"$\nu=10$"]):
        Z = np.full((len(ys), len(xs)), np.nan)
        for i, y in enumerate(ys):
            for j, x in enumerate(xs):
                if x * x + y * y >= 1: Z[i, j] = triangle_sum_box(A_tau(x + 1j * y), (nu,) * 3, 40)
        Z = np.log10(np.clip(Z / np.nanmin(Z) - 1, 1e-5, None))
        im = ax.pcolormesh(xs, ys, Z, cmap=cmap, vmin=-5, vmax=0, shading="auto", rasterized=True)
        th = np.linspace(np.pi / 3, np.pi / 2, 60); ax.plot(np.cos(th), np.sin(th), color=INK2, lw=0.6)
        ax.plot(0.5, np.sqrt(3) / 2, "o", ms=4, color="white", mec=INK, mew=0.6)
        ax.plot(0, 1, "s", ms=4, color="white", mec=INK, mew=0.6, clip_on=False)
        ax.grid(False); ax.set_title(lab, pad=2); ax.set_xlabel(r"Re $\tau$"); ax.set_xticks([0, 0.2, 0.4]); ax.set_xticklabels(["0", "0.2", "0.4"])
        ax.set_xlim(0, 0.5); ax.set_ylim(0.85, 1.6)
    axs[0].set_ylabel(r"Im $\tau$")
    cb = fig.colorbar(im, ax=axs, fraction=0.025, pad=0.015); cb.outline.set_linewidth(0.5)
    cb.set_label(r"$\log_{10}(T_\nu/\min T_\nu-1)$"); cb.ax.tick_params(width=0.5)
    save(fig, "fig_2d_landscape")

# ---------------------------------------------------------------- Fig. 5 unequal exponents
def fig_unequal():
    rows = json.load(open(os.path.join(E, "e3_unequal_exponents.json")))
    keys = ["hex", "square", "rect"]
    nus = sorted({r["nu"] for r in rows}); mus = sorted({r["mu"] for r in rows})
    Z = np.full((len(mus), len(nus)), np.nan)
    for r in rows: Z[mus.index(r["mu"]), nus.index(r["nu"])] = keys.index(r["phase"])
    fig, ax = plt.subplots(figsize=(COL, 2.7))
    ax.grid(False)
    ax.pcolormesh(nus, mus, Z, cmap=ListedColormap([BLUE, ORANGE, AQUA]), vmin=-0.5, vmax=2.5, shading="nearest", rasterized=True)
    ax.plot([2.5, 12], [2.5, 12], color=INK, lw=0.6, ls=(0, (2, 2)))
    ax.text(2.6, 11.3, "hexagonal", color="white", fontsize=7.5)
    ax.text(6.4, 6.9, "square", color="white", fontsize=7.5)
    ax.text(9.25, 10.3, "rect.", color="white", fontsize=7.5)
    ax.text(4.75, 5.55, r"$\mu=\nu$", color=INK, fontsize=7, rotation=45, ha="center", va="bottom")
    ax.set_xlabel(r"$\nu$  (exponent of $|x|$ and $|y|$)"); ax.set_ylabel(r"$\mu$  (exponent of $|x-y|$)")
    ax.set_xlim(2.25, 12.25); ax.set_ylim(2.25, 12.25); ax.set_aspect("equal")
    save(fig, "fig_e3_phase_map")

# ---------------------------------------------------------------- Fig. 6 mixed energies
def fig_mixed():
    d = json.load(open(os.path.join(E, "e8_lambda_c_dense.json")))
    x2 = [r["nu"] for r in d["2d"] if r["lam_c"] is not None]; y2 = [r["lam_c"] for r in d["2d"] if r["lam_c"] is not None]
    x3 = [r["nu"] for r in d["3d"]]; y3 = [r["lam_c"] for r in d["3d"]]
    fig, ax = plt.subplots(figsize=(COL, 2.3))
    ax.plot(x3, y3, "-o", color=BLUE, ms=3, mec="white", mew=0.5)
    ax.plot(x2, y2, "-s", color=ORANGE, ms=3, mec="white", mew=0.5)
    ax.axvline(NUSTAR, color=MUTED, lw=0.6, ls=(0, (2, 2)))
    ax.text(NUSTAR - 0.08, 12, r"$\nu^*$", color=INK2, fontsize=7, va="top", ha="right")
    ax.text(6.0, 0.42, r"$d=3$: FCC$\,\to\,$BCC", color=BLUE, fontsize=7, va="center", ha="left")
    ax.text(12.1, y2[-1] * 1.35, r"$d=2$: hex$\,\to\,$square", color=ORANGE, fontsize=7, va="bottom", ha="right")
    ax.set_yscale("log"); ax.set_xlim(3.3, 12.2); ax.set_ylim(0.04, 15)
    ax.set_xlabel(r"exponent $\nu$"); ax.set_ylabel(r"critical ratio $\lambda_c=c_3/c_2$")
    save(fig, "fig_e4_lambda_c")

# ---------------------------------------------------------------- Fig. asymptotics of the rectangular minimiser
def fig_asym():
    d = json.load(open(os.path.join(HERE, "..", "largenu", "asym_data.json")))
    yinf, kap = d["yinf"], d["kappa"]
    num = np.array(d["num"]); bal = np.array(d["bal"])
    fig, ax = plt.subplots(figsize=(COL, 2.3))
    ax.axhline(yinf, color=MUTED, lw=0.6, ls=(0, (2, 2)))
    ax.text(420, yinf + 0.004, r"$y_\infty$", color=INK2, fontsize=7, va="bottom", ha="right")
    nn = np.geomspace(9, 400, 200)
    ax.plot(nn, yinf - kap / nn, color=MUTED, lw=0.9, ls=(0, (4, 2)))
    ax.plot(bal[:, 0], bal[:, 1], color=AQUA, lw=1.4)
    ax.plot(num[:, 0], num[:, 1], "o", color=INK, ms=2.8, mec="white", mew=0.4, zorder=5)
    ax.set_xscale("log"); ax.set_xlim(8.5, 420); ax.set_ylim(1.02, 1.265)
    ax.set_xlabel(r"exponent $\nu$"); ax.set_ylabel(r"optimal aspect ratio $y_\nu$")
    ax.text(9.2, 1.236, "balance equation", color=AQUA, fontsize=7, ha="left", va="center")
    ax.text(40, 1.198, r"$y_\infty-\kappa/\nu$", color=INK2, fontsize=7, ha="left", va="top")
    ax.text(10.8, 1.058, "numerical minimisers", color=INK, fontsize=7, ha="left", va="center")
    ins = fig.add_axes([0.60, 0.25, 0.33, 0.30])
    ins.plot(num[:, 0], num[:, 0] * (yinf - num[:, 1]), "o", color=INK, ms=2.2, mec="white", mew=0.3)
    ins.plot(bal[:, 0], bal[:, 0] * (yinf - bal[:, 1]), color=AQUA, lw=1.0)
    ins.axhline(kap, color=MUTED, lw=0.6, ls=(0, (2, 2)))
    ins.set_xscale("log"); ins.set_xlim(15, 420); ins.set_ylim(0.9, 1.3)
    ins.set_title(r"$\nu(y_\infty-y_\nu)\to\kappa$", fontsize=6.5, pad=2)
    ins.tick_params(labelsize=5.5, length=2, which="both"); ins.grid(False)
    ins.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ins.set_xticks([20, 100]); ins.set_xticklabels(["20", "100"])
    save(fig, "fig_asym")

if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    fig1(); fig_energy(); fig_stability(); fig_landscape(); fig_unequal(); fig_mixed(); fig_asym()
    print("paper figures written")
