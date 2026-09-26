"""2D steep-decay asymptotics: minimiser of T_nu on the imaginary axis vs the prediction y_inf - kappa/nu."""
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.signal import fftconvolve

def T(x, y, nu, R=5):
    """direct sum over x, y in the cube [-R,R]^2 (tail negligible for nu >= 10), scaled by P*^nu"""
    ax = np.arange(-R, R + 1); m, n = [v.ravel() for v in np.meshgrid(ax, ax, indexing="ij")]
    keep = (m != 0) | (n != 0); m, n = m[keep], n[keep]
    X = (m + n * x) / np.sqrt(y); Y = n * np.sqrt(y)
    r = np.hypot(X, Y)
    d = np.hypot(X[:, None] - X[None, :], Y[:, None] - Y[None, :])
    np.fill_diagonal(d, np.inf)
    L = np.log(r)[:, None] + np.log(r)[None, :] + np.log(d) - np.log(1.4317344)
    return float(np.sum(np.exp(-nu * L)))

def T_old(x, y, nu, R=8):
    ax = np.arange(-R, R + 1); m, n = np.meshgrid(ax, ax, indexing="ij")
    X = (m + n * x) / np.sqrt(y); Y = n * np.sqrt(y)
    r = np.hypot(X, Y); r[R, R] = np.inf
    ax2 = np.arange(-2 * R, 2 * R + 1); m2, n2 = np.meshgrid(ax2, ax2, indexing="ij")
    r2 = np.hypot((m2 + n2 * x) / np.sqrt(y), n2 * np.sqrt(y)); r2[2 * R, 2 * R] = np.inf
    # log-scale to avoid underflow: factor out P*^-nu
    s = 1.431734
    f = (r / s ** (1 / 3)) ** (-nu); F = (r2 / s ** (1 / 3)) ** (-nu)
    c = fftconvolve(f, F[::-1, ::-1], mode="valid")  # c(x) = sum_y f(y) F(y - x)
    return float(np.sum(f * c))

a = ((1 + np.sqrt(17)) / 8) ** 0.25; yinf = a ** -2
F1p = 3 / a; F2p = (a - a ** -3) / (a ** 2 + a ** -2)
dz = np.log(4 * (-F2p) / F1p) / (F2p - F1p)          # a_nu - a* ~ dz/nu
kappa = 2 * a ** -3 * dz
print("y_inf = %.6f  kappa = %.6f" % (yinf, kappa))
for nu in [10, 15, 20, 30, 50, 80, 120, 200, 300]:
    res = minimize_scalar(lambda y: T(0.0, y, nu), bounds=(1.0, 1.3), method="bounded", options={"xatol": 1e-10})
    y = res.x
    h = 1e-4; d2x = (T(h, y, nu) - 2 * T(0, y, nu) + T(-h, y, nu)) / h ** 2
    print("nu=%5d  y_nu=%.7f  pred=%.7f  nu*(y_inf-y)=%.5f  d2T/dx2>0: %s" % (nu, y, yinf - kappa / nu, nu * (yinf - y), d2x > 0))
