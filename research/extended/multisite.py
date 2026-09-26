"""Triangle sums for crystals with a basis (non-Bravais), e.g. hcp.

Crystal C = union_j (L + s_j), L = A Z^3, n sites per cell.  Energy per atom
  T_nu(C) = (1/n) sum_i sum_{p, q in C, p, q, s_i pairwise distinct} (|p - s_i| |q - s_i| |p - q|)^-nu.
For fixed i and site pair (j, k): p = A m + s_j, q = A n + s_k, p - q = A(m - n) + s_j - s_k, so the inner
double sum is sum_m F_ij(m) (F_ik * G_jk)(m) with F_ij(m) = f(A m + s_j - s_i), G_jk(l) = f(A l + s_j - s_k)
(f(0) := 0), evaluated with the same cube truncation + FFT convolution + Richardson as lattice_sums.py.
"""
import numpy as np
from scipy.fft import rfftn, irfftn, next_fast_len

def _labels(R, d):
    ax = np.arange(-R, R + 1, dtype=float)
    return np.stack(np.meshgrid(*([ax] * d), indexing="ij"), -1)

def _f(A, lab, shift, nu):
    X = lab @ A.T + shift
    r = np.sqrt(np.sum(X * X, -1))
    with np.errstate(divide="ignore"):
        f = r ** (-nu)
    f[r < 1e-12] = 0.0
    return f

def triangle_sum_multisite_box(A, sites, nu, R):
    A = np.asarray(A, float); d = A.shape[0]; n = len(sites)
    lab = _labels(R, d); M = next_fast_len(3 * R + 2); shape = (M,) * d
    sl = tuple(slice(R, 3 * R + 1) for _ in range(d))
    tot = 0.0
    G = {(j, k): rfftn(_f(A, lab, sites[j] - sites[k], nu), s=shape) for j in range(n) for k in range(n)}
    for i in range(n):
        F = [_f(A, lab, sites[j] - sites[i], nu) for j in range(n)]
        Fh = [rfftn(Fj, s=shape) for Fj in F]
        for j in range(n):
            for k in range(n):
                c = irfftn(Fh[k] * G[(j, k)], s=shape)[sl]
                tot += float(np.sum(F[j] * c))
    return tot / n

def triangle_sum_multisite(A, sites, nu, Rs=(16, 24, 32, 40, 48)):
    d = np.asarray(A).shape[0]; p = 2 * nu - d
    vals = np.array([triangle_sum_multisite_box(A, sites, nu, R) for R in Rs])
    Rs = np.array(Rs, float)
    Mx = np.column_stack([np.ones_like(Rs)] + [Rs ** (-(p + j)) for j in range(len(Rs) - 1)])
    return float(np.linalg.solve(Mx, vals)[0])

def epstein_multisite(A, sites, nu, Rs=None):
    """pair energy per atom: (1/n) sum_i sum_{p != s_i} |p - s_i|^-nu  (no factor 1/2), via epsteinlib
    (exact Epstein zeta with shift; the cube sums converge too slowly for nu close to d)."""
    from epsteinlib import epstein_zeta
    A = np.asarray(A, float); n = len(sites); d = A.shape[0]
    return float(sum(epstein_zeta(nu, A, -(np.asarray(sites[j]) - np.asarray(sites[i])), np.zeros(d)).real
                     for i in range(n) for j in range(n)) / n)

def _epstein_multisite_cube(A, sites, nu, Rs=(24, 32, 48, 64)):
    """old cube-sum version, kept for reference (biased for nu - d small)"""
    A = np.asarray(A, float); d = A.shape[0]; n = len(sites); p = nu - d
    vals = []
    for R in Rs:
        lab = _labels(R, d)
        vals.append(sum(np.sum(_f(A, lab, sites[j] - sites[i], nu)) for i in range(n) for j in range(n)) / n)
    Rs = np.array(Rs, float)
    Mx = np.column_stack([np.ones_like(Rs)] + [Rs ** (-(p + j)) for j in range(len(Rs) - 1)])
    return float(np.linalg.solve(Mx, np.array(vals))[0])

# ---- structures, all normalised to one atom per unit volume -----------------
def normalise(A, sites):
    A = np.asarray(A, float); n = len(sites)
    s = (abs(np.linalg.det(A)) / n) ** (1 / A.shape[0])
    return A / s, [np.asarray(x, float) / s for x in sites]

def fcc_prim():   return normalise(np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]], float), [np.zeros(3)])
def bcc_prim():   return normalise(np.array([[-1, 1, 1], [1, -1, 1], [1, 1, -1]], float), [np.zeros(3)])
def fcc_cubic():
    return normalise(2 * np.eye(3), [np.zeros(3), np.array([0, 1, 1.]), np.array([1, 0, 1.]), np.array([1, 1, 0.])])
def hcp(c_over_a=np.sqrt(8 / 3)):
    a = 1.0; c = c_over_a * a
    A = np.array([[a, -a / 2, 0], [0, a * np.sqrt(3) / 2, 0], [0, 0, c]])
    return normalise(A, [np.zeros(3), A @ np.array([1 / 3, 2 / 3, 1 / 2])])
