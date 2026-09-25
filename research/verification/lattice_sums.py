"""Independent reference implementation of triangle lattice sums (no GZL code).

T(A; nu1, nu2, nu3) = sum_{x, y in Lambda, x != 0, y != 0, x != y}
                      |x|^-nu1 * |y|^-nu2 * |x - y|^-nu3,     Lambda = A Z^d.

This is the graph zeta function of the 3-cycle (triangle) with one vertex
pinned at the origin.  It is evaluated by direct summation over all integer
labels in the cube [-R, R]^d (all three vertices inside the cube), using an
FFT for the inner convolution, followed by Richardson extrapolation in R.
The truncation error of the cube sum is dominated by configurations with one
vertex far away and scales like R^(d - nu_i - nu_j).
"""
import numpy as np
from scipy.fft import rfftn, irfftn, next_fast_len


def _norms(A, R):
    d = A.shape[0]
    ax = np.arange(-R, R + 1, dtype=np.float64)
    grids = np.meshgrid(*([ax] * d), indexing="ij")
    labels = np.stack(grids, axis=-1)            # (2R+1,)*d + (d,)
    X = labels @ A.T                              # real-space vectors A m
    return np.sqrt(np.sum(X * X, axis=-1))


def _kernel(r, nu):
    with np.errstate(divide="ignore"):
        f = r ** (-float(nu))
    f[r == 0] = 0.0
    return f


def triangle_sum_box(A, nus, R):
    """Cube-truncated triangle sum with all vertex labels in [-R, R]^d."""
    A = np.asarray(A, dtype=np.float64)
    d = A.shape[0]
    r = _norms(A, R)
    f1, f2, f3 = (_kernel(r, nu) for nu in nus)
    n = 2 * R + 1
    M = next_fast_len(3 * R + 2)                 # >3R avoids aliasing onto the box
    shape = (M,) * d
    # c(x) = sum_y f2(y) f3(x - y), linear convolution, both supported in the box
    # place arrays with index 0 <-> label -R, then shift result by R
    c = irfftn(rfftn(f2, s=shape) * rfftn(f3, s=shape), s=shape)
    # conv index j corresponds to label (j - 2R); we need labels -R..R -> j = R..3R
    sl = tuple(slice(R, 3 * R + 1) for _ in range(d))
    c_box = c[sl]
    return float(np.sum(f1 * c_box))


def triangle_sum(A, nus, Rs=(32, 48, 64, 96), return_table=False):
    """Richardson-extrapolated triangle sum.

    Fits T(R) = T_inf + sum_k c_k R^-(p+k), k = 0..len(Rs)-2, with
    p = min_{i<j}(nu_i + nu_j) - d.
    """
    A = np.asarray(A, dtype=np.float64)
    d = A.shape[0]
    nus = [float(v) for v in nus]
    p = min(nus[0] + nus[1], nus[0] + nus[2], nus[1] + nus[2]) - d
    vals = np.array([triangle_sum_box(A, nus, R) for R in Rs])
    Rs_arr = np.array(Rs, dtype=np.float64)
    k = len(Rs) - 1
    M = np.column_stack([np.ones_like(Rs_arr)] + [Rs_arr ** (-(p + j)) for j in range(k)])
    coef = np.linalg.solve(M, vals)
    if return_table:
        return coef[0], vals
    return float(coef[0])


def triangle_sum_naive(A, nus, R):
    """O(N^2) pair loop over the same cube, for testing triangle_sum_box."""
    A = np.asarray(A, dtype=np.float64)
    d = A.shape[0]
    ax = np.arange(-R, R + 1)
    labels = np.stack(np.meshgrid(*([ax] * d), indexing="ij"), -1).reshape(-1, d)
    labels = labels[np.any(labels != 0, axis=1)]
    X = labels @ A.T
    r = np.linalg.norm(X, axis=1)
    D = np.linalg.norm(X[:, None, :] - X[None, :, :], axis=-1)
    np.fill_diagonal(D, np.inf)
    # same truncation as the FFT version: the label difference must lie in the cube too
    far = np.max(np.abs(labels[:, None, :] - labels[None, :, :]), axis=-1) > R
    D[far] = np.inf
    return float(np.sum(np.outer(r ** -nus[0], r ** -nus[1]) * D ** -nus[2]))


def epstein_direct(A, nu, Rs=(64, 128, 256)):
    """sum_{x != 0} |x|^-nu over Lambda = A Z^d by cube sums + Richardson (nu > d)."""
    A = np.asarray(A, dtype=np.float64)
    d = A.shape[0]
    p = nu - d
    vals = np.array([np.sum(_kernel(_norms(A, R), nu)) for R in Rs])
    Rs_arr = np.array(Rs, dtype=np.float64)
    k = len(Rs) - 1
    M = np.column_stack([np.ones_like(Rs_arr)] + [Rs_arr ** (-(p + j)) for j in range(k)])
    return float(np.linalg.solve(M, vals)[0])


# ---- lattices ---------------------------------------------------------------

def A_tau(tau):
    """2D lattice Z + tau Z scaled to unit covolume (columns = basis vectors)."""
    tau = complex(tau)
    return np.array([[1.0, tau.real], [0.0, tau.imag]]) / np.sqrt(tau.imag)


HEX = 0.5 + 0.5j * np.sqrt(3.0)
SQUARE = 1j

FCC = np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]], dtype=float)
BCC = np.array([[-1, 1, 1], [1, -1, 1], [1, 1, -1]], dtype=float)
SC = np.eye(3)


def unit_covolume(A):
    A = np.asarray(A, dtype=np.float64)
    return A / abs(np.linalg.det(A)) ** (1.0 / A.shape[0])
