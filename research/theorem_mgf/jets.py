"""Rigorous Taylor data for T_nu(tau) = sum' (|x||y||x-y|)^-nu on the unit-covolume lattice Lambda_tau.

Side function for a lattice vector v = (m, n):  phi_v(x, y) = l_v^-s,  s = nu/2,
    l_v = |m + n tau|^2 / y = ((m + n x)^2 + n^2 y^2) / y,  tau = x + i y.
Around (x0, y0) put A = a/y0, B = b/y0 (a, b = increments of x, y). With alpha = n y0/|z|, beta = u/|z|
(z = m + n tau0, u = m + n x0, alpha^2 + beta^2 = 1):
    l/l0 = (1 + 2 alpha beta A + alpha^2 A^2 + 2 alpha^2 B + alpha^2 B^2) / (1 + B)
         = 1 + [(2 alpha^2 - 1) B + 2 alpha beta A + alpha^2 (A^2 + B^2)] / (1 + B).
MAJORANT (coefficientwise, all |.| <= 1):  l/l0 - 1  <<  eps_hat(A,B) = (A + B + A^2 + B^2) / (1 - B).
Hence phi/phi0 << (1 - eps_hat)^-s and every triangle term t/t0 << G(A,B) = (1 - eps_hat)^-(3s),
uniformly in the lattice vectors.  This gives rigorous remainder and tail bounds.
"""
import numpy as np
from math import comb

def binom_gen(s, j):  # |binom(-s, j)| = binom(s + j - 1, j) for real s > 0
    out = 1.0
    for i in range(j): out *= (s + i) / (i + 1)
    return out

# ---- bivariate truncated power series: coefficient arrays c[i, j] for A^i B^j, i + j <= K ----
def smul(p, q, K):
    r = np.zeros(p.shape[:-2] + (K + 1, K + 1))
    for i in range(K + 1):
        for j in range(K + 1 - i):
            for k in range(i + 1):
                for l in range(j + 1):
                    r[..., i, j] += p[..., k, l] * q[..., i - k, j - l]
    return r

def majorant_eps(K):
    num = np.zeros((K + 1, K + 1)); num[1, 0] = num[0, 1] = 1; 
    if K >= 2: num[2, 0] = num[0, 2] = 1
    geo = np.zeros((K + 1, K + 1))
    for j in range(K + 1): geo[0, j] = 1.0            # 1/(1 - B)
    return smul(num, geo, K)

def majorant_power(eps, expo, K):
    """coefficients of (1 - eps)^-expo for a majorant series eps with zero constant term"""
    res = np.zeros((K + 1, K + 1)); res[0, 0] = 1.0
    pw = np.zeros((K + 1, K + 1)); pw[0, 0] = 1.0
    for j in range(1, K + 1):
        pw = smul(pw, eps, K)
        res += binom_gen(expo, j) * pw
    return res

def majorant_value(expo, Aabs, Babs):
    """G(A,B) = (1 - eps_hat)^-expo evaluated at nonnegative (A,B) (closed form), requires eps_hat < 1, B < 1"""
    e = (Aabs + Babs + Aabs**2 + Babs**2) / (1 - Babs)
    assert e < 1 and Babs < 1
    return (1 - e) ** (-expo)

def majorant_tail(expo, Aabs, Babs, K, coef=None):
    """sum_{i+j > K} Ghat_ij A^i B^j = G(A,B) - sum_{i+j<=K} Ghat_ij A^i B^j (all terms nonnegative)"""
    if coef is None: coef = majorant_power(majorant_eps(K), expo, K)
    low = sum(coef[i, j] * Aabs**i * Babs**j for i in range(K + 1) for j in range(K + 1 - i))
    return majorant_value(expo, Aabs, Babs) - low

# ---- exact side jets for all lattice vectors in a cube -------------------------------------
def side_jets(x0, y0, s, R, K):
    """coefficients (in A = a/y0, B = b/y0) of phi_v(x0 + a, y0 + b) for all v = (m, n) in [-R, R]^2.
    Returns array (2R+1, 2R+1, K+1, K+1) (index [m+R, n+R]); v = 0 gets 0."""
    ax = np.arange(-R, R + 1, dtype=float)
    m, n = np.meshgrid(ax, ax, indexing="ij")
    u = m + n * x0
    z2 = u * u + n * n * y0 * y0
    zero = z2 == 0
    z2s = np.where(zero, 1.0, z2)
    l0 = z2s / y0
    al2 = n * n * y0 * y0 / z2s                 # alpha^2
    ab = n * y0 * u / z2s                       # alpha * beta
    # eps = [(2 al2 - 1) B + 2 ab A + al2 (A^2 + B^2)] / (1 + B)
    num = np.zeros(m.shape + (K + 1, K + 1))
    num[..., 0, 1] = 2 * al2 - 1
    num[..., 1, 0] = 2 * ab
    if K >= 2:
        num[..., 2, 0] = al2; num[..., 0, 2] = al2
    inv = np.zeros((K + 1, K + 1))
    for j in range(K + 1): inv[0, j] = (-1.0) ** j
    eps = smul(num, np.broadcast_to(inv, num.shape).copy(), K)
    # (1 + eps)^-s = sum_j binom(-s, j) eps^j
    res = np.zeros_like(eps); res[..., 0, 0] = 1.0
    pw = np.zeros_like(eps); pw[..., 0, 0] = 1.0
    for j in range(1, K + 1):
        pw = smul(pw, eps, K)
        b = 1.0
        for i in range(j): b *= (-s - i) / (i + 1)
        res += b * pw
    phi0 = np.where(zero, 0.0, l0 ** (-s))
    return res * phi0[..., None, None]
