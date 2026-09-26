"""Rigorous Taylor coefficients of T_nu around tau0 = x0 + i y0 (variables A = a/y0, B = b/y0).

* side jets phi_v (v in a cube) are computed in interval arithmetic (iv.IV); phi0 = l0^-s via Arb;
* the triangle sums  sum_x J_a(x) sum_y J_b(y) J_c(x - y)  are evaluated with scipy's direct (non-FFT)
  convolution on the interval midpoints, and enclosed using
    |fl(S) - S_mid| <= gamma_N * sum |.|    (Higham, any summation order, N = number of roundings)
    |S_exact - S_mid| <= sum (|a|+ra)(|b|+rb)(|c|+rc) - sum |a||b||c|,
  where the absolute sums are evaluated in the same way and inflated by (1 + gamma_N);
* the truncation tail is bounded by B_R * Ghat_ij (universal majorant, see ../theorem_mgf/jets.py), with
  B_R = 2^(nu+2) Z * 8 smin^(-2nu) R^(2-2nu) / (2nu-2) and Z <= smin^-nu * Z_{Z^2}(nu), all in Arb.
"""
import numpy as np
from scipy.signal import convolve2d
from flint import arb
from iv import IV, up, down, pow_neg_arb, gamma

# ---------------------------------------------------------------- majorant (Arb)
def binom_gen_arb(s, j):
    out = arb(1)
    for i in range(j): out = out * (arb(s) + i) / (i + 1)
    return out

def majorant_coeffs(expo, K):
    """coefficients of (1 - eps_hat)^-expo, eps_hat = (A + B + A^2 + B^2)/(1 - B), as Arb (exact rationals up to rounding)"""
    def smul(p, q):
        r = [[arb(0)] * (K + 1) for _ in range(K + 1)]
        for i in range(K + 1):
            for j in range(K + 1 - i):
                acc = arb(0)
                for k in range(i + 1):
                    for l in range(j + 1): acc += p[k][l] * q[i - k][j - l]
                r[i][j] = acc
        return r
    num = [[arb(0)] * (K + 1) for _ in range(K + 1)]
    num[1][0] = num[0][1] = arb(1)
    if K >= 2: num[2][0] = num[0][2] = arb(1)
    geo = [[arb(1) if i == 0 else arb(0) for j in range(K + 1)] for i in range(K + 1)]
    eps = smul(num, geo)
    res = [[arb(1) if (i, j) == (0, 0) else arb(0) for j in range(K + 1)] for i in range(K + 1)]
    pw = [row[:] for row in res]
    for j in range(1, K + 1):
        pw = smul(pw, eps); b = binom_gen_arb(expo, j)
        for i1 in range(K + 1):
            for j1 in range(K + 1 - i1): res[i1][j1] += b * pw[i1][j1]
    return res

def majorant_value(expo, A, B):
    A = arb(A); B = arb(B)
    e = (A + B + A * A + B * B) / (1 - B)
    assert e.upper() < 1 and B.upper() < 1
    return (1 - e) ** (-arb(expo))

def majorant_tail(expo, A, B, K):
    """sum_{i+j>K} Ghat_ij A^i B^j (Arb upper bound)"""
    G = majorant_coeffs(expo, K)
    low = sum((G[i][j] * arb(A) ** i * arb(B) ** j for i in range(K + 1) for j in range(K + 1 - i)), arb(0))
    return majorant_value(expo, A, B) - low

# ---------------------------------------------------------------- side jets in interval arithmetic
def smul_iv(p, q, K):
    r = IV(np.zeros(p.shape), np.zeros(p.shape))
    for i in range(K + 1):
        for j in range(K + 1 - i):
            acc = IV(np.zeros(p.shape[:-2]), np.zeros(p.shape[:-2]))
            for k in range(i + 1):
                for l in range(j + 1):
                    acc = acc + p[..., k, l] * q[..., i - k, j - l]
            r[..., i, j] = acc
    return r

def side_jets(x0, y0, s, R, K):
    ax = np.arange(-R, R + 1, dtype=float)
    m, n = np.meshgrid(ax, ax, indexing="ij")
    X0, Y0 = IV.exact(np.full(m.shape, x0)), IV.exact(np.full(m.shape, y0))
    u = IV.exact(m) + IV.exact(n) * X0
    ny = IV.exact(n) * Y0
    z2 = u.sq() + ny.sq()
    zero = (m == 0) & (n == 0)
    z2.lo[zero] = 1.0; z2.hi[zero] = 1.0
    l0 = z2 / Y0
    al2 = ny.sq() / z2
    ab = (ny * u) / z2
    shp = m.shape + (K + 1, K + 1)
    num = IV(np.zeros(shp), np.zeros(shp))
    num[..., 0, 1] = al2 * 2.0 - 1.0
    num[..., 1, 0] = ab * 2.0
    if K >= 2: num[..., 2, 0] = al2; num[..., 0, 2] = al2
    inv = IV(np.zeros(shp), np.zeros(shp))
    for j in range(K + 1):
        inv.lo[..., 0, j] = (-1.0) ** j; inv.hi[..., 0, j] = (-1.0) ** j
    eps = smul_iv(num, inv, K)
    res = IV(np.zeros(shp), np.zeros(shp)); res.lo[..., 0, 0] = 1.0; res.hi[..., 0, 0] = 1.0
    pw = res.copy()
    for j in range(1, K + 1):
        pw = smul_iv(pw, eps, K)
        bj = arb(1)
        for i in range(j): bj = bj * (-arb(s) - i) / (i + 1)
        blo, bhi = down(float(bj.lower())), up(float(bj.upper()))
        res = res + pw * IV(np.full(shp, blo), np.full(shp, bhi))
    phi0 = pow_neg_arb(l0, s)
    phi0.lo[zero] = 0.0; phi0.hi[zero] = 0.0
    out = IV(np.zeros(shp), np.zeros(shp))
    for i in range(K + 1):
        for j in range(K + 1 - i):
            out[..., i, j] = res[..., i, j] * phi0
    return out

# ---------------------------------------------------------------- rigorous triangle-sum jets
_ZSQ = {}
def zsq_upper(nu, R=200):
    """rigorous upper bound of sum_{m in Z^2, m != 0} |m|^-nu"""
    if nu in _ZSQ: return _ZSQ[nu]
    ax = np.arange(-R, R + 1, dtype=float); M1, M2 = np.meshgrid(ax, ax, indexing="ij")
    r2 = IV.exact(M1 * M1 + M2 * M2); r2.lo[R, R] = r2.hi[R, R] = 1.0
    f = pow_neg_arb(r2, arb(nu) / 2); f.lo[R, R] = f.hi[R, R] = 0.0
    N = f.hi.size
    S = arb(float(np.sum(f.hi))) * (1 + arb(float(gamma(N)))) + 8 * arb(R) ** (2 - arb(nu)) / (arb(nu) - 2)
    _ZSQ[nu] = S; return S

def tail_bound(x0, y0, nu, R):
    x0, y0, nu = arb(x0), arb(y0), arb(nu)
    # A = [[1, x0], [0, y0]] / sqrt(y0);  s_min^2 = (t - sqrt(t^2 - 4)) / 2, t = trace(A^T A), det(A)=1
    t = (1 + x0 * x0 + y0 * y0) / y0
    smin2 = (t - (t * t - 4).sqrt()) / 2
    smin = smin2.sqrt()
    Z = smin ** (-nu) * zsq_upper(float(nu.mid()) if nu.rad() == 0 else nu)
    return 2 ** (nu + 2) * Z * 8 * smin ** (-2 * nu) * arb(R) ** (2 - 2 * nu) / (2 * nu - 2)

def T_jet(x0, y0, nu, R=40, K=2, return_trunc=False):
    """returns (coef, err) as lists of Arb balls: coefficient (i,j) of the Taylor expansion of T_nu in (A,B)"""
    s = arb(nu) / 2
    J = side_jets(x0, y0, s, 2 * R, K)
    Jc = J[R:3 * R + 1, R:3 * R + 1]
    idx = [(i, j) for i in range(K + 1) for j in range(K + 1 - i)]
    mid = {k: J.mid()[..., k[0], k[1]] for k in idx}; rad = {k: J.rad()[..., k[0], k[1]] for k in idx}
    midc = {k: v[R:3 * R + 1, R:3 * R + 1] for k, v in mid.items()}; radc = {k: v[R:3 * R + 1, R:3 * R + 1] for k, v in rad.items()}
    N = (2 * R + 1) ** 2 * 2 + 8
    g = float(gamma(N))
    coef = {k: arb(0) for k in idx}
    cache = {}
    for a in idx:
        for b in idx:
            for c in idx:
                i, j = a[0] + b[0] + c[0], a[1] + b[1] + c[1]
                if i + j > K: continue
                key = (b, c)
                if key not in cache:
                    cm = convolve2d(midc[b], mid[c], mode="valid")
                    ca = convolve2d(np.abs(midc[b]), np.abs(mid[c]), mode="valid")
                    cr = convolve2d(np.abs(midc[b]) + radc[b], np.abs(mid[c]) + rad[c], mode="valid")
                    cache[key] = (cm, ca, cr)
                cm, ca, cr = cache[key]
                S = float(np.sum(midc[a] * cm))
                Sabs = float(np.sum(np.abs(midc[a]) * ca))
                Srad = float(np.sum((np.abs(midc[a]) + radc[a]) * cr))
                # |S_exact - S| <= gamma*Sabs_exact + (Srad_exact - Sabs_exact);  Sabs_exact <= Sabs(1+g), Srad_exact <= Srad(1+g)
                err = arb(g) * arb(Sabs) * (1 + arb(g)) + (arb(Srad) * (1 + arb(g)) - arb(Sabs) * (1 - arb(g)))
                coef[(i, j)] += arb(S) + arb(0, 1) * err.upper()      # ball: S +- err  (arb(0, r) = [-r, r])
    BR = tail_bound(x0, y0, nu, R)
    G = majorant_coeffs(3 * s, K)
    out = {}
    for (i, j) in idx:
        out[(i, j)] = coef[(i, j)] + arb(0, 1) * (BR * G[i][j]).upper()
    if return_trunc:
        return out, BR, coef          # coef = enclosures of the cube-truncated sums (no tail)
    return out, BR
