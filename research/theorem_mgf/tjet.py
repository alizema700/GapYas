"""Rigorous Taylor coefficients of T_nu around tau0 = x0 + i y0 (in A = a/y0, B = b/y0):
enclosure = computed value +- err_ij, err_ij = (B_R + ROUND * T_R) * Ghat_ij, where
  B_R  = explicit tail bound of the cube-truncated triangle sum (x or y outside the cube),
  Ghat = coefficients of the universal majorant (1 - eps_hat)^-(3s),
  ROUND = 1e-10 >> (N1 + N2) * 2^-53 for sequential double sums (direct convolution, no FFT)."""
import numpy as np
from scipy.signal import convolve2d
from jets import side_jets, majorant_eps, majorant_power

ROUND = 1e-10

def tail_bound(x0, y0, nu, R, RZ=300):
    A = np.array([[1.0, x0], [0.0, y0]]) / np.sqrt(y0)
    sv = np.linalg.svd(A, compute_uv=False).min()
    ax = np.arange(-RZ, RZ + 1, dtype=float); M1, M2 = np.meshgrid(ax, ax, indexing="ij")
    r = np.sqrt((A[0, 0] * M1 + A[0, 1] * M2) ** 2 + (A[1, 1] * M2) ** 2); r[RZ, RZ] = np.inf
    Z = float(np.sum(r ** (-nu)))
    Zup = (Z + 8 * sv ** (-nu) * RZ ** (2 - nu) / (nu - 2)) * (1 + 1e-9)
    return 2 ** (nu + 2) * Zup * 8 * sv ** (-2 * nu) * R ** (2 - 2 * nu) / (2 * nu - 2)

def T_jet(x0, y0, nu, R=40, K=2):
    s = nu / 2
    J = side_jets(x0, y0, s, 2 * R, K)                  # on [-2R, 2R]^2
    Jc = J[R:3 * R + 1, R:3 * R + 1]                      # on [-R, R]^2
    idx = [(i, j) for i in range(K + 1) for j in range(K + 1 - i)]
    conv = {}
    coef = np.zeros((K + 1, K + 1))
    for a in idx:
        for b in idx:
            for c in idx:
                i, j = a[0] + b[0] + c[0], a[1] + b[1] + c[1]
                if i + j > K: continue
                key = (b, c)
                if key not in conv:
                    conv[key] = convolve2d(Jc[..., b[0], b[1]], J[..., c[0], c[1]], mode="valid")
                coef[i, j] += float(np.sum(Jc[..., a[0], a[1]] * conv[key]))
    BR = tail_bound(x0, y0, nu, R)
    Ghat = majorant_power(majorant_eps(K), 3 * s, K)
    err = (BR + ROUND * coef[0, 0]) * Ghat
    return coef, err, BR

if __name__ == "__main__":
    import time, sys
    sys.path.insert(0, "../verification")
    from lattice_sums import triangle_sum, A_tau
    for (x0, y0) in [(0.0, 1.0), (0.2, 1.1), (0.5, np.sqrt(3) / 2)]:
        t = time.time(); c, e, BR = T_jet(x0, y0, 4.0, R=40, K=2)
        ref = triangle_sum(A_tau(x0 + 1j * y0), (4, 4, 4), Rs=(48, 64, 96, 128))
        h = 1e-4
        d = (triangle_sum(A_tau(x0 + h + 1j * y0), (4,) * 3, Rs=(48, 64, 96, 128)) - triangle_sum(A_tau(x0 - h + 1j * y0), (4,) * 3, Rs=(48, 64, 96, 128))) / (2 * h) * y0
        print("tau=%.3f+%.3fi  T=%.12f (ref %.12f, err %.1e)  c10=%.8f (fd %.8f)  c01=%.6f c20=%.6f c11=%.6f c02=%.6f  tail %.1e  %.1fs"
              % (x0, y0, c[0, 0], ref, e[0, 0], c[1, 0], d, c[0, 1], c[2, 0], c[1, 1], c[0, 2], BR, time.time() - t))
