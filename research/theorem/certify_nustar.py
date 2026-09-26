"""Computer-assisted enclosure of the hex/square crossing nu^* (2D, three-body energy T_nu).

T = S_R + Rem_R with S_R = sum over labels x, y in the cube [-R,R]^2 (x - y unrestricted), computed by a
DIRECT (non-FFT) convolution, and 0 <= Rem_R <= B_R with the explicit tail bound
    B_R = 2^(nu+2) * Zup * 8 s^(-2nu) R^(2-2nu) / (2nu - 2),
where s = smallest singular value of the basis matrix and Zup >= sum_{x != 0} |x|^-nu
(itself bounded by a cube sum plus the analogous tail 8 s^-nu R^(2-nu)/(nu-2)).
Derivation: a term is missing only if x or y lies outside the cube; for fixed x,
sum_y |y|^-nu |x-y|^-nu <= 2 (|x|/2)^-nu Z, and |A m| >= s |m|_inf.
Floating point: all summands are positive; with sequential summation of N terms the relative error
is <= N*u (u = 2^-53), plus <= 20u per summand for the power evaluations; we use a relative
margin EPS = 1e-9, which exceeds these bounds by more than 100x for the sizes used here.
"""
import sys, json
import numpy as np
from scipy.signal import convolve2d
sys.path.insert(0, "../verification")
from lattice_sums import A_tau, HEX, SQUARE

EPS = 1e-9
def table(A, nu, R):
    ax = np.arange(-R, R + 1, dtype=float)
    M1, M2 = np.meshgrid(ax, ax, indexing="ij")
    X = A[0, 0] * M1 + A[0, 1] * M2; Y = A[1, 0] * M1 + A[1, 1] * M2
    r = np.sqrt(X * X + Y * Y)
    with np.errstate(divide="ignore"): f = r ** (-nu)
    f[R, R] = 0.0
    return f

def enclose(A, nu, R=100, RZ=400):
    s = np.linalg.svd(A, compute_uv=False).min()
    f = table(A, nu, R)
    F2 = table(A, nu, 2 * R)                     # f on [-2R, 2R]^2 for the difference vectors
    c = convolve2d(f, F2, mode="valid")          # c(x) = sum_y f(y) f(x - y) for x in the cube? (see check below)
    S = float(np.sum(f * c))
    Z = float(np.sum(table(A, nu, RZ)))
    Zup = (Z + 8 * s ** (-nu) * RZ ** (2 - nu) / (nu - 2)) * (1 + EPS)
    B = 2 ** (nu + 2) * Zup * 8 * s ** (-2 * nu) * R ** (2 - 2 * nu) / (2 * nu - 2)
    return S * (1 - EPS), S * (1 + EPS) + B, B

if __name__ == "__main__":
    # sanity check of the index convention against a brute-force loop on a tiny cube
    A = A_tau(0.21 + 1.37j); nu = 4.0; R = 4
    f = table(A, nu, R); F2 = table(A, nu, 2 * R); c = convolve2d(f, F2, mode="valid")
    lab = [(i, j) for i in range(-R, R + 1) for j in range(-R, R + 1)]
    brute = 0.0
    for (a, b) in lab:
        for (p, q) in lab:
            fx = f[a + R, b + R]; fy = f[p + R, q + R]; fd = F2[a - p + 2 * R, b - q + 2 * R]
            brute += fx * fy * fd
    print("index check: direct-convolution sum %.15g vs brute force %.15g" % (float(np.sum(f * c)), brute))
    out = {}
    for nu in (3.918364, 3.918366):
        lo_sq, hi_sq, Bs = enclose(A_tau(SQUARE), nu)
        lo_hx, hi_hx, Bh = enclose(A_tau(HEX), nu)
        dlo, dhi = lo_sq - hi_hx, hi_sq - lo_hx
        out[nu] = dict(square=[lo_sq, hi_sq], hex=[lo_hx, hi_hx], D=[dlo, dhi], tail_square=Bs, tail_hex=Bh)
        print("nu=%.6f  T(square) in [%.12f, %.12f]  T(hex) in [%.12f, %.12f]  D = T(sq)-T(hex) in [%.3e, %.3e]  (tails %.1e, %.1e)"
              % (nu, lo_sq, hi_sq, lo_hx, hi_hx, dlo, dhi, Bs, Bh))
    ok = out[3.918364]["D"][0] > 0 and out[3.918366]["D"][1] < 0
    print("CERTIFIED: hex and square exchange order at some nu in (3.918364, 3.918366):", ok)
    json.dump(out, open("certify_nustar.json", "w"), indent=1)
