"""V4b: C_{1,1,1} with the correct truncation ansatz.

With all three exponents at nu = d = 2 the cube-truncation error is
R^-2 (a log R + b) + R^-3 (...) + ...  (the inner sum over the near vertex
diverges logarithmically), so the plain power-law Richardson of V4 is biased.
"""
import numpy as np, mpmath as mp
from epsteinlib import epstein_zeta
from lattice_sums import triangle_sum_box

z3 = float(mp.zeta(3))
for tau in [1j, 0.5 + 0.5j * np.sqrt(3), 0.21 + 1.37j]:
    A = np.array([[1.0, tau.real], [0.0, tau.imag]])
    Rs = np.array([96, 128, 192, 256, 384, 512, 768, 1024], float)
    vals = np.array([triangle_sum_box(A, (2, 2, 2), int(R)) for R in Rs])
    L = np.log(Rs)
    for nb, cols in [("power law only", [Rs**-2, Rs**-3, Rs**-4]),
                     ("with log terms", [Rs**-2 * L, Rs**-2, Rs**-3 * L, Rs**-3, Rs**-4 * L, Rs**-4])]:
        M = np.column_stack([np.ones_like(Rs)] + cols)
        c = np.linalg.lstsq(M, vals, rcond=None)[0][0]
        C111 = (tau.imag / np.pi) ** 3 * c
        E3 = (tau.imag / np.pi) ** 3 * float(epstein_zeta(6, A, np.zeros(2), np.zeros(2)).real)
        print("tau=%-20s %-15s C111=%.15f  E3+z3=%.15f  rel %.1e" % (np.round(tau, 4), nb, C111, E3 + z3, abs(C111 - E3 - z3) / (E3 + z3)))
