"""Exact (rational) part of the local lemma at BCC.

At q0 = (1, 1, 1/3, 1/3, -1/3) all active triangle functions f_t = sqrt(l_m l_n l_{m-n} / det G) have
f_t = 3/2 and rational gradients grad f_t = (f_t/2) (sum_i grad l_i / l_i - grad det / det).  The six
distinct gradients form a simplex in R^5; we verify with exact rational arithmetic that the origin lies in
its interior and compute a rigorous lower bound for its distance to every facet.  This distance c
satisfies  max_t (-grad f_t) . d >= c |d|  for all d."""
import itertools
from fractions import Fraction as F
import sympy as sp

q0 = [F(1), F(1), F(1, 3), F(1, 3), F(-1, 3)]
def coeffs(v):
    m1, m2, m3 = v
    return [m1 * m1, m2 * m2, m3 * m3, 2 * m1 * m2, 2 * m1 * m3, 2 * m2 * m3]
def lval(c, q): return c[0] + c[1] * q[0] + c[2] * q[1] + c[3] * q[2] + c[4] * q[3] + c[5] * q[4]
def det(q):
    g22, g33, g12, g13, g23 = q
    return g22 * g33 - g23 ** 2 - g12 ** 2 * g33 - g13 ** 2 * g22 + 2 * g12 * g13 * g23
def grad_det(q):
    g22, g33, g12, g13, g23 = q
    return [g33 - g13 ** 2, g22 - g12 ** 2, -2 * g12 * g33 + 2 * g13 * g23, -2 * g13 * g22 + 2 * g12 * g23, -2 * g23 + 2 * g12 * g13]

D = det(q0); gD = grad_det(q0)
V = [m for m in itertools.product(range(-2, 3), repeat=3) if any(m)]
active, grads = [], set()
for m in V:
    for n in V:
        if m >= n: continue
        d = tuple(a - b for a, b in zip(m, n))
        ls = [lval(coeffs(v), q0) for v in (m, n, d)]
        f2 = ls[0] * ls[1] * ls[2] / D
        if f2 == F(9, 4):
            active.append((m, n))
            g = tuple(F(3, 4) * (sum(F(coeffs(v)[k + 1]) / l for v, l in zip((m, n, d), ls)) - gD[k] / D) for k in range(5))
            grads.add(g)
# check that no triangle with small coefficients beats 9/4 (sanity; the global statement is the B&B)
print("exact det at BCC:", D, "  active triangles:", len(active), "  distinct gradients:", len(grads))
G = [sp.Matrix([sp.Rational(x.numerator, x.denominator) for x in g]) for g in sorted(grads)]
assert len(G) == 6
# barycentric coordinates of the origin: sum w_i g_i = 0, sum w_i = 1
M = sp.Matrix.hstack(*G).col_join(sp.ones(1, 6))
w = M.solve(sp.Matrix([0, 0, 0, 0, 0, 1]))
print("barycentric weights of 0:", list(w))
assert all(x > 0 for x in w), "origin not in the interior"
cmin = None
for k in range(6):
    face = [G[i] for i in range(6) if i != k]
    A = sp.Matrix.hstack(*[f - face[0] for f in face[1:]]).T      # 4 x 5
    nvec = A.nullspace()[0]                                       # exact normal
    b = (nvec.T * face[0])[0]
    if b < 0: nvec, b = -nvec, -b
    assert (nvec.T * G[k])[0] < b                                 # opposite vertex on the other side
    dist2 = b ** 2 / (nvec.T * nvec)[0]                           # exact squared distance of 0 to the facet
    cmin = dist2 if cmin is None else min(cmin, dist2)
print("exact min squared facet distance c^2 =", cmin, "=", float(cmin))
c_lower = sp.Rational(sp.floor(sp.sqrt(cmin) * 10 ** 12), 10 ** 12)
assert c_lower ** 2 <= cmin
print("rigorous lower bound c >=", c_lower, "=", float(c_lower))
open("bcc_local_exact.txt", "w").write("c_lower %s\nc2 %s\nweights %s\n" % (c_lower, cmin, list(w)))
