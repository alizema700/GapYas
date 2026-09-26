"""Local lemma at BCC with interval arithmetic (mpmath.iv).

For the active triangles t at G0 = BCC (params q0 = (1,1,1/3,1/3,-1/3)), f_t(q) = sqrt(l1 l2 l3 / D).
On the box B_r = q0 + [-r, r]^5 we bound ||Hess f_t||_F <= M (interval arithmetic). With
c = min_{|d|=1} max_t (-grad f_t(q0) . d) > 0 (origin strictly inside conv{grad f_t}),
Taylor gives  min_t f_t(q0 + d) <= 3/2 - c|d| + M|d|^2/2 < 3/2  for 0 < |d| < 2c/M.
Hence Q < 3/2 on B_r minus {q0} as soon as sqrt(5) r < 2c/M.
"""
import itertools, json, sys
import numpy as np
from mpmath import iv, mpf
from scipy.spatial import ConvexHull
iv.dps = 30

q0 = [iv.mpf(1), iv.mpf(1), iv.mpf(1) / 3, iv.mpf(1) / 3, -iv.mpf(1) / 3]
def coeffs(v):
    m1, m2, m3 = v
    return [m1 * m1, m2 * m2, m3 * m3, 2 * m1 * m2, 2 * m1 * m3, 2 * m2 * m3]

# active triangles at BCC (exact rational check done in float then re-verified below)
V = [m for m in itertools.product(range(-2, 3), repeat=3) if any(m)]
def lval(c, q): return c[0] + c[1] * q[0] + c[2] * q[1] + c[3] * q[2] + c[4] * q[3] + c[5] * q[4]
qf = [1.0, 1.0, 1 / 3, 1 / 3, -1 / 3]
Df = 16 / 27
act = []
for m in V:
    for n in V:
        if m >= n: continue
        d = tuple(a - b for a, b in zip(m, n))
        val = (lval(coeffs(m), qf) * lval(coeffs(n), qf) * lval(coeffs(d), qf) / Df) ** 0.5
        if abs(val - 1.5) < 1e-9: act.append((m, n, d))
print("active unordered triangles at BCC:", len(act))

def D(q):  # determinant of the Gram matrix with g11 = 1
    g22, g33, g12, g13, g23 = q
    return g22 * g33 - g23 ** 2 - g12 ** 2 * g33 - g13 ** 2 * g22 + 2 * g12 * g13 * g23
def gradD(q):
    g22, g33, g12, g13, g23 = q
    return [g33 - g13 ** 2, g22 - g12 ** 2, -2 * g12 * g33 + 2 * g13 * g23, -2 * g13 * g22 + 2 * g12 * g23, -2 * g23 + 2 * g12 * g13]
def hessD(q):
    g22, g33, g12, g13, g23 = q
    z = 0 * g22
    return [[z, 1 + z, z, -2 * g13, z],
            [1 + z, z, -2 * g12, z, z],
            [z, -2 * g12, -2 * g33, 2 * g23, 2 * g13],
            [-2 * g13, z, 2 * g23, -2 * g22, 2 * g12],
            [z, z, 2 * g13, 2 * g12, -2 + z]]

def f_grad_hess(tri, q):
    ls = [lval(coeffs(v), q) for v in tri]
    gl = [coeffs(v)[1:] for v in tri]
    Dq = D(q); gD = gradD(q); HD = hessD(q)
    f = iv.sqrt(ls[0] * ls[1] * ls[2] / Dq) if isinstance(Dq, type(iv.mpf(1))) else (ls[0] * ls[1] * ls[2] / Dq) ** 0.5
    gF = [sum(gl[i][a] / ls[i] for i in range(3)) / 2 - gD[a] / Dq / 2 for a in range(5)]
    HF = [[-sum(gl[i][a] * gl[i][b] / ls[i] ** 2 for i in range(3)) / 2 - HD[a][b] / Dq / 2 + gD[a] * gD[b] / Dq ** 2 / 2
           for b in range(5)] for a in range(5)]
    Hf = [[f * (HF[a][b] + gF[a] * gF[b]) for b in range(5)] for a in range(5)]
    gf = [f * gF[a] for a in range(5)]
    return f, gf, Hf

# gradients at q0 (point intervals) and c
grads = []
for tri in act:
    f, g, _ = f_grad_hess(tri, q0)
    assert f.a <= mpf(3) / 2 <= f.b and f.delta < mpf(10) ** -20
    grads.append([float(x.mid) for x in g])
G = np.unique(np.round(np.array(grads), 12), axis=0)
hull = ConvexHull(G)
c = float(np.min(-hull.equations[:, -1]))
print("distinct gradients:", len(G), " c =", c)

def M_of(r):
    box = [x + iv.mpf([-r, r]) for x in q0]
    M = 0.0
    for tri in act:
        _, _, H = f_grad_hess(tri, box)
        fro = sum(max(abs(float(h.a)), abs(float(h.b))) ** 2 for row in H for h in row) ** 0.5
        M = max(M, fro)
    return M

res = []
for r in [0.002, 0.004, 0.006, 0.007, 0.008, 0.01, 0.015, 0.02]:
    M = M_of(r)
    ok = np.sqrt(5) * r < 2 * c / M
    res.append(dict(r=r, M=M, radius=2 * c / M, ok=bool(ok)))
    print("r=%.4f  M<=%.4f  2c/M=%.5f  sqrt5*r=%.5f  lemma holds on box: %s" % (r, M, 2 * c / M, np.sqrt(5) * r, ok))
json.dump(dict(c=c, n_active=len(act), results=res), open("local_lemma.json", "w"), indent=1)
