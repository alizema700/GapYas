"""Exact (rational) check of (H3) for symmetric max-min lattices.
For a lattice whose automorphism group acts absolutely irreducibly, the sum of the gradients over the (Aut-invariant)
active set is an invariant trace-free form, hence 0; so 0 is a strictly positive combination, and (H3) is equivalent
to the active gradients spanning the tangent space {dG : tr(G^-1 dG) = 0}. We verify the rank exactly."""
import itertools, sympy as sp
from fractions import Fraction as F
def run(G, k, cr, label):
    d = len(G); G = sp.Matrix(G); Ginv = G.inv()
    C = [sp.Matrix(v) for v in itertools.product(range(-cr, cr + 1), repeat=d) if any(v)]
    nrm = [(v.T * G * v)[0] for v in C]
    order = sorted(range(len(C)), key=lambda i: nrm[i])[:150]
    best, sets = None, []
    for combo in itertools.combinations(order, k - 1):
        pts = [sp.zeros(d, 1)] + [C[i] for i in combo]
        p = 1
        for a, b in itertools.combinations(range(k), 2):
            v = pts[a] - pts[b]; p *= (v.T * G * v)[0]
            if best is not None and p > best: break
        else:
            if best is None or p < best: best, sets = p, [pts]
            elif p == best: sets.append(pts)
    iu = [(i, j) for i in range(d) for j in range(i, d)]
    def grad(pts):
        g = []
        for (i, j) in iu:
            s = 0
            for a, b in itertools.combinations(range(k), 2):
                v = pts[a] - pts[b]; q = (v.T * G * v)[0]
                s += (v[i] * v[j] * (1 if i == j else 2)) / q
            g.append(s)
        return g
    ndir = [Ginv[i, j] * (1 if i == j else 2) for (i, j) in iu]     # d/dG_ij of log det
    A = sp.Matrix([grad(s) for s in sets])
    # tangent space basis: kernel of ndir; rank of A restricted = rank of [A; ndir] - 1
    rk = sp.Matrix.vstack(A, sp.Matrix([ndir])).rank() - 1
    tot = sp.zeros(1, len(iu))
    for r in range(A.rows): tot += A.row(r)
    # component of the sum along the tangent space: sum must be a multiple of ndir
    par = sp.Matrix.vstack(tot, sp.Matrix([ndir])).rank() == 1
    print("%s: min squared product = %s, active sets = %d, rank in tangent space = %d of %d, sum of gradients normal to tangent space: %s"
          % (label, best, len(sets), rk, len(iu) - 1, par))
run([[2, -1, 0, 0], [-1, 2, -1, -1], [0, -1, 2, 0], [0, -1, 0, 2]], 3, 2, "D4 (Cartan Gram), three-body")
run([[2, 1], [1, 2]], 4, 2, "A2 hexagonal, four-body")
run([[3, -1, -1], [-1, 3, -1], [-1, -1, 3]], 3, 1, "BCC (Gram of <(1,1,-1),(1,-1,1),(-1,1,1)>), three-body")
