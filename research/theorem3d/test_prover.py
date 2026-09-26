"""Sanity tests for the branch-and-bound prover.
(1) soundness: for random boxes and random points inside, the true Q(G) (brute force over many lattice
    triangles) never exceeds the box upper bound;
(2) the bound is attained: at BCC the true Q is 3/2;
(3) the prover must FAIL when asked to show Q < 1.49 (since Q(BCC) = 1.5 > 1.49)."""
import itertools, numpy as np, bnb

def Q_true(q, R=3):
    g22, g33, g12, g13, g23 = q
    G = np.array([[1, g12, g13], [g12, g22, g23], [g13, g23, g33]])
    L = np.linalg.cholesky(G).T                     # basis matrix with Gram G (columns)
    m = np.array([v for v in itertools.product(range(-R, R + 1), repeat=3) if any(v)], float)
    X = m @ L.T; r = np.linalg.norm(X, axis=1)
    D = np.linalg.norm(X[:, None] - X[None], axis=-1); np.fill_diagonal(D, np.inf)
    return np.min(np.outer(r, r) * D) / np.sqrt(np.linalg.det(G))

rng = np.random.default_rng(1); worst = np.inf; n = 0
for _ in range(400):
    c = np.array([rng.uniform(1, 2.6), rng.uniform(1, 3), rng.uniform(0, .5), rng.uniform(0, .5), rng.uniform(-1, 1)])
    c[1] = max(c[1], c[0]); c[4] = np.clip(c[4], -c[0] / 2, c[0] / 2)
    hw = rng.uniform(1e-3, 5e-2, 5)
    lo, hi = c - hw, c + hw; lo[2:4] = np.maximum(lo[2:4], 0)
    if np.linalg.eigvalsh(np.array([[1, c[2], c[3]], [c[2], c[0], c[4]], [c[3], c[4], c[1]]])).min() <= 0.05: continue
    ub = bnb.upper_bound(lo[None], hi[None])[0]
    for _ in range(3):
        q = lo + rng.uniform(0, 1, 5) * (hi - lo)
        G = np.array([[1, q[2], q[3]], [q[2], q[0], q[4]], [q[3], q[4], q[1]]])
        if np.linalg.eigvalsh(G).min() <= 0: continue
        Q = Q_true(q); n += 1
        worst = min(worst, ub - Q)
        assert Q <= ub + 1e-12, (q, Q, ub)
print("(1) soundness: %d samples, min(ub - Q_true) = %.3e (must be >= 0)" % (n, worst))
print("(2) Q_true(BCC) =", Q_true(np.array([1, 1, 1/3, 1/3, -1/3])))
bnb.TARGET = 1.49
ok, nproc = bnb.run([1, 1, 1/3, 1/3, -1/3], [0.007] * 5, max_boxes=3 * 10**6, log=lambda *a: None)
print("(3) prover with target 1.49 (must fail): success =", ok, "after", nproc, "boxes")
