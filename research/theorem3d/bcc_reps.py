"""Enumerate the Gram parameters (g22,g33,g12,g13,g23) of BCC in the search region R."""
import itertools, numpy as np
B = np.array([[-1, 1, 1], [1, -1, 1], [1, 1, -1]], float) / 2.0      # BCC primitive vectors (columns)
reps = set()
rng = range(-2, 3)
for M in itertools.product(rng, repeat=9):
    U = np.array(M, float).reshape(3, 3)
    if round(abs(np.linalg.det(U))) != 1: continue
    A = B @ U; G = A.T @ A; G = G / G[0, 0]
    g22, g33, g12, g13, g23 = G[1, 1], G[2, 2], G[0, 1], G[0, 2], G[1, 2]
    lam1 = min(np.linalg.norm(B @ np.array(v)) for v in itertools.product(range(-2, 3), repeat=3) if any(v)) ** 2
    if abs(A.T[0] @ A.T[0] - lam1) > 1e-9: continue
    if not (1 - 1e-9 <= g22 <= g33 + 1e-9 and 0 <= g12 <= .5 + 1e-9 and 0 <= g13 <= .5 + 1e-9 and abs(g23) <= g22 / 2 + 1e-9): continue
    ok = all(1 + g22 + 2 * (e1 * e2 * g12 + e1 * g13 + e2 * g23) >= -1e-9 for e1 in (1, -1) for e2 in (1, -1))
    if ok: reps.add(tuple(np.round([g22, g33, g12, g13, g23], 9)))
for r in sorted(reps): print(r)
np.save("bcc_reps.npy", np.array(sorted(reps)))
