"""Active sets and gradient geometry at the max-min lattices (2D: Lambda_*, 3D: BCC).
h_t = log(pi_t / P) for ordered pairs t = (m, n) of distinct nonzero lattice labels; T_nu = P^-nu sum exp(-nu h_t).
Checks hypothesis (A3): 0 lies in the interior of conv{grad h_t : t active}; reports the margin
c = min_{|V|=1} max_t (-grad h_t . V) (distance from 0 to the hull boundary), and the inactive gap gamma."""
import itertools, numpy as np
from scipy.optimize import linprog

def hull_margin(G):
    # distance from 0 to boundary of conv(G) = min over facets; compute as min_V max_t -g.V via LP over many V? use exact: 
    # margin = max r s.t. ball(0,r) in hull  <=>  for all unit V: max_t g.V >= r ; approximate by sampling + LP check
    d = G.shape[1]; best = np.inf
    rng = np.random.default_rng(0)
    V = rng.normal(size=(200000, d)); V /= np.linalg.norm(V, axis=1)[:, None]
    return float(np.min(np.max(V @ G.T, axis=1)))

# ---------- 2D ----------
a = ((1 + np.sqrt(17)) / 8) ** 0.25; y0 = a ** -2; P2 = 2 * a ** 3
def h2(x, y, R=4):
    ax = range(-R, R + 1); lab = [(m, n) for m in ax for n in ax if (m, n) != (0, 0)]
    Z = np.array([((m + n * x) / np.sqrt(y), n * np.sqrt(y)) for m, n in lab])
    out = {}
    for i, j in itertools.permutations(range(len(lab)), 2):
        p = np.linalg.norm(Z[i]) * np.linalg.norm(Z[j]) * np.linalg.norm(Z[i] - Z[j])
        out[(lab[i], lab[j])] = np.log(p / P2)
    return out
H0 = h2(0.0, y0); vals = np.array(sorted(H0.values()))
act = [t for t, v in H0.items() if v < 1e-9]
gap = min(v for v in H0.values() if v > 1e-9)
eps = 1e-6
Hx = h2(eps, y0); Hy = h2(0.0, y0 + eps)
G = np.array([[(Hx[t] - H0[t]) / eps, (Hy[t] - H0[t]) / eps] for t in act])
print("2D: active ordered pairs:", len(act), " inactive gap gamma = %.5f" % gap)
print("    distinct gradients (d/dx, d/dy):", np.unique(np.round(G, 6), axis=0))
print("    hull margin c = %.5f" % hull_margin(G))

# ---------- 3D (BCC, Gram chart q = (g22, g33, g12, g13, g23), g11 = 1) ----------
q0 = np.array([1, 1, 1/3, 1/3, -1/3])
def gram(q): g22, g33, g12, g13, g23 = q; return np.array([[1, g12, g13], [g12, g22, g23], [g13, g23, g33]])
labs = [np.array(v) for v in itertools.product(range(-2, 3), repeat=3) if any(v)]
def h3(q):
    Gm = gram(q); det = np.linalg.det(Gm); out = {}
    for i, j in itertools.permutations(range(len(labs)), 2):
        m, n = labs[i], labs[j]; d = m - n
        out[(i, j)] = 0.5 * np.log((m @ Gm @ m) * (n @ Gm @ n) * (d @ Gm @ d) / det) - np.log(1.5)
    return out
H0 = h3(q0); act = [t for t, v in H0.items() if v < 1e-9]
gap = min(v for v in H0.values() if v > 1e-9)
G = []
for t in act:
    g = []
    for k in range(5):
        e = np.zeros(5); e[k] = 1e-6; g.append((h3_t := None) or 0)
    G.append(g)
# vectorised gradient
Gr = np.zeros((len(act), 5))
for k in range(5):
    e = np.zeros(5); e[k] = 1e-6; Hk = h3(q0 + e); Hm = h3(q0 - e)
    Gr[:, k] = [(Hk[t] - Hm[t]) / 2e-6 for t in act]
print("3D: active ordered pairs:", len(act), " inactive gap gamma = %.5f" % gap)
print("    distinct gradients:", len(np.unique(np.round(Gr, 6), axis=0)), " sum of gradients:", np.round(Gr.sum(0), 8))
print("    hull margin c = %.5f" % hull_margin(Gr))
