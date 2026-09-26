"""Local analysis at BCC for Q(G) = P(G)/sqrt(det G), Gram parametrisation with g11 = 1:
params q = (g22, g33, g12, g13, g23).  Q <= min_t f_t, f_t(G) = sqrt(prod_i m_i^T G m_i)/sqrt(det G)."""
import itertools, numpy as np
from scipy.spatial import ConvexHull

def gram(q):
    g22, g33, g12, g13, g23 = q
    return np.array([[1, g12, g13], [g12, g22, g23], [g13, g23, g33]], float)

K = 2
VECS = np.array([m for m in itertools.product(range(-K, K + 1), repeat=3) if any(m)])
def triangles(G, nbest=None, tol=None):
    """all triangles (0, m, n) with m, n in VECS, m != n; returns list of (value, m, n)"""
    L2 = np.einsum("ij,jk,ik->i", VECS, G, VECS)
    D = VECS[:, None, :] - VECS[None, :, :]
    L3 = np.einsum("abj,jk,abk->ab", D, G, D)
    val = np.sqrt(np.outer(L2, L2) * L3) / np.sqrt(np.linalg.det(G))
    np.fill_diagonal(val, np.inf)
    idx = np.argwhere(val <= val.min() * (1 + tol)) if tol else None
    return val, idx

q_bcc = np.array([1.0, 1.0, 1/3, 1/3, -1/3])
G = gram(q_bcc)
val, idx = triangles(G, tol=1e-9)
print("Q(BCC) =", val.min(), " #active ordered triangles:", len(idx))

def f_t(q, m, n):
    G = gram(q); d = m - n
    return np.sqrt((m @ G @ m) * (n @ G @ n) * (d @ G @ d)) / np.sqrt(np.linalg.det(G))
grads = []
for a, b in idx:
    m, n = VECS[a], VECS[b]; h = 1e-6
    gr = np.array([(f_t(q_bcc + h * e, m, n) - f_t(q_bcc - h * e, m, n)) / (2 * h) for e in np.eye(5)])
    grads.append(gr)
grads = np.unique(np.round(np.array(grads), 10), axis=0)
print("distinct active gradients:", len(grads)); print(np.round(grads, 5))
# is 0 in the interior of conv{grad f_t}?  (then every direction decreases some active f_t at first order)
try:
    hull = ConvexHull(grads)
    dist = -hull.equations[:, -1]           # offsets: n.x + b <= 0 inside; distance of origin to facet = -b
    print("origin inside hull:", np.all(dist > 0), " min distance to facets c =", dist.min())
except Exception as ex:
    print("hull failed:", ex)
