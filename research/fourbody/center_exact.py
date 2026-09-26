"""Exact value of Q at the centre, by exhaustive enumeration with pruning in integer arithmetic.
Squared distances are values of an integral Gram form; the squared product of all mutual distances of the
k points {0, p_1, ..., p_{k-1}} is compared with Q*^2 (scaled to the Gram normalisation).
Pruning: points are added in order of increasing norm; each further point multiplies the squared product by at
least (its squared norm) * 1^(others) >= current squared norm, since all squared distances are >= lambda_1^2."""
import itertools, math
from fractions import Fraction as F
def run(G, k, target2, label, R=12):
    q = lambda v: G[0][0] * v[0] ** 2 + 2 * G[0][1] * v[0] * v[1] + G[1][1] * v[1] ** 2
    lam2 = min(q(v) for v in itertools.product(range(-3, 4), repeat=2) if any(v))
    pts = sorted([v for v in itertools.product(range(-R, R + 1), repeat=2) if any(v)], key=q)
    best = [None, 0]
    def rec(chosen, start, prod):
        if len(chosen) == k - 1:
            if best[0] is None or prod < best[0]: best[0], best[1] = prod, 1
            elif prod == best[0]: best[1] += 1
            return
        rem = k - 1 - len(chosen)
        for i in range(start, len(pts)):
            p = pts[i]
            # lower bound for the final squared product if p and further points (norm >= q(p)) are added
            lb = prod * (q(p) * lam2 ** len(chosen)) * (q(p) * lam2 ** (len(chosen) + 1)) ** (rem - 1) if rem > 1 else prod * q(p) * lam2 ** len(chosen)
            if lb > target2: break
            f = q(p)
            for c in chosen: f *= q((p[0] - c[0], p[1] - c[1]))
            if prod * f <= target2: rec(chosen + [p], i + 1, prod * f)
    rec([], 0, 1)
    # completeness: a set containing p has squared product >= q(p) * lam2^(k(k-1)/2 - 1); all points with
    # q(p) <= target2 / lam2^(...) must be in the list (q >= (3/4) max(|m|,|n|)^2 * min eigen-ish, checked crudely)
    bound = target2 / lam2 ** (k * (k - 1) // 2 - 1)
    outside_min = min(q(v) for v in itertools.product(range(-R - 1, R + 2), repeat=2) if max(abs(v[0]), abs(v[1])) == R + 1)
    assert outside_min > bound, ("enumeration radius too small", outside_min, bound)
    print("%s: min squared product = %s (target %s), attained by %d sets" % (label, best[0], target2, best[1]))
# hexagonal: Gram [[2,1],[1,2]] has covolume^2 = 3; unit covolume scaling of squared lengths by 1/sqrt(3)
# four points: 6 squared distances -> factor 3^(-3) ; Q*^2 = 64/9  ->  target2 = 64/9 * 27 = 192
if len(__import__("sys").argv) == 1: run([[2, 1], [1, 2]], 4, 192, "hexagonal, k=4")
# square Z^2: unit covolume; five points, Q*^2 = 160
if len(__import__("sys").argv) == 1: run([[1, 0], [0, 1]], 5, 160, "square, k=5")
if __name__ == "__main__" and len(__import__("sys").argv) > 1 and __import__("sys").argv[1] == "k6":
    import math
    t = math.ceil(106.0 ** 2 * 3 ** 7.5)
    run([[2, 1], [1, 2]], 6, t, "hexagonal, k=6 (Gram units; Q^2 = value / 3^7.5)", R=46)
