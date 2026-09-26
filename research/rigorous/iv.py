"""Vectorised interval arithmetic on numpy float64 arrays with outward rounding.

Soundness: IEEE-754 binary64 +, -, *, / and sqrt are correctly rounded (round to nearest), so the
computed result r of an operation on floats satisfies |r - exact| <= ulp(r)/2.  Moving r one ulp
outwards with nextafter therefore gives a valid lower/upper bound of the exact result.  Interval
operations evaluate the relevant endpoint combinations and widen the result in this way.
Transcendental functions (x^-s for real s) are evaluated with Arb (python-flint) and rounded outwards.
"""
import numpy as np
from flint import arb, ctx
ctx.prec = 128

INF = np.inf
def down(x): return np.nextafter(x, -INF)
def up(x): return np.nextafter(x, INF)

class IV:
    __slots__ = ("lo", "hi")
    def __init__(self, lo, hi=None):
        lo = np.asarray(lo, dtype=np.float64)
        self.lo = lo; self.hi = lo.copy() if hi is None else np.asarray(hi, dtype=np.float64)
    @staticmethod
    def exact(x):  # a float (array) that is represented exactly
        x = np.asarray(x, dtype=np.float64); return IV(x, x.copy())
    def __getitem__(self, k): return IV(self.lo[k], self.hi[k])
    def __setitem__(self, k, v): self.lo[k] = v.lo; self.hi[k] = v.hi
    @property
    def shape(self): return self.lo.shape
    def copy(self): return IV(self.lo.copy(), self.hi.copy())
    def __add__(self, o):
        o = o if isinstance(o, IV) else IV.exact(o)
        return IV(down(self.lo + o.lo), up(self.hi + o.hi))
    __radd__ = __add__
    def __neg__(self): return IV(-self.hi, -self.lo)
    def __sub__(self, o):
        o = o if isinstance(o, IV) else IV.exact(o)
        return IV(down(self.lo - o.hi), up(self.hi - o.lo))
    def __rsub__(self, o): return IV.exact(o) - self
    def __mul__(self, o):
        o = o if isinstance(o, IV) else IV.exact(o)
        c = [self.lo * o.lo, self.lo * o.hi, self.hi * o.lo, self.hi * o.hi]
        lo = np.minimum(np.minimum(c[0], c[1]), np.minimum(c[2], c[3]))
        hi = np.maximum(np.maximum(c[0], c[1]), np.maximum(c[2], c[3]))
        return IV(down(lo), up(hi))
    __rmul__ = __mul__
    def recip_pos(self):
        assert np.all(self.lo > 0), "recip_pos needs positive intervals"
        return IV(down(1.0 / self.hi), up(1.0 / self.lo))
    def __truediv__(self, o):
        o = o if isinstance(o, IV) else IV.exact(o)
        return self * o.recip_pos()
    def sqrt_pos(self):
        assert np.all(self.lo >= 0)
        return IV(down(np.sqrt(self.lo)), up(np.sqrt(self.hi)))
    def sq(self):
        a, b = self.lo * self.lo, self.hi * self.hi
        lo = np.where((self.lo <= 0) & (self.hi >= 0), 0.0, np.minimum(a, b))
        return IV(np.maximum(down(lo), 0.0) if np.all(lo >= 0) else down(lo), up(np.maximum(a, b)))
    def mag(self): return np.maximum(np.abs(self.lo), np.abs(self.hi))
    def mid(self): return 0.5 * self.lo + 0.5 * self.hi
    def rad(self): return up(np.maximum(self.hi - self.mid(), self.mid() - self.lo))
    def __repr__(self): return "IV(%s, %s)" % (self.lo, self.hi)

def arb_to_iv(a):
    """outward float enclosure of an arb ball"""
    lo = float(a.lower()); hi = float(a.upper())
    return down(lo), up(hi)

def pow_neg_arb(x_iv, s):
    """elementwise enclosure of x^-s for positive intervals x (Arb, then outward rounding)"""
    s = arb(s) if not isinstance(s, arb) else s
    lo = np.empty(x_iv.shape); hi = np.empty(x_iv.shape)
    flo, fhi = x_iv.lo.ravel(), x_iv.hi.ravel()
    olo, ohi = lo.ravel(), hi.ravel()
    for k in range(flo.size):
        a = arb(float(fhi[k])) ** (-s)     # decreasing in x: lower bound from the upper endpoint
        b = arb(float(flo[k])) ** (-s)
        olo[k] = arb_to_iv(a)[0]; ohi[k] = arb_to_iv(b)[1]
    return IV(lo, hi)

# ---- rigorous bound for sums/convolutions of floating-point data ---------------------------------
U = 2.0 ** -53
def gamma(n):
    """gamma_n = n u / (1 - n u) (Higham): |fl(sum_{i<=n+1} p_i) - sum p_i| <= gamma_n sum |p_i| for any
    summation order, and each product adds one more relative rounding."""
    nu = n * U
    assert nu < 0.01
    return up(nu / (1 - nu))
