"""Exact 'two-family balance' equation for the 2D minimiser on the imaginary axis:
   6 H1'(y) e^{-nu H1(y)} + 24 H2'(y) e^{-nu H2(y)} = 0,
   H1 = log(2 a^3), H2 = (1/2) log(a^2 + a^-2), a = y^(-1/2)  (collinear triples along u, right triangles).
Compare its root with the true minimiser and extract the 1/nu and 1/nu^2 coefficients."""
import mpmath as mp
mp.mp.dps = 50
def H1(y): a = y ** -0.5; return mp.log(2 * a ** 3)
def H2(y): a = y ** -0.5; return mp.log(a ** 2 + a ** -2) / 2
def root(nu):
    g = lambda y: mp.log(-4 * mp.diff(H2, y) / mp.diff(H1, y)) - nu * (H2(y) - H1(y))
    return mp.findroot(g, mp.mpf('1.24'))
a = ((1 + mp.sqrt(17)) / 8) ** mp.mpf(0.25); yinf = a ** -2
for nu in [10, 15, 20, 30, 50, 80, 120, 200, 300, 1000, 10000]:
    y = root(nu); print("nu=%6d  y_bal=%.10f  nu(yinf-y)=%.8f" % (nu, y, nu * (yinf - y)))
# asymptotic coefficients by Richardson on large nu
f = lambda nu: nu * (yinf - root(nu))
k1 = f(10**6); lam = (f(10**5) - f(10**6)) / (mp.mpf(1)/10**5 - mp.mpf(1)/10**6)
print("kappa ~", mp.nstr(k1 - lam / 10**6, 12), " lambda ~", mp.nstr(lam, 8))
