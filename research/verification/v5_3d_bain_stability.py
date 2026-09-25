"""V5: Bain path fcc <-> bcc and 5x5 Hessian (all lattice deformations) at BCC and FCC."""
import json
import numpy as np
from lattice_sums import triangle_sum, FCC, BCC, unit_covolume
from v3_3d_lattices import lll

def T(A, nu):
    return triangle_sum(lll(unit_covolume(A)), (nu,) * 3, Rs=(16, 24, 32, 40, 48))

def bct(r):   # conventional a=1, c=r ; r=1 -> BCC, r=sqrt(2) -> FCC
    return 0.5 * np.array([[-1, 1, 1], [1, -1, 1], [r, r, -r]], float)

def deform(A, e):
    # symmetric volume-preserving strain with 5 parameters
    E = np.array([[e[0], e[2], e[3]], [e[2], e[1], e[4]], [e[3], e[4], -e[0] - e[1]]])
    w, V = np.linalg.eigh(E)
    return V @ np.diag(np.exp(w)) @ V.T @ A

def hessian(A, nu, h=2e-3):
    f0 = T(A, nu); H = np.zeros((5, 5)); I = np.eye(5)
    for i in range(5):
        H[i, i] = (T(deform(A, h * I[i]), nu) - 2 * f0 + T(deform(A, -h * I[i]), nu)) / h**2
        for j in range(i):
            H[i, j] = H[j, i] = (T(deform(A, h * (I[i] + I[j])), nu) - T(deform(A, h * (I[i] - I[j])), nu)
                                 - T(deform(A, h * (-I[i] + I[j])), nu) + T(deform(A, -h * (I[i] + I[j])), nu)) / (4 * h * h)
    return np.linalg.eigvalsh(H)

def main():
    out = {"bain": {}, "hessian": {}}
    rs = np.linspace(0.8, 1.6, 33)
    for nu in [3.5, 4.5, 6.0, 9.0]:
        vals = [T(bct(r), nu) for r in rs]
        out["bain"][nu] = list(zip(rs.tolist(), vals))
        i = int(np.argmin(vals))
        print("nu=%.1f Bain path: min at c/a=%.3f (BCC c/a=1, FCC c/a=1.414) T=%.10f; value at FCC %.10f" % (nu, rs[i], vals[i], T(bct(np.sqrt(2)), nu)), flush=True)
    for nu in [3.5, 4.5, 6.0, 9.0]:
        hb = hessian(unit_covolume(BCC), nu); hf = hessian(unit_covolume(FCC), nu)
        out["hessian"][nu] = dict(BCC=hb.tolist(), FCC=hf.tolist())
        print("nu=%.1f Hessian eigenvalues BCC %s | FCC %s" % (nu, np.round(hb, 4), np.round(hf, 4)), flush=True)
    json.dump(out, open("results/v5_bain_stability.json", "w"), indent=1, default=float)


if __name__ == '__main__':
    main()
