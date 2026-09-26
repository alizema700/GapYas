# Lean 4 formalisation of Theorem 1 (planar max–min triangle problem)

Machine-checked with Lean 4 (`lean-toolchain`) and Mathlib.  Build:

```
curl -sSfL https://raw.githubusercontent.com/leanprover/elan/master/elan-init.sh | sh -s -- -y
lake exe cache get      # downloads the pre-built Mathlib
lake build              # checks TriangleLattice/Basic.lean
lake env lean Check.lean   # prints the axioms used by the main theorems
```

Main statements in `TriangleLattice/Basic.lean`:

| Lean name | content |
|---|---|
| `upper_bound` | Part A: for every reduced basis `(a,0), (p,1/a)` with `0 ≤ p ≤ a/2`, `a⁴ ≤ 4/3`: `min(4a⁶, Φ(a,p)) ≤ 4s*³` (squared version of `P ≤ P*`) |
| `upper_bound_eq` | equality forces `a² = s*` and `p = 0` (the rectangle `Λ*`) |
| `rect_lower` | Part B: in the rectangular lattice with parameter `s ∈ (1/√3, 1)` every triangle through `0` with distinct vertices has squared product `≥ min(4s³, s + 1/s)` |
| `rect_lower_star` | for `s = s*` (`4s*⁴ = s*² + 1`) this lower bound equals `4s*³ = P*²` |

`#print axioms` reports only `propext`, `Classical.choice`, `Quot.sound` (no `sorry`).

Not formalised (classical, cited in the paper): existence of a Lagrange–Gauss reduced basis, Hermite's
bound `a⁴ ≤ 4/3`, translation invariance of the triple product, and the compactness argument for Theorem 2.
