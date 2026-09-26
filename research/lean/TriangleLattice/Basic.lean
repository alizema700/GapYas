import Mathlib

/-!
# A max–min triangle problem for planar lattices (Theorem 1 of the paper, formalised core)

For a planar lattice `Λ` of unit covolume let `P(Λ)` be the smallest product of the three mutual
distances of three distinct lattice points.  Theorem 1 states `P(Λ) ≤ P*` with equality exactly for the
rectangular lattice `Λ* = a* ℤ × a*⁻¹ ℤ`, where `4 a*⁸ = a*⁴ + 1`.

This file formalises the two mathematical parts of the proof; all statements are about squared
quantities, so no square roots of the distances are needed.

* **Part A (upper bound, `upper_bound` and `upper_bound_eq`).**  For a reduced basis `u = (a, 0)`,
  `v = (p, a⁻¹)` with `0 ≤ p ≤ a/2` and `a⁴ ≤ 4/3`, the collinear triple `(0, u, 2u)` has squared product
  `4a⁶` and the triple `(0, u, v)` has squared product `Φ(a,p) = a² (p² + a⁻²) ((a-p)² + a⁻²)`.
  We prove `min (4a⁶) Φ ≤ 4 s*³` where `s* = a*²`, with equality only if `a² = s*` and `p = 0`.
* **Part B (the extremal lattice, `rect_lower`).**  In the rectangular lattice with side lengths
  `√s` and `1/√s` (squared norm `D s m n = s m² + n²/s`), every triangle with vertices
  `0, (m₁,n₁), (m₂,n₂)` (pairwise distinct) has squared side-length product
  `≥ min (4 s³) (s + 1/s)`, for every `s` with `1/3 < s² < 1`.  For `s = s*` both terms equal
  `4 s*³`, the squared value of `P*`.

Not formalised (classical and cited in the paper): the existence of a reduced basis
(Lagrange–Gauss reduction), Hermite's bound `a⁴ ≤ 4/3`, translation invariance of the triple product,
and the compactness argument deducing Theorem 2 from Theorem 1.
-/

namespace TriangleLattice

open Real

/-! ## Part B: the rectangular lattice -/

/-- squared length of the lattice vector `(m √s, n / √s)` -/
noncomputable def D (s : ℝ) (m n : ℤ) : ℝ := s * (m : ℝ) ^ 2 + (n : ℝ) ^ 2 / s

/-- the symmetric bilinear form belonging to `D` -/
noncomputable def Bf (s : ℝ) (m₁ n₁ m₂ n₂ : ℤ) : ℝ :=
  s * (m₁ : ℝ) * m₂ + (n₁ : ℝ) * n₂ / s

/-- Lagrange's identity: `D(w₁) D(w₂) = B(w₁,w₂)² + det(w₁,w₂)²`. -/
lemma lagrange (s : ℝ) (hs : 0 < s) (m₁ n₁ m₂ n₂ : ℤ) :
    D s m₁ n₁ * D s m₂ n₂ = (Bf s m₁ n₁ m₂ n₂) ^ 2 + ((m₁ : ℝ) * n₂ - (n₁ : ℝ) * m₂) ^ 2 := by
  unfold D Bf
  field_simp
  ring

lemma D_nonneg (s : ℝ) (hs : 0 < s) (m n : ℤ) : 0 ≤ D s m n := by
  unfold D; positivity

/-- a nonzero integer has square at least one -/
lemma one_le_sq_int (n : ℤ) (h : n ≠ 0) : (1 : ℝ) ≤ (n : ℝ) ^ 2 := by
  have h1 : (1 : ℤ) ≤ |n| := Int.one_le_abs h
  have h2 : (1 : ℝ) ≤ |(n : ℝ)| := by exact_mod_cast h1
  nlinarith [sq_abs (n : ℝ), abs_nonneg (n : ℝ)]

/-- every nonzero lattice vector has squared length at least `s` (when `s ≤ 1`) -/
lemma D_ge (s : ℝ) (hs : 0 < s) (hs1 : s ≤ 1) (m n : ℤ) (h : (m, n) ≠ (0, 0)) : s ≤ D s m n := by
  unfold D
  by_cases hm : m = 0
  · subst hm
    have hn : n ≠ 0 := by intro hn; exact h (by simp [hn])
    have : (1 : ℝ) ≤ (n : ℝ) ^ 2 := one_le_sq_int n hn
    have h1 : s ≤ 1 / s := by rw [le_div_iff₀ hs]; nlinarith
    have h2 : 1 / s ≤ (n : ℝ) ^ 2 / s := by gcongr
    have := le_trans h1 h2
    simpa using this
  · have : (1 : ℝ) ≤ (m : ℝ) ^ 2 := one_le_sq_int m hm
    have : 0 ≤ (n : ℝ) ^ 2 / s := by positivity
    nlinarith

/-- the only lattice vectors shorter than `√(s + 1/s)` are `±(1,0)` and `±(0,1)` -/
lemma short_vectors (s : ℝ) (hs : 0 < s) (hs1 : s < 1) (hs3 : 1 < 3 * s ^ 2) (m n : ℤ)
    (h0 : (m, n) ≠ (0, 0)) (h : D s m n < s + 1 / s) :
    (m = 0 ∧ (n = 1 ∨ n = -1)) ∨ (n = 0 ∧ (m = 1 ∨ m = -1)) := by
  unfold D at h
  by_cases hn : n = 0
  · subst hn
    right
    refine ⟨rfl, ?_⟩
    have hm0 : m ≠ 0 := by intro hm; exact h0 (by simp [hm])
    -- s m² < s + 1/s < 4 s  (as 1/s < 3 s)
    have h3 : s⁻¹ < 3 * s := by rw [inv_lt_iff_one_lt_mul₀ hs]; nlinarith
    have hlt : (m : ℝ) ^ 2 < 4 := by
      have : s * (m : ℝ) ^ 2 < 4 * s := by simp at h; linarith
      nlinarith
    have : m ^ 2 < 4 := by exact_mod_cast hlt
    have : -2 < m ∧ m < 2 := by constructor <;> nlinarith
    omega
  · left
    have hn1 : (1 : ℝ) ≤ (n : ℝ) ^ 2 := one_le_sq_int n hn
    have hmsq : 0 ≤ (m : ℝ) ^ 2 := by positivity
    have h1s : 1 / s ≤ (n : ℝ) ^ 2 / s := by gcongr
    have hm : m = 0 := by
      by_contra hm
      have : (1 : ℝ) ≤ (m : ℝ) ^ 2 := one_le_sq_int m hm
      nlinarith
    refine ⟨hm, ?_⟩
    subst hm
    have : (n : ℝ) ^ 2 / s < s + 1 / s := by simpa using h
    have hlt : (n : ℝ) ^ 2 < s ^ 2 + 1 := by
      have := (div_lt_iff₀ hs).1 this
      have e : (s + 1 / s) * s = s ^ 2 + 1 := by field_simp
      linarith
    have : (n : ℝ) ^ 2 < 2 := by nlinarith
    have : n ^ 2 < 2 := by exact_mod_cast this
    have : -2 < n ∧ n < 2 := by constructor <;> nlinarith
    omega

/-- collinear case: if `B² = x y` and `z = x + y - 2B` (three parallel vectors `w₁, w₂, w₂ - w₁`),
then `x y z ≥ 4 m³` for every common lower bound `m ≥ 0`. -/
lemma collinear_bound (x y z B m : ℝ) (hm : 0 ≤ m) (hx : m ≤ x) (hy : m ≤ y) (hz : m ≤ z)
    (hB : B ^ 2 = x * y) (hzB : z = x + y - 2 * B) : 4 * m ^ 3 ≤ x * y * z := by
  have hx0 : 0 ≤ x := le_trans hm hx
  have hy0 : 0 ≤ y := le_trans hm hy
  set α := Real.sqrt x with hα
  set β := Real.sqrt y with hβ
  have hαx : α ^ 2 = x := Real.sq_sqrt hx0
  have hβy : β ^ 2 = y := Real.sq_sqrt hy0
  have hα0 : 0 ≤ α := Real.sqrt_nonneg _
  have hβ0 : 0 ≤ β := Real.sqrt_nonneg _
  set μ := Real.sqrt m with hμ
  have hμm : μ ^ 2 = m := Real.sq_sqrt hm
  have hμ0 : 0 ≤ μ := Real.sqrt_nonneg _
  have hμα : μ ≤ α := Real.sqrt_le_sqrt hx
  have hμβ : μ ≤ β := Real.sqrt_le_sqrt hy
  -- B = ± α β
  have hB' : B = α * β ∨ B = -(α * β) := by
    have : B ^ 2 = (α * β) ^ 2 := by rw [hB, mul_pow, hαx, hβy]
    exact sq_eq_sq_iff_eq_or_eq_neg.mp this
  rcases hB' with h | h
  · -- z = (α - β)², and μ ≤ |α - β|
    have hz' : z = (α - β) ^ 2 := by rw [hzB, h, ← hαx, ← hβy]; ring
    have hμd : μ ^ 2 ≤ (α - β) ^ 2 := by rw [hμm, ← hz']; exact hz
    have hd : μ ≤ |α - β| := by
      have := Real.sqrt_le_sqrt hμd
      rwa [Real.sqrt_sq hμ0, Real.sqrt_sq_eq_abs] at this
    -- the larger of α, β is at least 2 μ
    have hmax : 2 * μ ≤ max α β := by
      rcases le_total α β with hab | hab
      · rw [abs_of_nonpos (by linarith)] at hd; rw [max_eq_right hab]; linarith
      · rw [abs_of_nonneg (by linarith)] at hd; rw [max_eq_left hab]; linarith
    have key : 4 * μ ^ 6 ≤ α ^ 2 * β ^ 2 * (α - β) ^ 2 := by
      rcases le_total α β with hab | hab
      · rw [max_eq_right hab] at hmax
        have h1 : μ ^ 2 ≤ α ^ 2 := by nlinarith
        have h2 : 4 * μ ^ 2 ≤ β ^ 2 := by nlinarith
        have := mul_le_mul (mul_le_mul h1 h2 (by positivity) (by positivity)) hμd (by positivity) (by positivity)
        nlinarith
      · rw [max_eq_left hab] at hmax
        have h1 : μ ^ 2 ≤ β ^ 2 := by nlinarith
        have h2 : 4 * μ ^ 2 ≤ α ^ 2 := by nlinarith
        have := mul_le_mul (mul_le_mul h2 h1 (by positivity) (by positivity)) hμd (by positivity) (by positivity)
        nlinarith
    calc 4 * m ^ 3 = 4 * μ ^ 6 := by rw [← hμm]; ring
      _ ≤ α ^ 2 * β ^ 2 * (α - β) ^ 2 := key
      _ = x * y * z := by rw [hαx, hβy, hz']
  · -- z = (α + β)² ≥ 4 μ²
    have hz' : z = (α + β) ^ 2 := by rw [hzB, h, ← hαx, ← hβy]; ring
    have h1 : μ ^ 2 ≤ α ^ 2 := by nlinarith
    have h2 : μ ^ 2 ≤ β ^ 2 := by nlinarith
    have h3 : 4 * μ ^ 2 ≤ (α + β) ^ 2 := by nlinarith
    have := mul_le_mul (mul_le_mul h1 h2 (by positivity) (by positivity)) h3 (by positivity) (by positivity)
    calc 4 * m ^ 3 = μ ^ 2 * μ ^ 2 * (4 * μ ^ 2) := by rw [← hμm]; ring
      _ ≤ α ^ 2 * β ^ 2 * (α + β) ^ 2 := this
      _ = x * y * z := by rw [hαx, hβy, hz']

/-- **Part B.** In the rectangular lattice every triangle `0, w₁, w₂` with pairwise distinct
vertices has squared side-length product at least `min (4 s³) (s + 1/s)`. -/
theorem rect_lower (s : ℝ) (hs : 0 < s) (hs1 : s < 1) (hs3 : 1 < 3 * s ^ 2)
    (m₁ n₁ m₂ n₂ : ℤ) (h1 : (m₁, n₁) ≠ (0, 0)) (h2 : (m₂, n₂) ≠ (0, 0))
    (h12 : (m₁, n₁) ≠ (m₂, n₂)) :
    min (4 * s ^ 3) (s + 1 / s) ≤ D s m₁ n₁ * D s m₂ n₂ * D s (m₂ - m₁) (n₂ - n₁) := by
  have h3 : (m₂ - m₁, n₂ - n₁) ≠ (0, 0) := by
    intro h; apply h12; simp only [Prod.mk.injEq] at h ⊢; constructor <;> omega
  have d1 := D_nonneg s hs m₁ n₁
  have d2 := D_nonneg s hs m₂ n₂
  have d3 := D_nonneg s hs (m₂ - m₁) (n₂ - n₁)
  set det : ℤ := m₁ * n₂ - n₁ * m₂ with hdet
  by_cases hcol : det = 0
  · -- collinear: all three vectors are parallel
    have hdR : (m₁ : ℝ) * n₂ - (n₁ : ℝ) * m₂ = 0 := by
      have : ((m₁ * n₂ - n₁ * m₂ : ℤ) : ℝ) = 0 := by exact_mod_cast hcol
      push_cast at this; linarith
    have hL := lagrange s hs m₁ n₁ m₂ n₂
    rw [hdR] at hL
    have hz : D s (m₂ - m₁) (n₂ - n₁) = D s m₁ n₁ + D s m₂ n₂ - 2 * Bf s m₁ n₁ m₂ n₂ := by
      unfold D Bf; push_cast; ring
    set m := min (D s m₁ n₁) (min (D s m₂ n₂) (D s (m₂ - m₁) (n₂ - n₁)))
    have hmin_s : s ≤ m := le_min (D_ge s hs hs1.le _ _ h1) (le_min (D_ge s hs hs1.le _ _ h2) (D_ge s hs hs1.le _ _ h3))
    have hcb := collinear_bound (D s m₁ n₁) (D s m₂ n₂) (D s (m₂ - m₁) (n₂ - n₁)) (Bf s m₁ n₁ m₂ n₂) m
      (by linarith) (min_le_left _ _) (le_trans (min_le_right _ _) (min_le_left _ _))
      (le_trans (min_le_right _ _) (min_le_right _ _)) (by linarith [hL]) hz
    have : 4 * s ^ 3 ≤ 4 * m ^ 3 := by gcongr
    exact le_trans (min_le_left _ _) (le_trans this hcb)
  · -- non-collinear: every pair of sides has |det| ≥ 1, and some side is long
    have hd1 : (1 : ℝ) ≤ ((m₁ : ℝ) * n₂ - (n₁ : ℝ) * m₂) ^ 2 := by
      have := one_le_sq_int det hcol
      simpa [hdet] using this
    -- the same determinant for the pairs (w₁, w₂ - w₁) and (w₂, w₂ - w₁)
    have p12 : 1 ≤ D s m₁ n₁ * D s m₂ n₂ := by
      rw [lagrange s hs]; nlinarith [sq_nonneg (Bf s m₁ n₁ m₂ n₂)]
    have p13 : 1 ≤ D s m₁ n₁ * D s (m₂ - m₁) (n₂ - n₁) := by
      rw [lagrange s hs]
      have : ((m₁ : ℝ) * ((n₂ - n₁ : ℤ) : ℝ) - (n₁ : ℝ) * ((m₂ - m₁ : ℤ) : ℝ)) ^ 2 = ((m₁ : ℝ) * n₂ - (n₁ : ℝ) * m₂) ^ 2 := by
        push_cast; ring
      nlinarith [sq_nonneg (Bf s m₁ n₁ (m₂ - m₁) (n₂ - n₁))]
    have p23 : 1 ≤ D s m₂ n₂ * D s (m₂ - m₁) (n₂ - n₁) := by
      rw [lagrange s hs]
      have : ((m₂ : ℝ) * ((n₂ - n₁ : ℤ) : ℝ) - (n₂ : ℝ) * ((m₂ - m₁ : ℤ) : ℝ)) ^ 2 = ((m₁ : ℝ) * n₂ - (n₁ : ℝ) * m₂) ^ 2 := by
        push_cast; ring
      nlinarith [sq_nonneg (Bf s m₂ n₂ (m₂ - m₁) (n₂ - n₁))]
    -- some side has D ≥ s + 1/s
    have hlong : s + 1 / s ≤ D s m₁ n₁ ∨ s + 1 / s ≤ D s m₂ n₂ ∨ s + 1 / s ≤ D s (m₂ - m₁) (n₂ - n₁) := by
      by_contra hcon
      push_neg at hcon
      obtain ⟨a1, a2, a3⟩ := hcon
      have c1 := short_vectors s hs hs1 hs3 _ _ h1 a1
      have c2 := short_vectors s hs hs1 hs3 _ _ h2 a2
      have c3 := short_vectors s hs hs1 hs3 _ _ h3 a3
      apply hcol
      rcases c1 with ⟨e1, f1 | f1⟩ | ⟨e1, f1 | f1⟩ <;> rcases c2 with ⟨e2, f2 | f2⟩ | ⟨e2, f2 | f2⟩ <;>
        simp only [hdet, e1, f1, e2, f2] at c3 ⊢ <;> omega
    rcases hlong with h | h | h
    · have : s + 1 / s ≤ D s m₁ n₁ * (D s m₂ n₂ * D s (m₂ - m₁) (n₂ - n₁)) := by nlinarith
      exact le_trans (min_le_right _ _) (by linarith [mul_assoc (D s m₁ n₁) (D s m₂ n₂) (D s (m₂ - m₁) (n₂ - n₁))])
    · have : s + 1 / s ≤ D s m₂ n₂ * (D s m₁ n₁ * D s (m₂ - m₁) (n₂ - n₁)) := by nlinarith
      exact le_trans (min_le_right _ _) (by nlinarith)
    · have : s + 1 / s ≤ D s (m₂ - m₁) (n₂ - n₁) * (D s m₁ n₁ * D s m₂ n₂) := by nlinarith
      exact le_trans (min_le_right _ _) (by nlinarith)

/-! ## Part A: the upper bound over reduced bases -/

/-- squared side-length product of the triangle `(0, u, v)` for `u = (a,0)`, `v = (p, 1/a)` -/
noncomputable def Φ (a p : ℝ) : ℝ := a ^ 2 * (p ^ 2 + 1 / a ^ 2) * ((a - p) ^ 2 + 1 / a ^ 2)

/-- the key identity: `Φ = a² + a⁻² + a² m (m - 2/a²)` with `m = p (a - p)` -/
lemma Φ_eq (a p : ℝ) (ha : 0 < a) :
    Φ a p = a ^ 2 + 1 / a ^ 2 + a ^ 2 * (p * (a - p)) * (p * (a - p) - 2 / a ^ 2) := by
  unfold Φ; field_simp; ring

lemma Φ_le (a p : ℝ) (ha : 0 < a) (ha4 : a ^ 4 ≤ 4 / 3) (hp0 : 0 ≤ p) (hp : p ≤ a / 2) :
    Φ a p ≤ a ^ 2 + 1 / a ^ 2 := by
  rw [Φ_eq a p ha]
  have hm0 : 0 ≤ p * (a - p) := by nlinarith
  have hm1 : p * (a - p) ≤ a ^ 2 / 4 := by nlinarith [sq_nonneg (a - 2 * p)]
  have ha2 : 0 < a ^ 2 := by positivity
  have hk : a ^ 2 / 4 < 2 / a ^ 2 := by
    rw [div_lt_div_iff₀ (by norm_num) ha2]; nlinarith
  have : p * (a - p) - 2 / a ^ 2 < 0 := by linarith
  have : a ^ 2 * (p * (a - p)) * (p * (a - p) - 2 / a ^ 2) ≤ 0 :=
    mul_nonpos_of_nonneg_of_nonpos (by positivity) this.le
  linarith

lemma Φ_eq_iff (a p : ℝ) (ha : 0 < a) (ha4 : a ^ 4 ≤ 4 / 3) (hp0 : 0 ≤ p) (hp : p ≤ a / 2)
    (h : Φ a p = a ^ 2 + 1 / a ^ 2) : p = 0 := by
  rw [Φ_eq a p ha] at h
  have ha2 : 0 < a ^ 2 := by positivity
  have hm1 : p * (a - p) ≤ a ^ 2 / 4 := by nlinarith [sq_nonneg (a - 2 * p)]
  have hk : a ^ 2 / 4 < 2 / a ^ 2 := by
    rw [div_lt_div_iff₀ (by norm_num) ha2]; nlinarith
  have hneg : p * (a - p) - 2 / a ^ 2 < 0 := by linarith
  have : a ^ 2 * (p * (a - p)) * (p * (a - p) - 2 / a ^ 2) = 0 := by linarith
  rcases mul_eq_zero.1 this with h1 | h1
  · rcases mul_eq_zero.1 h1 with h2 | h2
    · exact absurd h2 ha2.ne'
    · rcases mul_eq_zero.1 h2 with h3 | h3
      · exact h3
      · nlinarith
  · exact absurd h1 hneg.ne

/-- facts about `u* = s*²` with `4 u*² = u* + 1` -/
lemma ustar_lt (u : ℝ) (hu : 0 < u) (h : 4 * u ^ 2 = u + 1) : u < 3 / 4 := by nlinarith

/-- **Part A.** Let `s* > 0` with `4 s*⁴ = s*² + 1` (i.e. `s* = a*²`).  For every reduced basis
`(a, 0), (p, 1/a)` with `0 ≤ p ≤ a/2` and `a⁴ ≤ 4/3`:
`min (4 a⁶) (Φ a p) ≤ 4 s*³`. -/
theorem upper_bound (S a p : ℝ) (hS : 0 < S) (hSeq : 4 * S ^ 4 = S ^ 2 + 1)
    (ha : 0 < a) (ha4 : a ^ 4 ≤ 4 / 3) (hp0 : 0 ≤ p) (hp : p ≤ a / 2) :
    min (4 * a ^ 6) (Φ a p) ≤ 4 * S ^ 3 := by
  have hS1 : S < 1 := by nlinarith [ustar_lt (S ^ 2) (by positivity) (by nlinarith)]
  have hP : 4 * S ^ 3 = S + 1 / S := by field_simp; nlinarith
  set t := a ^ 2 with ht
  have ht0 : 0 < t := by positivity
  have hΦ := Φ_le a p ha ha4 hp0 hp
  rcases le_or_gt t S with h | h
  · -- a² ≤ s*: the collinear triple
    have : t ^ 3 ≤ S ^ 3 := by gcongr
    have e : a ^ 6 = t ^ 3 := by rw [ht]; ring
    exact le_trans (min_le_left _ _) (by rw [e]; linarith)
  · -- a² > s*: the right triangle, and t + 1/t < s + 1/s
    refine le_trans (min_le_right _ _) (le_trans hΦ ?_)
    rw [hP]
    have hu := ustar_lt (S ^ 2) (by positivity) (by nlinarith)
    rcases le_or_gt t 1 with h1 | h1
    · -- s* < t ≤ 1 : t + 1/t is decreasing on (0,1]
      have : t + 1 / t - (S + 1 / S) = (t - S) * (1 - 1 / (t * S)) := by field_simp; ring
      have h2 : 1 / (t * S) ≥ 1 := by
        rw [ge_iff_le, le_div_iff₀ (by positivity)]; nlinarith
      nlinarith
    · -- 1 < t with t² ≤ 4/3 : (t + 1/t)² ≤ 49/12 < (s + 1/s)²  (as s*² < 3/4)
      have ht2 : t ^ 2 ≤ 4 / 3 := by rw [ht]; nlinarith
      have hsq1 : (t + 1 / t) ^ 2 ≤ 49 / 12 := by
        have e : (t + 1 / t) ^ 2 = t ^ 2 + 2 + 1 / t ^ 2 := by field_simp; ring
        rw [e]
        have : 1 / t ^ 2 ≥ 3 / 4 := by rw [ge_iff_le, le_div_iff₀ (by positivity)]; nlinarith
        -- f(q) = q + 1/q on q = t² ∈ (1, 4/3] is at most 4/3 + 3/4
        have hq : t ^ 2 + 1 / t ^ 2 ≤ 4 / 3 + 3 / 4 := by
          have e2 : t ^ 2 + 1 / t ^ 2 - (4 / 3 + 3 / 4) = (t ^ 2 - 4 / 3) * (1 - 3 / (4 * t ^ 2)) := by
            field_simp; ring
          have : 1 - 3 / (4 * t ^ 2) ≥ 0 := by
            rw [ge_iff_le, sub_nonneg, div_le_one (by positivity)]; nlinarith
          nlinarith
        linarith
      have hsq2 : 49 / 12 < (S + 1 / S) ^ 2 := by
        have e : (S + 1 / S) ^ 2 = S ^ 2 + 2 + 1 / S ^ 2 := by field_simp; ring
        rw [e]
        have : 25 / 12 < S ^ 2 + 1 / S ^ 2 := by
          have hS2 : 0 < S ^ 2 := by positivity
          rw [← sub_pos]
          have e2 : S ^ 2 + 1 / S ^ 2 - 25 / 12 = (12 * (S ^ 2) ^ 2 - 25 * S ^ 2 + 12) / (12 * S ^ 2) := by
            field_simp; ring
          rw [e2]; apply div_pos _ (by positivity); nlinarith
        linarith
      have h0 : 0 ≤ t + 1 / t := by positivity
      have h0' : 0 ≤ S + 1 / S := by positivity
      nlinarith [sq_nonneg (t + 1 / t - (S + 1 / S))]

/-- **Part A, equality.** Equality forces `a² = s*` and `p = 0` (the rectangle `Λ*`). -/
theorem upper_bound_eq (S a p : ℝ) (hS : 0 < S) (hSeq : 4 * S ^ 4 = S ^ 2 + 1)
    (ha : 0 < a) (ha4 : a ^ 4 ≤ 4 / 3) (hp0 : 0 ≤ p) (hp : p ≤ a / 2)
    (h : min (4 * a ^ 6) (Φ a p) = 4 * S ^ 3) : a ^ 2 = S ∧ p = 0 := by
  have hS1 : S < 1 := by nlinarith [ustar_lt (S ^ 2) (by positivity) (by nlinarith)]
  have hP : 4 * S ^ 3 = S + 1 / S := by field_simp; nlinarith
  have hΦ := Φ_le a p ha ha4 hp0 hp
  have hge1 : 4 * S ^ 3 ≤ 4 * a ^ 6 := by rw [← h]; exact min_le_left _ _
  have hge2 : 4 * S ^ 3 ≤ Φ a p := by rw [← h]; exact min_le_right _ _
  -- a² ≥ S from the collinear term
  have hta : S ≤ a ^ 2 := by
    by_contra hlt; push_neg at hlt
    have : (a ^ 2) ^ 3 < S ^ 3 := by gcongr
    nlinarith
  -- a² ≤ S from the right triangle (strict decrease/increase argument of `upper_bound`)
  have hts : a ^ 2 = S := by
    rcases eq_or_lt_of_le hta with e | hlt
    · exact e.symm
    · exfalso
      -- rerun the strict part of `upper_bound`: Φ ≤ a² + a⁻² < S + 1/S
      have hstrict : a ^ 2 + 1 / a ^ 2 < S + 1 / S := by
        set t := a ^ 2 with ht
        have ht0 : 0 < t := by positivity
        have hu := ustar_lt (S ^ 2) (by positivity) (by nlinarith)
        rcases le_or_gt t 1 with h1 | h1
        · have e : t + 1 / t - (S + 1 / S) = (t - S) * (1 - 1 / (t * S)) := by field_simp; ring
          have h2 : 1 / (t * S) > 1 := by
            rw [gt_iff_lt, lt_div_iff₀ (by positivity)]; nlinarith
          nlinarith
        · have ht2 : t ^ 2 ≤ 4 / 3 := by rw [ht]; nlinarith
          have hsq1 : (t + 1 / t) ^ 2 ≤ 49 / 12 := by
            have e : (t + 1 / t) ^ 2 = t ^ 2 + 2 + 1 / t ^ 2 := by field_simp; ring
            rw [e]
            have hq : t ^ 2 + 1 / t ^ 2 ≤ 4 / 3 + 3 / 4 := by
              have e2 : t ^ 2 + 1 / t ^ 2 - (4 / 3 + 3 / 4) = (t ^ 2 - 4 / 3) * (1 - 3 / (4 * t ^ 2)) := by
                field_simp; ring
              have : 1 - 3 / (4 * t ^ 2) ≥ 0 := by
                rw [ge_iff_le, sub_nonneg, div_le_one (by positivity)]; nlinarith
              nlinarith
            linarith
          have hsq2 : 49 / 12 < (S + 1 / S) ^ 2 := by
            have e : (S + 1 / S) ^ 2 = S ^ 2 + 2 + 1 / S ^ 2 := by field_simp; ring
            rw [e]
            have : 25 / 12 < S ^ 2 + 1 / S ^ 2 := by
              have hS2 : 0 < S ^ 2 := by positivity
              rw [← sub_pos]
              have e2 : S ^ 2 + 1 / S ^ 2 - 25 / 12 = (12 * (S ^ 2) ^ 2 - 25 * S ^ 2 + 12) / (12 * S ^ 2) := by
                field_simp; ring
              rw [e2]; apply div_pos _ (by positivity); nlinarith
            linarith
          nlinarith [sq_nonneg (t + 1 / t - (S + 1 / S)), (by positivity : (0:ℝ) ≤ t + 1 / t),
            (by positivity : (0:ℝ) ≤ S + 1 / S)]
      linarith
  refine ⟨hts, Φ_eq_iff a p ha ha4 hp0 hp ?_⟩
  have : a ^ 2 + 1 / a ^ 2 = 4 * S ^ 3 := by rw [hP, hts]
  linarith

/-- For `s = s*` both terms of `rect_lower` coincide with `4 s*³`, so the extremal lattice attains
`P*²` (the upper bound of Part A). -/
theorem rect_lower_star (S : ℝ) (hS : 0 < S) (hSeq : 4 * S ^ 4 = S ^ 2 + 1)
    (m₁ n₁ m₂ n₂ : ℤ) (h1 : (m₁, n₁) ≠ (0, 0)) (h2 : (m₂, n₂) ≠ (0, 0))
    (h12 : (m₁, n₁) ≠ (m₂, n₂)) :
    4 * S ^ 3 ≤ D S m₁ n₁ * D S m₂ n₂ * D S (m₂ - m₁) (n₂ - n₁) := by
  have hu := ustar_lt (S ^ 2) (by positivity) (by nlinarith)
  have hS1 : S < 1 := by nlinarith
  have hS3 : 1 < 3 * S ^ 2 := by nlinarith
  have hP : 4 * S ^ 3 = S + 1 / S := by field_simp; nlinarith
  have := rect_lower S hS hS1 hS3 m₁ n₁ m₂ n₂ h1 h2 h12
  rwa [← hP, min_self] at this

end TriangleLattice
