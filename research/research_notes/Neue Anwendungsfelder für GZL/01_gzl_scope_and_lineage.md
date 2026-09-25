# GZL (Graph Zeta Library): exact scope, limits, roadmap and research lineage

Status of evidence: everything marked "verified" below was checked against primary material on 2026-09-25: the v1.0.0 wheel source (local copy of `gzl-1.0.0-py3-none-any.whl`, identical module layout to https://github.com/graph-zeta/gzl), the GitHub README / DOCUMENTATION.md / CHANGELOG.md / AGENTS.md / CITATION.cff on `main`, the arXiv HTML full texts of 2609.18918v1 and 2609.18761v1, the arXiv author feed of A. A. Buchheit, arXiv author searches, and the Semantic Scholar API. "Local run" means a small check executed with gzl 1.0.0 installed from PyPI. Where a GitHub `blob/main` link is given for a source file, the content was read from the v1.0.0 wheel; the link points to the matching file in the repository.

---

## Q1. What is the exact definition of the graph zeta function ζ_G(k; ν)? (distinct sites? pin? momentum? x = 0? which d?)

### Takeaway
ζ_G is an **unconstrained ("softcore") sum** over the placements of all vertices but one pinned vertex on an infinite Bravais lattice Λ = A·Z^d. It is a product of per-edge kernels K_e(x_{e+} − x_{e−}) times a single phase e^{−2πi(x_t − x_s)·k} between two terminals s and t. The power-law kernel is regularized as K_ν(0) = 0, so only the two endpoints of the *same edge* are kept off one site. Non-adjacent vertices may coincide. The physics paper handles the physical hardcore (distinct-site) constraint with a separate combinatorial "softcore mapping" that is already baked into the corpus prefactors. The library officially supports d = 1, 2, 3. Larger d runs for some blocks but is neither advertised nor calibrated.

### Cited Findings
- **Definition (README, verified):** G = (V, E, ν, s, t) is a finite multigraph with two distinguished terminals s, t, per-edge exponents ν_e > 0 and oriented edges. Λ = A·Z^d, and the columns of A ∈ R^{d×d} are the primitive vectors. ζ_G(k) = Σ_{x_v ∈ Λ, v ≠ p} e^{−2πi (x_t − x_s)·k} Π_{e∈E} K_{ν_e}(x_{e+} − x_{e−}). The pin p is arbitrary, and by translation invariance x_p = 0. "ζ_G is the lattice Fourier transform of the two-point correlator" with respect to the s–t displacement. A sufficient condition for absolute convergence is Re ν_e > d on every edge, or summability of K_e — [GZL README, "Definition of the graph zeta function"](https://github.com/graph-zeta/gzl#quick-start)
- **Regularization at x = 0 (verified):** K_ν(x) = |x|^{−ν} for x ≠ 0 and K_ν(0) = 0, "which sets the coincident-point self-interaction to zero." The general kernel K(x) = a(x) + Σ_j b_j K_{ν_j}(x) has a compact part a for which "a(0) is allowed and simply weights configurations whose two endpoints coincide" — [README](https://github.com/graph-zeta/gzl#quick-start); [gzl/interaction.py module docstring](https://github.com/graph-zeta/gzl/blob/main/gzl/interaction.py)
- **Vertices are NOT required to occupy distinct sites (verified):** the source says that on a torus that "cannot hold" a block, "its vertices cannot be placed on them with the two ends of every edge on distinct sites … The kernel vanishes at the origin, so every term of the sum vanishes." The only exclusion is therefore per edge, through K(0) = 0. The shipped corpora are named `tfim_softcore_corpus_0qp/1qp.npz` — [gzl/frontend.py, `NPointsRequiredError` / `_torus_holds`](https://github.com/graph-zeta/gzl/blob/main/gzl/frontend.py); [gzl/data/README.md](https://github.com/graph-zeta/gzl/tree/main/gzl/data)
- **Hardcore → softcore mapping (physics paper):** the linked-cluster embedding sums are "hardcore" sums, where no two vertices may occupy the same site. Duft et al. develop "an exact softcore mapping that removes explicit summation constraints … and expresses the perturbative series coefficients in terms of unconstrained graph lattice sums." The mapping is "purely combinatorial and independent of the choice of the interaction kernel" and gives "a linear map between the set of hardcore prefactors v_G and softcore prefactors v_G^(s) at each perturbative order" — [arXiv:2609.18761, Sec. III and VIII](https://arxiv.org/html/2609.18761)
- **No hardcore-mapping tool ships in GZL v1.0.0:** a grep for "hardcore" over the v1.0.0 package source finds nothing. The corpora ship already mapped, as softcore prefactors a_r(G) — [local source check; gzl/data/README.md](https://github.com/graph-zeta/gzl/tree/main/gzl/data)
- **Mathematical-paper definition (more general than the library):** a *graph lattice sum* 𝒵_{Λ,G}(k) is defined for any absolutely summable kernels K_e ∈ ℓ¹(Λ), with k ∈ BZ = A^{−T}T^d and an arbitrary pinning node p. A *graph zeta function* is the special case K_e = a_e + Σ_j b_{ej} K_{ν_ej}, where a_e: Λ → C is even and compactly supported and b_{ej} ∈ C. Summability corresponds to min Re ν_{ej} > d. The single-edge case (bridge) gives the Epstein zeta function Z_{Λ,ν}(k). On Z at k = 0 it is 2ζ(ν), twice the Riemann zeta function — [arXiv:2609.18918, Def. 2.3, 2.7](https://arxiv.org/html/2609.18918)
- **Elementary properties (paper, Lemma 2.4):** the sum is absolutely convergent with |𝒵| ≤ Π_{e∈T}‖K_e‖_1 Π_{e∉T}‖K_e‖_∞ for any spanning tree T. It is independent of the pin. For s = t it is k-independent ("vacuum", 0qp). Remark 2.5: ℓ¹ can be relaxed to ℓ^∞ on the chords if some spanning tree has summable kernels — [arXiv:2609.18918, §2.4–2.5](https://arxiv.org/html/2609.18918)
- **Edge merging (Hadamard):** parallel edges merge into one edge whose kernel is the pointwise product of the parallel kernels, with reflection for opposite orientation. For power laws K_ν·K_ν' = K_{ν+ν'} — [arXiv:2609.18918, Lemma 2.6](https://arxiv.org/html/2609.18918); [README step 1](https://github.com/graph-zeta/gzl#quick-start)
- **Momentum convention (verified):** `momentum` is given in fractional reciprocal-lattice coordinates β ∈ [0,1)^d, with β = A^T k. A physical wavevector maps as `momentum = A.T @ k_phys / (2π)`. Only the blocks on the s–t path of the block-cut tree (the "spine") carry k. Off-spine blocks contribute ζ_B(0) — [DOCUMENTATION.md, "Finite external momentum"](https://github.com/graph-zeta/gzl/blob/main/DOCUMENTATION.md)
- **Realness:** for Bravais lattices with real even kernels, ζ_G(k) = ζ_G(−k) ∈ R, and the library returns floats. An internal check raises instead of discarding an imaginary part "above the round-off", which it attributes to "a basis of several sites per cell … or a complex kernel such as a Peierls phase" — [gzl/_real.py docstring](https://github.com/graph-zeta/gzl/blob/main/gzl/_real.py); [evaluate_graph docstring](https://github.com/graph-zeta/gzl/blob/main/gzl/frontend.py)
- **Dimensions:** "Supported are general 1D, 2D, and 3D lattices", and the CHANGELOG says "1D, 2D and 3D Bravais lattices." The named lattices are exactly `chain`, `square`, `triangular` and `cubic`, but any non-singular square matrix A is accepted — [README](https://github.com/graph-zeta/gzl); [CHANGELOG 1.0.0](https://github.com/graph-zeta/gzl/blob/main/CHANGELOG.md); [gzl/_lattices.py](https://github.com/graph-zeta/gzl/blob/main/gzl/_lattices.py)
- **d ≥ 4 is not calibrated but partly runs:** `zeta_circle` covers d = 1, 2, 3 only, and "evaluate_graph sends a cycle in d ≥ 4 to the algebra or the tensor instead." Dense blocks at d ≥ 4 go to the real-space box. The calibration tables `_DENSE_ENGINE_BY_D`, `_SPLIT_SP_N_POINTS` and `_CORE_KAPPA` contain only d = 1, 2, 3, and a missing d means "not calibrated here" — [DOCUMENTATION.md API table](https://github.com/graph-zeta/gzl/blob/main/DOCUMENTATION.md); [AGENTS.md](https://github.com/graph-zeta/gzl/blob/main/AGENTS.md); [gzl/frontend.py](https://github.com/graph-zeta/gzl/blob/main/gzl/frontend.py)
- **Local run (gzl 1.0.0), for information only:** a bridge on Z^4 with ν = 6 returned 14.8298 and a bridge on Z^5 with ν = 7 returned 18.8361, both Epstein-zeta closed forms. A diamond graph on Z^4 at ν = 6 with n_points = 6 returned 42.63 with no accuracy guarantee. A triangle on Z^4 needed n_points > 0 because no closed-form cycle exists in d = 4 — local check, [gzl 1.0.0 on PyPI](https://pypi.org/project/gzl/)
- **The corpus is dimension-agnostic:** "The corpus stores combinatorial data alone and is therefore dimension-agnostic … one corpus drives every (Λ, d, ν)" — [README, "Evaluating a corpus"](https://github.com/graph-zeta/gzl#quick-start)

### Inferences
- Anyone whose target sum has a hardcore/distinct-site constraint, such as a cluster expansion or a Mayer-type sum with excluded volume, must first produce a softcore reformulation. The kernel-independent mapping is described in 2609.18761 Sec. III but is not shipped as code in GZL 1.0.0.
- The object is a *two-point* Fourier transform with one external momentum. Sums needing several independent external momenta, such as 2qp scattering or vertex functions, fall outside the definition as implemented (see Q2).
- d = 4+ bridges are exact, because EpsteinLib handles general d. Anything else at d ≥ 4 is uncalibrated and should be treated as unsupported for publication-grade numbers.

### Gaps
- Whether the hardcore→softcore mapping code exists in some unreleased repository of the Schmidt group could not be verified.
- The arXiv HTML of 2609.18761 Sec. III shows the mapping only by example (4-cycle). The general combinatorial rule was not fully extracted here.

---

## Q2. What exactly is supported and what is not? (non-Bravais, anisotropy, complex/odd/oscillatory/Yukawa kernels, finite lattices, slabs, mixed exponents, ν ≤ d)

### Takeaway
GZL 1.0.0 supports **infinite, single-site-per-cell (Bravais) lattices in d = 1–3**, with **real, even kernels** of the form (compact, possibly anisotropic, table a) + (finite sum of *isotropic* power laws b_j|x|^{−ν_j}). Mixed exponents and kernels per edge are allowed. The library supports 0qp (vacuum) and 1qp (one free terminal) objects and exponents ν > d, except for bridges, which accept analytic continuation.

It does **not** support:
- multi-site bases;
- anisotropic power-law tails, such as angle-dependent dipolar kernels;
- complex, odd or Peierls-phase kernels;
- genuinely infinite-range non-power-law tails (oscillatory RKKY, exponentially decaying Yukawa) except as truncated compact tables;
- finite or open lattices and physical slab geometries (`gzl.slab` is a numerical engine, not a slab geometry);
- two or more free terminals;
- non-bridge blocks with ν ≤ d.

### Cited Findings
**Lattices**
- The mathematical paper restricts to "a periodic point set Λ = A·Z^d … a monoatomic or Bravais lattice" — [arXiv:2609.18918, §2.1](https://arxiv.org/html/2609.18918)
- Multi-atomic lattices are future work. The method paper says "The method will also be adapted to multi-atomic lattices", and the physics paper lists "extension to multi-atomic lattices" under "Future developments" — [2609.18918 §10](https://arxiv.org/html/2609.18918); [2609.18761 §VIII](https://arxiv.org/html/2609.18761)
- The code treats "a basis of several sites per cell, whose sublattice phases survive inversion" as a violated hypothesis and raises if an imaginary part appears — [gzl/_real.py](https://github.com/graph-zeta/gzl/blob/main/gzl/_real.py)
- A local run with a non-square A (shape 1×2) raised a ValueError ("A must be a square (d, d) matrix"). A must be square and non-singular, so there are no embedded lower-dimensional lattices in higher-dimensional space — local check; [evaluate_graph docstring](https://github.com/graph-zeta/gzl/blob/main/gzl/frontend.py)

**Kernel class**
- `Interaction`: K(x) = a(x) + Σ_j b_j K_{ν_j}(x). Here a is "a real, EVEN, compactly supported function on the Bravais lattice", stored as a finite table over integer labels. The compact part "may be anisotropic" because it is a function of the displacement vector. The power-law part is the isotropic |x|^{−ν} — [gzl/interaction.py](https://github.com/graph-zeta/gzl/blob/main/gzl/interaction.py); [DOCUMENTATION.md "General interactions"](https://github.com/graph-zeta/gzl/blob/main/DOCUMENTATION.md)
- **Complex kernels are refused:** "K(x) must be real-valued … An Interaction is a real kernel." A local run with `exp(1j x)` raised a ValueError — [gzl/interaction.py `_refuse_complex`](https://github.com/graph-zeta/gzl/blob/main/gzl/interaction.py); local check
- **Odd kernels are refused:** "Compact parts must be even. The router discards edge orientation … and the finite-momentum engines use a real cos weight." A local run with an asymmetric table raised "the compact part must be even" — [gzl/interaction.py](https://github.com/graph-zeta/gzl/blob/main/gzl/interaction.py); local check
- The math paper's Def. 2.7 allows complex a_e and b_ej, and Def. 2.3 allows general ℓ¹ kernels. **The library is narrower than the theory** — [2609.18918 Def. 2.3, 2.7](https://arxiv.org/html/2609.18918)
- **Anisotropic power laws are NOT in v1.0.0.** For KTmSe₂ in 3D, "we used an experimental code extension to anisotropic dipolar sums" at order 5 only. "Higher orders are in principle also accessible for the angle-dependent dipolar kernel, but a stable and efficient implementation in the Graph Zeta Library is still under development" — [2609.18761 §VII.2](https://arxiv.org/html/2609.18761)
- The algebra's `GraphZeta` representation "may change when anisotropic Epstein zeta functions are included", and the anisotropic extension is therefore on the roadmap — [DOCUMENTATION.md "API stability / Provisional"](https://github.com/graph-zeta/gzl/blob/main/DOCUMENTATION.md)
- The physics paper states that the core elements "apply without modification to any absolutely summable, translationally invariant and even interaction kernel and only particular block evaluations must be adapted accordingly." This is a claim about the method, not about the v1.0.0 code — [2609.18761 §VIII](https://arxiv.org/html/2609.18761)
- `Interaction.from_function(K, A, radius, …)` samples a user function on |A m| ≤ radius, origin included, and stores it as the compact table. Per-edge sequences mixing floats and `Interaction`s are accepted. Also available: `from_shells`, `from_table`, `nearest_neighbour` and products `K1*K2`, `K**m` — [DOCUMENTATION.md API table](https://github.com/graph-zeta/gzl/blob/main/DOCUMENTATION.md)
- A block with a compact part runs on a torus of at least n_v·R + 2 points per dimension, where R is the support radius. `InteractionSupportError` is raised when an engine window cannot hold the support — [DOCUMENTATION.md "General interactions"](https://github.com/graph-zeta/gzl/blob/main/DOCUMENTATION.md)

**Mixed exponents per edge (verified supported)**
- `nu` may be "a per-edge float array … or a per-edge sequence mixing floats and Interactions." A local run of a triangle on the square lattice with ν = (3, 4, 5) returned 8.1027 — [evaluate_graph docstring](https://github.com/graph-zeta/gzl/blob/main/gzl/frontend.py); local check
- ν = ∞ is the nearest-neighbour indicator and gives exact integer embedding counts. The NN triangle on the triangular lattice returned 12.0 in a local run — [DOCUMENTATION.md evaluate_corpus example](https://github.com/graph-zeta/gzl/blob/main/DOCUMENTATION.md); local check

**ν ≤ d**
- "The lattice sum converges absolutely if ν_e > d on every edge, and for other exponents ζ_G is defined by meromorphic continuation. `evaluate_graph` returns this continuation for bridges … and raises `UnsupportedLatticeSumError` at its pole ν = d … Every other block with an exponent ν ≤ d after the Hadamard merge raises `UnsupportedLatticeSumError` as well. This is a limitation of the implementation and not a statement about convergence … no evaluator of GZL is validated there" — [DOCUMENTATION.md](https://github.com/graph-zeta/gzl/blob/main/DOCUMENTATION.md)
- In a local run, a triangle at ν = 0.9 on the chain raised `UnsupportedLatticeSumError`, and a bridge at ν = 0.5 on the chain returned the continuation −2.9207 (= 2ζ(0.5)) — local check
- Accuracy near ν → d is a design goal of the algebra: the algebra "makes accurate evaluation possible even for exponents close to the spatial dimension d." Benchmarks go down to σ = ν − d = 0.5, and MC comparisons use σ ≥ 1/2 because "for σ = 0.1 … Monte Carlo does not provide faithful results anymore" — [README](https://github.com/graph-zeta/gzl); [2609.18918 §8.3](https://arxiv.org/html/2609.18918)
- The algebra handles the Γ-pole family ν = d + 2n automatically, by compression or by shifting 10⁻⁶ off the pole — [DOCUMENTATION.md](https://github.com/graph-zeta/gzl/blob/main/DOCUMENTATION.md); [2609.18918 §6.6–6.7](https://arxiv.org/html/2609.18918)

**Terminals, sectors and momentum**
- "Two or more free terminals" raise `UnsupportedRequestError` (a `NotImplementedError`). So does a 1qp graph whose s–t path runs through a block that mixes ν = ∞ with finite exponents; the workaround is `Interaction.nearest_neighbour(A)` — [DOCUMENTATION.md](https://github.com/graph-zeta/gzl/blob/main/DOCUMENTATION.md); local check (`terminal=[1,2]` raised)
- Shipped corpora: the LRTFIM only, 0qp to order 13 (8,403 graphs) and 1qp to order 11 (22,677 graphs), plus a 31,080-graph topology snapshot — [gzl/data/README.md](https://github.com/graph-zeta/gzl/tree/main/gzl/data)

**Finite, open and slab geometries**
- All sums run over the infinite lattice Λ. The n_points torus is a numerical truncation, and the documented numerical values are "the graph zeta function … with … n_points grid points per dimension" carrying a truncation error that decreases with n_points — [DOCUMENTATION.md "Numerical values"](https://github.com/graph-zeta/gzl/blob/main/DOCUMENTATION.md)
- **`gzl.slab` is NOT a slab or film geometry.** It is "a dense core contracted with a second pin … It evaluates the vacuum (k = 0) value of one dense block … not a general-purpose engine." It is the same torus truncation associated differently, trading peak memory — [gzl/slab.py docstring](https://github.com/graph-zeta/gzl/blob/main/gzl/slab.py); [DOCUMENTATION.md API table](https://github.com/graph-zeta/gzl/blob/main/DOCUMENTATION.md)
- A MultiGraph with an isolated node is refused "because its infinite-lattice sum diverges" — [DOCUMENTATION.md API table](https://github.com/graph-zeta/gzl/blob/main/DOCUMENTATION.md)

**Engines and cost model (what is computed how)**
- Router for each biconnected block after the Hadamard merge — [DOCUMENTATION.md](https://github.com/graph-zeta/gzl/blob/main/DOCUMENTATION.md):
  - bridge → closed-form `epstein_zeta(ν, A, 0, k)`;
  - simple cycle at k = 0 → `zeta_circle`, a Brillouin-zone integral of products of Epstein zetas, d ≤ 3;
  - SP-reducible block (tw ≤ 2) with σ < 1.49 → semi-analytic algebra (Epstein-zeta singular part plus Fourier part on the n_points grid);
  - everything else → torus tensor network (hybrid FFT and dense core), with optional 3-point Richardson extrapolation.
- Tensor-network cost: O(N^tw) operations and O(N^{tw−1}) memory with N = n^d. If every maximal step is convolutional, this drops to O(N^{tw−1} log N) — [2609.18918, Thm 7.5](https://arxiv.org/html/2609.18918)
- "At d = 3 the edge-difference table caps practical use at n_points ≤ 16" for the treewidth-agnostic tensor path — [DOCUMENTATION.md "What this gives you"](https://github.com/graph-zeta/gzl/blob/main/DOCUMENTATION.md)

### Inferences
- **RKKY (cos(2k_F r)/r³):** only representable in v1.0.0 by truncating the oscillatory part into a finite compact table, optionally plus an isotropic power-law tail such as a non-oscillating envelope. The oscillating tail beyond the table radius cannot be written as Σ b_j|x|^{−ν_j}. A faithful RKKY treatment is therefore not supported and would be a methodological extension. A caution from a local run: `from_function` includes the origin, so a user function singular at x = 0 must set a(0) = 0 explicitly. A naive r^{−3} function produced ~10²⁷ from the x = 0 sample.
- **Yukawa / exponential decay:** truncation to a compact table is legitimate when the decay length is short relative to the table radius. Truncation error is the user's responsibility; there is no built-in tail estimate. The ℓ¹ theory of 2609.18918 covers it, but no dedicated evaluator exists.
- **Non-Bravais lattices** (honeycomb, kagome, pyrochlore, real multi-sublattice compounds) are out of scope for v1.0.0. The only workaround is a Bravais sublattice or a Bravais superlattice with the model rewritten; for genuinely multi-sublattice problems this generally fails without new theory.
- **Angle-dependent dipolar interactions** (LiHoF₄, Fe₈, Mn₁₂, RE(OH)₃, as listed by the authors) are explicitly named by the authors as a target. That strongly suggests they will pursue it themselves, so it is a crowded direction for novelty.
- 0qp and 1qp observables only: ground-state energy and 1qp dispersion. Spectral weights, 2qp or bound-state sectors, and finite temperature are not available (see Q4).

### Gaps
- No quantitative accuracy statement was found for compact kernels with a large support radius R, where the grid is lifted to n_v·R + 2. How cost grows with R was not benchmarked in the papers read.
- It was not checked whether complex-valued b_j (allowed in theory) could be emulated by linearity: the sum is multilinear in edge kernels, so real and imaginary parts could be split per edge at a cost of 2^|E| terms. This is an inference that was not verified with the library.

---

## Q3. What accuracy and runtime claims and benchmarks are in arXiv:2609.18918 and 2609.18761?

### Takeaway
The headline claim is that full LRTFIM series take **minutes on a laptop instead of ~10⁴–10⁵ core-hours** of cluster Monte Carlo, with a full-BZ grid for the cost of one momentum. Across 998 coefficients (114 series) the results agree with published MC within statistical errors (max 3.1 σ_MC, median 0.12 σ_MC). Series-parallel blocks cost linear time in |V|, with algebraic error decay ~n^{−(d+σ+2)}.

### Cited Findings
- Abstract of the method paper: the method "reduces the evaluation time for state-of-the art series expansions from tenthousands of core-hours to minutes on a single core with controlled numerical precision". A 3D dispersion on a 4096-point momentum grid takes "within ten minutes on a single core" — [arXiv:2609.18918 abstract](https://arxiv.org/abs/2609.18918)
- **MC comparison (Table 2):** 114 perturbative series with 998 coefficients on the chain, square, triangular and cubic lattices, 0qp and ferro/antiferro 1qp. Excluding the analytic orders and 17 coefficients of two underconverged reference top orders, the maximum error is 3.1 σ_MC and the median 0.12 σ_MC. Including all data, 89% are within 1 σ_MC and 97% within 2 σ_MC. Grid sizes were n = 2048 (d = 1), 48 (d = 2) and 10 (d = 3). The complete set was obtained "within 5 minutes on 8 cores of an Apple M1 Max", against about 24 h on 72 cores × ~20 seeds (~3×10⁴ core-hours) per MC series and momentum — [2609.18918 §8.3](https://arxiv.org/html/2609.18918)
- Orders r ≤ 3 are fully analytic (bridges and cycles, machine precision), so the deviations there measure the MC error alone. Evidence is given that the top-order MC coefficients in Fey's thesis for the square and triangular lattices are "underconverged": two MC runs differ by 25 σ_MC at order 10, and the uncertainties jump by a factor of 75 — [2609.18918 §8.3](https://arxiv.org/html/2609.18918)
- **Series-parallel accuracy:** cycles with |V| = 4, 8, 12 on Z and the triangular lattice were compared against the exact many-body (circle) zeta at n = 64 as a function of σ. The algebra is used for σ < 3/2, and Fourier discretization for σ > 3/2 to avoid Γ-pole cancellation. Error grows "only linearly with the number of operations" for iterated series/parallel compositions — [2609.18918 §8.1, Figs. 2–3](https://arxiv.org/html/2609.18918)
- **Error scaling:** for uniform ν = d + σ, the error of a series-parallel block scales as n^{−(d+σ+2)}, uniformly in k (Remark 6.7). This is confirmed numerically in Fig. 5 for σ = 0.5 and 1 — [2609.18918 §6.7, §8.1](https://arxiv.org/html/2609.18918)
- **Runtime scaling:** "the scaling should reduce from N^{|V|−1} to |V| N log N for all series-parallel graphs". This is confirmed as linear in |V| on an Apple M1 Max core at n = 64 — [2609.18918 §8.1, Fig. 4](https://arxiv.org/html/2609.18918)
- **tw ≥ 3 validation:** the box-truncated tensor network matches an independent explicit nested sum for K₄ and K₅ on Z and Z² at σ = 0.5, 1, 2, for L up to 28 (2.1×10⁸ terms), with maximum deviation 5.3×10⁻¹⁵. The torus tensor network for K₄ is benchmarked with and without Richardson extrapolation (Fig. 6). "An analytic treatment of the leading order truncation tails … would remove the need for extrapolation. This is, however, nontrivial … and will therefore be the topic of future work" — [2609.18918 §8.2](https://arxiv.org/html/2609.18918)
- **Physics result (method paper §9):** the 3D LRTFIM on the cubic lattice at σ = 1, antiferromagnetic λ = −0.03, order 10 (7,136 graphs, up to 30-dimensional sums) on a 16³ grid took about 10 min on a single M1 Max core. It shows a non-analytic |k|^σ dispersion at Γ ("anomalous dispersion"), traced analytically to the Epstein-zeta bridges — [2609.18918 §9, Figs. 8–9](https://arxiv.org/html/2609.18918)
- **Physics paper runtimes:** "we obtain the 0qp and 1qp series coefficients for the three-dimensional system in around 200 seconds" on a single Apple M4 core, at σ = 0.5 and orders up to 13/11 (Fig. 11) — [2609.18761 §VI.2](https://arxiv.org/html/2609.18761)
- **Physics outputs:** λ_c and zν from DlogPadé of the 1qp gap at order 11, for ferro- and antiferromagnetic LRTFIM on four lattices, with dense σ sampling "not previously achievable", consistent with earlier MC work — [2609.18761 §VI.4, Figs. 14–15](https://arxiv.org/html/2609.18761)
- **KTmSe₂:** three 2D models were compared at order 11 against INS data (Zheng et al., PRB 108, 054435 (2023)): J₁–J₂ with literature parameters, J₁–J₂ refit, and J₁ + dipolar. "The J₁-dipolar model provides the best overall description." The 3D ABC-stack dipolar model reached order 5 only, via the experimental anisotropic extension. The 3D fit gives J₁ ≃ 0.040 meV and h ≃ 0.534 meV, and the RMS deviation of the 2D J₁–J₂ model is about twice that of the dipolar models — [2609.18761 §VII](https://arxiv.org/html/2609.18761)
- **README example:** 1qp cubic, ν = 4, order 9, n = 8 runs "within seconds on a single core" and agrees with Fey's thesis Table F.10 within statistical errors. The README figure is an AFM dispersion at order 10 on a 16³ grid in about 10 min on one core — [README](https://github.com/graph-zeta/gzl)
- **Documented numerical fine print:** the `accuracy="floor"` corpus mode reports median relative truncation of 1.4e-05–3.3e-05 for the p = 7.5 class and 2.3e-03–5.6e-03 for the p = 4 class that "dominates order ≥ 10" at n = 8. Hybrid and tensor engines differ by up to 2·10⁻⁴ of the coefficient scale at order 11 on the chain at ν = 3, n = 64 — [evaluate_graph docstring](https://github.com/graph-zeta/gzl/blob/main/gzl/frontend.py); [DOCUMENTATION.md "Hybrid engine"](https://github.com/graph-zeta/gzl/blob/main/DOCUMENTATION.md)

### Inferences
- "Numerically exact" in the physics paper means free of sampling error, with algebraic convergence in n_points. At high orders in d = 3 the tw ≥ 3 blocks carry relative errors of about 10⁻³ at the default n = 8 (docstring numbers), which is the weakest link. Any application needing high precision at high order in 3D should budget for larger n or wait for the promised analytic tail treatment of tensor blocks.

### Gaps
- Exact per-order error tables of 2609.18761 (Fig. 10 relative deviations) were not transcribed.
- Memory limits for high-treewidth cores in d = 3 were not given as a single number in the papers. The code refuses with `HybridCoreTooLargeError` or `SlabTooLargeError` and falls back to the box.

---

## Q4. What open problems and outlook do the authors state, and what is on the roadmap?

### Takeaway
The stated roadmap covers anisotropic (angle-dependent) Epstein-zeta kernels, multi-atomic lattices, more models (Heisenberg named first), frustrated systems, other observables (spectral weights, higher quasiparticle sectors), finite temperature, correlated disorder, and an analytic singularity treatment for the tensor-network blocks. There are no GitHub issues yet, and issue creation is currently restricted.

### Cited Findings
- **Method paper outlook:** "Future work will be dedicated to accelerating numerical convergence of the algebra further by including anisotropic generalized zeta functions … The method will also be adapted to multi-atomic lattices … This will include frustrated systems, where Monte Carlo methods can be ineffective. Finally, extensions to other models, other perturbative approaches, to finite temperature [Burkard, Schneider, Sbierski 2026], as well as to systems with correlated disorder are worth exploring" — [2609.18918 §10](https://arxiv.org/html/2609.18918)
- **Physics paper outlook:** "Future developments include the extension to multi-atomic lattices, broader classes of interaction kernels, convergence improvements based on Epstein zeta derivatives, and refinements of the tensor-network architecture, including a singularity treatment comparable to the zeta algebra … We also aim to apply the method to additional observables, such as spectral weights and higher quasiparticle sectors, to other quantum models and lattice architectures, and finite temperatures" — [2609.18761 §VIII](https://arxiv.org/html/2609.18761)
- **Target platforms named by the authors:** Rydberg arrays with van-der-Waals LRTFIM, where dynamical structure factors are now measurable. Also 3D dipolar TFIM materials LiHoF₄, Fe₈, Mn₁₂ acetate and RE(OH)₃, whose "full dipolar coupling is angle dependent" and which require "the appropriate anisotropic lattice Fourier transform while adapting the semi-analytic zeta algebra" — [2609.18761 §VIII](https://arxiv.org/html/2609.18761)
- **Tensor-network roadmap:** "The full Γ-prefactor multiplication / convolution algebra … extends to N-terminal tensors and would give analytically-exact contractions on this path. That extension is the natural follow-up work" — [DOCUMENTATION.md "What this gives you"](https://github.com/graph-zeta/gzl/blob/main/DOCUMENTATION.md)
- **README:** "The method is general and additional default models, such as the Heisenberg model, will be included as time progresses" — [README](https://github.com/graph-zeta/gzl)
- **Anisotropic support is already in the dependency chain:** arXiv:2609.28282 (Buchheit, Busse, 23 Sep 2026), "Computation of anisotropic singular sums from high-order derivatives of Epstein zeta functions", is cited in the README as EpsteinLib's "recent anisotropic extension". The physics paper notes that anisotropic lattice Fourier transforms "are already implemented in EpsteinLib" — [README references](https://github.com/graph-zeta/gzl); [arXiv listing](https://arxiv.org/a/buchheit_a_1); [2609.18761 §VIII](https://arxiv.org/html/2609.18761)
- **Provisional API** (may change in a minor release): `Interaction`, the whole semi-analytic algebra (`GraphZeta` etc.; "representation … may change when anisotropic Epstein zeta functions are included"), the graph constructors and the engines. Stable: `evaluate_graph`, `evaluate_corpus`, `compute_series_coefficients`, `zeta_circle`, and the corpus file format — [DOCUMENTATION.md "API stability"](https://github.com/graph-zeta/gzl/blob/main/DOCUMENTATION.md)
- **CHANGELOG** contains only 1.0.0 (2026-09-25). Release candidates 1.0.0rc2 exist, as GitHub release pages and PRs #1 and #2 by andreasbuchheit — [CHANGELOG](https://github.com/graph-zeta/gzl/blob/main/CHANGELOG.md); [GitHub release rc2](https://github.com/graph-zeta/gzl/releases/tag/v1.0.0rc2)
- **GitHub issues:** 0 issues, 1 star and 0 forks at time of fetch. "Issue creation is restricted" — [github.com/graph-zeta/gzl/issues](https://github.com/graph-zeta/gzl/issues?q=is%3Aissue)
- **Licensing and contributions:** AGPL-3.0-or-later. Contributions are licensed MIT and AGPL, with a DCO sign-off required. Development used Claude Code — [README](https://github.com/graph-zeta/gzl)
- **Funding:** Klaus-Tschira Stiftung grant 00.025.2025 — [README Acknowledgements](https://github.com/graph-zeta/gzl)

### Inferences
- The authors have publicly claimed these directions: anisotropic dipolar, multi-atomic and frustrated lattices, Heisenberg/XXZ corpora, spectral weights, higher-qp sectors, finite temperature (with Sbierski-type high-T expansions), and correlated disorder. Candidate "new application fields" in these directions face high scoop risk from the Buchheit–Schmidt collaboration.
- Less claimed or unmentioned directions: classical statistical-mechanics cluster or virial expansions on lattices, lattice Feynman diagrams, crystal-lattice many-body energies beyond circle graphs (hinted at via the Lennard-Jones/ATM chemistry reference), and using the tw-based block algebra outside physics. These are candidates for novelty checks but need their own verification.

### Gaps
- No public ROADMAP.md exists (404). There were no discussions or issues to mine.
- A Heisenberg corpus timeline was not given.

---

## Q5. Prior-work lineage (Buchheit / Keßler / Busse zeta-function line, and Schmidt-group long-range linked-cluster line)

### Takeaway
GZL fuses two lines. The first is Buchheit and co-workers' numerical analysis of singular lattice sums: the singular Euler–Maclaurin expansion, EpsteinLib, the many-body (circle-graph) Epstein zeta method, anisotropic and periodic-BC zeta expansions. The second is the Schmidt group's pCUT/linked-cluster expansions for long-range quantum magnets, whose embedding sums were previously evaluated by Monte Carlo (Fey & Schmidt 2016 onwards). The direct precursor is the circle-zeta integral representation in Buchheit & Busse, "Epstein zeta method for many-body lattice sums".

### Cited Findings
**Buchheit line (arXiv IDs verified from the arXiv author feed)**
- A. A. Buchheit, T. Keßler, "Singular Euler–Maclaurin expansion", arXiv:2003.12422, published as J. Sci. Comput. 90, 53 (2022) ("On the efficient computation of large scale singular sums with applications to long-range forces in crystal lattices") — [arXiv feed](https://arxiv.org/a/buchheit_a_1); [GZL README refs](https://github.com/graph-zeta/gzl)
- A. A. Buchheit, T. Keßler, "Singular Euler–Maclaurin expansion on multidimensional lattices", arXiv:2102.10941, Nonlinearity 35, 3706 (2022) — [arXiv feed](https://arxiv.org/a/buchheit_a_1)
- A. A. Buchheit, T. Keßler, P. K. Schuhmacher, B. Fauseweh, "Exact continuum representation of long-range interacting systems and emerging exotic phases in unconventional superconductors", arXiv:2201.11101, Phys. Rev. Res. 5, 043065 (2023) — [arXiv feed](https://arxiv.org/a/buchheit_a_1)
- A. A. Buchheit, T. Keßler, K. Serkh, "On the computation of lattice sums without translational invariance", arXiv:2403.03213, Math. Comp. 94, 2533 (2025) — [arXiv feed](https://arxiv.org/a/buchheit_a_1)
- A. A. Buchheit, J. K. Busse, R. Gutendorf, "Computation and properties of the Epstein zeta function with applications to quantum systems", arXiv:2412.16317, IMA J. Numer. Anal. drag057 (2026). This is the EpsteinLib paper, and EpsteinLib is GZL's core dependency (`epsteinlib>=0.5`; PyPI latest 0.6.2) — [arXiv feed](https://arxiv.org/a/buchheit_a_1); [gzl METADATA](https://pypi.org/project/gzl/); [PyPI epsteinlib](https://pypi.org/project/epsteinlib/)
- A. A. Buchheit, J. K. Busse, "Epstein zeta method for many-body lattice sums", arXiv:2504.11989, Numer. Math. (2026). It gives the circle zeta function as a BZ integral of Epstein zetas, with cost linear in the number of nodes, and is the direct precursor ("[8, Theorem 2.6]" in 2609.18918) — [arXiv feed](https://arxiv.org/a/buchheit_a_1); [2609.18918 §1.1–1.2](https://arxiv.org/html/2609.18918)
- A. A. Buchheit, J. K. Busse, T. Keßler, F. N. Rybakov, "Zeta expansion for long-range interactions under periodic boundary conditions with applications to micromagnetics", arXiv:2509.26274, J. Comput. Phys. 559, 114885 (2026). It is anisotropic-kernel groundwork cited for the future anisotropic algebra — [Busse arXiv feed](https://arxiv.org/a/busse_j_1); [2609.18918 §10](https://arxiv.org/html/2609.18918)
- A. A. Buchheit, J. K. Busse, "Computation of anisotropic singular sums from high-order derivatives of Epstein zeta functions", arXiv:2609.28282 (23 Sep 2026) — [arXiv feed](https://arxiv.org/a/buchheit_a_1)
- A. Robles-Navarro, S. Cooper, A. A. Buchheit, J. K. Busse, A. Burrows, O. Smits, P. Schwerdtfeger, "Exact lattice summations for Lennard-Jones potentials coupled to a three-body Axilrod–Teller–Muto term applied to cuboidal phase transitions", J. Chem. Phys. 163, 094104 (2025). **No arXiv ID found**, and it is not in Buchheit's arXiv feed. It is cited as a crystal-stability application of circle-graph sums — [GZL README refs](https://github.com/graph-zeta/gzl); [2609.18918 §1](https://arxiv.org/html/2609.18918)
- Other Buchheit-coauthored works seen as citing papers: D. Haink, A. A. Buchheit, C. Weitenberg, B. Fauseweh, "Nonlocal edge mode hybridization in the long-range interacting Kitaev chain", arXiv:2509.26447. A. A. Buchheit, T. Keßler, S. Rjasanow, "Numerical Solution of the BCS Equation for Unconventional Superconductors", arXiv:2602.15911 — [Semantic Scholar API citations](https://api.semanticscholar.org/graph/v1/paper/arXiv:2403.03213)

**Schmidt-group long-range linked-cluster line**
- K. Coester, K. P. Schmidt, "Optimizing linked-cluster expansions by white graphs", PRE 92, 022118 (2015). This is the white-graph expansion GZL's corpora rely on. arXiv ID not verified here — [GZL README refs](https://github.com/graph-zeta/gzl)
- S. Fey, K. P. Schmidt, "Critical behavior of quantum magnets with long-range interactions in the thermodynamic limit", arXiv:1606.05111, PRB 94, 075156 (2016). It introduced linked-cluster expansions for long-range models (53 citations on Semantic Scholar) — [arXiv search](https://arxiv.org/abs/1606.05111); [Semantic Scholar](https://api.semanticscholar.org/graph/v1/paper/arXiv:1606.05111)
- S. Fey, S. C. Kapfer, K. P. Schmidt, "Quantum criticality of two-dimensional quantum magnets with long-range interactions", arXiv:1802.06684, PRL 122, 017203 (2019). This introduced MC evaluation of the embedding sums, and its supplementary tables ship as GZL reference data — [GZL NOTICE](https://github.com/graph-zeta/gzl/blob/main/NOTICE); [arXiv:1802.06684](https://arxiv.org/abs/1802.06684)
- S. Fey, PhD thesis, FAU Erlangen-Nürnberg (2020), "Investigation of zero-temperature transverse-field Ising models with long-range interactions" (Appendix F MC tables) — [open.fau.de/handle/openfau/13818](https://open.fau.de/handle/openfau/13818)
- P. Adelhardt, J. A. Koziol, A. Schellenberger, K. P. Schmidt, long-range anisotropic XY chain in a transverse field, arXiv:2007.16128, PRB 102, 174424 (2020) — [arXiv search](https://arxiv.org/abs/2007.16128)
- A. Langheld, J. A. Koziol, P. Adelhardt, S. C. Kapfer, K. P. Schmidt, "Scaling at quantum phase transitions above the upper critical dimension", arXiv:2203.08081, SciPost Phys. 13, 088 (2022), with Zenodo data 10.5281/zenodo.6645107 — [arXiv:2203.08081](https://arxiv.org/abs/2203.08081); [GZL NOTICE](https://github.com/graph-zeta/gzl/blob/main/NOTICE)
- P. Adelhardt, K. P. Schmidt, "Continuously varying critical exponents in long-range quantum spin ladders", arXiv:2209.01182, SciPost Phys. 15, 087 (2023) — [arXiv:2209.01182](https://arxiv.org/abs/2209.01182)
- J. A. Koziol, G. Morigi, K. P. Schmidt, "Quantum phases of hardcore bosons with repulsive dipolar density-density interactions on two-dimensional lattices", arXiv:2311.10632. It cites Buchheit–Keßler, so the two lines were already in contact in 2023 — [arXiv:2311.10632](https://arxiv.org/abs/2311.10632); [Semantic Scholar citations of 2102.10941](https://api.semanticscholar.org/graph/v1/paper/arXiv:2102.10941)
- J. A. Koziol, M. Mühlhauser, K. P. Schmidt, "Order-by-disorder and long-range interactions in the antiferromagnetic TFIM on the triangular lattice — a perturbative point of view", arXiv:2402.10584, Results in Physics 61, 107794 (2024) — [arXiv:2402.10584](https://arxiv.org/abs/2402.10584)
- P. Adelhardt, J. A. Koziol, A. Langheld, K. P. Schmidt, "Monte Carlo based techniques for quantum magnets with long-range interactions", arXiv:2403.00421, Entropy 26, 401 (2024). This is the MC-embedding review that GZL replaces — [arXiv:2403.00421](https://arxiv.org/abs/2403.00421)
- P. Adelhardt, A. Duft, K. P. Schmidt, "Quantum-critical and dynamical properties of the XXZ bilayer with long-range interactions", arXiv:2408.13145, PRB 111, 024409 (2025) — [arXiv:2408.13145](https://arxiv.org/abs/2408.13145)
- Related recent Schmidt-group long-range works: J. A. Koziol, A. Langheld, K. P. Schmidt, "Melting of devil's staircases in the long-range Dicke-Ising model", arXiv:2503.02734. J. A. Koziol, K. P. Schmidt, "Quantum annealing for lattice models with competing long-range interactions", arXiv:2511.08336. P. Adelhardt, S. R. Muleady, K. P. Schmidt, A. V. Gorshkov, long-range spin-one Heisenberg chain, arXiv:2604.12754. J. A. Koziol, "Quantum criticality of the ferromagnetic Dicke-Ising model", arXiv:2605.27484 — [arXiv author search Koziol](https://arxiv.org/a/koziol_j_1); [arXiv search Adelhardt](https://arxiv.org/abs/2604.12754)
- pCUT foundations cited: Knetter & Uhrig, EPJB 13, 209 (2000); Knetter, Schmidt & Uhrig, J. Phys. A 36, 7889 (2003); Hörmann & Schmidt, "Projective cluster-additive transformation", SciPost Phys. 15, 097 (2023). arXiv IDs were not verified — [2609.18918 references](https://arxiv.org/html/2609.18918)
- Graph-algorithmic ingredients cited in 2609.18918: bucket elimination (Dechter 1999), the generalized distributive law (Aji–McEliece 2000), tensor-network contraction (Markov–Shi 2008), series-parallel/treewidth-2 equivalence (Duffin 1965; Bodlaender 1998) and Hopcroft–Tarjan block-cut (1973) — [2609.18918 references](https://arxiv.org/html/2609.18918)

### Inferences
- Novelty in GZL lies in (i) block factorization of general graph lattice sums with momentum routed along the spine, (ii) the semi-analytic Epstein-plus-Fourier algebra for series-parallel graphs, and (iii) the combination with tensor-network bucket elimination for tw > 2. The single-cycle case (circle zeta) is prior work in 2504.11989.
- A related older idea is coordinate-space evaluation of lattice Feynman diagrams (e.g. hep-lat/9706014, which appeared in the web search). The GZL papers do not cite it. It may be relevant for a novelty check of "graph lattice sums" in lattice field theory, but its relation to GZL was not analysed here.

### Gaps
- The arXiv IDs of Coester & Schmidt (2015), Knetter/Uhrig/Schmidt pCUT papers and Hörmann & Schmidt (2023) were not verified in this pass.
- The Robles-Navarro et al. JCP 2025 paper appears not to be on arXiv. This was not confirmed.

---

## Q6. Who already cites or uses these works? Is GZL used outside the authors' group?

### Takeaway
As of 2026-09-25 there is **no evidence of any use of GZL outside the Buchheit/Schmidt collaboration**, which is expected nine days after the arXiv posting and on release day. The method paper has one citation (the companion paper) and the physics paper none. For the precursor Epstein-zeta papers, the only external citing group identified is the Sbierski/Lesanovsky Rydberg-thermometry paper (arXiv:2604.22743). The Schmidt group's MC-embedding papers are cited more broadly in long-range quantum-magnet literature (for example Sandvik; Meng; Del Maestro/Melko).

### Cited Findings
- Semantic Scholar citation counts on 2026-09-25 — [Semantic Scholar batch API](https://api.semanticscholar.org/graph/v1/paper/batch):

  | Paper | Citations |
  |---|---|
  | 2609.18918 | 1 (= 2609.18761) |
  | 2609.18761 | 0 |
  | 2609.28282 | 0 |
  | 2504.11989 | 4 |
  | 2412.16317 | 2 |
  | 2403.03213 | 8 |
  | 2201.11101 | 12 |
  | 2003.12422 | 6 |
  | 2102.10941 | 7 |
  | 2509.26274 | 2 |
  | 2403.00421 (MC review) | 23 |
  | 1606.05111 | 53 |
  | 1802.06684 | 39 |

- **Citers of 2504.11989 (Epstein zeta method for many-body sums):** 2609.18761; E. Fitzner, I. Lesanovsky, B. Sbierski, "Thermometry for a kagome lattice dipolar Rydberg simulator", arXiv:2604.22743 (**external group**); and two Buchheit self-citations — [Semantic Scholar citations](https://api.semanticscholar.org/graph/v1/paper/arXiv:2504.11989/citations)
- **Citers of 2412.16317 (EpsteinLib):** only 2609.18761 and 2504.11989, both from the group — [Semantic Scholar](https://api.semanticscholar.org/graph/v1/paper/arXiv:2412.16317/citations)
- **Citers of 2403.03213:** among others, Koziol–Langheld–Schmidt 2503.02734 and Koziol–Morigi–Schmidt 2311.10632 (the Schmidt group), and P. Schwerdtfeger & D. J. Wales, "100 Years of the Lennard-Jones Potential" (2024, chemistry) — [Semantic Scholar](https://api.semanticscholar.org/graph/v1/paper/arXiv:2403.03213/citations)
- **Citers of 2201.11101:** includes external papers such as Diessel, Diehl, Defenu, Rosch, Chiocchetta, "Generalized Higgs mechanism in long-range-interacting quantum systems" (arXiv:2208.10487), and Jefremovas et al., "The role of magnetic dipolar interactions in skyrmion lattices" (arXiv:2407.00539) — [Semantic Scholar](https://api.semanticscholar.org/graph/v1/paper/arXiv:2201.11101/citations)
- **Citers of the MC-embedding review 2403.00421 outside the Schmidt group:** S. Yang, G. Schumm, A. Sandvik, "Dynamic structure factor of a spin-1/2 Heisenberg chain with long-range interactions" (arXiv:2412.15168); J. Zhao, N. Laflorencie, Z.-Y. Meng (arXiv:2406.02685); Z.-J. Fan, C. Zhang, Y. Deng, "Clock factorized QMC for long-range interacting systems" (arXiv:2305.14082); B. Schneider, R. Burkard, B. Olmos, I. Lesanovsky, B. Sbierski (arXiv:2407.18156); M. S. Moeed, C. Pennaforti, A. Del Maestro, R. Melko (arXiv:2601.20058); Sfairopoulos et al. (arXiv:2503.19572); Yi et al. (arXiv:2502.04165) — [Semantic Scholar](https://api.semanticscholar.org/graph/v1/paper/arXiv:2403.00421/citations)
- **Repository uptake:** 1 star, 0 forks, 0 issues. PyPI shows the single release 1.0.0, uploaded 2026-09-25 13:40 UTC — [GitHub issues page](https://github.com/graph-zeta/gzl/issues?q=is%3Aissue); [PyPI gzl JSON](https://pypi.org/pypi/gzl/json)
- **Web search** for "graph zeta" / GZL found only the repository, its releases and PRs, and the arXiv paper. No third-party usage turned up — [WebSearch results: GitHub graph-zeta/gzl, arXiv 2609.18918](https://github.com/graph-zeta/gzl)

### Inferences
- The Sbierski group (Burkard, Schneider, Sbierski; Fitzner, Lesanovsky, Sbierski) is the closest external adopter of the Buchheit zeta tools and is explicitly named in both GZL papers' finite-temperature outlook. It is the most likely first external user and a likely competitor for finite-temperature or high-T graph-expansion applications.
- The long-range QMC and entanglement community (Sandvik, Meng, Deng, Del Maestro/Melko) cites the MC-embedding approach and is a natural audience for benchmarking against GZL series, but none has used GZL yet.

### Gaps
- Google Scholar was not queried; Semantic Scholar may lag by days. A re-check is advisable before any novelty claim.
- GitHub code search for `import gzl` in third-party repositories could not be run, because GitHub API access for this repository was not enabled in the session and the HTML issue scrape returned nothing. This is a gap, though unlikely to change the conclusion given release day.
