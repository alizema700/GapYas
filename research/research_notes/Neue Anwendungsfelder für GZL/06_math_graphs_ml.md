# GZL transfer opportunities in pure/applied mathematics (number theory, lattice sums, modular graph functions, lattice energy optimisation, graph invariants) and ML/signal processing

Scope note: all arXiv IDs below were checked by fetching their arxiv.org abstract pages on 2026-09-25, except where marked "(unverified)". Items marked "OWN COMPUTATION" come from experiments I ran with `pip install gzl` (GZL 1.0.0) and `epsteinlib` 0.6.2 in a scratch venv. They are not literature results. The scripts are summarised inline so they can be reproduced.

## Q1. Which mathematical objects are exactly graph zeta functions (directly or after Fourier transform)? Check against the GZL definition

### Takeaway
GZL's object is ζ_{Λ,G}(k) = Σ_{x_v∈Λ, v≠p} e^{-2πi(x_t−x_s)·k} Π_e K_e(x_{e+}−x_{e−}), with K_e = a_e + Σ_j b_ej|x|^{-ν_ej}, K(0)=0 for the power-law part, and min Re ν > d required. Four families match this definition exactly:
- Epstein zeta functions: one edge.
- Many-body/"circle" zeta functions: cycles.
- 1D symmetrised Mordell–Tornheim sums: cycles on Z.
- The whole class of planar modular graph functions (MGFs) from string theory, evaluated on the momentum lattice Λ_τ = Z+τZ: an MGF of a planar graph Γ equals the graph zeta of the planar dual Γ*, with ν_e = 2a_e.

I confirmed the MGF dictionary numerically at arbitrary τ against Zagier's identity and against the C_{2,2,1} identity. It holds, with two limits. First, most string-theory MGFs have edges with a_e=1 (ν=d=2), which is outside GZL's stated domain. They can only be reached by extrapolating ν→d. Second, modular graph *forms* with a_e≠b_e are anisotropic and are not supported.

### Cited Findings
**GZL definition and scope (from the method paper)**
- GZL definition (Def. 2.3/2.7 of the method paper): "graph lattice sum" Z_{Λ,G}(k) = Σ_{x_v∈Λ, v≠p} e^{−2πi δx_{(s,t)}·k} Π_{e∈E} K_e(δx_e) on a two-terminal multigraph with pinned vertex p. It is called a "graph zeta function" when K_e(x) = a_e(x) + Σ_j b_ej K_{ν_ej}(x), with a_e even and compactly supported, K_ν(x)=|x|^{−ν} and K_ν(0)=0. "The summability condition K_e∈ℓ¹(Λ) then corresponds to min_{e,j} Re(ν_ej) > d." — [arXiv:2609.18918](https://arxiv.org/abs/2609.18918)
- Special cases stated in the paper:
  - A single bridge on Z gives 2ζ(ν). A bridge on Λ gives the Epstein zeta Z_{Λ,ν}(k).
  - Multi-edges merge via Hadamard products (K_ν K_ν' = K_{ν+ν'}).
  - Trees factorise into products of Epstein zetas, so for single power laws ζ_G(k) = Z_{Λ,μ}(k) with μ = Σ_e ν_e.

  — [arXiv:2609.18918](https://arxiv.org/abs/2609.18918)
- The "circle zeta function" is ζ^{(n)}(ν) = Σ'_{x_1..x_{n−1}∈Λ} Π_j |x_j − x_{j−1}|^{−ν_j} (x_0 = x_n = 0). It has the integral representation ζ^{(n)} = V_Λ ∫_BZ Π_i Z_{Λ,ν_i}(k) dk, valid for Re ν_j ≥ d, with meromorphic continuation via a Hadamard integral. It comes from the earlier many-body paper. — [arXiv:2609.18918](https://arxiv.org/abs/2609.18918); [arXiv:2504.11989](https://arxiv.org/abs/2504.11989)
- Algorithmic structure:
  - Block factorisation. Bridges and cycles are analytic blocks.
  - Series-parallel (tw ≤ 2) blocks go through a "zeta algebra": multiplication and convolution on the BZ, with the singular part written as Epstein zeta functions and the rest as a rapidly decaying Fourier series.
  - tw > 2 blocks use tensor-network bucket elimination.
  - An FFT gives the full k-grid.

  — [arXiv:2609.18918](https://arxiv.org/abs/2609.18918)
- Precision reported in the paper:
  - Bridges and cycles are "computed to full precision". In practice this means double precision.
  - The tensor network with box truncation matched an independent nested sum to a largest relative deviation of 5.3×10⁻¹⁵.
  - The zeta-algebra error scales as E_rel(n) ~ n^{−(d+σ+2)}, with σ = ν−d.
  - Stated future work: anisotropic generalised zeta functions and multi-atomic lattices.

  — [arXiv:2609.18918](https://arxiv.org/abs/2609.18918)
- The installed API enforces the domain. `gzl.zeta_circle(nu_vec, A)` docstring: "must satisfy nu_i > d … UnsupportedLatticeSumError: An exponent is <= d"; d = 1, 2, 3 only. `evaluate_graph` accepts arbitrary lattice matrices, per-edge ν, general `Interaction` kernels and a `momentum` argument. — [gzl on PyPI](https://pypi.org/project/gzl/) (GZL 1.0.0, AGPL-3.0, author A. A. Buchheit; docstrings read from the installed package)
- EpsteinLib (same group) is a C library with a Python wrapper, `epsteinlib` 0.6.2. It exposes `epstein_zeta`, `epstein_zeta_reg` and anisotropic variants (`epstein_zeta_aniso`). — [epsteinlib on PyPI](https://pypi.org/project/epsteinlib/); method paper [arXiv:2412.16317](https://arxiv.org/abs/2412.16317)
- New companion paper (2026-09-23): anisotropic Epstein zeta functions are built from high-order k-derivatives of isotropic ones. This is the natural route to non-isotropic edge kernels. — [arXiv:2609.28282](https://arxiv.org/abs/2609.28282)

**Modular graph functions**
- MGFs are "non-holomorphic modular functions associated with Feynman graphs for a conformal scalar field theory on a two-dimensional torus". The same paper introduces single-valued elliptic multiple polylogarithms, "associated with Feynman graphs with vanishing external momenta at all but two vertices", which depend on a point ζ on the elliptic curve. — [arXiv:1512.06779](https://arxiv.org/abs/1512.06779) (anchor, 2015)
- In the tetrahedral MGF family, discrete momenta p_i = m_i + n_i τ are summed over integers, "with restrictions ensuring that all propagators … are non-zero". This matches GZL's convention K(0)=0 on every edge. — [arXiv:1706.01889](https://arxiv.org/abs/1706.01889) (anchor, 2017; search-result summary of the PDF)

**1D analogue: symmetrised Mordell–Tornheim zeta**
- Dobrowolski (2026) evaluates the symmetrised Mordell–Tornheim zeta ζ̄_n(w_1..w_n) = Σ_{a_i∈Z*, a_1+…+a_n=0} 1/|a_1^{w_1}…a_n^{w_n}|. He proves ζ̄_n(1,…,1) = B_n(f'(0),…,f^{(n)}(0)) with f(x) = ln binom(−2x, −x). — [arXiv:2603.20550](https://arxiv.org/abs/2603.20550)
- Nearby 2025–26 work on Mordell–Tornheim sums: Kronecker-limit-type formulas and special values [arXiv:2510.10093](https://arxiv.org/abs/2510.10093), [arXiv:2501.01380](https://arxiv.org/abs/2501.01380); integral analogues [arXiv:2409.19980](https://arxiv.org/abs/2409.19980).
- Kalinin's telescoping over Conway topographs gives "arithmetic proofs for modular graph function identities" and covers Mordell–Tornheim series:
  - Σ_{(x,y)∈SL_+(2,Z)} 1/(‖x‖²‖y‖²‖x+y‖²) = π/4.
  - A short proof of Zagier's identity D_{1,1,1} = 2E(z,3) + π³ζ(3).

  — [arXiv:2510.02082](https://arxiv.org/abs/2510.02082); [arXiv:2410.10884](https://arxiv.org/abs/2410.10884)

**Many-body lattice energies**
- The ATM three-body lattice energy was converted "into an integral involving products of Epstein zeta functions", evaluated "to machine precision within minutes", including meromorphic continuations. — [arXiv:2504.07338](https://arxiv.org/abs/2504.07338); [arXiv:2504.11989](https://arxiv.org/abs/2504.11989)

**OWN COMPUTATION: MGF dictionary check** ([gzl 1.0.0](https://pypi.org/project/gzl/), [epsteinlib 0.6.2](https://pypi.org/project/epsteinlib/))

Setup:
- Λ_τ = A·Z² with A = [[1, Re τ], [0, Im τ]], so that |A(m,n)| = |m+nτ|.
- C_{a,b,c}(τ) := (τ₂/π)^{a+b+c} · `gzl.zeta_circle([2a,2b,2c], A)`.
- E_s(τ) := (τ₂/π)^s · `epstein_zeta(2s, A, 0, 0)`.
- Test points τ ∈ {i, e^{iπ/3}, 0.21+1.37i, −0.4+0.95i}.

Results:
- C_{2,2,1+ε/2} → (2/5)E_5 + ζ(5)/30. The deviation is linear in ε: ≈ 0.9ε relative, i.e. 9×10⁻⁷ at ε=10⁻⁶.
- After a 3-point Richardson extrapolation in ε: relative error 1.1–1.2×10⁻⁸ at all τ tested.
- C_{1+ε,1+ε,1+ε} → E_3 + ζ(3) (Zagier), relative error ≈ 9×10⁻⁷ after extrapolation.
- Each circle evaluation takes ~0.1 s on one core.
- Tetrahedral graph K₄ (tw=3, tensor network) with all ν=4 (MGF weight 12): converges as n_points = 32/48/64 → 2.29641310909 / 2.29641310926 / 2.29641310927 at τ=i, 0.1–0.3 s each.

### Inferences
**Planar MGF ↔ graph zeta dictionary (my derivation, consistent with the numerics above).**
- An MGF is C_Γ(τ) = Σ_{p_e∈Λ_τ\{0}} Π_v δ(Σ_{e∋v} ±p_e) Π_e (τ₂/(π|p_e|²))^{a_e}.
- The momentum-conservation lattice is the cycle lattice of Γ. It can be written as "differences of vertex variables" p_e = x_{f(e)} − x_{f'(e)} exactly when Γ is planar. Then the vertices are the faces of Γ, i.e. the vertices of the dual Γ*. This is Whitney's criterion: the cycle matroid of Γ is graphic iff Γ is planar.
- Hence, for planar Γ: C_Γ(τ) = (τ₂/π)^{Σa_e} · ζ_{Λ_τ, Γ*}(k=0) with ν_e = 2a_e.
- Examples:
  - The dihedral/"banana" MGFs C_{a_1…a_r} (2 vertices, r edges) have the r-cycle as dual, so they are GZL circle zetas.
  - The trihedral MGFs C_{a,b,c; d,e,f; …} are graph zetas of the dual (a series-parallel theta-type graph).
  - The tetrahedral MGF (K₄, self-dual) is ζ_{K₄} (tw=3, tensor-network engine).
- Non-planar MGFs (first appearing at higher loop order, e.g. K_{3,3} topologies) are NOT graph zetas in GZL's sense.
- Consistency check: GZL's BZ-integral representation ∫ Π Z_{Λ,ν}(k) is exactly the MGF "position-space" representation ∫_{torus} Π G_a(z) d²z/τ₂. This holds because the torus Green function G_a(z|τ) = Σ_{p≠0} (τ₂/π|p|²)^a e^{2πi⟨p,z⟩} is an Epstein zeta at wave-vector z.

**GZL's momentum k and elliptic MGFs.** GZL's external momentum k ∈ BZ (the dual torus) plays the role of the torus point z/ζ in elliptic MGFs. The single-edge case is exact: Z_{Λ_τ,2a}(k) is the Kronecker–Eisenstein / Green-function building block. For planar graphs, the two-terminal phase e^{−2πi(x_t−x_s)·k} corresponds to a phase on the Γ-edges cut by a dual path, i.e. momentum inserted at two vertices. That is the 1512.06779 "all but two vertices" setting. Plausible, but I did not verify it beyond the dihedral case. The GZL FFT feature, which returns the whole k-grid in one run, would then give an eMGF on a full grid of z.

**1D graph zetas are Mordell–Tornheim/MZV objects.**
- The GZL circle zeta on Λ=Z is exactly Dobrowolski's ζ̄_n (a_j = x_j − x_{j−1}, Σa_j = 0), with ν_j = w_j.
- Splitting sign regions turns a 1D triangle zeta into symmetrised Tornheim sums T(a,b,c). For integer arguments these reduce to MZVs.
- So d=1 is well understood theoretically. GZL is useful there only as a validation bench, and at the edge w_j = 1 = d it is outside GZL's strict domain.

**Three-body lattice energies are triangle graph zetas.**
- A pure three-body Riesz energy E₃(Λ) = Σ'_{x,y} |x|^{−ν}|y|^{−ν}|x−y|^{−ν} is the triangle (3-cycle) graph zeta.
- In 2D at unit covolume it equals π^{3s}·C_{s,s,s}(τ) with s = ν/2. Minimising a 2D three-body Riesz lattice energy is therefore the same problem as minimising the dihedral MGF C_{s,s,s} over the fundamental domain.
- ATM: 1 + 3cosθ₁cosθ₂cosθ₃ is a polynomial in r_ij² over r_ij² products. The ATM energy is therefore a finite linear combination of triangle graph zetas with per-edge exponents, some of them < d, which need the meromorphic continuation from 2504.11989.

### Gaps
- Could not verify whether GZL's general `evaluate_graph` (not only `zeta_circle`) supports the ν=d borderline via analytic continuation. The documentation/API says exponents ≤ d are rejected. The ε→0 extrapolation works but tops out at ~10⁻⁸ relative, probably because of ε log ε terms. This was not analysed.
- The exact mapping of GZL two-terminal phases to the elliptic MGFs of 1512.06779 and 2511.15883 is a derivation sketch, not checked against a published eMGF value.
- Modular graph *forms* (a_e ≠ b_e, kernels p^{−a} p̄^{−b}) would need the anisotropic Epstein machinery of 2609.28282, which GZL does not yet expose for graphs.

## Q2. Which recent (2025–2026) works state open problems or need high-precision numerics that GZL could supply?

### Takeaway
The strongest fit is lattice-energy optimisation with many-body (graph-structured) power-law interactions:
- 2D/3D three- and four-body Riesz/ATM energies over all lattices.
- The 2026 claimed proof of the Sarnak–Strömbergsson conjecture (FCC minimises the 3D Epstein zeta) invites many-body analogues.
- The Buchheit group itself has only explored the cubic Bain path.

The second-best fit is modular graph functions. GZL can evaluate planar MGFs (including tw=3 tetrahedral ones) at arbitrary τ to 10⁻¹⁰–10⁻¹⁵, including non-integer exponents. However, the MGF community already has exact tools (iterated-Eisenstein-integral representations, Mathematica packages, telescoping proofs) that give arbitrarily many digits. GZL would be a cross-check or an extension (non-integer s, eMGFs on grids, d=3 analogues), not a missing capability.

PSLQ-style discovery is limited by double precision.

### Cited Findings
**Lattice energy optimisation**
- Luo & Wei (2026-09-15) claim a proof of the Sarnak–Strömbergsson conjecture:
  - Among unit-covolume 3D lattices, the theta function Θ(α,L) is minimised by FCC for α>1 and by BCC for α<1.
  - The Epstein zeta E(L,s) is minimised by FCC for all s>3/2.
  - Very recent and unrefereed.

  — [arXiv:2609.17356](https://arxiv.org/abs/2609.17356)
- Luo & Wei classify minimisers of ratios and differences of theta and Epstein zeta functions. The hexagonal lattice is pivotal. They cite applications to crystallisation (Bétermin). — [arXiv:2605.07580](https://arxiv.org/abs/2605.07580)
- Deng & Luo prove hexagonal optimality for Σ|P|⁴e^{−πα|P|²} for α ≥ 3/2, "partially answer[ing] some open questions proposed by Bétermin". — [arXiv:2411.17199](https://arxiv.org/abs/2411.17199)
- Other signs-of-derivatives work: [arXiv:2501.01265](https://arxiv.org/abs/2501.01265)
- Cooper & Schwerdtfeger: the Epstein zeta of a Conway–Sloane lattice family has a local minimum at BCC. — [arXiv:2501.05746](https://arxiv.org/abs/2501.05746)
- LJ + ATM three-body work by the GZL authors' circle (Robles-Navarro, Cooper, Buchheit, Busse, Burrows, Smits, Schwerdtfeger):
  - Scope: Bain-type cuboidal transformations only (fcc–mcc–bcc–acc).
  - The ATM cohesive energy has an extremum at bcc, numerically a minimum for repulsive three-body forces along the Bain path.
  - Strong repulsive ATM can make bcc favourable over fcc, "however… the bcc phase remains susceptible to further cuboidal distortions".

  — [arXiv:2504.07338](https://arxiv.org/abs/2504.07338)
- Buchheit & Busse: the n-body method cuts the ATM lattice sum runtime "from weeks to minutes". They find an fcc→bcc transition with increasing ATM coupling at finite pressure. — [arXiv:2504.11989](https://arxiv.org/abs/2504.11989)
- Bétermin & Furlanetto: numerical minimisation of LJ and of Epstein zeta for p-norms shows "a new and unexpected phase transition for the minimizers with respect to p". — [arXiv:2407.20762](https://arxiv.org/abs/2407.20762)
- Bétermin, Šamaj & Travěnec: general theory of structural transitions (square→rectangular, first/second order, tricritical point) for pair potentials on 2D rectangular lattices. — [arXiv:2312.01395](https://arxiv.org/abs/2312.01395)
- Pair-potential ground-state anchors: [arXiv:2107.14020](https://arxiv.org/abs/2107.14020) (3D Riesz/LJ), [arXiv:2104.09795](https://arxiv.org/abs/2104.09795) (computer-assisted triangular optimality)
- A web search for 2025 work on three-body lattice energy minimisation over all 2D/3D lattices found none. — [search result set incl. 2312.01395, 2107.14020, 2104.09795](https://arxiv.org/abs/2312.01395)

**Modular graph functions / forms (2025–26)**
- Claasen & Doroudiani convert lattice-sum MGFs into iterated Eisenstein integrals, "with a Mathematica package that implements all modular graph form topologies up to four vertices". Used for the α'⁸ ζ₃ζ₅ amplitude coefficient. — [arXiv:2502.05531](https://arxiv.org/abs/2502.05531); earlier seventh-order amplitude: [arXiv:2412.04381](https://arxiv.org/abs/2412.04381)
- Schlotterer, Sohnle & Tao solve the differential equations of elliptic MGFs via equivariant iterated integrals. They construct single-valued elliptic polylogarithms. — [arXiv:2511.15883](https://arxiv.org/abs/2511.15883)
- Fedosova & Klinger-Logan prove convolution identities for complex-indexed divisor sums, which explain non-critical L-values appearing in MGFs. The relevant objects have complex exponents. — [arXiv:2512.21413](https://arxiv.org/abs/2512.21413)
- Dorigoni, Green & Wen write two-loop MGFs and N=4 SYM integrated correlators as 4D lattice sums from theta lifts, with L-values of cusp forms in their Fourier modes. — [arXiv:2501.03996](https://arxiv.org/abs/2501.03996)
- Kleinschmidt, Porkert & Schlotterer: genus-one motivic coaction from zeta generators. — [arXiv:2508.02800](https://arxiv.org/abs/2508.02800)
- Zeta-generator canonicalisation. — [arXiv:2406.05099](https://arxiv.org/abs/2406.05099)
- The frontier has moved to higher genus (flat connections, polylogarithms). There, lattice sums over a single Bravais lattice no longer describe the objects:
  - [arXiv:2608.29169](https://arxiv.org/abs/2608.29169)
  - [arXiv:2607.05656](https://arxiv.org/abs/2607.05656)
  - [arXiv:2602.09108](https://arxiv.org/abs/2602.09108)
  - [arXiv:2501.07640](https://arxiv.org/abs/2501.07640)
- Anchor: D'Hoker & Kaidi obtained and proved "all identities at weight six and all dihedral identities at weight seven". They extended to a hierarchy of inhomogeneous Laplace equations. — [arXiv:1608.04393](https://arxiv.org/abs/1608.04393)
- Anchor: tetrahedral MGFs, inhomogeneous Laplace equations from representation theory. — [arXiv:1706.01889](https://arxiv.org/abs/1706.01889)

**Lattice sums, Euler/Tornheim sums, multiple Eisenstein series**
- Bailey, McPhedran & Salvy give a closed-form algorithm for Euler sums Σ R(k)H_k. "Computation of Euler sums directly to very high precision enables us to rigorously check" formulas. This is the high-precision-numerics culture GZL would need to match. — [arXiv:2604.02384](https://arxiv.org/abs/2604.02384); earlier [arXiv:2311.06294](https://arxiv.org/abs/2311.06294)
- He & Hu: closed-form boundary term for conditionally convergent Madelung direct sums in arbitrary triclinic lattices. — [arXiv:2608.10041](https://arxiv.org/abs/2608.10041)
- Samart's conjecture n₄(81) = 40·L'(g₇,0) status report: the L-value side is proven via lattice sums over the ring of integers of Q(√−7). The Mahler-measure side is obstructed. — [arXiv:2608.02265](https://arxiv.org/abs/2608.02265)
- Multiple Eisenstein series are very active (2026):
  - [arXiv:2601.13626](https://arxiv.org/abs/2601.13626) (symmetric MES; every modular form is a combination of symmetric triple Eisenstein series)
  - [arXiv:2609.03777](https://arxiv.org/abs/2609.03777)
  - [arXiv:2602.08176](https://arxiv.org/abs/2602.08176)

  These are ordered/holomorphic sums (Σ over m_1 ≻ … ≻ m_r of (m_iτ+n_i)^{−k}), not products over graph edges of |x|^{−ν}.

**Feynman periods**
- Portner: canonical graph integrals evaluate to single-valued MZVs, and "every single-valued multiple zeta value occurs as a rational linear combination of Feynman periods". — [arXiv:2607.25595](https://arxiv.org/abs/2607.25595)
- Schnetz: five-twist identity for Feynman periods [arXiv:2505.02578](https://arxiv.org/abs/2505.02578); φ³ at six loops [arXiv:2505.15485](https://arxiv.org/abs/2505.15485); graphical functions and HyperFORM [arXiv:2604.25739](https://arxiv.org/abs/2604.25739), [arXiv:2511.19992](https://arxiv.org/abs/2511.19992)
- These are continuum integrals with propagator exponents below d (e.g. 1/|x|² in d=4).

**ML / signal processing**
- Long-range machine-learning interatomic potentials need reciprocal-space/Ewald sums of pairwise inverse-power kernels. Examples:
  - PSWF-LR, "exponent-aware" treatment of arbitrary 1/r^p channels. — [arXiv:2606.06617](https://arxiv.org/abs/2606.06617)
  - Latent Ewald, "long-range electrostatics … is easier than we thought". — [arXiv:2512.18029](https://arxiv.org/abs/2512.18029)
  - [arXiv:2507.14302](https://arxiv.org/abs/2507.14302)
  - Differentiable PME. — [arXiv:2606.01598](https://arxiv.org/abs/2606.01598)
  - [arXiv:2412.03281](https://arxiv.org/abs/2412.03281)

  All of these are single-edge (two-body) lattice sums in a periodic cell.
- Gerken et al. compute finite-width neural tangent kernels from Feynman diagrams. The graph sums there are not lattice sums. — [arXiv:2508.11522](https://arxiv.org/abs/2508.11522)

### Inferences
**(A) Many-body lattice optimisation: strong bridge, GZL-native.** Each item is exactly a graph zeta with ν > d, the lattice matrix is a free parameter, and evaluation is cheap:
- 3-body Riesz energies = triangle.
- 4-body = K₄ / 4-cycle.
- ATM = combination of triangles.
- Axilrod–Teller–Muto plus higher-order dispersion (4-body dipole terms) = cycles.

The 3D search space is 5-dimensional (lattices modulo rotation and scale), which is small enough for global optimisation. Natural open questions:
- (i) Many-body Sarnak–Strömbergsson: which 3D lattice minimises the triangle zeta for each ν? FCC vs BCC vs others.
- (ii) Is there a first-order hexagonal→square transition in 2D (my numerics below say yes)?
- (iii) Stability of BCC under general (non-Bain) distortions with ATM. 2504.07338 left this open ("susceptible to further cuboidal distortions").

**(B) MGFs: real but modest bridge.**
- GZL could supply double-precision values of planar MGFs at any τ, which is useful for numerical cross-checks of conjectured identities or Laplace equations away from the cusp.
- It could also supply values of generalised MGFs with non-integer or complex exponents (C_{s1,s2,s3} for real s_i > 1). Fedosova–Klinger-Logan-type complex-index objects and Dorigoni–Green "generalised Eisenstein series" sit near these.
- Honest caveat: iterated-Eisenstein-integral representations (2502.05531) already give many digits at any τ, and all weight ≤ 6 and dihedral weight-7 identities are proven (1608.04393). GZL does not unlock a blocked problem there.
- Where it could matter:
  - Higher-loop planar MGFs beyond 4 vertices, where no package exists (2502.05531 covers "up to four vertices").
  - Genus-one objects with a 3D analogue: lattice sums on Λ ⊂ R³, e.g. MGF-like sums for T³ compactifications. Speculative.

**(C) PSLQ / closed-form discovery: weak with GZL as is.**
- GZL/EpsteinLib are double precision (paper: machine precision at best; zeta-algebra errors ~10⁻⁸–10⁻¹³).
- PSLQ needs roughly (basis size × digits per coefficient) digits. At 15 digits, only tiny bases with small integer coefficients are reliable.
- A multiprecision re-implementation of the circle-zeta BZ integral (mpmath/Arb) would be needed for serious PSLQ work on e.g. C_{a,b,c}(i) or 3-body lattice constants.

**(D) Feynman periods: speculative bridge.**
- For a primitive graph, the lattice graph zeta ζ_{Z^d,G}(ν) should have a pole at the "overall" exponent ν|E| = d(|V|−1). A lattice sum of a homogeneous function of degree −(|V|−1)d diverges logarithmically, and its residue is the projective (sphere) integral, i.e. the position-space Feynman period.
- The lattice has no UV subdivergences because K(0)=0. This is analogous to the Epstein zeta residue at ν=d.
- But this regime has per-edge ν < d (φ⁴: ν=2 in d=4). That is outside GZL's domain and needs analytic continuation GZL does not provide. Periods are computed far better by graphical functions/HyperInt. I would not pursue this.

**(E) Random walks with power-law jumps: exact but narrow bridge.**
- The return probability after n steps of a lattice Lévy flight with p(x) ∝ |x|^{−ν} is P_n(0) = ζ^{(n)}_Λ(ν,…,ν) / Z_{Λ,ν}(0)^n. This is exactly a circle zeta.
- The Green function/escape probability is a resolvent ∫_BZ 1/(1−p̂). That needs EpsteinLib plus quadrature, not graph machinery.

**(F) ML: weak.**
- Long-range ML potentials need two-body inverse-power lattice sums (EpsteinLib-level, not graph-level).
- Graph-structured multi-body kernels (ACE/cluster expansions) mostly factorise into stars/trees. On a perfect lattice these reduce to products of lattice Fourier transforms, which is trivial for GZL.
- I found no ML architecture that needs Σ over lattice placements of a non-tree graph with power-law edges.
- Possible niche: GZL as ground truth or feature generator for long-range spin-model datasets in ML-for-quantum-materials. That is the physics application, not a new field.

### Gaps
- Did not find any 2025–26 paper explicitly posing "global lattice optimum for 3-body Riesz/ATM energy among all 2D/3D lattices" as an open problem. The novelty of that question is inferred from absence in my searches (arXiv abstract search, web search). Bétermin's "open questions" list (cited by 2411.17199) was not fetched.
- Did not verify whether any group has published high-precision MGF values at interior τ for 5+ vertex planar graphs.
- The Luo–Wei Sarnak–Strömbergsson proof (2609.17356) is 10 days old. Its acceptance is unknown.
- A. Basu's 2025–26 arXiv output could not be listed (the arXiv author search returned no 2025–26 hits for "Basu, Anirban"). McPhedran's 2025–26 lattice-sum work beyond Euler sums was not found.

## Q3. What would be genuinely new results? (with pilot numerics)

### Takeaway
The most credible "first" is a systematic map of optimal lattices for pure many-body power-law energies in 2D and 3D.

OWN COMPUTATION: pilot numerics already show two things.
- In 2D, the three-body Riesz energy is minimised by the hexagonal lattice for ν ≲ 3.918. There is then a first-order transition to the square lattice. Using the MGF dictionary this is a statement about the minimum of C_{s,s,s}(τ).
- In 3D, BCC beats FCC for the triangle energy at all ν tested. Contrast Luo–Wei: FCC is optimal for two-body Epstein zeta.

Second-tier firsts:
- Double-precision evaluation of higher-vertex planar MGFs, eMGF grids and non-integer-exponent MGFs at arbitrary τ.
- A "many-body Sarnak–Strömbergsson" conjecture backed by GZL scans.

### Cited Findings
Experimental setup: GZL 1.0.0 `zeta_circle`, unit-covolume lattices, A(τ) = [[1, Re τ], [0, Im τ]]/√(Im τ), equal exponents ν on all three edges ([gzl on PyPI](https://pypi.org/project/gzl/)).

**OWN COMPUTATION, 2D three-body Riesz energy E₃(τ;ν) = ζ_{C₃}(ν,ν,ν):**
- Coarse grid over the fundamental domain (x ∈ [0, 0.5], |τ| ≥ 1, y ≤ 2) plus Nelder–Mead from 6 starts.

| ν | global min at | hex value | square value |
|---|---|---|---|
| 2.5 | hex | 20.6364817949 | 21.0272603756 |
| 3.0 | hex | 13.3917941102 | 13.6528937154 |
| 3.5 | hex | 9.7723761021 | 9.8932825030 |
| 3.75 | hex | 8.5685550844 | 8.6172422133 |
| 4.0 | square | 7.6061707526 | 7.5826699011 |
| 5.0 | square | — | 4.8380268467 |
| 6.0 | square | 3.7457712638 | 3.2499743124 |

- Hex = square crossing at ν* ≈ 3.91836 (Brent root).
- Along the boundary arc |τ|=1 at ν* the energy rises from 7.898512 (hex, 60°) to a maximum of ≈7.9096 (~72°) and falls back to 7.898512 (square, 90°). This is a barrier, so the transition is first order, with both lattices local minima.
- Rectangular deformations of the square raise the energy (y = 1.05: 7.922704).
- Brute-force check: a direct double sum over |x|,|y| < 12 on Z² at ν=6 gives 3.2499743122 vs GZL 3.2499743124.
- Mechanism (my estimate): at large ν the energy is dominated by the smallest triangles. Hexagonal has 12 ordered equilateral nearest-neighbour triangles (≈12·(2/√3)^{−3ν/2}). Square's smallest triangles are (1,1,√2) (8·2^{−ν/2}). So repulsive three-body terms penalise equilateral triangles.

**OWN COMPUTATION, 3D triangle energy at unit covolume:**

| ν | FCC | BCC | SC |
|---|---|---|---|
| 3.5 | 38.1993 | 38.0744 | 41.9368 |
| 4.5 | 17.9624 | 17.7406 | 21.0981 |
| 6.0 | 8.0434 | 7.7102 | 10.2962 |
| 9.0 | 2.3133 | 1.9785 | 3.2723 |

- About 28 s per set of three lattices on one core.
- BCC < FCC throughout. This is consistent in direction with the ATM-on-Bain-path finding of [arXiv:2504.07338](https://arxiv.org/abs/2504.07338), and opposite to the two-body Epstein result (FCC optimal for all s > 3/2) of [arXiv:2609.17356](https://arxiv.org/abs/2609.17356).

**OWN COMPUTATION, MGF corollaries:**
- Zagier's identity (C_{1,1,1} = E_3 + ζ(3)) and C_{2,2,1} = (2/5)E_5 + ζ(5)/30 were reproduced numerically (Q1).
- E_s(τ) is minimised at the hexagonal point for all s>1 (classical Rankin/Cassels/Ennola/Diananda/Montgomery). This is from background knowledge, not fetched in this session.
- Therefore the 2D three-body energies with exponents (2,2,2) and (4,4,2) are exactly minimised by the hexagonal lattice. The (2,2,2) case sits at the borderline ν=d, where the sum still converges.
- Attribution of the C_{2,2,1} identity to D'Hoker–Green–Vanhove, arXiv:1502.06698 (unverified). The identity itself is confirmed numerically to 1.2×10⁻⁸.

### Inferences
**Candidate result 1: "Optimal lattices for many-body power-law energies".** A systematic study would cover:
- 2D and 3D; triangle, 4-cycle and K₄ (4-body) energies.
- Pure and mixed with two-body LJ/Riesz; ATM decomposed into triangle zetas.
- Global optimisation over the 2- or 5-dimensional lattice moduli space.
- Phase diagrams in ν and in coupling ratio.

Why it fits:
- GZL is the only tool I know that evaluates these sums to near machine precision in seconds, for ν arbitrarily close to d.
- The 2D hex→square first-order transition at ν* ≈ 3.918, and 3D BCC < FCC, are concrete hooks.
- Proving the 2D statement for integer s could use the MGF identity/Laplace-equation machinery. For s=1 and (2,2,1) exponents it reduces to Epstein minimisation, which is already proven.

This is the clearest bridge from GZL into math.MG / math-ph / chemistry. Novelty is unverified: related work exists in Bétermin's group and the Buchheit/Schwerdtfeger group, and the latter is GZL's own circle.

**Candidate result 2: "C_{s,s,s}(τ) and generalised dihedral MGFs as functions of real s".** The MGF dictionary gives a physics-free reformulation: the minimum of the two-loop MGF C_{s,s,s} over the fundamental domain moves from the hexagonal to the square point at s* ≈ 1.959. That is a statement about modular functions which, as far as I found, has not been posed. Worth checking against Laplace-equation structure: for integer s, C_{s,s,s} satisfies inhomogeneous Laplace equations with Eisenstein-series sources.

**Candidate result 3 (weaker).** A double-precision atlas of planar MGFs with ≥5 vertices, and of eMGFs on a full z-grid, at interior τ. This would provide cross-checks for iterated-integral packages (2502.05531) and for eMGF constructions (2511.15883). Value to that community is moderate, because they have exact methods.

**Candidate result 4 (needs extra engineering).** A multiprecision circle-zeta/graph-zeta evaluator (mpmath/Arb port of the BZ-integral + Epstein representation) to run PSLQ on:
- lattice-sum constants such as C_{a,b,c}(i), C_{a,b,c}(e^{iπ/3}) and 3-body Madelung-like constants of FCC/BCC;
- values against bases of products of ζ(s)β(s), L-values of CM forms and Kronecker limit constants.

Not possible with GZL at double precision.

**Not recommended:** Feynman periods (the domain is ν<d, and better tools exist), ML (only two-body sums are needed there, i.e. EpsteinLib), and 1D Tornheim/MZV sums (already exact).

### Gaps
- The 2D scan covered only the fundamental domain with y ≤ 2–3 and six local starts. Strongly rectangular or rhombic minima at larger y are unlikely, given energies rise with y, but were not excluded. The 3D study compared only FCC/BCC/SC; no global 5-parameter search was done.
- Precision of the 2D/3D values: the brute-force check agrees to ~6×10⁻¹¹ (truncation-limited). GZL's internal error was not separately quantified for 3D.
- s* ≈ 1.959 is a numerical root, not verified against any analytic result.
- Whether Bétermin et al. or Luo–Wei have studied three-body lattice energies in 2025–26 was not confirmed. The web search found none, but the conference and journal literature was not searched.

## Q4. How robust are the bridges? (limits of GZL that affect every transfer)

### Takeaway
Four constraints decide which transfers are real:
- Summability ν > d per edge. Borderline ν = d is reachable only by extrapolation to ~10⁻⁶–10⁻⁸.
- Isotropic |x|^{−ν} kernels only. Anisotropic Epstein functions now exist in 2609.28282 but are not yet in GZL graphs.
- d ≤ 3 Bravais lattices; multi-atomic lattices are future work.
- Double precision.

Within these, many-body lattice energies and planar MGFs with a_e ≥ 2 (or non-integer a_e > 1) are fully in scope. Madelung/Coulomb sums (ν<d), Feynman periods (ν<d), non-planar MGFs, modular graph forms with a≠b, and high-digit PSLQ are out of scope or weak.

### Cited Findings
- Summability requirement min Re ν > d. Integral representation of circle zetas valid for Re ν ≥ d, with meromorphic continuation via Hadamard integral. Future work: anisotropic zeta functions and multi-atomic lattices. — [arXiv:2609.18918](https://arxiv.org/abs/2609.18918)
- `zeta_circle` raises UnsupportedLatticeSumError for exponents ≤ d; d ∈ {1, 2, 3}. — [gzl on PyPI](https://pypi.org/project/gzl/) (installed docstring)
- Anisotropic Epstein zeta functions from k-derivatives, including removal of Rayleigh–Wood singularities. — [arXiv:2609.28282](https://arxiv.org/abs/2609.28282)
- Epstein zeta's meromorphic continuation and ν<d handling (needed for Coulomb/Madelung and ATM terms with small exponents) are available at the single-edge and circle level. — [arXiv:2412.16317](https://arxiv.org/abs/2412.16317); [arXiv:2504.07338](https://arxiv.org/abs/2504.07338)
- Conditionally convergent Madelung sums are an active, separate numerical topic (shape-dependent boundary terms). — [arXiv:2608.10041](https://arxiv.org/abs/2608.10041)

### Inferences
- Madelung-type sums (Coulomb, ν=1, d=3) belong to EpsteinLib, not GZL.
- MGFs with unit-weight edges (most string amplitudes) are borderline. Clean support needs GZL's analytic continuation extended from circles to general graphs. The Hadamard-integral technique already exists for circles.
- Adding anisotropic kernels would unlock modular graph *forms* (holomorphic/antiholomorphic exponents) and dipolar many-body lattice energies.
- Adding multiprecision would unlock PSLQ.

### Gaps
- No documentation was found on whether GZL exposes the circle-zeta Hadamard continuation (ν < d) through its Python API. The public `zeta_circle` rejects ν ≤ d.
- EpsteinLib precision beyond double was not checked. I assume double; the C library signature returns `complex`.
