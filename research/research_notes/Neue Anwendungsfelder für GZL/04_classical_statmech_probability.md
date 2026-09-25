# GZL transfer opportunities in classical statistical physics, probability and spreading processes with power-law kernels on lattices (focus 2025–2026)

Scope note: every arXiv ID below was checked against its arxiv.org abstract page or the arXiv search listing on 2026-09-25. The two GZL papers exist: arXiv:2609.18918 (Buchheit & Rupp, method) and arXiv:2609.18761 (Duft, Adelhardt, Koziol, Buchheit & Schmidt, application). Anything I could not verify is marked "UNVERIFIED" or listed under Gaps. arXiv searches (full-text search sorted by date) that returned nothing relevant from 2025–2026 are reported as gaps, because a null result is informative too.

Facts about GZL that matter for the assessment (taken from the papers themselves):
- In GZL a graph lattice sum pins one vertex and sums the others over the lattice. Power-law kernels use the convention K_ν(x)=|x|^(−ν) for x≠0 and K_ν(0)=0 (Def. 2.7). The method is limited to absolutely summable kernels K_e ∈ ℓ¹(Λ), which means ν>d (Lemma 2.5). Momentum enters as a phase on a source–sink displacement. Blocks with tw≤2 use Epstein-zeta algebra; blocks with tw>2 use tensor-network bucket elimination; one FFT gives the whole momentum grid. — [arXiv:2609.18918 HTML](https://arxiv.org/html/2609.18918), [abs](https://arxiv.org/abs/2609.18918)
- In the quantum linked-cluster application, the embedding sums are hard-core ("hardcore constraint … explicitly prohibits vertices from occupying the same lattice site"). They are mapped to unconstrained "soft-core" sums by summing over valid contractions, where two non-adjacent vertices are merged. The result is deterministic coefficients "within minutes on standard desktop hardware", replacing cluster-scale Monte Carlo integration. — [arXiv:2609.18761 HTML](https://arxiv.org/html/2609.18761), [abs](https://arxiv.org/abs/2609.18761)
- Reported scale: order 13 with 8,403 graphs and order 11 with 22,677 graphs. A 3D dispersion on 4096 momentum points takes about 10 minutes on one core. — [arXiv:2609.18918 HTML](https://arxiv.org/html/2609.18918)

## Q1: Which quantities in recent (2025–2026) classical/probabilistic work are lattice sums of graph-structured products of power laws? (formulas)

### Takeaway
The cleanest matches are the **coefficients of high-temperature (small-β) expansions of classical long-range spin, Potts/percolation and self-avoiding-walk models**. At every order these are sums over lattice placements of graph vertices of products of J(x)=|x|^{−(d+σ)} (and powers of it). They are the same objects GZL already evaluates for quantum linked-cluster expansions. Nobody in the 2025–26 literature I found computes them. That literature is dominated by rigorous non-perturbative proofs (Hutchcroft, Bäumler, Liu/Slade) and large-scale Monte Carlo (Deng group). Lace-expansion diagrams (bubble, triangle) and k-point tree formulas are graph lattice sums, but their edge kernels are Green functions ~|x|^{−(d−α)}. Those are *not* summable power laws, so they fall outside GZL's ℓ¹ domain in real space.

### Cited Findings

**Long-range percolation (LRP), rigorous side**
- Model: points x,y∈Z^d are joined with probability 1−exp(−β‖x−y‖^{−d−α}). Eight conjectured critical regimes exist, with the mean-field/low-dimension transition at d=min{6,3α} and a long-range/short-range crossover α_c(d). Hutchcroft's three-paper series (Aug 2025) develops a non-perturbative real-space RG. — [Hutchcroft, Critical LRP I, arXiv:2508.18807](https://arxiv.org/abs/2508.18807)
- In the LR low-dimensional regime (d/3<α<α_c(d)) the results are η=2−α, γ=(2−η)ν, δ=(d+α)/(d−α) and d_f=(d+α)/2. The paper says computing α_c(d) for 2<d<6 "appears to be beyond the scope of current techniques". — [Critical LRP II, arXiv:2508.18808](https://arxiv.org/abs/2508.18808)
- At d=3α<6 the paper gives P_βc(|K|≥n) ~ C (log n)^{1/4}/√n and a critical three-point function ≍ sqrt(∏ ‖x_i−x_j‖^{−d+α} / log(1+min‖·‖)). The constant C is not computed. — [Critical LRP III, arXiv:2508.18809](https://arxiv.org/abs/2508.18809)
- A survey of the dimension dependence and the eight regimes. — [Hutchcroft, arXiv:2510.03951](https://arxiv.org/abs/2510.03951)
- The triangle condition (finiteness of ∇=Σ_{x,y} τ(0,x)τ(x,y)τ(y,0), where τ is the two-point function) is proved non-perturbatively in d=1,2,3 for 0<α<d/3, using the critical two-point bound τ(x) ≍ ‖x‖^{−d+α} for α<1. — [Hutchcroft, arXiv:2404.07276](https://arxiv.org/abs/2404.07276) (2024 anchor)
- In the spread-out, lace-expansion setting, long-range percolation, Ising and SAW with couplings ~|x|^{−(d+α)} above the upper critical dimension have two-point functions that coincide exactly with a random-walk two-point function up to a constant. The lace-expansion convergence proof is given "assuming diagrammatic estimates". — [Y. Liu, arXiv:2502.12104](https://arxiv.org/abs/2502.12104)
- The high-dimensional critical k-point function is T(x_1..x_k) ~ V^{k−2}A^{2k−3} Σ_{trees T} Σ_{Φ: V(T)→Z^d, Φ(i)=x_i} ∏_{edges} G(Φ(u),Φ(v)). This is literally a lattice sum over placements of the internal vertices of cubic trees, with G the lattice Green function. The constants A and V are not explicit. — [Blanc-Renaudie & Hutchcroft, arXiv:2607.22387](https://arxiv.org/abs/2607.22387) (nearest-neighbour/high-d, not long-range)
- Other 2025–26 LRP papers are purely qualitative or rigorous and contain no numerically computed lattice sums: strict monotonicity of critical points [arXiv:2510.26314](https://arxiv.org/abs/2510.26314); a discontinuous transition on the hierarchical lattice [arXiv:2512.12124](https://arxiv.org/abs/2512.12124); the truncated one-arm exponent [arXiv:2601.07808](https://arxiv.org/abs/2601.07808); recurrence of inhomogeneous LRP [arXiv:2608.25201](https://arxiv.org/abs/2608.25201); the LR contact-process truncation property [arXiv:2604.15753](https://arxiv.org/abs/2604.15753).

**Long-range models, numerical/field-theory side (Deng/Fan/Xiao/Chen group, very active 2025–26)**
- 2D bond percolation with p_ij = ρ·C(σ,L)/r_ij^{2+σ}, normalized by C(σ,L)Σ_{j≠i} r_ij^{−(2+σ)} = 4 (minimum-image distance, periodic boundaries). Event-based Monte Carlo reaches L=16384. Reported ρ_c(σ=1)=0.307591(4) and ρ_c(σ=1/2)=0.26263(2). η deviates from 2−σ near σ≈3/2. — [Liu, Xiao, Fan, Deng, arXiv:2608.20750](https://arxiv.org/abs/2608.20750); details from the [HTML version](https://arxiv.org/html/2608.20750)
- 2D LR Ising Monte Carlo up to L=8192 (FK critical polynomial, Binder ratio) puts the LR→SR crossover at σ*=2, against Sak's criterion. — [Xiao, Liu, Fan, Deng, arXiv:2512.04805](https://arxiv.org/abs/2512.04805)
- Proposed (d,σ) universality diagrams for percolation, O(n) and FK-Ising, motivated by Lévy flights. Lévy flights are used only heuristically; no Green functions or bubble diagrams are computed. — [Xiao, Fan, Deng, arXiv:2512.02948](https://arxiv.org/abs/2512.02948); [HTML](https://arxiv.org/html/2512.02948)
- 1D long-range SAW with jump weight r^{−(d+σ)}: Monte Carlo estimates of the critical fugacity z_c, the Binder ratio, ν and η. The LR/SR boundary is at σ=1. — [Zhou et al., arXiv:2609.23346](https://arxiv.org/abs/2609.23346)
- Long-range loop-erased random walks with P(r)~|r|^{−(d+σ)} in d=1–5, studied by Monte Carlo. — [Xiao, Pan, Fan, Deng, arXiv:2603.27992](https://arxiv.org/abs/2603.27992)
- Field theory: a two-loop 4−ε expansion for LR O(n) [arXiv:2602.07818](https://arxiv.org/abs/2602.07818); a one-loop 6−ε expansion for LR percolation and Lee–Yang [arXiv:2608.15120](https://arxiv.org/abs/2608.15120); a two-loop analysis of LR quantum O(n) [arXiv:2606.22407](https://arxiv.org/abs/2606.22407). These are continuum momentum integrals and involve no lattice sums.
- 1D LR Ising: FRG in the local-potential approximation compared with Dyson's hierarchical model; ν is benchmarked against Monte Carlo and expansions. — [Pagni, Giachetti, Trombettoni, Defenu, arXiv:2510.02458](https://arxiv.org/abs/2510.02458)

**Cluster expansions of LR Ising (rigorous)**
- A convergent low-temperature cluster expansion of the 1D LR Ising model with J(r)=r^{−α}, α∈(1,2]. — [Bissacot & Corsini, arXiv:2602.12447](https://arxiv.org/abs/2602.12447)
- Multidimensional LR Ising with J/|x−y|^α: a cluster expansion over multiscale contours. — [Affonso et al., arXiv:2508.15666](https://arxiv.org/abs/2508.15666)
- These are existence and decay results; they contain no numerical coefficients.

**Long-range spin glasses**
- A zero-temperature spin glass in a field on a 1D long-range model, studied with a loop expansion in the Bethe M-layer formalism. — [Angelini, Palazzi, Parisi, Rizzo, arXiv:2602.03488](https://arxiv.org/abs/2602.03488)
- A phase transition on the Nishimori line for the 1D LR spin glass with J(r)~r^{−α}, 1<α<3/2. The case 3/2≤α≤2 remains open. — [Okuyama & Ohzeki, arXiv:2604.07130](https://arxiv.org/abs/2604.07130)

**Long-range random walks / Lévy flights on lattices**
- Pólya-type recurrence for "fractional random walks" generated by L^{α/2}: transient for d>α. The paper gives closed forms for the 1D Green matrix and ever-passage probabilities with Riesz-potential decay. — [Michelitsch et al., arXiv:1707.05843](https://arxiv.org/abs/1707.05843) (2017 anchor)
- Effective-medium theory for lattice random walks with power-law rates. — [Thiel & Sokolov, arXiv:1604.06621](https://arxiv.org/abs/1604.06621) (anchor)
- A closed form for Pólya's return probability (nearest neighbour, d≥3) via a Lauricella F_C. — [Gaunt, Nadarajah, Pogány, arXiv:2311.11326](https://arxiv.org/abs/2311.11326)
- 2025–26 LR random-walk papers are rigorous (homogenization with critical jump index |x−y|^{−d−2}; circulant LR walks on Z_q^d) and do not compute lattice constants. — [arXiv:2604.21162](https://arxiv.org/abs/2604.21162), [arXiv:2510.22554](https://arxiv.org/abs/2510.22554)

**Long-range voter / coarsening / discrete Gaussian chain**
- Voter imitation probability P(r) ∝ r^{−α}, and the LR Ising/voter coarsening correspondence. — [arXiv:2603.14165](https://arxiv.org/abs/2603.14165), [arXiv:2510.15410](https://arxiv.org/abs/2510.15410)
- Occupation-time fluctuations of the LR voter model use a local CLT for the LR random walk. — [arXiv:2509.17518](https://arxiv.org/abs/2509.17518)
- The disordered LR discrete Gaussian chain is rigorous. — [arXiv:2609.18786](https://arxiv.org/abs/2609.18786)

**Lattice energies (Riesz/Epstein), noted briefly**
- The 3D Sarnak–Strömbergsson conjecture is proved: FCC minimizes E(L,s)=Σ|v|^{−2s} among unit-covolume lattices for s>3/2. — [Luo & Wei, arXiv:2609.17356](https://arxiv.org/abs/2609.17356)
- Epstein-zeta toolchain: EpsteinLib [arXiv:2412.16317](https://arxiv.org/abs/2412.16317); n-body (Axilrod–Teller–Muto) lattice sums, "weeks to minutes" [arXiv:2504.11989](https://arxiv.org/abs/2504.11989); anisotropic Epstein zeta and the lattice Euler–Maclaurin correction [arXiv:2609.28282](https://arxiv.org/abs/2609.28282).

**Spreading processes (epidemics/dispersal)**
- Power-law distance-decaying infection between locations appears in metapopulation SIR models, for example [arXiv:2007.08002](https://arxiv.org/abs/2007.08002) and [arXiv:2110.13956](https://arxiv.org/abs/2110.13956) (not 2025–26). Heavy-tailed dispersal kernels from stopped subdiffusive fBM: [arXiv:2606.21681](https://arxiv.org/html/2606.21681) (seen via web search; abstract not independently re-checked).

### Inferences
- **High-temperature (small-β) linked-cluster expansion of classical LR models — strong bridge.** In the standard free-embedding linked-cluster expansion for classical spins (vertex cumulants times free lattice sums), the coefficient of β^n in χ(k) or ln Z is Σ_graphs (combinatorial/vertex factor) × Σ_{x_2..x_m ∈ Λ} ∏_{edges e} J(x_{i_e}−x_{j_e})^{m_e} e^{ik·(x_s−x_t)}, with J(x)=|x|^{−(d+σ)}, J(0)=0 and m_e the edge multiplicity (J^{m} is again a power law with ν=m(d+σ)>d). This is exactly a graph zeta function with GZL's K(0)=0 convention. The vertex factors (Ising, O(n), φ⁴, Potts cumulants) are model-specific and user-supplied. (Inference; the classical free-graph linked-cluster expansion is textbook, Wortis/Englert/Lüscher–Weisz. I did not verify a 2025–26 source for it.)
- **LR percolation / q→1 Potts series in β.** Here 1−exp(−β J(x)) = Σ_m (−1)^{m+1} β^m J(x)^m/m! (formula from [2508.18807](https://arxiv.org/abs/2508.18807)). Every β-order is a finite sum of power laws, which fits GZL's K=a+Σ b_j|x|^{−ν_j} per edge. The FK/subgraph expansion needs distinct sites (hard-core). That is handled by the contraction mapping already used in [2609.18761](https://arxiv.org/html/2609.18761), at combinatorial cost.
- **LR SAW generating function** c_n = Σ_{distinct x_0..x_n} ∏ D(x_i−x_{i−1}) is a hard-core path-graph sum. Inclusion–exclusion over coincidences produces contracted graphs (cycles and multigraphs), and their number grows quickly (Bell-number-like). Only moderate orders (maybe n≈12–16) look plausible. This is weaker than the spin case.
- **1D LR spin-glass high-temperature series** (quenched average): the Edwards–Anderson susceptibility series involves lattice sums of products of J_ij² ∝ r^{−2σ} over even graphs. This is also a graph zeta function. It is plausible but speculative; I found no recent series work.
- **Lace-expansion bubble/triangle and tree k-point sums:** graph-shaped, but the edge kernels are Green functions. At criticality G(x) ≍ |x|^{−d+α}, which is not in ℓ¹ (ν=d−α<d), so it violates GZL's summability requirement. Off criticality G is not a finite power-law sum. The workable GZL route is indirect (see Q3).

### Gaps
- No 2025–26 paper computing high-temperature, low-temperature or virial series coefficients for classical LR Ising/O(n)/Potts/percolation on lattices turned up. arXiv search for "high-temperature series long-range" and "linked-cluster classical power-law" returned only quantum works ([arXiv:2408.13145](https://arxiv.org/abs/2408.13145), [2609.18761](https://arxiv.org/abs/2609.18761)). A web search surfaced a 1990 J. Phys. A paper "The Ising model with long-range ferromagnetic interactions" ([link](http://fulltext.calis.edu.cn/iop/0305-4470/23/11/036/jav23i11p2157.pdf)). Its authors, method and whether it contains series are UNVERIFIED. Older LR-Ising series works (e.g. Nagle–Bonner 1970, Glumac–Uzelac ~1989) are recalled from memory only and are UNVERIFIED.
- No 2025–26 epidemic/ecology/neural-field paper on arXiv with lattice power-law kernel sums was found (searches for "dispersal kernel fat-tailed lattice", "long-distance dispersal power-law", "neural field power-law", "Levy epidemic lattice" returned nothing relevant).

## Q2: How are these quantities computed today, and what precision/size limits are reported?

### Takeaway
Today these quantities are handled in three ways. (a) Rigorous work proves bounds or asymptotics up to constants and needs no numbers, except the computer-assisted nearest-neighbour lace expansion. (b) Monte Carlo on finite tori, 2D up to L=8192–16384, gives critical points to about 5–6 digits. (c) Continuum ε-expansions to two loops. Graph lattice sums with power laws on Bravais lattices are not being evaluated deterministically anywhere in this community. In the quantum community the same sums were done by Monte Carlo integration until GZL.

### Cited Findings
- Rigorous LRP results are "up-to-constants" estimates: pointwise τ(x) ≍ ‖x‖^{−d+α} [arXiv:2404.07276](https://arxiv.org/abs/2404.07276), and the LR-LD k-point function "up-to-constants" [arXiv:2508.18808](https://arxiv.org/abs/2508.18808). The amplitudes (C in the (log n)^{1/4}/√n tail; A and V in the tree formula) are not computed ([2508.18809](https://arxiv.org/abs/2508.18809), [2607.22387](https://arxiv.org/abs/2607.22387)).
- The lace-expansion convergence argument for long-range models is conditional on "diagrammatic estimates" and works in the spread-out setting (large range L), not the plain L=1 model. — [arXiv:2502.12104](https://arxiv.org/abs/2502.12104)
- Precedent for numerics in the lace expansion: the nearest-neighbour percolation mean-field proof for d≥11 is computer-assisted. It uses "rigorous numerical upper bounds on various simple random walk integrals" (Hara–Slade 1992) together with Mathematica notebooks, and gives p_c bounds "at most 1.306% off in d=11". — [Fitzner & van der Hofstad, arXiv:1506.07977](https://arxiv.org/abs/1506.07977); related [arXiv:1506.07969](https://arxiv.org/abs/1506.07969), [arXiv:1905.02785](https://arxiv.org/abs/1905.02785). No long-range analogue was found.
- Monte Carlo scale: 2D LRP up to L=16384 with ρ_c to about 6 significant digits, e.g. 0.307591(4) [arXiv:2608.20750](https://arxiv.org/html/2608.20750); 2D LR Ising up to L=8192 [arXiv:2512.04805](https://arxiv.org/abs/2512.04805); 1D LR SAW "high-precision" z_c by irreversible Monte Carlo [arXiv:2609.23346](https://arxiv.org/abs/2609.23346). These simulations truncate interactions via the torus (minimum image), and the kernel normalization is a finite lattice sum ([2608.20750 HTML](https://arxiv.org/html/2608.20750)).
- The Luijten–Blöte cluster algorithm samples from the cumulative bond-probability distribution, costing O(N log N) for integrable power laws. Periodic images or Ewald-type sums are used to reduce finite-size effects. — [arXiv:1611.05659](https://arxiv.org/abs/1611.05659) (as summarized by web search; anchor)
- In the quantum analogue, before GZL the long-range linked-cluster lattice sums were done by classical Monte Carlo integration inside the white-graph embedding scheme. — [Adelhardt, Koziol, Langheld, Schmidt review, arXiv:2403.00421](https://arxiv.org/abs/2403.00421). GZL replaces this with deterministic evaluation in minutes. — [arXiv:2609.18761](https://arxiv.org/abs/2609.18761)
- The ε-expansions for LR O(n), percolation and quantum O(n) are two-loop or one-loop continuum calculations. — [2602.07818](https://arxiv.org/abs/2602.07818), [2608.15120](https://arxiv.org/abs/2608.15120), [2606.22407](https://arxiv.org/abs/2606.22407)
- Lévy-flight usage in the universality-diagram papers is heuristic; no Green functions or bubble diagrams are computed. — [arXiv:2512.02948 HTML](https://arxiv.org/html/2512.02948)

### Inferences
- **The main competitor to a GZL-based series is large-scale Monte Carlo.** Monte Carlo is strong on universal exponents and weak near marginal points with logarithmic corrections (σ→2 crossover, d=3α). Series would give non-universal quantities (β_c, z_c, ρ_c, amplitudes) deterministically and the whole σ-range in one scan. They would give exponents only after series analysis (ratio or Padé methods). Near crossovers, series analysis is also hard, so this is not a guaranteed win.
- Nobody computes lattice-exact D̂(k)=Σ_x D(x)e^{ik·x} for pure power-law steps as an input to diagrammatic numerics, because no long-range computer-assisted lace expansion exists yet.

### Gaps
- I did not find the exact vertex-exclusion conventions or order limits of any classical LR series expansion, because none was found.
- Reported Monte Carlo precision for 3D LR Ising critical temperatures (e.g. Luijten–Blöte 2002) was not re-checked this session.

## Q3: Would GZL apply directly, and what is missing?

### Takeaway
GZL applies **directly** to high-temperature and small-β series of classical LR models on Bravais lattices in d=1,2,3. Those series need isotropic |x|^{−(d+σ)} kernels with σ>0 (so ν>d, which is ℓ¹), K(0)=0, per-edge powers J^m, momentum-resolved χ(k), and free or hard-core embeddings via contractions. GZL applies **indirectly** to lace expansion and Lévy-flight lattice Green functions: it supplies D̂(k) on a full Brillouin-zone grid (a one-edge graph with momentum), which then feeds a singular k-integral that GZL does not do. It does **not** apply to hierarchical lattices, non-Bravais or random environments, anisotropic dipolar kernels (only via the separate anisotropic Epstein work), critical Green-function kernels (not in ℓ¹), or computer-assisted proofs that need rigorous interval bounds.

### Cited Findings
- GZL supports K_ν(0)=0, the ℓ¹ restriction (ν>d), one pinned vertex, a momentum phase on the source–sink displacement, and FFT for the full grid. — [arXiv:2609.18918](https://arxiv.org/html/2609.18918)
- Hard-core sums are converted to soft-core sums via contraction of non-adjacent vertices. — [arXiv:2609.18761](https://arxiv.org/html/2609.18761)
- Epstein zeta functions with a phase are evaluated to full precision for arbitrary real ν, with singularity removal, which "facilitates the evaluation of integrals". — [EpsteinLib, arXiv:2412.16317](https://arxiv.org/abs/2412.16317)
- Anisotropic Epstein zeta functions (dipolar-type kernels) and removal of Rayleigh–Wood singularities are covered by separate new work, not by GZL's stated scope. — [arXiv:2609.28282](https://arxiv.org/abs/2609.28282)
- Many recent LRP results live on the hierarchical lattice, which is not Bravais. — [arXiv:2512.12124](https://arxiv.org/abs/2512.12124), [arXiv:2509.09589](https://arxiv.org/abs/2509.09589)
- Spread-out models (the lace-expansion setting) use D(x) with range L≫1 and power-law tails. — [arXiv:2502.12104](https://arxiv.org/abs/2502.12104)

### Inferences
Assessment by area:

| Area | Quantity | GZL fit | Verdict |
|---|---|---|---|
| HT/β-series, LR Ising/O(n)/φ⁴ (d=1,2,3) | Σ_graphs w_G · Σ_x ∏_e J(x_e)^{m_e} e^{ik·Δx} | Direct: free embeddings, K(0)=0, per-edge powers, χ(k) grid | **Strong** |
| LR percolation / Potts q→1, β-series | subgraph expansion with v_e=e^{βJ_e}−1 or p_e expansions | Direct after hard-core→soft-core contractions | **Moderate–strong** (combinatorics heavy) |
| LR SAW / LR LERW series | hard-core path sums | Contractions blow up | **Moderate–weak** |
| 1D/2D LR spin glass HT series | even graphs with J² ∝ r^{−2σ} | Direct (disorder-averaged kernels are deterministic power laws) | **Moderate**, speculative |
| Lace-expansion bubble/triangle, long-range L=1 | ∫ D̂(k)^a/(1−D̂(k))^b dk | Indirect: GZL/EpsteinLib gives D̂ on the grid; the singular BZ integral is extra; floating-point is not rigorous | **Moderate** as numerics, **weak** for proofs |
| k-point tree sums (2607.22387) | Σ ∏ G over cubic trees | G not a power law or ℓ¹ at criticality | **Weak** |
| Lévy-flight return probability / Pólya constants / escape probability | 1−1/G(0), G(0)=∫dk/(1−D̂(k)) | Indirect as above; known singular part of D̂ ~ |k|^σ | **Moderate**, low impact |
| Discrete fractional Laplacian with power-law kernel | symbol Σ_x(1−cos k·x)|x|^{−d−2s} | One-edge graph: EpsteinLib suffices, GZL overkill | **Weak** (for GZL specifically) |
| Riesz/Epstein lattice energies | E(L,s)=Σ|v|^{−2s} | One-edge: EpsteinLib; multi-body: Buchheit–Busse 2504.11989 | **Weak / already covered** |
| Epidemic/dispersal moment closures on grids | e.g. Σ_{x,y}K(x)K(y)K(x−y) (triangle "clustering" term) | tw=2, trivial for GZL; the problem is small-scale and usually continuum | **Weak** |
| Hierarchical-lattice LRP, random environments | – | Not Bravais / not translation invariant | **Does not apply** |

- Momentum resolution matters specifically for long-range models. For σ<2 the second moment Σ|x|²G(x) diverges, so the standard second-moment correlation length from HT series is not defined. One instead needs the nonanalytic |k|^σ coefficient of χ(k)^{−1}. GZL's full-grid χ(k) plus EpsteinLib's analytically known singular part of power-law lattice sums (Epstein zeta with phase ~|k|^{ν−d}) make that extraction natural. (Inference.)
- "Precise even for ν close to d" matters because σ→0 (J ∝ r^{−d−σ} barely summable) is exactly where truncated direct sums fail. (Inference, based on the GZL/EpsteinLib claims.)
- Missing pieces for classical use: (i) graph generation and vertex-weight (cumulant) libraries for classical models. Users must bring these, and graph counts at order 15–20 could be about 10^5–10^7 (order-of-magnitude guess, UNVERIFIED). (ii) Series-analysis tooling. (iii) Rigorous error enclosures if proofs are the goal. (iv) Lattices with a basis (e.g. honeycomb, kagome, pyrochlore) are excluded.

### Gaps
- I could not confirm whether GZL's public API exposes the hard-core→soft-core contraction step or whether it lives in the application code of 2609.18761.
- The typical treewidth of high-order classical linked-cluster graphs was not checked (1PI graphs at order ~15–20 may reach tw 4–6, which raises tensor-network cost).

## Q4: What would be genuinely new results, and how honest is each bridge?

### Takeaway
The most credible new results are deterministic high-order (β^10–β^16 or more) high-temperature series for classical long-range Ising/O(n)/Potts and long-range percolation on Z^d (d=1,2,3) and other Bravais lattices. They would give critical couplings β_c(σ) and ρ_c(σ) as smooth functions of σ, benchmarked against the Deng group's Monte Carlo (ρ_c=0.307591(4) at σ=1 in 2D). They would also give momentum-resolved χ(k) with the |k|^σ amplitude. A second, more modest line is lattice-exact Lévy-flight and random-walk diagram constants (return probabilities, bubble and triangle values of the long-range random walk). These could feed a future computer-assisted long-range lace expansion, but they are not proofs themselves. "High-precision triangle diagram for LR percolation at s near d" is **not** a well-posed direct GZL computation: the critical percolation triangle involves the unknown τ, and GZL cannot use τ as a kernel.

### Cited Findings
- Benchmarks available for series checks: ρ_c(σ=1)=0.307591(4) and ρ_c(σ=1/2)=0.26263(2) for 2D LRP [arXiv:2608.20750](https://arxiv.org/html/2608.20750); the 2D LR Ising critical points at L≤8192 [arXiv:2512.04805](https://arxiv.org/abs/2512.04805); 1D LR SAW z_c(σ) [arXiv:2609.23346](https://arxiv.org/abs/2609.23346).
- Open questions where independent deterministic data would be welcome: the location of the LR/SR crossover (σ*=2 versus Sak's 2−η_SR) [arXiv:2512.04805](https://arxiv.org/abs/2512.04805), [arXiv:2602.07818](https://arxiv.org/abs/2602.07818), [arXiv:2608.15120](https://arxiv.org/abs/2608.15120); the non-Gaussian variation of ν(σ) in 2D LRP [arXiv:2608.20750](https://arxiv.org/abs/2608.20750); α_c(d) for 2<d<6, which is "beyond the scope of current techniques" [arXiv:2508.18808](https://arxiv.org/abs/2508.18808); the 1D LR spin-glass transition for 3/2≤α≤2 [arXiv:2604.07130](https://arxiv.org/abs/2604.07130).
- The quantum precedent shows the pipeline works at 10^4 graphs and orders 11–13 in minutes per coefficient set. — [arXiv:2609.18918](https://arxiv.org/html/2609.18918), [arXiv:2609.18761](https://arxiv.org/abs/2609.18761)
- The computer-assisted lace-expansion precedent (nearest neighbour, d≥11) relies on numerically bounded random-walk integrals. — [arXiv:1506.07977](https://arxiv.org/abs/1506.07977)

### Inferences
Ranked opportunities, most to least credible:

1. **High-temperature series for classical LR O(n)/Ising (d=1,2,3). The bridge holds.** It is the classical twin of the quantum linked-cluster expansion GZL was built for, and the lattice-embedding factor is exactly a graph zeta function (free sums, J(0)=0). New outputs: β_c(σ) on dense σ-grids, including σ→0 (ν→d, where direct summation fails); χ(k) over the whole Brillouin zone; the |k|^σ amplitude; and ratio or Padé estimates of γ(σ) for comparison with ε-expansions (2602.07818) and FRG (2510.02458). The 3D LR Ising and 2D LR XY/Heisenberg models, where Monte Carlo is expensive, are natural targets. Risk: series order may be limited by graph counts, and near σ≈2 logarithmic corrections make series analysis unreliable.
2. **Small-β series for LR percolation (1−exp(−βJ) expanded) and q-state LR Potts. The bridge holds with extra combinatorics** (hard-core contractions). New: series estimates of β_c(α,d) and ρ_c(σ), checked against 2608.20750, plus mean-cluster-size amplitudes. These are non-universal numbers that the rigorous programme (Hutchcroft) never provides.
3. **1D/2D LR spin-glass high-temperature series (Edwards–Anderson susceptibility). Plausible but speculative.** It could address the open 3/2≤α≤2 window (2604.07130) numerically. The disorder average makes edge kernels J² ∝ r^{−2σ}, which is still isotropic and summable for 2σ>d.
4. **Lévy-flight lattice constants: return and escape probabilities, and bubble/triangle values of the LR random walk with D(x) ∝ |x|^{−d−α} on Z^d or other Bravais lattices. The bridge is weak-to-moderate.** GZL (or just EpsteinLib) gives D̂(k) exactly; the singular BZ integral needs separate treatment. The value lies in inputs for a hypothetical plain (non-spread-out) long-range lace expansion. Honest caveat: proofs need rigorous interval bounds that GZL does not provide, and GZL's graph machinery is barely used because these are one-edge objects.
5. **Finite-size-scaling plateau constants for long-range models above d_c** (Liu–Park–Slade theory, [arXiv:2412.08814](https://arxiv.org/abs/2412.08814); RG for LR |φ|⁴ in d≥4, [arXiv:2511.03495](https://arxiv.org/abs/2511.03495)). This needs torus sums Σ_{k≠0} 1/(1−D̂(k)) at k∈(2π/L)Z^d. GZL's grid gives the periodized D̂ exactly. It is a one-edge object again, so modest; it depends on whether the torus model uses periodized versus minimum-image kernels (the Deng group uses minimum image).
6. **Spreading processes (SIR, contact process, voter, dispersal) on grids with power-law kernels. The bridge is weak.** The needed sums are low-order (single-edge normalizations, occasionally triangle-type closure terms), cheaply done by FFT convolution or in the continuum, and I found no 2025–26 lattice work where these sums are a bottleneck. A long-range contact-process series expansion (analogous to Dickman–Jensen series for the nearest-neighbour contact process; UNVERIFIED anchor) would be the only nontrivial use and is speculative.
7. **Riesz/Coulomb lattice energies and universally optimal lattices. The bridge is weak / already covered** by EpsteinLib and the n-body Epstein method ([2412.16317](https://arxiv.org/abs/2412.16317), [2504.11989](https://arxiv.org/abs/2504.11989)). Current results such as Luo–Wei ([2609.17356](https://arxiv.org/abs/2609.17356)) are analytic proofs.

A forced bridge to avoid in the final report: the claim that GZL could "compute the triangle diagram of long-range percolation near s≈d". The percolation triangle uses the unknown critical two-point function τ(x) ≍ |x|^{−d+α}, which is not a user-supplied kernel and is not in ℓ¹. Only the random-walk (mean-field) proxy triangle ∫D̂³/(1−D̂)³ is computable, and only indirectly.

### Gaps
- It is not established whether classical LR high-temperature series have ever been pushed beyond low order (older 1970s–1990s work is UNVERIFIED here). A literature check (e.g. Nagle–Bonner; Glumac–Uzelac; Monroe; Luijten's thesis) is needed before claiming "first high-order series".
- No quantitative estimate exists of the achievable series order for classical LR models with GZL (graph counts, treewidth distribution, runtime per graph).
- The practical accuracy of ratio or Padé analysis for LR models, where corrections to scaling come from ω ~ (2−σ) or logarithms, is unknown and could limit exponent estimates even with exact coefficients.
