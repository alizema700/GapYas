# Quantum simulators and AMO platforms with power-law interactions: where GZL could produce new results (2025–2026)

How I checked: I read every arXiv ID below on arxiv.org/abs (title, authors, dates and abstract) on 2026-09-25. For the papers I read in full HTML (2609.18761, 2603.20372, 2608.07178, 2311.11726, 2604.12754, 2606.13499, 2605.07685, 2601.20058), the details come from the full text. Items marked **[unverified]** were not read beyond the title or are known only from memory. Under each question, "Cited Findings" lists facts from the sources. "Inferences" gives my own judgement about how GZL applies.

Reference points for GZL itself:
- The physics paper, arXiv:2609.18761 (Duft, Adelhardt, Koziol, Buchheit, Schmidt, 16 Sep 2026), computes the long-range transverse-field Ising model (LRTFIM) on four lattices: chain, square, triangular and cubic. It covers ferromagnetic (FM) and antiferromagnetic (AFM) couplings and gives 0qp and 1qp series. Critical points and zν come from DlogPadé extrapolation of the 1qp gap at maximal order 11 (Fig. 14/15). — [arXiv:2609.18761](https://arxiv.org/abs/2609.18761)
- The method paper is arXiv:2609.18918 (Buchheit, Rupp). — [arXiv:2609.18918](https://arxiv.org/abs/2609.18918)

---

## Q1: Which 2025–2026 experiments measure dispersions, gaps, spectra, critical points or phase diagrams with power-law interactions on Bravais lattices? What theory do they compare against, and where is that theory limited?

### Takeaway
Rydberg tweezer arrays dominate this area. The two main settings are the van der Waals (vdW, 1/r^6) transverse-field Ising model (TFIM) on chains, square and triangular lattices, and the dipolar (1/r^3) XY model on square and kagome arrays. Trapped ions supply tunable-α Ising chains and 2D crystals. Momentum-resolved spectra now exist in several experiments: quench spectroscopy, modulation spectroscopy of CFT spectra, and dynamical structure factors. The theory they are compared against is finite-size (MPS, DMRG, tensor networks, tVMC, exact diagonalization), linear spin-wave (LSW) theory, or QMC. None of these experiments is compared against a converged thermodynamic-limit, full-Brillouin-zone calculation that keeps the complete power-law tail.

### Cited Findings

**Rydberg arrays, vdW Ising (the GZL target that exists today)**
- The Pasqal/Browaeys 256-qubit processor ran the 2D TFIM on a 16×16 square array of 87Rb atoms with vdW interactions U_ij = U(R/r_ij)^6 and J = U/4. The transverse field is tuned through the atom spacing R. The critical point used, h*_t/J ≈ 2.6 for L = 16, comes from QMC-SSE thermal matching. Tensor networks lose control at late times (tJ ≳ 1 for L = 16). The paper contains no excitation spectra or gaps. — [arXiv:2608.07178](https://arxiv.org/abs/2608.07178)
- Earlier Pasqal work benchmarks classical dynamics for the same 2D TFIM (MPS, TTN, 2D TN with belief propagation, NQS), published in Phys. Rev. Research 8, 023311 (2026). — [arXiv:2511.19340](https://arxiv.org/abs/2511.19340)
- A "one-to-one" simulation of the frustrated triangular magnet TmMgGaO4 used 256 qubits. Atoms sit on a triangular rhombus with R₁ = 9 µm and U_ij = C6/r^6. The material model is J1–J2 Ising with J2 ≈ 0.05 J1, transverse field Δx ≈ 1.08 J1 and longitudinal field Δz. The paper says only that the Rydberg interactions have "almost the same neighbour dependence". Comparisons are to DMRG, QMC and MPS/TDVP (D = 600). The MPS runs took about two weeks, versus about one day on the QPU. — [arXiv:2603.20372](https://arxiv.org/abs/2603.20372)
- Endres group, Caltech (Nature 657, 98 (2026) as cited in 2609.18761): modulation spectroscopy resolved the finite-size spectra of a Rydberg chain tuned to Ising and tricritical-Ising critical points. The same group measured the dynamical structure factor at criticality. — [arXiv:2601.16275](https://arxiv.org/abs/2601.16275)
- The GZL paper itself cites this as a direct comparison target: "As dynamical structure factors are now experimentally accessible in these platforms [Sun et al.], the momentum-resolved excitation spectra and, in future extensions, spectral weights produced by the graph zeta method can be compared directly with measurements." — [arXiv:2609.18761](https://arxiv.org/html/2609.18761)
- Meson/E8 spectroscopy in Rydberg chains and ladders (Pasqal, Vovrosh et al.) measured mass spectra at the Ising critical point with a longitudinal perturbation. — [arXiv:2506.21299](https://arxiv.org/abs/2506.21299)
- Samajdar group plus QuEra: a Lieb-lattice Rydberg array (multi-site basis) showed density-wave phases and "good agreement between theory and experiment". — [arXiv:2508.05737](https://arxiv.org/abs/2508.05737)
- A Sierpinski-fractal Rydberg array (88Sr, vdW TFIM) is not a Bravais lattice. — [arXiv:2509.03514](https://arxiv.org/abs/2509.03514)

**Rydberg arrays, dipolar XY (1/r^3)**
- Chen et al., Science 389, 483 (2025), an older anchor paper: quench spectroscopy on a 10×10 square array with a = 15 µm. The quantization axis is perpendicular to the plane "to ensure isotropic dipolar interactions", giving H = −(J/2)Σ(a³/r³)(σxσx+σyσy). The measured dispersion was compared to LSW (with modified PBC), tVMC and exact diagonalization on 4×4. For the ferromagnet the agreement was good. For the antiferromagnet: "the AFM data are in rather poor agreement with LSW theory, which does not predict any damping". The authors attribute this to one-to-three magnon decay. — [arXiv:2311.11726](https://arxiv.org/abs/2311.11726)
- A Dirac-spin-liquid candidate was studied with 114 dipolar Rydberg atoms on a kagome array (multi-site basis). — [arXiv:2602.14323](https://arxiv.org/abs/2602.14323)
- High-temperature-expansion thermometry of that experiment gives T = 0.55 J and S/N = 0.67 ln 2, published in PRB 114, L171112 (2026). — [arXiv:2604.22743](https://arxiv.org/abs/2604.22743)
- A Tomonaga–Luttinger liquid was observed in a dipolar XY Rydberg chain (PRX 15, 031021 (2025)). The dipolar tail renormalizes the TLL parameters relative to a nearest-neighbour model. — [arXiv:2501.08179](https://arxiv.org/abs/2501.08179)
- A dipolar Rydberg chain with 14 atoms showed weak integrability breaking. — [arXiv:2602.02251](https://arxiv.org/abs/2602.02251)
- A doped antiferromagnet with dipolar tunnelings was realized on a square array (Nature 644, 889 (2025)) as a t–J–V model. — [arXiv:2501.08233](https://arxiv.org/abs/2501.08233)
- STM-like spectral functions were measured for single holes and magnetic polarons. — [arXiv:2602.17600](https://arxiv.org/abs/2602.17600)
- "Kinetically-induced bound states in a frustrated Rydberg tweezer array" **[unverified, title only]**. — [arXiv:2510.17183](https://arxiv.org/abs/2510.17183)

**Trapped ions**
- Bounded-error quantum simulation (Innsbruck/Zoller, PRX 16, 031037 (2026)) used long-range Ising interactions with up to 51 ions. Hamiltonian learning infers the actual J_ij rather than assuming a power law. — [arXiv:2511.23392](https://arxiv.org/abs/2511.23392)
- Bubble nucleation was observed in a mixed-field Ising chain with tunable, time-dependent interactions (Monroe group). — [arXiv:2505.09607](https://arxiv.org/abs/2505.09607)
- The field-theory counterpart treats false-vacuum decay in a 1/r^α mixed-field Ising chain (σ = α − 1). — [arXiv:2607.03274](https://arxiv.org/abs/2607.03274)
- A 2D long-range XY model with 200 Yb+ ions (Duan group) has couplings engineered through one phonon mode, J_ij ∝ b_ik b_jk. This is **not a power law**. No dispersion was measured, and the comparison was to mean-field and Holstein–Primakoff theory. — [arXiv:2606.13499](https://arxiv.org/html/2606.13499)
- Topological spin textures were prepared in a 2D Penning-trap crystal of more than 150 ions. — [arXiv:2604.13872](https://arxiv.org/abs/2604.13872)
- A proposal for Rydberg ions in a Penning trap gives MHz dipolar interactions in 2D (PRX Quantum 7, 033064 (2026)). — [arXiv:2601.01626](https://arxiv.org/abs/2601.01626)
- A 20-qubit trapped-ion QPU ran a numerical linked-cluster expansion with a quantum algorithm (NLCE+QA, with K. P. Schmidt) to get thermodynamic-limit ground-state energies and 1qp dispersions. This was **only for the nearest-neighbour TFIM** (chain, ladder, and the chain with a longitudinal field). — [arXiv:2605.28599](https://arxiv.org/abs/2605.28599)

**Polar molecules and magnetic atoms**
- Erbium quantum gas microscope (Greiner/Grusdt): a mixed-dimensional XXZ model and Z2 lattice gauge theory in coupled chains, with stripe melting observed. — [arXiv:2509.16200](https://arxiv.org/abs/2509.16200)
- Chromium (spin-3) itinerant magnetism in a 3D lattice, PRL 136, 103401 (2026). — [arXiv:2501.11402](https://arxiv.org/abs/2501.11402)
- Theory for microwave-dressed molecules with a sign-changing anisotropic tail (x²−y²)/r^5 on the square lattice, using worm QMC. — [arXiv:2605.25019](https://arxiv.org/abs/2605.25019)
- In 2025–2026 I found **no** molecule or magnetic-atom experiment that measures a full dispersion or critical point on a Bravais lattice with an isotropic tail.

**Non-AMO comparison point**
- Google measured 2D XY magnon spectra and lifetimes across the Brillouin zone on up to 97 superconducting qubits. MPS "become inaccurate" away from small systems or low temperature. This is a nearest-neighbour coupler graph, not a power law. — [arXiv:2607.13301](https://arxiv.org/abs/2607.13301)

### Inferences
- The Rydberg vdW Ising Hamiltonian H = Σ U_ij n_i n_j + (Ω/2)Σσx − δΣn maps onto the shipped AFM LRTFIM with ν = 6 **only on the particle–hole-symmetric line**. There the longitudinal field vanishes: δ = Σ_j U_ij / 2, a single Epstein-zeta sum. Covering the full (Ω, δ) phase diagram needs a longitudinal field in the pCUT corpus, which is not shipped.
- Pasqal writes the model as a TFIM with J = U/4 and FM/PM phases (2608.07178). How they handle the sign or sublattice mapping for repulsive vdW needs checking in the full text **[unverified detail]**.
- The TmMgGaO4 mapping (2603.20372) is a clear GZL target. On the triangular lattice the vdW tail gives J2/J1 = (1/√3)^6 = 1/27 ≈ 0.037 and J3/J1 = 1/64, while the material model has J2 ≈ 0.05 J1 and no J3. GZL can quantify exactly how the infinite 1/r^6 tail and the compact correction a(x) (for example, adding back ΔJ2) shift the gap, the dispersion and the transverse-field critical point in the thermodynamic limit. This is the same logic as the KTmSe₂ comparison in 2609.18761, applied to a simulator–material pair.
- Trapped-ion J_ij are not pure power laws: they come from finite, inhomogeneous crystals and are learned by Hamiltonian learning. GZL's infinite-Bravais-lattice idealization is therefore a reference model, not a one-to-one prediction for ion experiments. Rydberg and molecule tweezer arrays are much closer to "ideal lattice plus exact power law".

### Gaps
- I found no 2025–2026 experiment reporting a full-Brillouin-zone 1qp dispersion of a vdW Rydberg TFIM in 2D. Only the 1D structure factor (2601.16275) and quench or dynamics data exist. Absence from my search is not proof that none exists.
- I found no 2025–2026 3D quantum-simulator experiment with an isotropic power law on a cubic lattice.
- The arXiv ID of the earlier trapped-ion quasiparticle-spectroscopy anchor (Jurcevic et al., PRL 2015) was not checked **[unverified]**.

---

## Q2: Which recent theory papers use series, pCUT, linked-cluster, spin-wave or perturbative methods with truncated or Monte-Carlo-evaluated lattice sums, and what accuracy limits would GZL remove?

### Takeaway
The pCUT+MC line from the Erlangen group is the direct predecessor. Its stated bottleneck is statistically noisy, HPC-scale MC embedding that must be rerun for every α, k and anisotropy, and 2609.18761 removes exactly that. Most other 2025–2026 theory still truncates the tail (DMRG with R_max = W/2), approximates it with sums of exponentials (iMPS), or uses LSW. Several of these papers state that such approximations bias the results.

### Cited Findings
- 2609.18761 on pCUT+MC: "MC simulations converge slowly, with statistical errors scaling as ∼N_steps^(-1/2)… Simulations must be repeated for each set of parameters, including the decay exponent, the momentum, and additional model parameters, such as interaction anisotropies. In practice, this restricts the method to a relatively limited set of accessible data points." The paper lists earlier pCUT+MC applications to Ising, XY and XXZ models (anisotropic XY chain, PRB 102, 174424 (2020); XXZ bilayer, PRB 111, 024409 (2025)) and to spin-1 models. — [arXiv:2609.18761](https://arxiv.org/html/2609.18761)
- The long-range spin-1 Heisenberg chain with single-ion anisotropy (Adelhardt, Muleady, Schmidt, Gorshkov; accepted in PRR) uses pCUT+MC to **order 10 for the gap** and order 9 for spectral weights, together with MPS. It finds continuously varying critical exponents versus α. Open problems it states: "the precise properties of the quantum phase transition along the Gaussian line as well as the mechanism underlying the continuously varying critical exponents". It names trapped ions, ultracold atoms and Rydberg arrays as platforms. — [arXiv:2604.12754](https://arxiv.org/html/2604.12754)
- Koziol and Schmidt treat long-range Ising models on triangular and kagome lattices in the thermodynamic limit with a unit-cell quantum-annealing scheme, for "lattice problems with resummable long-range interactions". — [arXiv:2511.08336](https://arxiv.org/abs/2511.08336)
- Devil's staircases in the long-range Dicke–Ising model on square and triangular lattices, treated with unit-cell mean-field theory plus wormhole QMC. — [arXiv:2503.02734](https://arxiv.org/abs/2503.02734)
- Dipolar XY DMRG on nine Archimedean lattices: "We choose to simply truncate, so that only sites with distance r_ij < R_max = W/2 are coupled". LSW overestimates susceptibility by about 10–30% and is off in stiffness by up to about 40% compared with DMRG. — [arXiv:2605.07685](https://arxiv.org/html/2605.07685)
- Quench spectroscopy for the 2D XY model with power-law interactions in a staggered field, treated analytically with LSW and a low-energy phase theory (Brechtelsbauer, Büchler). — [arXiv:2608.16631](https://arxiv.org/abs/2608.16631)
- iMPS without a surrogate: finite-pole sum-of-exponentials surrogates "can bias a critical diagnostic away from criticality", while summing the exact algebraic tail through a polylog of the transfer matrix does not. — [arXiv:2606.20522](https://arxiv.org/abs/2606.20522)
- Rigorous bound: truncating a two-body r^−p tail at range R changes local ground-state observables by O(R^−(p−d)) for p > 2d, and the bound is optimal. — [arXiv:2608.15576](https://arxiv.org/abs/2608.15576)
- Tensor-network VMC with O(N) local energy for all-to-all interactions reproduced the 10×10 dipolar XY Rydberg adiabatic protocol, which was "previously beyond the reach of classical simulation". — [arXiv:2607.05178](https://arxiv.org/abs/2607.05178)
- 2D LRTFIM quench dynamics were studied with neural quantum states plus an effective theory. Long-lived oscillations come from magnon bound states formed by attractive long-range magnon interactions. — [arXiv:2512.09037](https://arxiv.org/abs/2512.09037)
- A universal squeezing dynamical transition in power-law bilayer XXZ models on square, triangular and honeycomb lattices, with Bogoliubov scaling a_Z* ∝ L^(2/(α−d)) for α > d+2. — [arXiv:2605.13969](https://arxiv.org/abs/2605.13969)
- NLCE+QA on a trapped-ion QPU for nearest-neighbour TFIM dispersions. — [arXiv:2605.28599](https://arxiv.org/abs/2605.28599)

### Inferences
- For every pCUT+MC model already set up in Erlangen, GZL turns the published sparse α grids into dense, deterministic scans. This covers the LR XY/XXZ chain, the XXZ bilayer and the spin-1 Heisenberg chain. The graph-zeta step is a drop-in replacement, but the white-graph corpora for those models must be produced (not shipped).
- The bilayer and ladder kernels have the form (r² + d²)^(−α/2) between dimers. This is not a pure power law of the in-plane lattice vector. I think it can still be represented with the allowed kernel form K = a(x) + Σ_j b_j|x|^(−ν_j): an asymptotic series Σ b_j |x|^(−α−2j) for the tail plus a compact correction a(x) for the near field. This is an inference that needs a convergence check. If it works, dimer and bilayer models (squeezing bilayers, XXZ bilayers) become accessible without a multi-site basis.
- The DMRG truncation to R_max = W/2 in 2605.07685 and the LSW mismatches suggest a clean benchmark. GZL-powered series for the square-lattice dipolar XXZ/XY model, expanded from a gapped limit (Ising-anisotropic or staggered-field), would give thermodynamic-limit numbers with the complete tail.
- Speculative: the loop corrections of 1/S spin-wave theory are Brillouin-zone convolutions of J(k). GZL's Epstein-zeta product/convolution algebra might speed these up for AFM magnon decay (2311.11726). That would be a new use of the library beyond linked-cluster expansions.

### Gaps
- I did not find which exact orders and error bars the earlier pCUT+MC XY/XXZ papers reached; they are referenced in 2609.18761 but I did not read them.
- It is unclear whether GZL's per-graph external momentum can take extra staggered phases e^{iπ·x} on individual edges. That feature would matter for staggered-field XY models, and I did not check it.

---

## Q3: Are there open disputes where precise high-order series could decide?

### Takeaway
Yes, several concrete ones, but most need corpora beyond the shipped LRTFIM. The one open question that is fully inside today's GZL is the nature of the triangular-lattice AFM LRTFIM transition at small σ (first-order or continuous), and the dense-σ map of crossover boundaries. The biggest disputes about α_c in continuous-symmetry chains need XY, XXZ and spin-1 corpora.

### Cited Findings
- 2609.18761 on the triangular AFM LRTFIM: the model "remains in that universality class [3d XY] down to σ ≈ 1… For smaller σ there is evidence for a potential first-order transition towards a stripe-ordered low-field phase". It cites Humeniuk (PRB 93, 104412 (2016), QMC), Koziol et al. (PRB 100, 144411 (2019)), Fey et al. (PRL 122, 017203 (2019)) and Saadatmand et al. (PRB 97, 155116 (2018), iDMRG). The same paper notes that rounding at regime boundaries is "typically attributed to finite perturbative orders" and that zν is overestimated in the nearest-neighbour limit. — [arXiv:2609.18761](https://arxiv.org/html/2609.18761)
- LR spin-1/2 XY chain: "LSW predicts a phase boundary of αc = 3.0". SSE-QMC finds a minimum in the normalized superfluid density at αc ≈ 2.7, "in line with perturbative field theory predictions". Perturbative RG and DMRG also suggest αc slightly below 3. — [arXiv:2601.20058](https://arxiv.org/html/2601.20058)
- LR spin-1 Heisenberg chain (staggered power law): QMC in the split-spin representation gives αc = 2.49(1) for the Haldane-to-Néel transition, which is nonconformal with z ≠ 1. Version 2 "improved location of phase transition point". — [arXiv:2604.20831](https://arxiv.org/abs/2604.20831)
- The pCUT+MC plus MPS study of the same family finds continuously varying exponents with a strong boundary-condition dependence and leaves the mechanism open. — [arXiv:2604.12754](https://arxiv.org/abs/2604.12754)
- Classical 2D LR Heisenberg: MC up to L = 8192 places the LR–SR crossover at σ* = 2. — [arXiv:2512.01956](https://arxiv.org/abs/2512.01956)
- BKT persists for all LR exponents in the 2D XY model, "distinct from previous results", with explicit relevance to Rydberg experiments. — [arXiv:2511.07305](https://arxiv.org/abs/2511.07305)
- EnLRO/ReLRO crossover at σ ≈ 1.575 (1D) and σ ≈ 3.2 (2D) in classical LR XY/Heisenberg models. — [arXiv:2507.10018](https://arxiv.org/abs/2507.10018)
- The dipolar XY model on the triangular lattice has competing coplanar, stripe and possible spin-liquid phases whose "relative stability is sensitive to the long-range couplings". — [arXiv:2605.07685](https://arxiv.org/abs/2605.07685)

### Inferences
- **Triangular AFM LRTFIM, available today.** Deterministic 1qp gaps and 0qp energies on a dense σ grid can locate where DlogPadé stops showing a clean 3D-XY pole. They can also track the momentum of the gap minimum over the full BZ, where a jump signals a competing stripe instability. A gap-based series cannot prove a first-order transition by itself. It would need ground-state-energy comparisons with the stripe phase from the ordered side, which is another expansion and another corpus. For ν = 6 (Rydberg), the physical point σ = ν − d = 4 lies deep in the nearest-neighbour-like regime. The dispute therefore matters for trapped-ion-like α, not for Rydberg vdW.
- **α_c of the LR XY chain (LSW 3.0 vs about 2.7) and of the LR spin-1 chain (QMC 2.49(1) vs pCUT+MC/MPS).** Noise-free high-order series with dense α sampling could settle the location and the continuously varying exponents. This needs XXZ/XY and spin-1 corpora, which would be a GZL extension.
- **Crossovers between universality regimes in the Ising case (d = 1, 2, 3).** 2609.18761 already shows dense σ scans (Fig. 14/15). A dedicated high-resolution study around σ_uc = 2d/3 (with multiplicative logs) and the σ* boundary, including d = 3 on the cubic lattice, is in scope today. It is also a natural benchmark table for trapped-ion and 3D simulators.
- **Deconfinement, Higgs or Lifshitz points** with power-law tails: I found no 2025–2026 source that frames such a dispute in a GZL-compatible model. (2509.16200 studies Z2 LGT confinement with dipolar erbium, but in a mixed-dimensional XXZ model.)

### Gaps
- I found no explicit 2025–2026 paper that disputes 3D quantum Ising criticality with dipolar or vdW tails on cubic lattices.
- I did not check the order and value that 2604.12754 gives for αc of the Haldane-to-Néel transition, so I could not say whether it disagrees with the QMC value 2.49(1).

---

## Q4: Which applications need only today's GZL, and which need roadmap features?

### Takeaway
Today's GZL covers Bravais lattices, isotropic power laws plus a compact anisotropic part, and the LRTFIM 0qp/1qp corpora. That is enough for Rydberg vdW Ising arrays on chain, square and triangular lattices on the symmetric-detuning line, for 1D and 2D ion-like tunable-α Ising models, for 3D cubic benchmarks, and for "simulator vs material" comparisons in the style of KTmSe₂ or TmMgGaO4. Most of the popular recent Rydberg geometries need a multi-site basis: kagome, Lieb, ruby, breathing kagome and Archimedean tilings. Molecule and magnetic-atom setups with in-plane or 3D dipoles need the anisotropic kernel. The popular XY, XXZ and t–J models need new corpora.

### Cited Findings
- 2609.18761 limits Bravais lattices to the LRTFIM with isotropic power laws. It says "The pCUT formalism can be extended to multi-site unit cells and higher qp sectors as well". Its future developments are "extension to multi-atomic lattices, broader classes of interaction kernels… spectral weights and higher quasiparticle sectors… finite temperatures". — [arXiv:2609.18761](https://arxiv.org/html/2609.18761)
- The 3D KTmSe₂ model with angle-dependent dipolar coupling was computed only to order 5 with "an experimental code extension to anisotropic dipolar sums… a stable and efficient implementation in the Graph Zeta Library is still under development". Named 3D dipolar targets are LiHoF₄, Fe₈, Mn₁₂ acetate and RE(OH)₃. — [arXiv:2609.18761](https://arxiv.org/html/2609.18761)
- A new numerical method for anisotropic Epstein zeta functions reaches machine precision "for any lattice, power-law decay exponent and anisotropy order". It is the mathematical basis for the anisotropic roadmap. — [arXiv:2609.28282](https://arxiv.org/abs/2609.28282)
- Rydberg dipolar XY with the quantization axis perpendicular to a 2D array has **isotropic** 1/r^3 couplings (Browaeys square array). — [arXiv:2311.11726](https://arxiv.org/html/2311.11726)
- Multi-site geometries in current experiments and proposals include kagome with 114 atoms ([arXiv:2602.14323](https://arxiv.org/abs/2602.14323)), breathing kagome with a dipolar chiral spin liquid ([arXiv:2603.25784](https://arxiv.org/abs/2603.25784)), the Lieb lattice ([arXiv:2508.05737](https://arxiv.org/abs/2508.05737)), Archimedean lattices ([arXiv:2605.07685](https://arxiv.org/abs/2605.07685); [arXiv:2505.01513](https://arxiv.org/abs/2505.01513)), and ruby with spin lakes ([arXiv:2512.09040](https://arxiv.org/abs/2512.09040), title only **[unverified content]**).

### Inferences

| Application | Today | Needs roadmap or extension |
|---|---|---|
| vdW (ν=6) AFM LRTFIM on chain, square, triangular (Rydberg, symmetric line) | yes | longitudinal field / detuning for the full (Ω, δ) diagram |
| Tunable-α FM/AFM LRTFIM chain (ion-chain reference), 2D triangular (Penning-like) | yes | real ion J_ij are not pure power laws |
| LRTFIM on cubic (3D) with isotropic tail, crossovers | yes | — |
| Triangular Ising with J1 + vdW or dipolar tail vs material (TmMgGaO4, KTmSe₂) | yes (2D, isotropic in-plane dipolar) | 3D stacked dipolar needs anisotropy |
| Isotropic dipolar XY/XXZ on square/triangular (perpendicular quantization axis) | kernel yes | XY/XXZ corpora and a gapped expansion limit |
| LR spin-1 Heisenberg chain, LR XY chain α_c | kernel yes | spin-1 and XXZ corpora |
| Kagome, Lieb, ruby, honeycomb Rydberg arrays | no | multi-site basis |
| Molecules or magnetic atoms with in-plane/3D dipoles, LiHoF₄ | no | anisotropic 1−3cos²θ kernel |
| Magnon bound states (2512.09037), spectral weights / DSF (2601.16275) | no | 2qp sector, observables |
| Bilayer or dimer models with (r²+d²)^(−α/2) | maybe (see Q2 inference) | dimer corpus plus kernel-series check |

### Gaps
- I did not verify which corpora beyond the LRTFIM the GZL GitHub repository ships (for example, whether any XXZ corpus exists in a development branch).

---

## Q5: For each candidate, what exactly would be the new result?

### Takeaway
There are four strong, near-term candidates that use today's features, plus three high-leverage extension projects. The best "new number vs experiment" project is the thermodynamic-limit, full-BZ spectrum and critical line of the vdW-tailed AFM TFIM on the lattices used by the 256-qubit Pasqal machines and the Caltech chain. That result would serve as a quasi-exact reference for those machines and for the TmMgGaO4 mapping.

### Cited Findings (anchors per candidate)
- C1: [arXiv:2608.07178](https://arxiv.org/abs/2608.07178), [arXiv:2511.19340](https://arxiv.org/abs/2511.19340), [arXiv:2603.20372](https://arxiv.org/abs/2603.20372)
- C2: [arXiv:2601.16275](https://arxiv.org/abs/2601.16275), [arXiv:2506.21299](https://arxiv.org/abs/2506.21299)
- C3: [arXiv:2609.18761](https://arxiv.org/html/2609.18761), [arXiv:2601.01626](https://arxiv.org/abs/2601.01626)
- C4: [arXiv:2605.28599](https://arxiv.org/abs/2605.28599), [arXiv:2511.23392](https://arxiv.org/abs/2511.23392)
- C5: [arXiv:2311.11726](https://arxiv.org/abs/2311.11726), [arXiv:2605.07685](https://arxiv.org/abs/2605.07685), [arXiv:2608.16631](https://arxiv.org/abs/2608.16631)
- C6: [arXiv:2601.20058](https://arxiv.org/abs/2601.20058), [arXiv:2604.20831](https://arxiv.org/abs/2604.20831), [arXiv:2604.12754](https://arxiv.org/abs/2604.12754)
- C7: [arXiv:2512.09037](https://arxiv.org/abs/2512.09037)

### Inferences (candidate list; the new results are my own proposals)

1. **C1: Thermodynamic-limit reference for 2D Rydberg vdW Ising arrays (square and triangular). Today's features.**
   - New result: the 1qp dispersion over the full BZ, the gap and the critical point λ_c for K(x) = U(a/|x|)^6 on the square and triangular lattices, with DlogPadé zν.
   - Scan the compact part a(x) to model deviations from pure C6, or to reproduce the TmMgGaO4 J2 ≈ 0.05 J1 against the vdW 1/27.
   - Compare λ_c with the QMC-based estimate h*/J ≈ 2.6 (L = 16) used by Pasqal, and give a quantified "tail correction" to the material mapping in 2603.20372.
   - Novelty: no thermodynamic-limit spectral reference with the full vdW tail exists for these experiments; they use finite-size DMRG, QMC and MPS.
   - Caveat: this holds only on the Z2-symmetric detuning line.

2. **C2: 1D Rydberg chain at criticality (Caltech CFT-spectra experiment). Today's features for the Ising line.**
   - New result: the vdW-tailed AFM LRTFIM chain dispersion and λ_c at high precision. These give the non-universal velocity and gap needed to convert measured finite-size CFT ratios and the dynamical structure factor into absolute energies.
   - Spectral weights (the DSF itself) are on the roadmap.
   - The tricritical-Ising and E8 (longitudinal-field) settings need new corpora.

3. **C3: A 3D and dimension-crossover benchmark atlas. Today's features.**
   - New result: dense-σ tables of λ_c and zν on cubic (plus square, triangular, chain) for FM and AFM. Crossover boundaries σ_uc = 2d/3 and the NN limit, with published deterministic coefficients.
   - Use as benchmark data for trapped-ion and Penning (2D triangular) simulators and for classical methods such as TN-VMC (2607.05178) and iMPS (2606.20522).
   - Novelty: the previous MC data were sparse; 2609.18761 shows this is feasible but frames it as a method demonstration.

4. **C4: Long-range extension of NLCE+QA benchmarks. Today's features.**
   - The trapped-ion NLCE+QA paper computed only nearest-neighbour TFIM dispersions.
   - New result: exact long-range reference dispersions (α = 1–3 chain) against which hybrid quantum–classical linked-cluster schemes, or Hamiltonian-learned ion models, can be scored.
   - Cheap, and it involves the same group (Schmidt).

5. **C5: Dipolar XY/XXZ on the square lattice (Browaeys quench spectroscopy). Needs an XXZ/XY corpus.**
   - New result: beyond-LSW, thermodynamic-limit magnon dispersions with the complete 1/r^3 tail. One route is an expansion from a gapped limit (Ising-anisotropy or staggered field, as in Büchler's quench-spectroscopy setting); the other is GZL-accelerated 1/S loop sums.
   - These could address the AFM damping and nonlinearity puzzle in 2311.11726 and the 10–40% LSW vs DMRG discrepancies in 2605.07685.

6. **C6: α_c and continuously varying exponents in LR XY and spin-1 Heisenberg chains. Needs XXZ and spin-1 corpora.**
   - Existing values: LSW 3.0 vs QMC/DMRG/RG ≈ 2.7 for the LR XY chain; QMC αc = 2.49(1) for the spin-1 Haldane–Néel transition.
   - New result: deterministic, densely sampled high-order series settling these values.
   - The pCUT+MC infrastructure for these models already exists in Erlangen.

7. **C7: Magnon bound states in 2D LRTFIM (Naik & Heyl). Needs the 2qp sector.**
   - New result: the exact two-magnon bound-state spectrum over the BZ versus α. This would give a quantitative prediction for the persistent oscillations seen in NQS and in Rydberg or ion quench experiments.

Lower priority or blocked:
- Kagome, breathing-kagome and Lieb Rydberg spin liquids (multi-site).
- Anisotropic molecule stripes (2605.25019) and LiHoF₄ (anisotropic kernel).
- Cavity-mediated Dicke–Ising (2503.02734, 2605.27484): the photon mode is all-to-all, not a lattice power law, so GZL helps only with the LR Ising part.
- Engineered-mode ion couplings (2606.13499): not a power law.

### Gaps
- I did not check the GZL repository for the exact API support of a longitudinal field, extra observables (spectral weights) or 2qp sectors. The table assumes the limits stated in the task brief and in 2609.18761.
- Whether the experimental groups named here (Pasqal/Browaeys, Endres, Monroe, Roos/Blatt, Duan) plan follow-ups that would publish full-BZ Rydberg dispersions is unknown.
