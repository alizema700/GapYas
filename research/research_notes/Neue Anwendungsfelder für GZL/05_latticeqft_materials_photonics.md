# GZL transfer opportunities: lattice field theory, electrostatics/materials, photonics/metamaterials

Scope note: I verified every arXiv ID below by fetching its arxiv.org abstract page, or by finding it in an arxiv.org search listing (marked "listing"). The dates are v1 submission dates where the abs page gave one. Any claim not backed by a fetched source is labelled as an inference or listed under Gaps. GZL facts come from the GZL GitHub README and the two GZL papers.

Baseline GZL facts used throughout (all verified):
- GZL is "Graph Zeta Library for the efficient and precise computation of high-dimensional lattice sums occurring in linked-cluster expansions". Kernels are |x|^(-ν) for x≠0, with K(0)=0, so **coincident vertices are excluded**. Kernels can be finite sums of power laws plus a "real, even, compactly supported short-range part a". Each edge can carry its own kernel. It works in d=1,2,3 on general Bravais lattices A∈R^(d×d), and "Re(ν_e) > d on every edge" is stated as a *sufficient* convergence condition. The full BZ is obtained via FFT "at the cost of a single momentum evaluation". Main API: `evaluate_graph`, `evaluate_corpus`, `compute_series_coefficients`. The only shipped model so far is the LR transverse-field Ising model; the roadmap mentions the Heisenberg model. — [GZL GitHub](https://github.com/graph-zeta/gzl)
- The method paper covers blocks of treewidth ≤2, which it computes at linear cost with Epstein-zeta algebra, and uses tensor-network bucket elimination for tw>2. It cuts series-expansion cost "from tenthousands of core-hours to minutes". — [arXiv:2609.18918](https://arxiv.org/abs/2609.18918) (Buchheit & Rupp, 16 Sep 2026)
- The physics application reproduces MC data for the LR-TFIM in 1d/2d/3d and models KTmSe2 with dipolar interactions. — [arXiv:2609.18761](https://arxiv.org/abs/2609.18761) (Duft, Adelhardt, Koziol, Buchheit, Schmidt, 16 Sep 2026)

---

## Q1. Which recent works need graph-structured power-law lattice sums, and how are they computed now?

### Takeaway
The clearest real consumer outside linked-cluster expansions is **many-body (≥3-body) dispersion lattice energies of crystals**: the Axilrod–Teller–Muto (ATM) triangle sum, its strain derivatives, and higher n-body terms. The Buchheit group has **already solved the plain ATM lattice sum** (2025, "weeks to minutes"), so GZL's added value there is k-dependence, per-edge kernels and higher graphs, not the ATM energy itself. In lattice QFT (Lüscher zeta functions, QED_L finite-volume coefficients) and photonics (coupled-dipole arrays), the sums in current use are mostly **single lattice sums**, or kernels outside GZL's class (oscillatory, pole-type, d=4). There the graph machinery is not the bottleneck.

### Cited Findings

**A. Many-body dispersion / cohesive energies (materials, chem-ph)**
- Buchheit & Busse, "Epstein zeta method for many-body lattice sums" (v1 16 Apr 2025, revised Jun 2026). They write n-body lattice sums as "singular integrals over products of Epstein zeta functions". For ATM in 3D this cuts runtime "from weeks to minutes", and the cost grows "only linearly with n". They report "full precision for exponents greater than the system dimension". In an LJ+ATM model they find an fcc→bcc transition at finite pressure as the ATM coupling grows. — [arXiv:2504.11989](https://arxiv.org/abs/2504.11989)
- Robles-Navarro, Cooper, Buchheit, Busse, Burrows, Smits, Schwerdtfeger (9 Apr 2025) study Bain-path cuboidal lattices (fcc, mcc, bcc, acc; **all Bravais**) with a general (n,m) LJ potential plus repulsive ATM. They evaluate the three-body sums "and their meromorphic continuations to machine precision within minutes on a standard laptop" by turning the sum into an integral over products of Epstein zeta functions. The ATM energy has an extremum at bcc, and strong repulsive ATM can make bcc favourable for soft LJ. Because of the "wrong short-range behavior of the triple-dipole ATM model", the LJ repulsion must satisfy n>9. — [arXiv:2504.07338](https://arxiv.org/abs/2504.07338)
- Related Schwerdtfeger/Cooper cuboidal-lattice-sum series (two-body level): [arXiv:2501.05746](https://arxiv.org/abs/2501.05746) (Cooper & Schwerdtfeger, "A Minimum Property for Cuboidal Lattice Sums", 2025); [arXiv:2502.09828](https://arxiv.org/abs/2502.09828) (hcp→cuboidal phase-transition mechanism, Feb 2025); [arXiv:2406.09635](https://arxiv.org/abs/2406.09635) (lattice instabilities hcp→cuboidal, revised Oct 2025). Older anchors: [arXiv:2105.08922](https://arxiv.org/abs/2105.08922) (cuboidal lattice sums), [arXiv:2012.05413](https://arxiv.org/abs/2012.05413) (analytical LJ vibrational effects in cubic and hcp lattices, Einstein model), [arXiv:2107.11380](https://arxiv.org/abs/2107.11380) (bcc instability from exact lattice sums). — all listing
- In DFT dispersion corrections for layered materials, Rumson, Bryenton & Johnson (8 Apr 2026) show that adding ATM three-body terms to XDM(BJ)/XDM(Z) improves exfoliation energies on the LM26 benchmark (graphite, h-BN, PbO, TMDs). They report "the best performance achieved on LM26 using semi-local functionals to date" relative to RPA. — [arXiv:2604.06539](https://arxiv.org/abs/2604.06539). The companion XDM(Z) molecular-crystal work is [arXiv:2602.04172](https://arxiv.org/abs/2602.04172) (listing). *The lattice-summation technique they use for ATM (real-space cutoff or otherwise) is not given in the abstract; see Gaps.*
- Many-body dispersion (MBD, Tkatchenko): in libMBD the periodic dipole tensor is evaluated "with Ewald summation" with real- and reciprocal-space cutoffs, and the energy is a Brillouin-zone sum. — [libMBD arXiv:2308.03140](https://arxiv.org/pdf/2308.03140). A 2016 reciprocal-space MBD implementation also expresses the energy as BZ-sampled contributions — [J. Phys.: Condens. Matter 28, 045201](https://iopscience.iop.org/article/10.1088/0953-8984/28/4/045201). Recent MBD work in 2026 is about ML surrogates and force decompositions, not lattice-sum precision — [arXiv:2602.22086](https://arxiv.org/abs/2602.22086) (MBD-ML, Feb 2026); [arXiv:2603.28518](https://arxiv.org/abs/2603.28518) (structured force reformulation, Mar 2026).

**B. Madelung / electrostatic lattice sums (single sums, very active in 2025–26)**
- He & Hu (10 Aug 2026) give a closed-form boundary term for conditionally convergent direct Madelung sums in arbitrary triclinic lattices. The residual finite-size correction decays as (2p+1)^(-2). — [arXiv:2608.10041](https://arxiv.org/abs/2608.10041)
- Zhao, He & Hu (Dec 2025) study the Madelung problem for finite crystals. The leading finite-size correction is [24r^4−40(x^4+y^4+z^4)]/[9√3(2p+1)^2] for cubic crystals, which is accurate already at 3^3 cells. — [arXiv:2512.18138](https://arxiv.org/abs/2512.18138)
- Calara & Miller (Mar 2025, revised Aug 2026) use axial-multipole real-space summation with r^(-13) far-field decay. They get the NaCl Madelung constant to 13 decimals and validate up to d≤6. — [arXiv:2503.00977](https://arxiv.org/abs/2503.00977)

**C. Periodic long-range sums in magnetism (Buchheit group, single-sum/Epstein level)**
- Zeta expansion for power-law interactions of bodies with periodic images in micromagnetics (Sep 2025). Common practice is to truncate image sums "introducing uncontrolled errors"; the new method reaches machine precision. — [arXiv:2509.26274](https://arxiv.org/abs/2509.26274)
- Anisotropic Epstein zeta functions (23 Sep 2026) are built from wave-vector derivatives of isotropic lattice sums and reach machine precision "for any lattice, power-law decay exponent and anisotropy order". Rayleigh–Wood singularities are removed analytically, with the dipole interaction as the target. — [arXiv:2609.28282](https://arxiv.org/abs/2609.28282). **This is the building block the GZL roadmap needs for anisotropic dipolar edges.**
- Epstein zeta function computation and properties (EpsteinLib background). — [arXiv:2412.16317](https://arxiv.org/abs/2412.16317)
- Classical dipolar Heisenberg models on the 11 Archimedean and 8 Laves lattices (Aug 2026). The paper computes the Fourier interaction matrix A(k), Luttinger–Tisza wave vectors and linear spin waves, and gives canting angles "to thirty digits". Truncation tests show the canting-angle value "is set by the long-range tail". — [arXiv:2609.02934](https://arxiv.org/abs/2609.02934). These are **non-Bravais lattices with anisotropic dipolar kernels**, and at linear spin-wave level they are single sums.
- DMRG study of dipolar XY models on Archimedean lattices (polar molecules, Rydberg atoms), compared against linear spin-wave theory. — [arXiv:2605.07685](https://arxiv.org/abs/2605.07685)

**D. Lattice QFT**
- Modified Lüscher zeta function in the presence of a long-range force (Bubna, Hammer, Hoid, Pang, Rusetsky, Wu; Jul 2025) uses "an efficient numerical algorithm" and discusses regularization and renormalization and cutoff problems. — [arXiv:2507.18399](https://arxiv.org/abs/2507.18399)
- Discretization effects built into Lüscher's formalism give "modified versions of the finite-volume zeta functions" with two sets of angular-momentum indices. — [arXiv:2408.07062](https://arxiv.org/abs/2408.07062) (Hansen & Peterken, 2024, anchor)
- The three-neutron quantization condition implementation (Schaaf & Sharpe, Dec 2025) — [arXiv:2512.24508](https://arxiv.org/abs/2512.24508) (listing). Staggered ππ Lüscher formalism — [arXiv:2603.28967](https://arxiv.org/abs/2603.28967) (listing).
- Finite-volume QED: the QED_r action (Di Carlo, Hansen, Hermansson-Truedsson, Portelli; Jan 2025) removes kinematics-independent O(1/L^3) effects of QED_L and pushes the leading contamination to O(1/L^4). — [arXiv:2501.07936](https://arxiv.org/abs/2501.07936). Structure-dependent O(1/L^3) effects — [arXiv:2310.13358](https://arxiv.org/abs/2310.13358) (listing).
- Power-law finite-volume coefficients are single sums C_j = L^(-3) Σ_{n≠0} |2πn/L|^(-j). For point charges, **the two- and three-loop O(α^2), O(α^3) double and triple mode sums factorize into products of single sums** (e.g. Σ^(2,1) ∝ C_3 C_2, Σ^(3,1) ∝ C_2[(C_3)^2 + 2C_2C_4]), and the power-law FV effects cancel after renormalization. — [Matzelle & Tiburzi, arXiv:1702.01296](https://arxiv.org/html/1702.01296) (anchor, 2017)
- Finite-volume formalism for electroweak loop integrals gives power-law FV correction formulas for two-particle intermediate states. — [arXiv:2407.16930](https://arxiv.org/abs/2407.16930) (Tuo & Feng, 2024, revised 2025)
- Position-space photon propagator in QED corrections to lattice QCD (Erb, Meyer, Ottnad; Mar 2026). The bottleneck is the "volume-squared sum over the endpoints of the photon propagator". They factorize it through integral representations (Fourier, or 5D-propagator autoconvolution) and test on a 48^3×128 ensemble. — [arXiv:2603.13086](https://arxiv.org/abs/2603.13086)
- Euclidean coordinate-space perturbation theory with a single mass scale (Schröder & Meyer, Nov 2025) is a continuum method (Gegenbauer expansion), intended "at finite temperature and/or in finite volume" and for auxiliary PT in lattice-QCD hadronic precision observables. — [arXiv:2511.17349](https://arxiv.org/abs/2511.17349)
- 3D lattice perturbation theory anchor: the "Four-loop plaquette in 3d with a mass regulator" (Torrero, Di Renzo, Miccio, Laine, Schröder, 2005) appeared in an arXiv listing. *I did not capture its arXiv ID; do not cite it without checking.*

**E. Photonics / atomic arrays**
- Single-photon band structure of atomic lattices in 1D/2D/3D (Xie & Schotland, Dec 2025). The review of prior work says the lattice sums are "complex-valued, long-ranged … often only conditionally convergent" and are handled by Poisson summation with regularization cutoffs, because the on-site self-energy (divergent single-atom Lamb shift) has to be regularized. — [arXiv:2512.24596](https://arxiv.org/abs/2512.24596) (abstract verified; the method quote is from the search snippet of the PDF)
- Exceptional topology in 2D subwavelength atomic arrays under square→triangular deformation (Sep 2026), with "singular radiative dipolar couplings" near the light cone. The search snippet notes the dyadic-Green's-function sum is "computationally demanding" and done with Ewald's method. — [arXiv:2609.17313](https://arxiv.org/abs/2609.17313)
- Collective spectroscopy of a 1D atomic array (Browaeys/Ferrier-Barbut, 2024/25) — [arXiv:2412.02541](https://arxiv.org/abs/2412.02541) (listing). Finite trapping and dark states in emitter arrays — [arXiv:2502.09851](https://arxiv.org/abs/2502.09851) (listing).

### Inferences
- In materials, the demand is real and recent: ATM changes phase stability (2504.07338) and benchmark accuracy (2604.06539). But the method-level bottleneck for the *ATM energy on a Bravais lattice* is already closed by 2504.11989/2504.07338, which come from the same author as GZL. GZL can be a unifying and generalizing tool here, not a first solution.
- Madelung constants are single (shifted) sums, and 2025–26 papers already reach 13+ digits by elementary means. GZL's graph machinery adds nothing; EpsteinLib is the right tool.
- In lattice QFT, the finite-volume sums that practitioners actually use (C_j, Z_lm) are single sums, and multi-loop point-charge sums factorize (1702.01296). So there is no graph bottleneck to remove there.

### Gaps
- I could not confirm how Johnson's group (2604.06539) sums ATM in periodic layered solids (cutoff radius, precision). Layered materials have multi-atom bases anyway.
- I found no 2025–26 paper that computes **four-body dispersion** lattice sums for crystals. arXiv searches for "four-body dispersion crystal" returned nothing relevant.
- I found no 2025–26 arXiv paper that computes ATM three-body contributions to **phonon dispersions** or **elastic constants (Cauchy-relation violation)** with controlled lattice sums. The searches "Axilrod-Teller phonon" and "Cauchy relation three-body" returned nothing relevant. This is either a gap in the literature or in the search.
- The arXiv API (export.arxiv.org/api) returned HTTP 406 through the proxy, so the searches used the arxiv.org web search. Its keyword recall is limited, and some 2025–26 papers may have been missed.

---

## Q2. Which of these map onto GZL today, and which need extensions?

### Takeaway
**Maps today:** isotropic many-body dispersion energies on Bravais lattices (ATM and its strain derivatives, per-edge mixtures of multipole orders, Laplacian-type Einstein-model and phonon second-moment quantities), and dipolar or van der Waals quantum spin models whose couplings are isotropic in the lattice plane.
**Needs extensions:** anything with a multi-atom basis (hcp, NaCl, graphite, TMDs, kagome and other Archimedean lattices), full phonon dynamical matrices and dipolar tensors (anisotropic, which 2609.28282 is preparing), oscillatory Green's functions e^{ik0 r}/r (photonics), pole kernels 1/(p²−q²) and vertex-dependent weights (Lüscher-type sums), and d=4 (lattice QCD).

### Cited Findings
- GZL kernels are functions of vertex differences only, with isotropic power laws plus compact, possibly anisotropic, short-range parts. Coincident vertices are excluded (K(0)=0), the lattices are Bravais, d≤3, and "Re(ν_e) > d on every edge" is the stated sufficient convergence condition. — [GZL GitHub](https://github.com/graph-zeta/gzl)
- Anisotropic Epstein zeta functions for any lattice, exponent and anisotropy order now exist (Sep 2026) and target the dipole interaction. — [arXiv:2609.28282](https://arxiv.org/abs/2609.28282)
- The ATM-type lattice sums needed "meromorphic continuations" in 2504.07338, which means some exponents fall outside the naive convergence region. — [arXiv:2504.07338](https://arxiv.org/abs/2504.07338)
- Coupled-dipole lattice sums in atomic arrays are complex, oscillatory, conditionally convergent and self-energy-singular. They are handled with Ewald or Poisson plus regularization. — [arXiv:2512.24596](https://arxiv.org/abs/2512.24596); [arXiv:2609.17313](https://arxiv.org/abs/2609.17313)
- MBD energies are computed from the Ewald dipole tensor plus BZ sampling, i.e. Tr ln(1−AT(k)) per k. — [libMBD arXiv:2308.03140](https://arxiv.org/pdf/2308.03140)
- Multi-loop finite-volume QED sums for point charges factorize into single sums. — [arXiv:1702.01296](https://arxiv.org/html/1702.01296)

### Inferences

Bridge-by-bridge assessment. Everything here is my own analysis built on the cited facts.

1. **ATM triangle energy on Bravais lattices: maps today (strong, but already solved).** By the law of cosines, 1 + 3cosγ1cosγ2cosγ3 = 1 + 3(a²+b²−c²)(a²+c²−b²)(b²+c²−a²)/(8a²b²c²). So E_ATM/C9 = Σ [(abc)^(-3) + (3/8)·P(a²,b²,c²)/(abc)^5], where P is a degree-3 polynomial in squared edge lengths. Each monomial a^{2i}b^{2j}c^{2k} (i+j+k=3) gives edge exponents ν = 5−2i ∈ {5,3,1,−1}. ATM is therefore a **finite sum of triangle graph zeta functions with isotropic power-law edges**, and the sum runs over distinct sites exactly as in GZL.
   - **Caveat:** several terms have ν ≤ d on some edge, including ν=−1, which is a *growing* edge factor. They violate GZL's stated sufficient condition and need the meromorphic continuation used in 2504.07338. Whether GZL's `evaluate_graph` accepts ν≤d or negative ν, and continues it correctly, is **unverified**.
   - This works only because every dot product in a *triangle* is expressible through its edge lengths.
2. **Strain derivatives (pressure, bulk modulus, elastic constants C11/C12/C44, Cauchy-relation violation from 3-body forces): maps today (strong).** GZL takes an arbitrary lattice matrix A, so strain derivatives follow from derivatives with respect to A. That can be done numerically, or analytically through the Epstein-zeta representation. This works for cubic, tetragonal (Bain path) and triclinic Bravais lattices. hcp (2 atoms per cell) is excluded, which matters because fcc vs hcp is the classic rare-gas question in the Schwerdtfeger literature.
3. **Laplacian-type vibrational quantities: maps today (medium).**
   - Einstein-model zero-point energy (the approach of 2012.05413, extended to ATM) needs Σ ∇_0² V. The Laplacian of an isotropic pair term is isotropic.
   - For the three-body term, ∇_0²V_ATM expands in edge lengths again, because the gradients' dot products within a triangle are edge-expressible. So it stays a sum of triangle graph zeta functions.
   - The same holds for the phonon **second spectral moment** Σ_s ω_s²(q) = Tr D(q) with full q-dependence: the traced dynamical matrix with a Bloch phase on one edge is exactly GZL's external-momentum mode.
4. **Full phonon dispersions and dynamical matrices with 3-body forces: needs the anisotropic extension.** The Hessian ∂²V/∂u_α∂u_β contains x_αx_β/|x|^(ν+4) tensors. This is roadmap-level (built on 2609.28282).
5. **Four-body and higher dispersion (e.g. the 4-atom dipole ring Tr(T12T23T34T41)): forced or needs an extension.** For a 4-cycle, dot products of non-adjacent edges (r12·r34) involve the diagonals r13 and r24. Those enter as *positive-power polynomial kernels* on non-edges, which gives a denser graph with growing kernels, and the tensor contractions are anisotropic. The Buchheit–Busse n-body method (2504.11989) claims linear-in-n cost for "a broad class of n-body lattice sums", so it may cover some of these already.
6. **MBD (full many-body dispersion): weak bridge.** On a one-atom Bravais lattice, the ring expansion Tr[(AT)^n] is diagonal in k and reduces to a 3×3 T(k) per k-point. That needs anisotropic *single* sums (2609.28282 or Ewald), not graph zeta functions. Real MBD systems also have multi-atom cells.
7. **Madelung and Coulomb single sums: no graph structure.** EpsteinLib covers these, and GZL adds nothing.
8. **Dipolar spin models: partial.**
   - Dipoles polarized perpendicular to a 2D Bravais plane (square, triangular) have isotropic 1/r³ couplings. Rydberg van der Waals 1/r⁶ couplings are isotropic in any geometry. Both work in GZL today, but they belong to the linked-cluster field, which another note covers.
   - In-plane or tilted dipoles (2605.25019 in the listing uses a (x²−y²)/r⁵ tail) and Archimedean or Laves lattices (2609.02934, 2605.07685) need the anisotropic and multi-site extensions.
   - Beyond linear spin-wave theory, the 1/S corrections are BZ sums of products of A(k) and are graph-structured in real space. They need both extensions.
9. **Photonics (coupled dipoles, collective Lamb shifts, plasmonic lattices): does not map.**
   - The kernel is the oscillatory, complex dyadic Green's tensor e^{ik0 r}(…)/r. GZL's momentum k is a Bloch phase e^{iq·x} on a real even kernel, not an outgoing-wave kernel.
   - The single-excitation band structure is a single sum anyway, which Ewald and Poisson methods handle.
   - A quasi-static near-field expansion (1/r³ leading term) discards the radiative light-cone physics that 2609.17313 and 2512.24596 study, so it is forced.
10. **Lüscher zeta Z_lm(1;q²): does not map.** The kernel is 1/(|n+d|²−q²), a pole kernel rather than a power law. It is a single sum with vertex- and angle-dependent weights Y_lm. The modified long-range version (2507.18399) nests sums over p and p′ with a Coulomb-like 1/|p−p′|² edge, but its vertex weights 1/(p²−q²) are not difference kernels. An expansion in q² (Z_00 as a series of Epstein zeta values) converges only for |q²| below the first free level, so it is forced.
11. **QED_L power-law coefficients: Epstein level, not GZL level.** The coefficients C_j need Epstein zeta values, some at ν=j≤d, i.e. continuation. Multi-loop point-charge sums factorize (1702.01296). GZL's graph machinery is not needed, but EpsteinLib is a direct fit.
12. **Position-space photon propagator in lattice QCD+QED (2603.13086): does not map.** The volume-squared sum couples numerical hadronic correlator data at the endpoints. Those are vertex functions from Monte Carlo, not analytic kernels.
13. **Lattice perturbation theory in coordinate space: conceptual match, blocked in practice.**
    - Vacuum or n-point diagrams on the lattice are Σ over vertex positions of products of lattice propagators G(x_i−x_j). That is exactly a graph lattice sum on Z^d with per-edge kernels.
    - The Lüscher–Weisz-style split "exact G on a finite region + asymptotic expansion outside" matches GZL's "compact a(x) + power-law tails" design.
    - **Blockers:**
      - (a) Lattice QCD is d=4, but GZL supports d≤3.
      - (b) The asymptotic expansion of the lattice propagator has cubic-anisotropic subleading terms, which needs the anisotropic extension.
      - (c) Massless 3D propagators ~1/(4π|x|) have ν=1<d, which gives IR-divergent products and needs a mass regulator (see the "mass regulator" anchor). Massive propagators decay exponentially, which is neither compact nor a power law.
      - (d) Beyond one loop, lattice vertices carry derivative couplings (gauge fields), which are anisotropic.
    - **Only a narrow feasible target:** multi-loop vacuum graphs of 3D *scalar* lattice theories (e.g. 3D φ⁴ or EFTs for electroweak and first-order phase-transition studies), with the propagator's isotropic leading tail plus an exact compact core. The anisotropic subleading terms would still be required for precision.

### Gaps
- I did not verify whether GZL accepts ν ≤ d or negative ν (growing edge factors), or composite kernels with polynomial growth. This decides whether the ATM decomposition (item 1) runs in GZL as it stands.
- I found no 2025–26 paper doing coordinate-space **3D** lattice perturbation theory where multi-loop lattice sums are the bottleneck. The field's recent activity (EQCD improvement: [arXiv:1911.02351](https://arxiv.org/abs/1911.02351), listing) predates 2025.
- I did not verify whether the Buchheit–Busse n-body method (2504.11989) already handles Bloch phases or 4-body tensor terms. Only the abstract was read.

---

## Q3. What would be genuinely new results, and how strong is each bridge?

### Takeaway
The most defensible new results are in **classical crystal physics with many-body dispersion on Bravais lattices**. They are all feasible with isotropic kernels:
- machine-precision **elastic constants and Cauchy-relation violations** from ATM (plus per-edge mixtures of multipole-order three-body terms) along full Bain and triclinic deformation paths;
- **Einstein-model and phonon-second-moment vibrational corrections including three-body forces**, with full q-dependence obtained through FFT.

Photonics and Lüscher-type lattice QFT are **not** credible near-term GZL targets. QED_L is an EpsteinLib target, not a GZL one. Coordinate-space lattice PT is a long-term target that needs d=4 and anisotropy.

### Cited Findings
- The ATM energy alone already runs in minutes to machine precision for cuboidal lattices, and the fcc/bcc stability results exist. — [arXiv:2504.07338](https://arxiv.org/abs/2504.07338); [arXiv:2504.11989](https://arxiv.org/abs/2504.11989)
- Schwerdtfeger's group derived analytical Einstein-model vibrational corrections (ZPE, anharmonicity, Grüneisen) for LJ in cubic and hcp lattices at the two-body level. — [arXiv:2012.05413](https://arxiv.org/abs/2012.05413) (listing)
- The Schwerdtfeger/Cooper line cites vibrational, temperature and pressure effects, as well as >2-body forces, as the factors deciding bcc stability. — [arXiv:2504.07338](https://arxiv.org/abs/2504.07338)
- ATM demonstrably changes computed exfoliation energies of layered materials. — [arXiv:2604.06539](https://arxiv.org/abs/2604.06539)
- Anisotropic Epstein zeta functions are now available as a building block. — [arXiv:2609.28282](https://arxiv.org/abs/2609.28282)
- GZL obtains a full BZ grid at the cost of one momentum evaluation. — [arXiv:2609.18918](https://arxiv.org/abs/2609.18918)

### Inferences
Ranked candidate projects, from strongest to weakest bridge. All are inferences.
1. **(Strong, feasible now if ν≤d continuation is supported)** *Three-body elastic constants and Cauchy violations for all cubic and Bain-path Bravais lattices.* Compute C_ij, bulk modulus and pressure from LJ(n,m)+ATM (plus optional dipole-dipole-quadrupole three-body terms as extra per-edge power laws) by differentiating graph zeta functions with respect to A.
   - This is new relative to 2504.07338, which reports energies only according to its abstract.
   - It gives machine-precision C12−C44 (Cauchy violation), a classic signature of many-body forces.
   - Risk: the Buchheit–Busse framework can likely do this too, so GZL would be "the general tool" rather than the unique enabler.
2. **(Medium–strong)** *Vibrational ZPE and the traced phonon spectrum Σ_s ω_s²(q) including ATM over the full BZ.* Laplacians of isotropic pair and triangle terms stay isotropic, so this runs with isotropic kernels plus a Bloch phase on one edge. It extends the two-body Einstein-model results (2012.05413) to three-body forces and feeds straight into the fcc/bcc/cuboidal stability question.
3. **(Medium, needs the anisotropic extension)** *Full phonon dispersions with three-body forces for Bravais crystals*, combining GZL's k-FFT with anisotropic Epstein zeta edges from 2609.28282. This would be genuinely new as a machine-precision tool, but it depends on the roadmap.
4. **(Medium–weak)** *Four-body dispersion lattice energies.* Probably novel, but the graph acquires diagonal positive-power kernels and anisotropic contractions. Check first whether 2504.11989's n-body class already covers it.
5. **(Weak, Epstein level)** *QED_L and QED_r power-law coefficients, Madelung constants, micromagnetic periodic images.* These are single sums; the right tools are EpsteinLib and 2509.26274. Do not sell them as GZL results.
6. **(Weak, long-term)** *Coordinate-space multi-loop lattice PT.* The structure matches, but it needs d=4, anisotropic tails, IR regulators and derivative vertices. A restricted 3D scalar version is the only near-term test case, and I found no 2025–26 demand signal for it.
7. **(Not a bridge)** *Photonics and atomic-array collective shifts, plasmonic lattice sums, Lüscher zeta functions.* The kernels are oscillatory or pole-type, the sums are single and have non-difference vertex weights. Calling these GZL applications would be forced.

Cross-cutting limitation: many headline materials (hcp rare-gas solids, NaCl-type ionic crystals, graphite, h-BN, TMDs, kagome and Archimedean magnets) are **non-Bravais**. Multi-site (sublattice-resolved) graph zeta functions are the single most valuable extension for the materials transfer, more than anisotropy.

### Gaps
- There is no quantitative benchmark comparing GZL with the Buchheit–Busse n-body Epstein integral method on ATM (runtime, precision). It is needed to judge the added value.
- I could not confirm whether anyone in 2025–26 has published machine-precision ATM elastic constants or three-body phonon moments. Absence from my limited arXiv keyword searches is weak evidence of novelty.
- I did not check journal-only literature (J. Chem. Phys., Phys. Rev. B) for rare-gas-crystal three-body lattice dynamics.
