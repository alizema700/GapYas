# Novelty check, author by author: three-body product energy T_ν, max-min triple product P(L), 3D BCC

Date: 2026-09-26. This note complements `07_novelty_math.md`. That note was organised by topic; this one goes **author by author**, as requested. It also adds a separate check for the max-min "triple product" problem (R2).

**Labels.** The labels below follow the current request. They differ from the labels in 07.
- **R1** — 2D, unit covolume. The energy is T_ν(L) = Σ_{x≠y∈L\{0}} (|x||y||x−y|)^{−ν}. Numerically, the minimiser is hexagonal for ν < 3.9184, square for 3.9184 < ν < 8.606, and rectangular for larger ν. Other names for the same object:
  - the radial part of Axilrod-Teller-Muto (ATM) at ν=3;
  - π^{3s} C_{s,s,s}(τ) with s = ν/2;
  - Buchheit-Busse ζ^{(3)}_L(ν,ν,ν).
- **R2** — Theorem. Define P(L) = min over triples of distinct lattice points of the product of the three mutual distances. Then P(L) ≤ P* = 2((1+√17)/8)^{3/4} ≈ 1.4313. Equality holds only for the rectangle with aspect ratio √((√17−1)/2) ≈ 1.2496. As a consequence, minimisers of T_ν converge to that rectangle as ν→∞. Also covered here: any max-min problem for the smallest lattice triangle (by area, perimeter or product).
  - Sanity check (mine). For a rectangle a×(1/a), the optimum equalises the collinear triple 2a³ with the right triangle √(a²+a⁻²). This gives a⁴ = (1+√17)/8, which reproduces P* and the aspect ratio. For comparison, hex has P = (2/√3)^{3/2} ≈ 1.2408.
- **R3** — 3D. Numerically, BCC minimises T_ν among Bravais lattices for ν ∈ [3.2, 12]. Conjecture: BCC maximises P, with P = 3/2.

**Method.**
- Author listings came from the arxiv.org/search HTML (author field, up to 200 results, abstracts shown).
- I parsed every listing and read the abstracts after a keyword filter (lattice, crystal, three-body, many-body, theta, Epstein, energy, minimiser, triangle, packing, modular, Axilrod). I also skimmed the full title lists by eye.
- Full texts were grepped (arXiv HTML or pdftotext) for 2504.07338, 2504.11989 (v1 and v2), and 1504.01153.

**Tools that did not work.** arxiv.org/a/<id> pages returned "Not Found" for the ids I guessed. The export API returned 406. Semantic Scholar returned 429 on every attempt. The arXiv listing search and WebSearch worked.

All arXiv IDs below were seen in arXiv search output or on abs pages on 2026-09-26. Journal-only items are marked *(no arXiv ID verified)*.

---

## 1. Laurent Bétermin
The arXiv listing shows 31 papers under both "Betermin" and "Bétermin".

| ID | Title (year) | Relevance |
|---|---|---|
| 1402.2751 | Minimization of energy per particle among Bravais lattices in R²: Lennard-Jones and Thomas-Fermi cases (2014) | Two-body only |
| 1502.03839 | 2D Theta Functions and Crystallization among Bravais Lattices (2015) | Two-body; the framework for T_ν-type problems in 2D |
| 1505.08047 | Sufficient Condition for a Compact Local Minimality of a Lattice (2015) | Two-body |
| 1607.08716 | Dimension reduction techniques for the minimization of theta functions on lattices (w/ Petrache, 2017) | Two-body |
| 1611.07798 | Local optimality of cubic lattices for interaction energies (2017) | Shows that SC, FCC and BCC are critical points by symmetry. This applies to any symmetric energy, including T_ν, so criticality of BCC is not news |
| 1611.07820 | Local variational study of 2d lattice energies and application to Lennard-Jones type interactions (2018) | Two-body. Says square and hex are always critical points, which again holds for T_ν |
| 1704.02887 | On Born's conjecture about optimal distribution of charges for an infinite ionic crystal (w/ Knüpfer, 2018) | Two-body |
| 1710.05581 | Optimal lattice configurations for interacting spatially extended particles (w/ Knüpfer, 2018) | Two-body |
| 1806.02233 | Optimal and non-optimal lattices for non-completely monotone interaction potentials (w/ Petrache, 2019) | Two-body; square can win for designed potentials |
| 1809.00473 | Minimal Soft Lattice Theta Functions (2019) | Two-body |
| 1901.08957 | Minimizing lattice structures for Morse potential energy in two and three dimensions (2019) | Two-body; 2D and 3D, including BCC and FCC |
| 1907.06105 | Crystallization to the square lattice for a two-body potential (w/ De Luca, Petrache, 2019) | Two-body; square optimum |
| 1908.01515 | On a lattice generalisation of the logarithm and a deformation of the Dedekind eta function (2020) | Two-body |
| 2004.04553 | On the optimality of the rock-salt structure among lattices with charge distributions (w/ Faulhuber, Knüpfer, 2020) | Two-body |
| 2007.15977 | Maximal Theta Functions - Universal Optimality of the Hexagonal Lattice for Madelung-Like Lattice Energies (w/ Faulhuber, 2023) | Two-body |
| 2008.00676 | Effect of periodic arrays of defects on lattice energy minimizers (2020) | Two-body |
| 2009.11503 | Stability of Z² configurations in 3D (w/ Friedrich, Stefanelli, 2021) | Two- plus three-body **angular**, finite configurations, angle-rigidity. Not a lattice sum |
| 2010.08264 | Theta functions and optimal lattices for a grid cells model (2021) | Two-body |
| 2101.05602 | Lattice ground states for Embedded-Atom Models in 2D and 3D (w/ Friedrich, Stefanelli, 2021) | **Many-body (EAM)**, but a function of pair sums. Triangular, square or orthorhombic optimum in 2D; FCC, BCC or SC in 3D. Different energy |
| 2103.10286 | On energy ground states among crystal lattice structures with prescribed bonds (2021) | Two-body; BCC and FCC among lattices with prescribed bonds |
| 2104.09795 | Optimality of the triangular lattice for Lennard-Jones type lattice energies: a computer-assisted method (2023) | Two-body |
| 2105.07922 | How well-conditioned can the eigenvalue problem be? (2021) | Not relevant |
| 2107.14020 | Three-dimensional lattice ground states for Riesz and Lennard-Jones type energies (w/ Šamaj, Travěnec, 2022) | Two-body Riesz in 3D. FCC wins for large s via a packing argument, the pair analogue of R2's ν→∞ mechanism. BCC wins only for small s |
| 2110.06008 | A variational principle for Gaussian lattice sums (w/ Faulhuber, Steinerberger, 2021) | Two-body |
| 2110.09368, 2307.06002 | Zeros of 2D Epstein zeta and lattice energies (w/ Šamaj, Travěnec) | Two-body |
| 2312.01395 | Structural transitions in interacting lattice systems (w/ Šamaj, Travěnec, 2023/25) | Two-body square→rectangle transitions at fixed density. The analogue of R1's square→rectangle step, but for pair potentials |
| 2407.20762 | On crystallization in the plane for pair potentials with an arbitrary norm (w/ Furlanetto, 2024/26) | Two-body |
| 2502.16639 | Equidistant versus bipartite ground states for 1D classical fluids (w/ Šamaj, Travěnec, 2025) | 1D, two-body |
| 1404.4485, 1804.05743 | Sphere log-energy; alternating chains | Not relevant |

**Verdict.** No paper with a triangle/product three-body lattice sum and no max-min triple problem, so R1, R2 and R3 are **not contained**. The closest are 2101.05602 (EAM many-body), 2312.01395 (pair square→rectangle) and 2107.14020 (large-exponent packing argument, the method template for R2).

## 2. Mircea Petrache
The listing shows 47 papers, most of them ML or geometric analysis.

Relevant:
- 1607.08716, 1806.02233, 1907.06105 (with Bétermin; above).
- 1908.09714, Crystallization for Coulomb and Riesz Interactions as a Consequence of the Cohn-Kumar Conjecture (w/ Serfaty).
- 1409.7534, 1609.03849, 1707.07664, 1706.06008 — Riesz/jellium next-order asymptotics.
- 2408.02136 — discrete energies with topological singularities.
- 2101.11977 — Wulff shapes in lattices.

**Verdict.** All of these are two-body (or OT/jellium). None of R1–R3 appears.

## 3. Markus Faulhuber
The listing shows 29 papers.

Relevant:
- 2007.15977, 2004.04553, 2110.06008 (with Bétermin; above).
- 1903.06856, An Extremal Property of the Hexagonal Lattice (w/ Steinerberger, 2019). Distances to a deep hole; one-point sums.
- 2306.16266, A note on energy minimization in dimension 2 (w/ Shafkulovska, Zlotnikov, 2023). Hex vs periodic rivals and the honeycomb; pair energies.
- 2509.24687, The polarization problem for the honeycomb structure (2025). Polarization; pair/one-point.
- 1608.01168, 1601.02972, 1709.06006, 1901.01218. Theta functions, rectangular tori, heat kernels.
- Gabor-frame papers: 2204.02917, 2403.10503, 2408.08975, 2502.09510, 2609.01296, 2609.02610.

**Verdict.** There are no three-body or multi-point lattice energies. **Not contained.**

## 4. Ulisse Stefanelli, Manuel Friedrich, Leonard Kreutz, Edoardo Mainini, Paolo Piovano, Bernd Schmidt
The listings show 87, 76, 25, 32, 16 and 76 papers respectively.

Relevant (crystallisation / three-body angular):
- 2009.11503 and 2101.05602 (above).
- 1807.00811, N^{3/4} law in the cubic lattice (Mainini, Piovano, Schmidt, Stefanelli).
- 1604.02077, The geometry of C_60 (Friedrich, Piovano, Stefanelli).
- 1802.05049, Graphene ground states (Friedrich, Stefanelli), and 1802.05053, Ripples in graphene. Two- plus three-body angular potentials.
- 1706.01494, 1909.12023 — carbon nanotubes / Cauchy-Born (Friedrich, Mainini, Piovano). Angular three-body.
- 1903.00331 and 1808.10675 — ionic compounds, square and hex finite crystallisation (Friedrich, Kreutz).
- 2209.14880 — finite crystallisation via stratification (Kreutz, Friedrich).
- 2509.05642 — Winterbottom shape (Friedrich, Kreutz, Stefanelli, 2025).
- 2006.01558 and 2604.19239 — polycrystals (Friedrich, Kreutz, Schmidt; Kreutz, Ziereis 2026).
- 2609.21508 — charge-dependent hard spheres (Kreutz, Ziereis, 2026).
- 1703.01981 — continuum limits of multi-body lattice energies (Braides, Kreutz).
- 2402.03281 — Winterbottom (Kreutz, Schmidt).
- 2003.01679 — Wulff fluctuations (Mainini, Schmidt).
- 1302.6513 — sticky disc N^{3/4} (Schmidt).
- Journal-only, cited in Blanc-Lewin *(no arXiv ID verified)*: Mainini, Piovano, Stefanelli, Finite crystallization in the square lattice, Nonlinearity 27 (2014); Mainini, Stefanelli, Crystallization in carbon nanostructures, CMP 328 (2014).

**Verdict.** These papers use short-range **angular** three-body terms that favour prescribed bond angles, or sticky-disc models, applied to finite-N crystallisation. None is a power-law product sum over a lattice, and none is a max-min triple problem. **Not contained.**

## 5. Senping Luo, Juncheng Wei
The listing shows 17 papers for Luo. I read the lattice-related papers by Wei.

| ID | Title | Relevance |
|---|---|---|
| 2004.13882 | On minima of sum of theta functions and Mueller-Ho Conjecture (Luo, Wei) | Two-species pair energy with a hex→rhombic→square→rectangular sequence. The phenomenological analogue of R1, but a different energy and a different mechanism |
| 2110.08728 | On universally optimal lattice phase transitions… (Luo, Wei, Zou) | Same as above; square is optimal on an interval |
| 1902.09611 | Non-hexagonal lattices from a two species interacting system (Luo, Ren, Wei) | Same |
| 2203.00264, 2212.10727, 2302.05042 | Differences of theta / Epstein; LJ exact; non-monotone | Two-body |
| 2312.02497 | Variational model with hexagonal to square phase transitions (Luo, Wei) | Engineered modular-invariant function with a direct hex→square jump (no rhombic phase), like R1's first-order jump. Not a lattice sum |
| 2411.17199, 2412.09201, 2501.01265 | Deng and Luo | Two-body (Gaussian-weighted), including a "hexagonal to skinny-rhombic" minimiser |
| 2605.07580 | On ratios of theta functions (Luo, Wei, 2026) | Two-body ratios |
| 2609.17356 | On Sarnak-Strömbergsson conjecture (Luo, Wei, Sep 2026) | 3D theta: FCC for α>1, BCC for α<1. Two-body. Contrast for R3 |
| 2208.00528 | The BCC lattice in a long range interaction system (Ren, Wei) | BCC for a diblock-copolymer-type functional. Not a three-body lattice sum |

**Verdict.** **Not contained.** The phenomenology is related (hex→square→rectangle sequences; BCC in 3D), but every energy is two-body or engineered.

## 6. Henry Cohn
The listing shows 64 papers.

Relevant:
- 1103.0485, Three-point bounds for energy minimization (w/ Woo). Three-point SDP bounds for **pair** energies in RP².
- 2206.15373, Three-point bounds for sphere packing (w/ de Laat, Salmon).
- 1902.05438, universal optimality of E8/Leech.
- 1603.09684, 0911.2169, 0811.1236 — Gaussian core / ground states.
- 1003.3053, 1603.05202 — surveys.
- 2609.03121 — stealthy configurations (2026).

**Verdict.** "Three-point" here means SDP bounds on two-point problems. **Not contained.**

## 7. Dmitriy Bilyk
The listing shows 34 papers.

Relevant:
- 2303.12283, Optimizers of three-point energies and nearly orthogonal sets (w/ Ferizović, Glazyrin, Matzke, Park, Vlasiuk, 2023). Three-point energies **and a three-point packing (max-min) problem on the sphere**. This is the closest conceptual relative of R2, but it is on S^{d−1}, with a different three-point function, and not about lattices.
- 2303.14258, Optimal Measures for Multivariate Geometric Potentials (2023). k-point potentials such as simplex volume; optimal measures on the sphere.
- 2510.25442, tensor-product energies on the torus (2025). Pair energies, permutation sets.
- 2604.14326, 2409.16508, 2302.13067, 1908.10354 — Riesz, polarization and sphere energies.

**Verdict.** k-point energies appear only on spheres and for measures, never for Bravais lattices. **Not contained.** 2303.12283 should be cited as the conceptual precursor for "three-point packing".

## 8. Blanc & Lewin
- 1504.01153, The Crystallization Conjecture: A Review (2015).
- I grepped the full text. Three-body terms are mentioned only as angular terms in finite crystallisation: Flatley-Theil 1407.0692 (FCC, 3D), and [78, 169, 170], i.e. E-Li, and Mainini-Piovano-Stefanelli / Mainini-Stefanelli.
- There is nothing on lattice sums of three-body power laws and no max-min triangle problem.

**Verdict.** **Not contained.** Lewin's other papers (104 listed) are quantum/Coulomb and not relevant.

## 9. Schwerdtfeger / Burrows / Cooper / Buchheit / Busse group
The listings show 24 papers for Schwerdtfeger, 5 for Burrows, 14 for Buchheit and 5 for Busse.

| ID | Title | Relevance |
|---|---|---|
| **2504.07338** | Exact lattice summations for LJ potentials coupled to a three-body ATM term applied to cuboidal phase transitions (Robles-Navarro, Cooper, Buchheit, Busse, Burrows, Smits, Schwerdtfeger; J. Chem. Phys. 163, 094104, 2025) | **Closest.** The full text defines f_r^{(3)} = (1/6)Σ'(|x||y||x−y|)^{−3}, which is exactly T_3/6. See the notes below the table |
| **2504.11989** | Epstein zeta method for many-body lattice sums (Buchheit, Busse; v2 June 2026, Numer. Math.) | Defines ζ^{(n)}_Λ(ν), so T_ν = ζ^{(3)}(ν,ν,ν). Benchmarks on Z², hex, FCC and BCC; n-body sums up to n=51 for sq/hex at ν=3; LJ+ATM fcc→bcc transition at finite pressure. **No minimisation over lattice space.** It says "ongoing investigation into the influence of many-body interactions on the stability of matter", which is the **scooping risk** |
| 2609.18918 | Graph lattice sums and graph zeta functions… (Buchheit, Rupp, Sep 2026) | Computational: graph zeta functions for linked-cluster expansions. No optimisation |
| 2609.28282, 2412.16317, 2509.26274, 2403.03213, 2102.10941, 2201.11101 | Epstein / anisotropic / singular sums | Computational, two-body |
| 2501.05746 | A Minimum Property for Cuboidal Lattice Sums (Cooper, Schwerdtfeger, 2025) | **Two-body** Epstein zeta has a local minimum at BCC along the Conway-Sloane cuboidal family. A pair-energy analogue of R3 along a one-parameter path |
| 2105.08922, 2107.11380, 2406.09635, 2502.09828, 2012.05413 | Cuboidal lattice sums; BCC instability (sticky-hard-sphere / LJ); hcp↔cuboidal paths | Two-body |

Notes on 2504.07338:
- 2D (Sec. III.2): only square and hex are evaluated, at **fixed nearest-neighbour distance**, with f_r(sq) = 2.27548 and f_r(hex) = 4.26383. They conclude hex is more destabilised. They also study a rectangular distortion γ for LJ+ATM, where the lattice collapses into chains.
- 3D: the Bain path fcc–mcc–bcc–acc. They show analytically that the ATM energy has an extremum at BCC, numerically a minimum along the Bain path. Their Appendix F (Theorem F.1) is for a **two-body** sum L(A;s) along that path.
- Everything is at ν=3, uses the full ATM (angular) energy, and is normalised by nearest-neighbour distance.
- Converting their f_r values to unit covolume gives T_3(hex) = 13.392 < T_3(sq) = 13.653. This is consistent with R1 at ν=3, but the paper does not state it.

**Verdict.**
- R1 is not contained. They have no ν-scan, no fixed-density comparison and no fundamental-domain search.
- R2 is not contained.
- R3 is **partially anticipated** in a weaker and different form: BCC is favoured and extremal along the Bain path for three-body ATM at ν=3, NN-normalised. There is no global Bravais search, no pure radial energy, no ν range, and no P-maximisation.

## 10. Keyword searches (discrete geometry / modular graph functions)

**Discrete geometry.** I searched arXiv abstracts and the web for: "product of distances" + lattice; "triangle" + lattice + max-min / maximal minimal area / perimeter; "Heilbronn" + lattice; "three-body" + Bravais / crystallization / lattice ground state; "Axilrod-Teller" + lattice; "many-body lattice energy minimizer".
- Nothing matched a max-min triple-product (or triangle-perimeter) problem over unit-covolume lattices.
- The area version is trivial: every non-degenerate lattice triangle has area ≥ covol/2, and collinear triples give area 0. This is why the product (or perimeter) is the non-trivial functional.
- Only the classical Pick / minimal-area facts came up, plus unrelated items (Laugesen et al. 1609.06172 and 1701.03217 on optimal stretching for lattice points; Leblé 2511.03353 on local universal optimality of hex for pair potentials; Zhou 2609.05622 on distance energies on flat tori).

**Modular graph functions.** A quoted search on "modular graph functions" (41 abstracts) and web searches found no paper on the minimum or extremum of C_{a,b,c} over the fundamental domain. This agrees with the full-text checks in 07.

---

## Overall verdicts

| Result | Verdict | Confidence |
|---|---|---|
| **R1** (T_ν over 2D unit-covolume lattices: hex below 3.9184, square up to 8.606, then rectangular) | **Not found in any listed author's work.** The object T_ν = ζ^{(3)}(ν,ν,ν) = π^{3s}C_{s,s,s} is known and computed (Buchheit-Busse 2504.11989; Robles-Navarro et al. 2504.07338 at ν=3, sq/hex only, NN-normalised). The phase sequence resembles two-body results (Luo-Wei 2004.13882 / 2110.08728, Luo-Wei 2312.02497, Bétermin-Šamaj-Travěnec 2312.01395), but for different energies | High, about 85% novel |
| **R2** (P(L) ≤ 2((1+√17)/8)^{3/4}, uniquely attained by the rectangle with aspect √((√17−1)/2); ν→∞ limit of T_ν minimisers) | **Not found.** No max-min "smallest triangle / triple-product" lattice problem turned up anywhere. The large-exponent reduction to a max-min problem is a standard technique (e.g. FCC for large s in Bétermin-Šamaj-Travěnec 2107.14020), so the method is known; the specific extremal problem and its solution appear new. Closest conceptual relative: three-point packing on the sphere (Bilyk et al. 2303.12283) | High, about 85-90% novel. The elementary nature leaves some risk that it sits in an old geometry-of-numbers note or problem column |
| **R3** (3D: BCC minimises T_ν for ν ∈ [3.2,12]; BCC maximises P = 3/2) | **Partially anticipated.** 2504.07338 shows that repulsive three-body ATM energy is extremal and minimal at BCC along the Bain path, with BCC favoured over FCC (ν=3, NN-normalised, angular ATM). Cooper-Schwerdtfeger 2501.05746 gives the two-body analogue. BCC criticality follows from symmetry (Bétermin 1611.07798). The global Bravais-lattice statement for the pure product energy over a ν range, and the P-maximisation conjecture, were not found | Medium: about 65% novel for the global statement and the P conjecture; the Bain-path aspect is anticipated |

**Scooping risk.** The Buchheit / Busse / Schwerdtfeger group is the main risk. They state an ongoing programme on "many-body interactions and stability of matter" and have the tools (ζ^{(n)}). Their most recent items (2609.18918, 2609.28282, Sep 2026) are still purely computational.

**Could not open or check:**
- arxiv.org/a/<id> author pages (the ids I tried were not found; the listing search was used instead).
- Semantic Scholar (HTTP 429).
- Journal-only papers by Mainini-Piovano-Stefanelli (2014), Mainini-Stefanelli (2014) and E-Li (2009). I know these only through the Blanc-Lewin reference list.
- Non-arXiv literature, in particular older geometry-of-numbers or problem-column sources on minimal triangle products. This could only be probed by web search.
