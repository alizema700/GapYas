# Novelty check: lattice minimisation of the three-body product ("triangle") energy T_ν and the MGF C_{s,s,s}

Date: 2026-09-25.

**How IDs were checked.** Every arXiv ID marked [v] was checked on 2026-09-25. For each one I fetched the arxiv.org/abs page (title, authors, date) or ran the arXiv listing search. Full texts were grepped (HTML or pdftotext) for the following papers: 2504.11989, 2504.07338, 2609.18918, 1502.06698, 1512.06779, 1608.04393, 1708.07998, 1902.04180, 1905.06217, 1706.01889, 2104.09916, 1603.00839, 1606.07084, 2208.07242, 1504.01153, 2101.05602, 2107.14020, 2110.08728. Items without [v] were not checked.

**Tools that did not work.**
- The export.arxiv.org API returned HTTP 406.
- Semantic Scholar returned 429 (rate limit).
- OpenAlex returned "daily budget exhausted".
- Google Scholar returned 403.

**Tools used instead.** The arxiv.org/search HTML (author, title and abstract fields), WebSearch, and full-text greps.

**Object under test.**
- Definition: T_ν(Λ) = Σ'_{x≠y} |x|^{-ν}|y|^{-ν}|x−y|^{-ν}.
- In Buchheit–Busse notation this is exactly the "three-body zeta function" ζ^{(3)}_Λ(ν,ν,ν).
- In Robles-Navarro et al. notation it is 6·f_r^{(3)} at ν=3, the radial part of the ATM sum.
- In 2D, T_{2s}(Λ_τ) = π^{3s} C_{s,s,s}(τ).

---

## 1. Closest prior work (directly relevant)

### 1a. Robles-Navarro, Cooper, Buchheit, Busse, Burrows, Smits, Schwerdtfeger
"Exact lattice summations for Lennard-Jones potentials coupled to a three-body Axilrod-Teller-Muto term applied to cuboidal phase transitions". arXiv:2504.07338 [v], April 2025; published in J. Chem. Phys. 163(9), 094104 (2025).

**What it does.**
- Treats the ATM three-body energy f_coh = f_r + f_a, with f_r = (1/6) ζ^{(3)}(3,3,3). This is our T_3/6. The angular part is f_a.
- 2D, Sec. III.2:
  - It evaluates the square and hexagonal lattices only, at **fixed nearest-neighbour distance R=1, not unit density**.
  - Values: f_r(sq) = 2.27548228589309, f_r(hex) = 4.263827935989311.
  - It concludes the hexagonal lattice is "more strongly destabilized" by three-body forces. That is an artefact of normalising by nearest-neighbour distance.
  - It studies a rectangular distortion γ, but only for LJ + ATM. There the ATM term turns attractive for chain-like lattices, so the lattice collapses into chains.
  - No optimisation over the fundamental domain, no rhombic lattices, no pure product energy, no scan in ν.
- 3D:
  - Studies the Bain path fcc ↔ mcc ↔ bcc ↔ acc.
  - Shows analytically that the ATM energy has an extremum at bcc; numerically it is a minimum along the Bain path for repulsive ATM.
  - Strong ATM coupling destabilises fcc in favour of bcc.
  - Table 4 lists f_r(A) along the path, again at nearest-neighbour distance 1.
  - All of this is at ν=3 and for the ATM (angular) potential. There is no global search over 3D lattices and no Hessian or saddle analysis of FCC.

**Consistency check (my conversion of their published numbers to unit covolume; not stated in the paper).**
- 2D, ν=3: T_3(hex) = 6·4.263828·(√3/2)^{4.5} = **13.39179** and T_3(Z²) = 6·2.275482 = **13.65289**. The ratio is 0.981, so hex is lower. This agrees with R1 (hex optimal for ν<ν*≈3.918).
- 3D, ν=3: T_3(BCC) ≈ **66.714** and T_3(FCC) ≈ **66.794**, so BCC is lower. This agrees with our 3D claim.
- Caveat for the 3D conversion: I assumed the Bain-path values in Table 4 are at nearest-neighbour distance 1, with V_bcc = (2/√3)³/2 and V_fcc = 1/√2.
- These are useful external reference values for testing our code at ν=3.

**Verdict.** Same building block (the ζ^{(3)}(3,3,3) product sum), but a different question: fixed nearest-neighbour distance, ATM angular factor, a handful of lattices, only ν=3.
- It does **not** contain R1 or R2.
- It partially anticipates the 3D statement "three-body repulsion favours BCC over FCC; BCC is a minimum along the Bain path", but for the ATM energy and at ν=3 only.

### 1b. Buchheit & Busse, "Epstein zeta method for many-body lattice sums"
arXiv:2504.11989 [v], April 2025, v2 June 2026; published in Numer. Math. (2026), doi 10.1007/s00211-026-01558-y.

**What it does.**
- Defines the n-body zeta function ζ^{(n)}_Λ(ν). This contains T_ν = ζ^{(3)}(ν,ν,ν) exactly.
- Computes it for Z², the hexagonal lattice, FCC and BCC as benchmarks.
- Fig. 9 shows ζ^{(n)}(3,…,3)/Z^n against n for square vs hex.
- Uses a skewed family A_κ only for conditioning tests.
- Sec. 5 gives the LJ + ATM fcc/bcc enthalpy phase diagram.

**What it does not do.** No minimisation over lattice space, no fundamental-domain study, no ν-dependent phase transition of the pure product energy.

**Competition risk.** Both this paper and the GZL paper (arXiv:2609.18918 [v], Buchheit & Rupp) state that there is "an ongoing investigation into the influence of many-body interactions on the stability of crystal lattices" (Buchheit/Schwerdtfeger group). This group has the tools to find R1–R3 quickly. **This is the main scooping risk.**

### 1c. Bétermin, Friedrich, Stefanelli, "Lattice ground states for Embedded-Atom Models in 2D and 3D"
arXiv:2101.05602 [v], 2021; Lett. Math. Phys. (2021).

**What it does.**
- A many-body energy, but of EAM type: E = Σ pair + F(Σ ρ(|x|)). This is a function of **pair** sums, not a triangle sum.
- 2D numerics: triangular, square and "orthorhombic" lattices are each optimal for different parameters. The square lattice is a local minimiser at some parameters, with first-order-like switching.
- 3D: FCC, BCC or SC is optimal depending on parameters.

**Verdict.** Related phenomenology (a many-body term can favour square or BCC), but a different energy. No triangle or product sums.

### 1d. Finite crystallisation with angular three-body terms
These use short-range two-body plus three-body **angular** potentials that favour particular bond angles. They are finite-N crystallisation proofs, not lattice sums of a power-law product.
- Mainini, Piovano, Stefanelli, "Finite crystallization in the square lattice", Nonlinearity 27 (2014). The arXiv ID was not verified.
- Mainini & Stefanelli, "Crystallization in carbon nanostructures", CMP 2014 (hexagonal/graphene). The arXiv ID was not verified.
- E & Li, "On the crystallization of 2D hexagonal lattices", CMP 2009. The arXiv ID was not verified.
- Flatley & Theil, "Face-centered cubic crystallization of atomistic configurations", arXiv:1407.0692 [v].
- Friedrich & Kreutz, "Finite crystallization and Wulff shape emergence for ionic compounds in the square lattice", arXiv:1903.00331 [v].
- Friedrich & Kreutz, "A proof of finite crystallization via stratification", arXiv:2209.14880 [v].
- Friedrich, Kreutz & Stefanelli, "Crystallization in the Winterbottom shape and sharp fluctuation laws", arXiv:2509.05642 [v].
- Mainini, Piovano, Schmidt & Stefanelli, "N^{3/4} law in the cubic lattice", arXiv:1807.00811 [v].
- The review by Blanc & Lewin, "The Crystallization Conjecture: A Review", arXiv:1504.01153 [v], mentions three-body terms only in this angular, finite-crystallisation sense.

**Verdict.** Different energy: angle-based, short range, the square lattice is built in by design. No overlap with R1–R3.

### 1e. Bilyk, Ferizović, Glazyrin, Matzke, Park, Vlasiuk
- "Optimizers of three-point energies and nearly orthogonal sets", arXiv:2303.12283 [v].
- "Optimal Measures for Multivariate Geometric Potentials", arXiv:2303.14258 [v].

These treat three-point and k-point energies **on spheres** (volume of a simplex, three-point p-frame energies). They are not lattice sums.

**Verdict.** Conceptually related ("multi-input Riesz"), but a different setting.

## 2. Pair-energy lattice results with similar phenomenology (for context and contrast)

- **Rankin 1953, Cassels 1959, Ennola 1964, Diananda 1964.** The hexagonal lattice uniquely minimises the Epstein zeta E_s for all s>0. This gives R3 at s=1, since C_{1,1,1} = E_3 + ζ(3) (Zagier; D'Hoker–Green–Vanhove, arXiv:1502.06698 [v]). This part of R3 is known and trivial.
- **Montgomery 1988** (theta functions) and **Cohn–Kumar–Miller–Radchenko–Viazovska**, arXiv:1902.05438 [v]: universal optimality for pair energies (E8, Leech). In 2D the hexagonal lattice is only conjecturally universally optimal. For completely monotone **pair** potentials the square lattice is always a saddle. For T_ν we find the square lattice becomes a **local minimum** for ν>3.8063 and the **global minimum** for ν*<ν<ν2. This is qualitatively impossible for any completely monotone pair energy. It is a genuine "three-body" effect (T_ν is not a function of pair sums).
- **Luo & Wei**, "On minima of sum of theta functions and Mueller–Ho conjecture", arXiv:2004.13882 [v]. **Luo, Wei & Zou**, arXiv:2110.08728 [v]. **Luo, Ren & Wei**, arXiv:1902.09611 [v].
  - For two-species or two-term **pair** energies, the minimiser moves hexagonal → rhombic → square → rectangular as a parameter varies. The square lattice is optimal on a closed interval.
  - This is the closest analogue to the R1+R2 sequence (hex → square → rectangular). It is a different energy and mechanism (competition of intertwined lattices).
  - Their sequence passes through **rhombic** lattices. We see a first-order jump hex → square along |τ|=1, with no rhombic phase.
- **Luo & Wei**, "On a variational model for the continuous mechanics exhibiting hexagonal to square phase transitions", arXiv:2312.02497 [v]. Constructs modular-invariant Landau-type functions with a direct hex–square transition that skips rhombic lattices. It is an engineered energy, not a lattice sum. It is the closest analogue to our first-order hex → square jump.
- **Bétermin, Šamaj & Travěnec**, "Structural transitions in interacting lattice systems", arXiv:2312.01395 [v].
  - Square → rectangular transitions for pair potentials (double Yukawa, Yukawa–Coulomb) at fixed density, driven by density.
  - Covers first- and second-order transitions and tricritical points.
  - This is the analogue of R2 (a continuous square → rectangle symmetry breaking), but for pair potentials, restricted to rectangular lattices, with density as the parameter.
- **Bétermin & Petrache**, arXiv:1806.02233 [v], and **Bétermin, De Luca & Petrache**, arXiv:1907.06105 [v]. Square lattice beats triangular for designed one-well pair potentials.
- **Bétermin, Šamaj & Travěnec**, "Three-dimensional lattice ground states for Riesz and Lennard-Jones type energies", arXiv:2107.14020 [v]. Also **Luo & Wei**, "On Sarnak–Strömbergsson conjecture", arXiv:2609.17356 [v], Sep 2026.
  - For **pair** Riesz/Epstein energies, FCC is optimal for s>3/2 and BCC for small exponents.
  - Our result is BCC for all ν ∈ [3.2, 12] for the three-body product energy, with FCC a saddle for ν ≥ 4.5. That is the **opposite** of the pair-Riesz behaviour in the same exponent range. This contrast is new, as far as I found.
  - Luo & Wei, "On ratios of theta functions", arXiv:2605.07580 [v], concerns theta/Epstein ratios. Hex plays the pivotal role. No three-body energies.
- **Cooper & Schwerdtfeger**, "A minimum property for cuboidal lattice sums", arXiv:2501.05746 [v]. Along the Conway–Sloane cuboidal family, the **pair** Epstein zeta has a local minimum at BCC. This is the pair analogue of "BCC is the minimum on the Bain path" and is relevant background for our Bain-path statement.
- **Batle**, arXiv:2609.02936 [v]. Dipolar (anisotropic pair) order over all 14 Bravais lattices. Not relevant beyond the "global search over Bravais space" method.

## 3. Modular graph functions: any extremal results?

I searched for "minimum", "extremum", "positivity", "lower bound", "special values at τ=i or ρ" of C_{a,b,c}, using full-text greps of the main MGF papers and web and arXiv searches. The papers checked:
- D'Hoker–Green–Vanhove, arXiv:1502.06698 [v]
- D'Hoker–Green–Gürdoğan–Vanhove, arXiv:1512.06779 [v]
- D'Hoker–Green, arXiv:1603.00839 [v]
- Basu, arXiv:1606.07084 [v]
- D'Hoker–Kaidi, arXiv:1608.04393 [v]
- Kleinschmidt–Verschinin, arXiv:1706.01889 [v]
- D'Hoker–Duke, arXiv:1708.07998 [v]
- D'Hoker–Kaidi, arXiv:1902.04180 [v]
- D'Hoker, arXiv:1905.06217 [v]
- Laplace-eigenvalue paper, arXiv:2104.09916 [v]
- D'Hoker–Kaidi lectures, arXiv:2208.07242 [v]

What the literature covers:
- Laplace equations, Fourier and Poincaré series, cusp asymptotics, and integrals over the fundamental domain.
- DGV note that the C_{a,b,c} formulas hold for complex a, b, c where convergent. So C_{s,s,s} with real s is a legitimate object in their framework.
- D'Hoker–Kaidi note that the summands of C_{a1,a2,a3} are positive for real a_i.

What it does not contain:
- **No paper found that locates the minimum of any C_{a,b,c} (other than those reducible to Eisenstein series) on the fundamental domain.**
- Nothing on its dependence on s, or on symmetry-breaking minima.

**Consequences of R3 at integer s.** These follow from s* = 1.959 and s2 = 4.303, and would be new statements about genuine (non-Eisenstein-reducible) MGFs:
- C_{2,2,2}, C_{3,3,3} and C_{4,4,4} have their global minimum at τ=i, not at ρ.
- C_{s,s,s} for s ≥ 5 is minimised at a non-elliptic point τ = iy with y>1.
- Worth checking independently, for example with the Fourier expansions of D'Hoker–Duke, arXiv:1708.07998.
- C_{2,1,1} = (2/5)E_4 + ζ(5)/30 (DGV) is still minimised at ρ. So "mixed" C_{a,b,c} can behave differently from C_{s,s,s}, which is an interesting side remark.

## 4. Tornheim, double-Epstein and multiple lattice sums

The Mordell–Tornheim(–Witten) literature is 1D (Borwein, Bailey et al.; e.g. arXiv:1205.0037, arXiv:2501.01380 — IDs from search results, not individually verified). It is about evaluation and MZV relations, not lattice optimisation.

I found no paper minimising a "double Epstein" or "multiple lattice sum" over the space of lattices.

## 5. Side remark on the convergence domain

Absolute convergence of ζ^{(3)} requires ν_i+ν_j > d and Σν_i > 2d (Robles-Navarro et al., Sec. II.5). For T_ν this means **ν > 2d/3**, not ν > d.
- R1 (hex optimal) might therefore extend down to ν ∈ (4/3, 2] in 2D.
- The constraint ν > d comes from the GZL implementation, not from the math.
- In MGF language, C_{s,s,s} converges for s > 2/3.
- C_{1,1,1} (s=1) lies below the stated d<ν range, since ν=2=d. It is still convergent. It is fine to quote it as a proven anchor point.

---

## Verdicts

| Claim | Verdict | Confidence |
|---|---|---|
| **R1** (hex optimal for d<ν<ν*=3.91836, square optimal on (ν*, ν2), first-order jump along \|τ\|=1, local-stability thresholds 3.8063 / 4.2784) | **No prior work found.** Related but different: Robles-Navarro et al. 2025 compare only sq vs hex at ν=3 and fixed nearest-neighbour distance (their numbers, converted, agree with our hex<sq at ν=3). Luo–Wei and Mueller–Ho give hex→…→square sequences for different pair energies. | High (~85%) that the exact statement for T_ν is new |
| **R2** (continuous square → rectangle transition at ν2≈8.606; aspect ratios; limit √((√17−1)/2)) | **No prior work found.** Related but different: Bétermin–Šamaj–Travěnec 2023 (square→rectangle for pair potentials, density-driven); Robles-Navarro et al. (rectangular distortion for LJ+ATM with chain collapse, different mechanism). | High (~85%) |
| **R3** (T_{2s} = π^{3s} C_{s,s,s}; minimiser of C_{s,s,s} jumps ρ→i at s*=1.959, rectangular for s>4.303) | Identity: essentially **known by construction**. It is the MGF lattice sum rescaled to unit covolume, implicit in DGV/DGGV and in the GZL MGF dictionary. s=1: **known** (Rankin/Cassels/Ennola/Diananda via C_{1,1,1}=E_3+ζ(3)). The extremal statements for s>1 (in particular C_{2,2,2}, C_{3,3,3}, C_{4,4,4} minimised at τ=i; C_{s,s,s}, s≥5, at τ=iy, y>1): **no prior work found**. | Medium-high (~75%). The MGF literature is large and a plot or remark could hide in a physics paper, but targeted full-text searches found nothing |
| **3D** (BCC global minimiser for ν∈[3.2,12]; FCC local min at 3.5, saddle for ν≥4.5; BCC minimum along the Bain path) | **Related, partially anticipated.** Robles-Navarro et al. 2025 show that repulsive three-body (ATM) energy has a minimum at BCC along the Bain path and favours BCC over FCC (ν=3, ATM, NN-normalised). Their f_r data converts to T_3(BCC) < T_3(FCC). Buchheit–Busse 2025 do fcc/bcc LJ+ATM phase diagrams. The global search over all 3D lattices, the ν range, and the FCC saddle result are **not found**. | Medium (~65%) for novelty of the global/saddle statements. The Bain-path part is largely anticipated |

**Main risk.** The Buchheit/Busse/Schwerdtfeger group states it has an ongoing programme on many-body interactions and crystal-lattice stability, and owns the tooling (ζ^{(n)}, GZL). Any 2D/3D optimal-lattice paper for ζ^{(3)} from them would directly overlap. None was found as of 2026-09-25. The latest items found are arXiv:2609.28282 and arXiv:2609.18918, which are computational papers.

**Suggested citations if we write this up:**
- arXiv:2504.07338 and arXiv:2504.11989 (closest; three-body zeta, ATM, Bain path)
- arXiv:2609.18918 (GZL)
- arXiv:2101.05602 (EAM many-body lattice ground states)
- arXiv:2004.13882, arXiv:2110.08728 and arXiv:2312.02497 (hex/square/rectangular transitions for pair or engineered energies)
- arXiv:2312.01395 (square→rectangle transitions)
- arXiv:2107.14020 and arXiv:2609.17356 (3D pair Riesz: FCC vs BCC)
- arXiv:2501.05746 (BCC minimum on the cuboidal family for Epstein zeta)
- arXiv:1502.06698, arXiv:1512.06779 and arXiv:1708.07998 (MGFs, C_{a,b,c} for complex a,b,c)
- arXiv:1504.01153 (review)
- arXiv:1902.05438 (universal optimality)
