# Novelty check (physics / chemistry / soft matter): pure three-body power-law lattice energies

Date: 2026-09-25.

**Object under test.** T_ν(Λ) = Σ'_{x,y∈Λ, x≠y} |x|^{-ν}|y|^{-ν}|x−y|^{-ν} over Bravais lattices Λ at unit density. This is a purely three-body, product-form, repulsive power law, i.e. the triangle graph zeta. It is identical to the "three-body zeta function" ζ_Λ^{(3)}(ν,ν,ν) of Robles-Navarro et al. / Buchheit–Busse, and to the isotropic (radial) part of ATM when ν=3.

**Claims under test.**
- (C1) 2D, hexagonal→square transition.
  - Hex is optimal for ν < 3.9184.
  - Square is optimal for 3.9184 < ν < 8.606. The hex→square transition is first order.
- (C2) 2D, square→rectangular transition.
  - The transition at ν ≈ 8.606 is continuous.
  - The aspect ratio then grows to about 1.2–1.25 at large ν.
- (C3) 3D, global search over Bravais lattices.
  - BCC is optimal for ν ∈ [3.2, 12].
  - FCC is a saddle for ν ≥ 4.5.
  - BCC is the minimum along the Bain path.

**Verification status.**
- arXiv IDs were verified by fetching arxiv.org abstract pages on 2026-09-25, unless marked otherwise.
- For the two key papers (arXiv:2504.07338, arXiv:2504.11989) the full PDF text was read (downloaded and converted with pdftotext). Quotes are taken from those texts.
- "Unverified" means I saw the item only in a search snippet or catalogue entry and did not read its content.

---

## 1. Verdict table

| Claim | Verdict | Confidence | Closest prior work |
|---|---|---|---|
| C1 (hex→square at ν≈3.918 for pure product-form three-body) | **No prior work found** | medium-high | 2504.07338 §III B (ATM incl. angular term, ν=3 only, SL vs HL at fixed nearest-neighbour distance) |
| C2 (continuous square→rectangular at ν≈8.606, aspect ratio ~1.2–1.25) | **No prior work found** for this energy. The mechanism and transition type are *related but different* from prior work | medium-high | 2504.07338 Fig. 6 (rectangular distortion of SL under LJ+ATM, driven by ATM turning *attractive* for collinear triples); 2312.01395 (general theory of square→rectangle transitions, pair potentials) |
| C3a (BCC < FCC for pure product-form three-body) | **Partly anticipated / related.** Published ν=3 data imply it (see §3), but nobody has stated it for the pure product form, at unit density, or for general ν | high that ν=3 data exist; high that general-ν statement is absent | 2504.07338 (Bain path, Table IV, Appendix H) |
| C3b (BCC is the global minimiser over *all* Bravais lattices, ν∈[3.2,12]) | **No prior work found** | medium-high | 2504.07338 only covers the 1-parameter Bain path and states "the bcc phase remains susceptible to further cuboidal distortions" (for LJ+ATM) |
| C3c (BCC stationary along the Bain path) | **Found in prior work** (general-ν proof) | high | 2504.07338 Appendix H proves ∂_A ζ^{(3)}_{Λ(A)}(ν⃗)=0 at A=1/2 for **any** ν⃗ |
| C3d (FCC a saddle for ν ≥ 4.5) | **No prior work found** | medium | none |
| Mechanism ("repulsive three-body penalises the many small equilateral triangles of hex/FCC") | **Qualitatively stated before** | high | 2504.07338; Sellin & Babaev 2013 (strong 3-body repulsion turns triangles into lines) |

---

## 2. Key prior work in detail

### 2.1 Robles-Navarro, Cooper, Buchheit, Busse, Burrows, Smits, Schwerdtfeger (2025). The closest match

**Bibliographic data**
- Title: "Exact lattice summations for Lennard-Jones potentials coupled to a three-body Axilrod-Teller-Muto term applied to cuboidal phase transitions".
- Identifiers: [arXiv:2504.07338](https://arxiv.org/abs/2504.07338), v1 9 Apr 2025. Published as J. Chem. Phys. 163, 094104 (2025), DOI [10.1063/5.0276677](https://doi.org/10.1063/5.0276677). The DOI was confirmed via Semantic Scholar and the reference list of 2504.11989.
- Correction to note 01: that note says "No arXiv ID found" for this paper. That is wrong; the ID is 2504.07338.

**What it contains that overlaps with us**

- **Defines exactly our object.** Eq. (26) defines the "three-body zeta function" ζ_Λ^{(3)}(ν⃗) = Σ' |x|^{-ν1}|y|^{-ν2}|y−x|^{-ν3}.
  - The ATM energy is split into a radial part and an angular part.
  - The radial part is f_r^{(3)} = (1/6) ζ_Λ^{(3)}(3,3,3), Eq. (27). This equals our T_3 up to the factor 1/6.
- **Computes 2D numbers for SL and HL at ν=3 only.** The lattices are normalised to nearest-neighbour distance R=1 (§III B, Eqs. 44, 51):
  - Square: f_r = 2.27548228589309, f_a = −1.50538863537637, f_coh = 0.7700936505167162.
  - Hex: f_r = 4.263827935991082, f_a = −2.3454945711432025, f_coh = 1.9183333648478795.
  - Their conclusion is the triangle-penalty mechanism: "the hexagonal lattice is more strongly destabilized by adding three-body interactions compared to the square lattice because f_coh,hex > f_coh,sq. This is due to the hexagonal lattice being a close-packed structure in 2D with the highest packing density and kissing number."
  - **Caveat.** This comparison is made at fixed nearest-neighbour distance, not at fixed density, and always with an LJ two-body part. They never state that SL beats HL at equal density for any exponent. My conversion to unit density at ν=3 (see §3) gives HL < SL. This agrees with our C1, where hex wins below ν≈3.918.
- **Rectangular distortion in 2D is a different mechanism.**
  - They parametrise SL→rectangular lattices by an aspect parameter γ and scan γ against the ATM coupling λ (Fig. 6).
  - The rectangular lattice with γ=2 (one side twice the other) becomes more stable than SL at large λ "because it is located at the region where the three-body potential becomes attractive, similar to the case of the linear chain". In other words, ATM's angular factor makes collinear triples attractive, and in the limit the distortion goes "into a set of linear chains".
  - Our pure product-form T_ν is strictly repulsive. Our continuous square→rectangle transition at ν≈8.6, with a modest aspect ratio of about 1.2, has a different origin and is a different transition. **Related but different.**
- **3D Bain path.**
  - Table IV (Appendix H) lists f_r(A), f_a(A) and f_coh(A) for A ∈ [0.1, 1] in steps of 1/90. They are normalised to fixed body-centre distance R_bc = 1, *not* fixed density.
  - Text: "the highest repulsive three-body energy occurs for the densely packed fcc lattice, while the energy minimum is reached for the bcc lattice within the studied parameter range". This refers to the *full* ATM f_coh at fixed R_bc.
  - Appendix H proves, for **any** ν⃗ ∈ C³, that ∂/∂A ζ^{(3)}_{Λ(A)}(ν⃗) = 0 at A = 1/2 (bcc).
    - That covers our T_ν stationarity at bcc along the Bain path for all ν.
    - The unit cell volume V(A) = 2√A/(A+1)^{3/2} is also stationary at A=1/2, so the stationarity carries over to fixed density (my check).
  - Abstract: "strong repulsive three-body interactions can destabilize the fcc phase and render bcc energetically favorable for soft LJ potentials. However, even in this scenario, the bcc phase remains susceptible to further cuboidal distortions."
    - This concerns LJ+ATM with R optimised, not the pure three-body energy.
    - It does NOT claim global optimality over all Bravais lattices.
  - Fig. 20 plots ζ^{(3)}_{SL} and ζ^{(3)}_{HL} for ν⃗ = (ν−2, ν, ν+2), not (ν,ν,ν). They do not interpret it as an energy comparison.
- **Scope limits relative to us.**
  - Exponent is ν=3 only (triple-dipole).
  - The angular factor is always included in the physics conclusions.
  - A two-body LJ part is always present.
  - Normalisation is fixed distance, not density.
  - 2D covers only SL, HL and the rectangular family.
  - 3D covers only the Bain path.

### 2.2 Buchheit & Busse (2025/2026). Method paper

**Bibliographic data:** "Epstein zeta method for many-body lattice sums", [arXiv:2504.11989](https://arxiv.org/abs/2504.11989), v1 16 Apr 2025, v2 12 Jun 2026. Numer. Math., DOI 10.1007/s00211-026-01558-y. The DOI appears on the Springer page found in search; I did not open the article itself.

**What it contains**
- It computes ζ^{(3)} and the n-body ζ^{(n)}(ν,…,ν) on Λ_sq and Λ_hex. These are *accuracy benchmarks only* (Figs. 3, 7, 8). They note that such sums "correspond to the isotropic part of the n-body cohesive energy obtained from a Drude model", citing Schwerdtfeger et al., Angew. Chem. 2016.
- The only physics application is §5: LJ + full ATM, fcc vs bcc, enthalpy at finite pressure. They find "for sufficiently large coupling strengths λ, three-body interactions can destabilize fcc and lead to the formation of a bcc phase".
- There is no lattice optimisation of the pure product-form sum and no 2D phase comparison. **Related (tool), not a prior result.**

### 2.3 Other related but different work

**2D, three-body or many-body repulsion and non-hexagonal order**
- **Sellin & Babaev**, "Stripe, gossamer, and glassy phases in systems with strong non-pairwise interactions", [arXiv:1308.2109](https://arxiv.org/abs/1308.2109), Phys. Rev. E 88, 042305 (2013).
  - Setup: 2D classical particles with a two-body potential (short-range repulsive, long-range attractive) plus a *repulsive Gaussian three-body term*, inspired by type-1.5 superconductor vortices. This is not a power law and not product-form.
  - Result: for strong three-body repulsion "the ground state of three particles … will be that of a straight line instead of a triangular configuration". Monte Carlo gives hexagonal → gossamer → stripe phases, not a square lattice.
  - Verdict: same "triangles are penalised" intuition, different energy and different outcome.
- **Edström**, "Three and Four-Body Intervortex Forces in the Ginzburg-Landau Models…", [arXiv:1209.4334](https://arxiv.org/abs/1209.4334), Physica C (2013), DOI 10.1016/j.physc.2013.01.020.
  - The three-body vortex force is repulsive in the type-1.5 case. The paper does not discuss lattice structure.
- **Wolf, Vagov, Shanenko, Axt, Aguiar**, "Vortex matter stabilized by many-body interactions", Phys. Rev. B 96, 144515 (2017). **Unverified:** seen only as an APS search hit.
  - It covers many-body vortex clusters in intertype superconductors, not square lattices from a product-form three-body term.
- **Bétermin, Friedrich, Stefanelli**, "Lattice ground states for Embedded-Atom Models in 2D and 3D", [arXiv:2101.05602](https://arxiv.org/abs/2101.05602), Lett. Math. Phys. (2021), DOI 10.1007/s11005-021-01446-6.
  - Setup: an EAM many-body energy, i.e. an embedding F of a star sum ρ = Σ_x r^{-s}, plus a pair term.
  - 2D: the square lattice becomes a *local* minimiser for some parameters. The authors call this "the first occurrence of such minimality among all possible lattices, without a density constraint". Triangular→orthorhombic transitions also occur.
  - 3D: simple cubic is favoured over FCC and BCC for small s.
  - Verdict: many-body, but star/tree-type rather than triangle-type. Related but different.
- **Mainini, Piovano, Stefanelli**, "Finite crystallization in the square lattice", Nonlinearity 27, 717 (2014), DOI 10.1088/0951-7715/27/4/717. arXiv ID not verified.
  - Two-body plus *angular* three-body terms that explicitly favour 90° bond angles. The square lattice is built in. Different.
- **Bétermin, De Luca, Petrache**, "Crystallization to the square lattice for a two-body potential", [arXiv:1907.06105](https://arxiv.org/abs/1907.06105).
  - A square ground state from a pair potential. Relevant only as a contrast: a square lattice does not require three-body forces.
- **Bétermin, Šamaj, Travěnec**, "Structural transitions in interacting lattice systems", [arXiv:2312.01395](https://arxiv.org/abs/2312.01395), Anal. Math. Phys. 14, 27 (2024).
  - A general theory of first- versus second-order square→rectangle transitions at fixed density, with tricritical points, for *pair* potentials (double Yukawa, Yukawa–Coulomb).
  - Useful framework for our C2 (continuous square→rectangle). It does not treat three-body energies.

**Colloids with three-body interactions**
- **Brunner, Dobnikar, von Grünberg, Bechinger**, "Direct measurement of three-body interactions amongst charged colloids", PRL 92, 078301 (2004), [arXiv:0801.3917](https://arxiv.org/abs/0801.3917).
- **Dobnikar et al.**, PRE 69, 031402 (2004), [arXiv:0801.3922](https://arxiv.org/abs/0801.3922).
- **Dobnikar et al.**, SPIE 5514 (2004), [arXiv:0801.3920](https://arxiv.org/abs/0801.3920).
- Findings of these three papers: the three-body term is **attractive**, and the four-body term is repulsive. There is no square-lattice stabilisation.
- **Hynninen, Dijkstra, van Roij**, "Effect of three-body interactions on the phase behavior of charge-stabilized colloidal suspensions", PRE 69, 061407 (2004). Seen via the APS abstract page; arXiv ID not verified.
  - Triplet *attractions* are found. At low salinity a dilute fluid coexists with an almost-close-packed fcc phase, and at intermediate salinity there is bcc–fcc coexistence.
  - Opposite sign to ours. Different.

**Ultracold gases and Rydberg systems**
- **Büchler, Micheli, Zoller**, "Three-body interactions with cold polar molecules", [cond-mat/0703688](https://arxiv.org/abs/cond-mat/0703688) (Nature Phys. 2007).
- Capogrosso-Sansone et al., [arXiv:0805.1408](https://arxiv.org/abs/0805.1408). Seen only as a search hit.
- Honeycomb follow-up, [arXiv:0911.0312](https://arxiv.org/abs/0911.0312). Seen only as a search hit.
- **Samajdar, Lukin, Walther**, "Three-body interactions in Rydberg lattices", [arXiv:2604.11870](https://arxiv.org/abs/2604.11870).
- All of these are lattice-gas occupation problems on a *fixed* optical or tweezer lattice. They do not select a Bravais geometry in the continuum. Not relevant to C1–C3.
- **Cooper, Rezayi, Simon**, [cond-mat/0505759](https://arxiv.org/abs/cond-mat/0505759), PRL 95, 200402 (2005).
  - Triangular→square→stripe vortex lattices appear, but they are driven by *two-body* dipolar interactions. This is an analogue only.

**Rare-gas and metal crystal-structure literature (fcc/hcp/bcc with three-body forces)**
- **Axilrod**, "Triple-Dipole Interaction. II. Cohesion in crystals of the rare gases", J. Chem. Phys. 19, 724 (1951). Content not read: citation only.
- **Jansen & Lombardi**, "Dreikörperkräfte und die Stabilität einfacher Kristallgitter", Z. Phys. 190, 161 (1966), DOI 10.1007/BF01327141. Metadata was verified via Crossref and Semantic Scholar; the abstract was not available. My presumption, not verified: it covers three-body *exchange* forces between rare-gas atoms and the fcc–hcp stability question.
- **Rossi & Danon**, "Stability of crystal structures of heavy rare gases" (1965), DOI 10.1016/0022-3697(65)90003-X. Metadata only.
- **Doran & Zucker** (1971), "Higher order multipole three-body van der Waals interactions and stability of rare gas solids". Seen only as a ResearchGate hit.
- The classic "rare-gas crystal structure problem" in these papers is fcc vs hcp. Jansen explained fcc through three-body *exchange* terms that favour large triangle angles. No paper in this list establishes that "pure product-form three-body repulsion selects bcc over all Bravais lattices".
- **Soft matter bcc/A15 selection by many-body forces:** Ziherl & Kamien, "Soap froths and crystal structures", PRL 85, 3528 (2000). Seen via the APS abstract page and PubMed 11030938.
  - Many-body (corona surface-area, Kelvin-problem) interactions favour bcc/A15 over fcc.
  - Qualitatively the same "many-body ⇒ bcc" message, but a different energy (Voronoi cell area, not a power law).
- **Metals:** Jerabek, Burrows, Schwerdtfeger, "A smooth bcc to fcc phase transition for metallic lithium", [arXiv:2211.06317](https://arxiv.org/abs/2211.06317). It is DFT-based and does not involve a model three-body term.

---

## 3. OWN CHECK: our claims against published numbers

This is my computation from published data, not a literature result.

**Method.** Table IV of arXiv:2504.07338 gives f_r(A) = (1/6) ζ^{(3)}_{Λ(A)}(3,3,3) on the Bain path with R_bc = 1. I rescaled to unit density with factor V(A)^3, using V(A) = 2√A/(A+1)^{3/2} and energy ∝ length^{-9}.

**3D results at ν=3, unit density** (T_3/6):

| A | lattice | f_r at unit density |
|---|---|---|
| 1/3 | acc/sc | 11.195884 |
| 0.4 | | 11.136627 |
| 1/2 | bcc | **11.118936** (minimum along the path) |
| ≈0.845 | | 11.133063 (maximum along the path) |
| 0.9889 | | 11.132412 |
| 1 | fcc | 11.132403 (a shallow *local minimum* along the path) |

- bcc is lower than fcc by 0.12%.
- This is consistent with our C3: BCC beats FCC, and FCC is still locally stable along Bain at ν=3, below our ν≈4.5 saddle onset.
- So a careful reader of 2504.07338 could have extracted "bcc < fcc for the pure radial ATM term at unit density at ν=3" from their data, although the paper never states it.
- For comparison, the *full* ATM (with angular part) at unit density is monotone along the path: bcc 6.742596 < fcc 6.782190, with fcc a local maximum along the path.

**2D results at ν=3, unit density.**
- SL: 2.275482.
- HL: 4.263828 × (√3/2)^{4.5} = 2.231966.
- HL < SL by 1.9%. This is consistent with C1 (hex optimal below ν≈3.918).
- No published number exists for ν ≠ 3, so the crossover at 3.9184 is not in prior work.

---

## 4. Search log (queries with no relevant hit)

**Web searches**
- "three-body interactions stabilize square lattice two-dimensional crystal": hits were lattice-gas Hubbard models and angular-potential crystallization.
- "three-body interaction induces square lattice supersolid droplet crystal": dipolar-droplet square lattices come from trap geometry or LHY terms, not from three-body product terms (e.g. [arXiv:2209.14373](https://arxiv.org/abs/2209.14373), seen only as a search hit).
- "rotating BEC three-body interaction vortex lattice triangular to square": only dipolar two-body mechanisms.
- "non-additive three-body repulsion … square crystal": found Sellin & Babaev 2013 and the three-body hard-core model [arXiv:1410.1454](https://arxiv.org/abs/1410.1454). The latter gives a dimer solid in 2D and a simple-hexagonal lattice in 3D; it is a hard cut-off, not a power law.
- "triple-dipole energy two-dimensional triangular lattice square lattice adsorbed monolayer": nothing relevant.
- "Axilrod-Teller energy two-dimensional lattice square hexagonal": only 2504.07338.
- "modular graph function C_{a,b,c} minimum fundamental domain": no paper found on minimising C_{s,s,s}(τ) over τ.
  - The 2D problem is equivalent to minimising the dihedral MGF C_{s,s,s}, with s = ν/2; see note 06.
  - The MGF literature is about Laplace equations, identities and integrals, not extrema.

**arXiv listing searches:** "three-body zeta function lattice", "Axilrod-Teller-Muto lattice sum", "three-body lattice sum Epstein". The only hits were 2504.07338, 2504.11989 and 2412.16317.

**Author feeds**
- P. Schwerdtfeger and A. Robles-Navarro (arXiv author search): no paper after 2504.07338 on three-body lattice optimisation.
- L. Bétermin: no three-body Riesz lattice paper.

**Semantic Scholar**
- Citations of 2504.07338 (6 found) and of 2504.11989 (4 found) include 2609.18761, 2604.22743, 2509.26274, 2412.16317 and a Francium DFT paper (DOI 10.1039/d5cp04548g).
- None of them treats three-body lattice optimisation.
- Coverage of journal-only citations may be incomplete.

---

## 5. Caveats and gaps

- I did not read the full text of the classic 1950s–1970s rare-gas three-body papers (Axilrod 1951, Jansen & Lombardi 1966, Rossi & Danon 1965, Doran & Zucker 1971). They may contain bcc or sc triple-dipole lattice sums, though for ATM with its angular factor and ν=3. A radial-only, variable-ν study there is unlikely but not excluded.
- I did not check journal-only literature on adsorbed monolayers with three-body forces (Bruch, Cole, Zaremba, *Physical Adsorption*).
- Semantic Scholar rate-limited several keyword queries (HTTP 429). Those gaps were covered by web search and arXiv listing searches instead.
- 2504.07338 Appendix H already gives the general-ν stationarity of bcc on the Bain path. Any write-up should cite it rather than claim it.
- Some mathematical-physics work on EAM or angular three-body terms reports "square lattice from many-body terms" (Bétermin–Friedrich–Stefanelli; Mainini–Piovano–Stefanelli). The novelty statement should therefore be phrased narrowly:
  - "for the pure product-form (triangle) power-law three-body energy, at fixed density, over all 2D Bravais lattices", and
  - "over all 3D Bravais lattices for ν ∈ [3.2, 12]".

---

## 6. Summary

- **C1 (hex→square at ν≈3.918): no prior work found.** Medium-high confidence.
  - The only direct precedent is 2504.07338, and only at ν=3. It uses nearest-neighbour normalisation, and its physics includes the ATM angular term and an LJ part.
  - It already states qualitatively that the hexagonal lattice is "more strongly destabilized" by three-body forces because it is close-packed.
- **C2 (continuous square→rectangle at ν≈8.6): no prior work found.** Medium-high confidence.
  - 2504.07338 does find rectangular distortions of SL under LJ+ATM, but through ATM becoming *attractive* for collinear triples, and it heads towards linear chains. That is a different mechanism.
  - Bétermin–Šamaj–Travěnec 2024 give a general pair-potential theory for this transition type.
- **C3 (BCC in 3D): partially anticipated.**
  - "Repulsive three-body forces favour bcc over fcc" is established for LJ+ATM (2504.07338; 2504.11989).
  - Bcc stationarity on the Bain path is proven there for arbitrary exponents (Appendix H).
  - Their Table IV implies bcc < fcc for the pure radial ν=3 term at unit density (my conversion). They never state it.
  - Not found anywhere:
    - global BCC optimality over *all* Bravais lattices;
    - the ν-range [3.2, 12];
    - FCC becoming a saddle at ν≥4.5;
    - any statement for pure product-form three-body at general ν.
  - Medium-high confidence for these gaps.
