
**Kurzfazit.** Die stärkste neue Idee lässt sich mit der heutigen GZL sofort starten: ein **Atlas optimaler Gitter für Mehrkörper-Potenzgesetzenergien**. Die Dreikörper-Riesz-Energie eines Gitters ist exakt eine Dreiecks-Graph-Zetafunktion. Eigene Rechnungen, inzwischen mit einer von GZL unabhängigen Referenzimplementierung nachgeprüft (Abschnitt 7), zeigen: In 2D springt das Optimum bei **ν\* = 3,9183649026 in einem Übergang erster Ordnung vom Hexagonal- zum Quadratgitter** und geht bei **ν₂ ≈ 8,606 kontinuierlich in ein Rechteckgitter** über. In 3D ist **BCC bei allen untersuchten ν das beste gefundene Gitter**, und FCC wird ab ν ≈ 3,752 instabil. Das ist die Gegenrichtung zum gerade beanspruchten Beweis, dass für Zweikörperenergien FCC optimal ist. Den größten wissenschaftlichen Hebel hat eine **deterministische Hochtemperaturreihe für klassische langreichweitige Ising-, O(n)- und Perkolationsmodelle**. Deren Einbettungssummen sind genau die Objekte, die GZL berechnet. Im laufenden Streit um die Sak-Grenze (σ\* = 2 oder 2 − η) wäre sie die erste Stimme, die nicht auf Monte Carlo beruht. Am nächsten am Experiment liegt eine **Simulator-Material-Zwillingsstudie**: Dispersion, Gap und kritischer Punkt im thermodynamischen Limes, mit vollem van-der-Waals- bzw. Dipolschwanz, für die 256-Qubit-Simulation von TmMgGaO₄ und für NaTmSe₂. Hier ist das Scoop-Risiko hoch, weil die GZL-Autoren Rydberg-Arrays selbst als Ziel nennen. Die überraschendste Brücke: **Planare modulare Graphfunktionen der Stringtheorie sind exakt GZL-Graph-Zetas des dualen Graphen.** Das ist an der Zagier-Identität und an C₂,₂,₁ numerisch auf 10⁻¹⁴ bis 10⁻¹⁵ bestätigt (Abschnitt 7).

Kennzeichnung im ganzen Bericht: **[belegt]** heißt, es steht in einer zitierten Quelle. **[abgeleitet]** ist eine Schlussfolgerung aus belegten Fakten. **[spekulativ]** muss erst geprüft werden. **[Pilot]** markiert eigene, **nicht begutachtete** Rechnungen eines Recherche-Agenten mit GZL 1.0.0 und epsteinlib 0.6.2. **[nachgeprüft]** heißt: mit zwei unabhängigen Verfahren übereinstimmend berechnet (Abschnitt 7, Skripte in `research/verification/`). Die zugehörigen Skripte liegen unter `/home/user/GapYas/research/pilots/` (`mgf_test.py`, `mgf_test2.py`, `opt_test.py`, `opt2.py`, `opt3.py`). `mgf_test2.py` lädt `mgf_test.py` über einen relativen Pfad und muss deshalb aus diesem Verzeichnis gestartet werden.

---

## Was GZL wirklich kann: eine Softcore-Gittersumme mit genau einem Impuls

### Fünf Sätze zur Klasse der Summen

**(1)** GZL berechnet für einen endlichen Multigraphen G mit zwei Terminals s, t die Summe ζ_G(k) = Σ_{x_v ∈ Λ, v ≠ p} e^{−2πi(x_t − x_s)·k} Π_e K_e(x_{e+} − x_{e−}). Dabei ist ein Knoten p gepinnt, Λ = A·ℤ^d ein Bravais-Gitter mit d = 1, 2, 3, und das Ergebnis ist die Gitter-Fourier-Transformierte eines Zweipunktkorrelators **[belegt]** ([README](https://github.com/graph-zeta/gzl), [arXiv:2609.18918](https://arxiv.org/abs/2609.18918)). **(2)** Die Summe ist **softcore**: Nur die beiden Endpunkte derselben Kante werden durch K_ν(0) = 0 getrennt. Nicht benachbarte Knoten dürfen auf demselben Platz liegen. Die physikalische Hardcore-Bedingung wird in der Physik-Arbeit durch eine kombinatorische „softcore mapping" in die Korpus-Vorfaktoren eingerechnet, die **nicht als Code in GZL 1.0.0 enthalten** ist **[belegt]** ([arXiv:2609.18761](https://arxiv.org/html/2609.18761), [README](https://github.com/graph-zeta/gzl)). **(3)** Die Kerne haben die Form K = a(x) + Σ_j b_j |x|^{−ν_j}. Der Anteil a ist reell, gerade, kompakt getragen und darf anisotrop sein, der Potenzgesetzanteil ist isotrop. Jede Kante darf einen eigenen Kern haben, und ν = ∞ liefert den Nächste-Nachbarn-Indikator mit exakten ganzzahligen Einbettungszahlen **[belegt]** ([DOCUMENTATION.md](https://github.com/graph-zeta/gzl/blob/main/DOCUMENTATION.md)). **(4)** Die Rechnung faktorisiert in Blöcke. Brücken sind geschlossene Epstein-Zetas. Einfache Zyklen bei k = 0 laufen über ein Brillouin-Zonen-Integral von Epstein-Produkten. Serien-parallele Blöcke laufen über eine semi-analytische Algebra mit Fehler ~ n^{−(d+σ+2)} und Kosten linear in |V|. Blöcke mit höherer Baumweite laufen über ein Torus-Tensornetz mit Kosten O(N^{tw}), N = n^d. Das ganze k-Gitter kostet per FFT so viel wie ein einzelner Punkt **[belegt]** ([arXiv:2609.18918](https://arxiv.org/html/2609.18918)). **(5)** GZL berechnet **nicht**: nicht-Bravais-Gitter, anisotrope Potenzgesetzschwänze, komplexe oder ungerade Kerne, oszillierende Schwänze (RKKY) oder Yukawa-Schwänze außer als abgeschnittene Tabelle, Summen mit zwei oder mehr freien Terminals, nicht-Brücken-Blöcke mit ν ≤ d und Ergebnisse jenseits doppelter Genauigkeit **[belegt]** ([DOCUMENTATION.md](https://github.com/graph-zeta/gzl/blob/main/DOCUMENTATION.md)).

### Grenzen und Korrekturen am Auftrag (Abschnitt 3)

| Punkt | Aussage im Auftrag | Tatsächlicher Stand v1.0.0 | Status |
|---|---|---|---|
| Summationsbereich | „Platzierungen der Knoten auf Λ" | Softcore: Nur Kantenendpunkte sind getrennt. Hardcore→Softcore-Abbildung ist nicht mitgeliefert. | **Korrektur/Präzisierung** |
| Exponenten nahe d | „präzise auch für ν nahe d" | Stimmt nur für ν > d. Brücken werden meromorph fortgesetzt. **Jeder andere Block mit ν ≤ d wirft `UnsupportedLatticeSumError`.** Benchmarks reichen bis σ = ν − d = 0,5. | **Einschränkung** |
| „analytisch" | Brücken und Zyklen analytisch | Brücken sind geschlossen. Zyklen nur bei k = 0 über `zeta_circle`, und nur für d ≤ 3. SP-Blöcke sind semi-analytisch mit algebraischem Fehler. | Präzisierung |
| Kosten tw ≥ 3 | „polynomial in der Gittergröße" | O(N^{tw}) Operationen und O(N^{tw−1}) Speicher. **In d = 3 ist praktisch n_points ≤ 16.** Bei n = 8 liegt der relative Fehler hoher Ordnungen um 10⁻³. | Präzisierung |
| `slab_zeta` | klingt nach Film-/Slab-Geometrie | **Ist keine Slab-Geometrie**, sondern eine speichersparende Engine für den Vakuumwert eines dichten Blocks | **Korrektur** |
| Kerne | reell, gerade | Stimmt. Die Mathe-Arbeit erlaubt komplexe ℓ¹-Kerne, **die Bibliothek ist enger als die Theorie**. | Präzisierung |
| Anisotropie | Roadmap | Stimmt. KTmSe₂ in 3D wurde mit einer experimentellen Erweiterung nur bis Ordnung 5 gerechnet. Die anisotrope Epstein-Zeta ist in EpsteinLib bereits vorhanden ([arXiv:2609.28282](https://arxiv.org/abs/2609.28282)). | bestätigt |
| Dimension | d = 1, 2, 3 | Offiziell ja. d ≥ 4 läuft teilweise (Brücken exakt), ist aber unkalibriert und nicht publikationsfähig. | Ergänzung |
| Gittermatrix | jede Matrix A | A muss quadratisch und regulär sein, also keine eingebetteten niedrigdimensionalen Gitter | Präzisierung |
| Impuls | e^{ik·(…)} | Genau **ein** externer Impuls, in fraktionalen reziproken Koordinaten. Zwei oder mehr freie Terminals werden abgelehnt, 2qp-Sektoren gibt es also nicht. | **Einschränkung** |
| Korpora | TFIM 0qp bis Ordnung 13, 1qp bis Ordnung 11 | Bestätigt: 8.403 bzw. 22.677 Graphen. Die Korpora sind dimensionsunabhängig. | bestätigt |

Die Leistungsbehauptungen halten der Prüfung stand **[belegt]**. 114 Reihen mit 998 Koeffizienten stimmen mit publizierten MC-Daten überein: **maximal 3,1 σ_MC, Median 0,12 σ_MC**. Der Satz wurde in **5 Minuten auf 8 M1-Kernen** gerechnet, gegenüber etwa **3×10⁴ Kernstunden pro MC-Reihe und Impuls**. Eine 3D-Dispersion der Ordnung 10 auf 16³ Impulspunkten dauert rund 10 Minuten auf einem Kern ([arXiv:2609.18918](https://arxiv.org/html/2609.18918)). Die Physik-Arbeit liefert 0qp- und 1qp-Reihen in 3D bis Ordnung 13 bzw. 11 in etwa 200 s ([arXiv:2609.18761](https://arxiv.org/html/2609.18761)).

Stand der Verbreitung am 25.09.2026: 1 Stern, 0 Forks, 0 Issues, keine Nutzung außerhalb der Gruppe Buchheit/Schmidt. **Das Feld ist also offen.** Die Gruppe Sbierski/Lesanovsky ist als einzige externe Nutzerin der Vorläufer-Werkzeuge nachweisbar und wäre die wahrscheinlichste erste Konkurrenz ([arXiv:2604.22743](https://arxiv.org/abs/2604.22743)).

### Scoop-Landkarte: was die GZL-Autoren selbst angekündigt haben

Die Autoren beanspruchen öffentlich folgende Richtungen **[belegt]** ([arXiv:2609.18918 §10](https://arxiv.org/html/2609.18918), [arXiv:2609.18761 §VIII](https://arxiv.org/html/2609.18761)):

- anisotrope Dipolkerne, konkret LiHoF₄, Fe₈, Mn₁₂ und RE(OH)₃;
- mehratomige Gitter und frustrierte Systeme;
- Heisenberg- und weitere Modelle;
- Spektralgewichte und höhere Quasiteilchensektoren;
- endliche Temperatur, ausdrücklich mit Burkard, Schneider und Sbierski;
- korrelierte Unordnung;
- Rydberg-vdW-Arrays im Vergleich mit dynamischen Strukturfaktoren;
- eine analytische Singularitätsbehandlung für die Tensornetz-Blöcke.

Wer dort ohne Kooperation arbeitet, läuft Gefahr, überholt zu werden. **Nicht genannt** sind:

- klassische Hochtemperatur- und Perkolationsreihen;
- Gitteroptimierung;
- modulare Graphfunktionen;
- Gitter-Störungstheorie.

Im Umfeld der Autoren liegen allerdings ATM-Energien entlang des Bain-Pfads ([arXiv:2504.07338](https://arxiv.org/abs/2504.07338), [arXiv:2504.11989](https://arxiv.org/abs/2504.11989)).

---

## Top-5-Projektideen: Zwei Transfer-Ideen führen, zwei Physik-Ideen sind experimentnah

| # | Idee | Spur | Neuheit | Passung | Hebel | Machbarkeit | Wow | Anschluss | Σ | Scoop |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | Optimale Gitter für Mehrkörper-Potenzgesetzenergien | B | 5 | 5 | 4 | 5 | 5 | 4 | **28** | mittel |
| 2 | Klassische LR-Hochtemperatur- und Perkolationsreihen gegen den Sak-Streit | B | 4 | 3 | 4 | 2 | 3 | 5 | **21** | mittel |
| 3 | Modulare Graphfunktionen als duale Graph-Zetas | B | 5 | 3 | 2 | 4 | 4 | 3 | **21** | gering |
| 4 | Simulator-Material-Zwillinge: Rydberg-vdW und Tm-Dreiecksmagnete | A | 3 | 4 | 4 | 4 | 5 | 4 | **24** | **hoch** |
| 5 | Quanten-Sak-Test und Log-Korrekturen im LRTFIM (d = 1–3) | A | 3 | 5 | 3 | 4 | 3 | 4 | **22** | mittel–hoch |

Die Reihenfolge gewichtet die Neuheit stärker als die bloße Punktsumme, weil der Auftrag „etwas völlig Neues" verlangt.

### Idee 1: Optimale Gitter für Mehrkörper-Potenzgesetzenergien („Mehrkörper-Sarnak–Strömbergsson")

**Pitch.** GZL als erstes Werkzeug, das die Frage „Welches Gitter minimiert eine Drei- oder Vierkörper-Potenzgesetzenergie?" über den *gesamten* Modulraum der 2D- und 3D-Gitter in Sekunden pro Punkt beantwortet. Die Pilotrechnungen deuten auf eine andere Antwort als im Zweikörperfall.

**Das Problem im Zielfeld.** Die Optimierung von Gitterenergien ist in der Mathematik gerade ein heißes Thema. Luo und Wei beanspruchen einen Beweis der Sarnak–Strömbergsson-Vermutung: FCC minimiert die 3D-Epstein-Zeta für alle s > 3/2. Der Beweis ist zehn Tage alt und noch nicht begutachtet **[belegt]** ([arXiv:2609.17356](https://arxiv.org/abs/2609.17356)). Dazu kommen Klassifikationen von Minimierern für Verhältnisse von Theta- und Epstein-Funktionen ([arXiv:2605.07580](https://arxiv.org/abs/2605.07580)) und überraschende Phasenübergänge der Minimierer bei p-Normen ([arXiv:2407.20762](https://arxiv.org/abs/2407.20762)). Die ganze Linie behandelt aber **Zweikörper-Energien**, also Summen über eine Kante. In der Kristallchemie ist dagegen belegt, dass Dreikörperterme (ATM) die Phasenstabilität verschieben. Starke repulsive ATM-Beiträge können BCC gegenüber FCC begünstigen, doch „the bcc phase remains susceptible to further cuboidal distortions". Untersucht wurde dort nur der Bain-Pfad, also eine Ein-Parameter-Kurve **[belegt]** ([arXiv:2504.07338](https://arxiv.org/abs/2504.07338), [arXiv:2504.11989](https://arxiv.org/abs/2504.11989)). Eine Websuche nach Dreikörper-Gitteroptimierung über alle 2D- und 3D-Gitter fand keine Arbeit aus 2025 **[belegt als Negativbefund]**. Der Grund für die Lücke ist schlicht Rechenaufwand: Eine Dreikörpersumme ist eine Doppelsumme über ein unendliches Gitter mit langsamem Abfall nahe ν = d.

**Abbildung auf GZL.** Die Dreikörper-Riesz-Energie E₃(Λ) = Σ'_{x,y} |x|^{−ν}|y|^{−ν}|x−y|^{−ν} ist **exakt der Dreiecksgraph (3-Zyklus)**. Vierkörperenergien sind der 4-Zyklus oder K₄ (tw = 3, Tensornetz). Das Gitter A ist der freie Parameter. In 2D ist der Modulraum die Fundamentaldomäne (2 Parameter), in 3D sind es 5 Parameter (Gitter modulo Rotation und Skala) **[abgeleitet]**.

```python
import numpy as np, gzl, epsteinlib
def A2(tau):                        # 2D-Gitter mit Einheitsvolumen, tau in der Fundamentaldomäne
    A = np.array([[1, tau.real], [0, tau.imag]]); return A/np.sqrt(tau.imag)
E2 = lambda A, nu: epsteinlib.epstein_zeta(nu, A, np.zeros(2), np.zeros(2)).real  # Zweikörper
E3 = lambda A, nu: gzl.zeta_circle([nu, nu, nu], A)                               # Dreieck
K4 = np.array([[0,1],[0,2],[0,3],[1,2],[1,3],[2,3]])
E4 = lambda A, nu: gzl.evaluate_graph(K4, np.full(6, nu), A, n_points=48)         # Vierkörper, tw=3
# Observable: argmin_A [c2*E2 + c3*E3 + c4*E4] über den Modulraum -> Phasendiagramm in (nu, c3/c2)
```

**Pilotergebnisse [Pilot, inzwischen nachgeprüft und ergänzt, siehe Abschnitt 7]** (`opt_test.py`, `opt2.py`, `opt3.py`; GZL `zeta_circle`, Einheitsvolumen, gleiche ν auf allen drei Kanten; [gzl auf PyPI](https://pypi.org/project/gzl/)):

In 2D ist das Hexagonalgitter für ν ≲ 3,918 global optimal, darüber das Quadratgitter. **Korrektur aus der Nachprüfung:** Das Quadrat bleibt nur bis ν₂ ≈ 8,606 optimal, danach gewinnt ein Rechteckgitter (Abschnitt 7.2).

| ν | Globales Minimum | Hexagonal | Quadrat |
|---|---|---|---|
| 2,5 | hex | 20,6364817949 | 21,0272603756 |
| 3,0 | hex | 13,3917941102 | 13,6528937154 |
| 3,5 | hex | 9,7723761021 | 9,8932825030 |
| 3,75 | hex | 8,5685550844 | 8,6172422133 |
| 4,0 | Quadrat | 7,6061707526 | 7,5826699011 |
| 6,0 | Quadrat | 3,7457712638 | 3,2499743124 |

Der Schnittpunkt liegt bei **ν\* ≈ 3,91836** (Brent-Nullstelle). Bei ν\* steigt die Energie entlang des Randbogens |τ| = 1 von 7,898512 (hex, 60°) auf etwa **7,9096 bei ~72°** und fällt wieder auf 7,898512 (Quadrat, 90°). Es gibt also eine Barriere, und der Übergang ist **erster Ordnung**. Ein direkter Brute-Force-Vergleich auf ℤ² bei ν = 6 ergibt 3,2499743122 gegenüber 3,2499743124 aus GZL. Die Erklärung (eigene Abschätzung): Bei großem ν dominieren die kleinsten Dreiecke, und repulsive Dreikörperterme bestrafen die zwölf gleichseitigen Nächste-Nachbarn-Dreiecke des Hexagonalgitters.

In 3D (nur FCC, BCC und SC verglichen, ~28 s pro Satz) liegt **BCC bei allen getesteten ν vor FCC**:

| ν | FCC | BCC | SC |
|---|---|---|---|
| 3,5 | 38,1993 | **38,0744** | 41,9368 |
| 4,5 | 17,9624 | **17,7406** | 21,0981 |
| 6,0 | 8,0434 | **7,7102** | 10,2962 |
| 9,0 | 2,3133 | **1,9785** | 3,2723 |

**Was neu herauskäme.** Die erste Karte der Gitteroptima für reine Mehrkörperenergien über den gesamten Modulraum, mit drei Bausteinen:

- Phasendiagrammen in ν und im Kopplungsverhältnis Zwei- zu Dreikörper;
- einer belastbar formulierten Vermutung „BCC minimiert die Dreiecks-Zeta in 3D" als Mehrkörper-Gegenstück zu Sarnak–Strömbergsson;
- dem ersten dokumentierten Hexagonal-Quadrat-Übergang erster Ordnung für eine reine Dreikörperenergie.

**[abgeleitet]** Für ganzzahlige Exponenten gibt es beweisbare Ankerpunkte. Die 2D-Energien mit Exponenten (2,2,2) und (4,4,2) reduzieren sich über die Identitäten C₁,₁,₁ = E₃ + ζ(3) und C₂,₂,₁ = (2/5)E₅ + ζ(5)/30 auf Epstein-Zetas. Diese werden bekanntlich vom Hexagonalgitter minimiert (Hintergrundwissen, in den Notizen nicht belegt). Das passt zur Pilotaussage „hex für kleine ν" und eröffnet einen Beweisweg über modulare Graphfunktionen (Idee 3).

**Erster Meilenstein (2–4 Wochen).** Die Piloten mit Konvergenzstudie reproduzieren: n_points variieren, `zeta_circle` gegen `evaluate_graph` prüfen, Brute-Force-Vergleiche bei mehreren ν. Die 2D-Suche auf große Im τ ausdehnen und ν\* samt Barriere bestätigen.

**Ausbaupfad.**

- **Monat 2–3:** globale 3D-Suche im 5-dimensionalen Modulraum, zum Beispiel über reduzierte Gram-Matrizen, auf einem ν-Gitter von 3,2 bis 12. Ergänzend 4-Zyklus und K₄.
- **Monat 4–6:** gemischte Energien aus Zweikörper-Riesz bzw. Lennard-Jones und Dreieck mit Phasendiagramm. Damit direkt die offene Frage aus [arXiv:2504.07338](https://arxiv.org/abs/2504.07338) angehen: Ist BCC gegen *nicht-kuboidale* Verzerrungen stabil?
- **Monat 6–12:** ATM über Zerlegung in Dreiecks-Zetas. Dafür ist eine **GZL-Erweiterung nötig**: meromorphe Fortsetzung für Dreiecksblöcke mit ν ≤ d, siehe Risiken. Außerdem Beweisversuche für ganzzahlige s in 2D über MGF-Laplace-Gleichungen.

**Vorwissen und Einstieg.** Gittergeometrie (Fundamentaldomäne, Reduktionstheorie), globale Optimierung, Grundlagen der Epstein-Zeta. Einstieg über [arXiv:2609.17356](https://arxiv.org/abs/2609.17356), [arXiv:2312.01395](https://arxiv.org/abs/2312.01395) (Strukturübergänge auf Rechteckgittern) und [arXiv:2504.11989](https://arxiv.org/abs/2504.11989).

**Risiken.**

1. Neuheit ist nur über das Fehlen von Treffern gesichert. Die Bétermin-Gruppe (offene Fragen, zitiert in [arXiv:2411.17199](https://arxiv.org/abs/2411.17199)) und der Schwerdtfeger-Buchheit-Kreis könnten Ähnliches in Arbeit haben.
2. Die 2D-Suche deckte nur Im τ ≤ 2–3 mit sechs Startpunkten ab. In 3D wurden nur drei Gitter verglichen.
3. ATM enthält pro Kante Exponenten ν ∈ {5, 3, 1, −1}, also ν ≤ d. **Das lehnt GZL 1.0.0 für Nicht-Brücken ab.** Die Fortsetzung existiert für Zyklen im Buchheit-Busse-Verfahren, ist in GZL aber nicht freigelegt.
4. Die Zahlen sind inzwischen unabhängig nachgerechnet (Abschnitt 7), aber nicht begutachtet. Die 3D-Suche mit 16 Starts ist kein Beweis des globalen Minimums.

**Scoop-Risiko: mittel.** Die GZL-Roadmap nennt Gitteroptimierung nicht, aber der Autorenkreis arbeitet an ATM-Energien auf dem Bain-Pfad. Eine frühe Kontaktaufnahme mit Buchheit und Schwerdtfeger könnte aus dem Risiko eine Kooperation machen.

**Schlüsselpapers.**

| arXiv | Titel / Inhalt (laut Notizen) | Jahr | Relevanz |
|---|---|---|---|
| [2609.17356](https://arxiv.org/abs/2609.17356) | Luo & Wei, beanspruchter Beweis der 3D-Sarnak–Strömbergsson-Vermutung | 2026 | Zweikörper-Referenz: FCC optimal, Gegenpol zum BCC-Pilot |
| [2504.07338](https://arxiv.org/abs/2504.07338) | Robles-Navarro et al., LJ + ATM auf kuboidalen (Bain-)Gittern | 2025 | ATM begünstigt BCC. Offene Stabilitätsfrage außerhalb des Bain-Pfads |
| [2504.11989](https://arxiv.org/abs/2504.11989) | Buchheit & Busse, „Epstein zeta method for many-body lattice sums" | 2025 | Methodischer Vorläufer, ATM „from weeks to minutes" |
| [2407.20762](https://arxiv.org/abs/2407.20762) | Bétermin & Furlanetto, Minimierer für p-Normen | 2024 | Präzedenz für unerwartete Minimierer-Übergänge |
| [2312.01395](https://arxiv.org/abs/2312.01395) | Bétermin, Šamaj, Travěnec, Strukturübergänge auf 2D-Rechteckgittern | 2023 | Theorie von Übergängen erster und zweiter Ordnung, Vorlage für die Analyse |
| [2605.07580](https://arxiv.org/abs/2605.07580) | Luo & Wei, Minimierer von Verhältnissen und Differenzen von Theta- und Epstein-Funktionen | 2026 | Aktuelle Technik, zentrale Rolle des Hexagonalgitters |

**Bewertung.**

| Kriterium | Punkte | Begründung |
|---|---|---|
| Neuheit | 5 | Kein Treffer für Mehrkörper-Optimierung über den ganzen Modulraum. Der Pilot liefert zwei neue qualitative Befunde. |
| Passung | 5 | Läuft mit `zeta_circle` und `evaluate_graph` in v1.0.0, ν > d. Nur ATM braucht eine Erweiterung. |
| Hebel | 4 | Sekunden pro Gitterpunkt statt Doppelsummen mit langsamer Konvergenz. Globale Suche wird erst dadurch praktikabel. |
| Machbarkeit | 5 | Die Piloten existieren, der erste Meilenstein ist in Wochen erreichbar. |
| Wow | 5 | Heatmap über die Fundamentaldomäne, springendes Optimum, BCC schlägt FCC. |
| Anschluss | 4 | math.MG, math-ph und Kristallchemie. Mit ATM auch materialrelevant. |

### Idee 2: Klassische langreichweitige Hochtemperatur- und Perkolationsreihen für den Sak-Streit

**Pitch.** Die klassische Schwester der Quanten-Linked-Cluster-Entwicklung, für die GZL gebaut wurde: deterministische Hochtemperaturreihen für LR-Ising, LR-O(n) und LR-Perkolation. Sie liefern β_c(σ), ρ_c(σ) und χ(k) über die ganze Brillouin-Zone und geben im Sak-Streit eine Stimme, die nicht auf Monte Carlo beruht.

**Das Problem im Zielfeld.** Wo die Grenze zwischen lang- und kurzreichweitiger Universalität liegt, ist 2025–2026 offen umstritten **[belegt]**:

- Die Deng-Gruppe findet mit MC bis L = 8192 **σ\* = 2** für 2D LR-Ising, XY, Heisenberg und Perkolation, gegen Saks σ\* = 2 − η_SR = 7/4 ([arXiv:2512.04805](https://arxiv.org/abs/2512.04805), [arXiv:2512.01956](https://arxiv.org/abs/2512.01956)).
- Ein Quanten-RG-Papier beansprucht σ\* = 2 auch für Quanten-O(n)-Modelle ([arXiv:2606.22407](https://arxiv.org/abs/2606.22407)).
- In der 2D-LR-Perkolation bis L = 16384 weicht η nahe σ ≈ 3/2 von 2 − σ ab. Dort sind ρ_c(σ = 1) = 0,307591(4) und ρ_c(σ = 1/2) = 0,26263(2) bestimmt ([arXiv:2608.20750](https://arxiv.org/abs/2608.20750)).
- Hutchcroft nennt die Berechnung von α_c(d) für 2 < d < 6 „beyond the scope of current techniques" ([arXiv:2508.18808](https://arxiv.org/abs/2508.18808)).

Heute gibt es drei Zugänge: rigorose Schranken ohne Konstanten, MC auf endlichen Tori mit Minimum-Image-Abschneidung und Kontinuums-ε-Entwicklungen ([arXiv:2602.07818](https://arxiv.org/abs/2602.07818), [arXiv:2608.15120](https://arxiv.org/abs/2608.15120)). **Eine Reihe mit exakten Gitterkoeffizienten für klassische LR-Modelle wurde für 2025–26 nicht gefunden** **[belegt als Negativbefund]**.

**Abbildung auf GZL [abgeleitet].** In der Free-Embedding-Linked-Cluster-Entwicklung ist der Koeffizient von βⁿ in χ(k) oder ln Z eine Summe Σ_G w_G · Σ_x Π_e J(x_e)^{m_e} e^{ik·Δx}, mit J(x) = |x|^{−(d+σ)}, J(0) = 0 und Kantenvielfachheit m_e. Da J^m wieder ein Potenzgesetz mit ν = m(d+σ) > d ist, ist das **exakt eine Graph-Zeta mit GZLs Konvention K(0) = 0**. Die Vertexfaktoren (Kumulanten für Ising, O(n), φ⁴) bringt der Nutzer mit. Für LR-Perkolation liefert die Entwicklung 1 − e^{−βJ} = Σ_m (−1)^{m+1} β^m J^m/m! pro Ordnung endliche Summen von Potenzgesetzen, also genau GZLs Kernklasse ([arXiv:2508.18807](https://arxiv.org/abs/2508.18807)). Die Subgraph-Entwicklung verlangt allerdings Hardcore-Einbettungen und damit die Kontraktionsabbildung aus [arXiv:2609.18761](https://arxiv.org/html/2609.18761), die neu implementiert werden muss.

```python
# Pseudocode – Korpus und Vertexgewichte sind Eigenbau
for n in range(1, n_max + 1):
    for G, w_G, mult in classical_LR_corpus(model="ising", order=n):   # Graphenenumeration + Kumulanten
        nu_e = [m * (d + sigma) for m in mult]                           # J^m bleibt ein Potenzgesetz
        chi_k[n] += w_G * gzl.evaluate_graph(G.edges, nu_e, A, terminal=(s, t), n_points=N)  # ganze BZ per FFT
# Analyse: Ratio/Padé für beta_c(sigma), gamma(sigma); nichtanalytischer |k|^sigma-Koeffizient von 1/chi(k)
```

**Was neu herauskäme.**

- β_c(σ) und ρ_c(σ) als glatte Funktionen auf dichten σ-Gittern, bis nahe σ → 0, wo direkte Summation versagt.
- γ(σ) aus Ratio- und Padé-Analyse zum Vergleich mit ε-Entwicklungen und FRG ([arXiv:2510.02458](https://arxiv.org/abs/2510.02458)).
- Die **Amplitude des |k|^σ-Terms** von χ(k)^{−1}. Für σ < 2 divergiert das zweite Moment, sodass die Standard-Korrelationslänge der Hochtemperaturreihe gar nicht definiert ist. Das volle k-Gitter aus GZL ist dafür das natürliche Werkzeug **[abgeleitet]**.
- Als spekulative Erweiterung: die Edwards-Anderson-Suszeptibilität der 1D-LR-Spinglas-Kette mit J² ∝ r^{−2σ}, für das offene Fenster 3/2 ≤ α ≤ 2 ([arXiv:2604.07130](https://arxiv.org/abs/2604.07130)).

**Erster Meilenstein (6–10 Wochen).**

1. **Literaturschranke:** Prüfen, ob alte LR-Reihen (Nagle–Bonner 1970, Glumac–Uzelac um 1989; beide in den Notizen *nicht verifiziert*) schon hohe Ordnungen erreichten.
2. Eine eigene Graphenenumeration mit Ising-Kumulanten schreiben und χ(k = 0) für 1D bis etwa Ordnung 8 berechnen.
3. Niedrige Ordnungen gegen direkte Summen prüfen.

**Ausbaupfad.**

- **Monat 3–6:** 2D auf dem Quadrat- und Dreiecksgitter mit dichtem σ-Gitter, Vergleich mit den kritischen Punkten der Deng-Gruppe.
- **Monat 6–9:** χ(k) über die ganze Brillouin-Zone und die |k|^σ-Amplitude.
- **Monat 9–15:** LR-Perkolation mit Hardcore-Kontraktionen, Vergleich mit ρ_c = 0,307591(4). Dabei muss die Normierungskonvention angeglichen werden: Unendliches Gitter und Epstein-Zeta stehen gegen Torus mit Minimum-Image.

**Vorwissen und Einstieg.** Klassische Linked-Cluster-Technik nach Wortis, Englert, Lüscher–Weisz (Lehrbuchstoff, in den Notizen ohne arXiv-Quelle), Graphenenumeration, Reihenanalyse. Einstieg: [arXiv:2403.00421](https://arxiv.org/abs/2403.00421) (MC-Einbettungs-Review, zeigt die Quanten-Pipeline), [arXiv:2512.04805](https://arxiv.org/abs/2512.04805), [arXiv:2608.20750](https://arxiv.org/html/2608.20750).

**Risiken.**

1. Die Graphenzahl bei Ordnung 15–20 wird grob auf 10⁵–10⁷ geschätzt (*nicht verifiziert*). Die Baumweite hoher 1PI-Graphen kann 4–6 erreichen, und das Tensornetz kostet dann O(N^{tw}).
2. Nahe σ ≈ 2 machen Log-Korrekturen die Reihenanalyse unzuverlässig. GZL entfernt das Koeffizientenrauschen, nicht den Extrapolationsfehler begrenzter Ordnung.
3. Der Hardcore-Code fehlt in GZL.
4. Alte Literatur könnte die Neuheit schmälern.

**Scoop-Risiko: mittel.** Klassische LR-Modelle nennen die Autoren nicht. „Endliche Temperatur" mit der Gruppe Sbierski steht aber auf ihrer Roadmap, und die Sbierski-Gruppe arbeitet bereits mit Hochtemperaturentwicklungen ([arXiv:2604.22743](https://arxiv.org/abs/2604.22743)).

**Schlüsselpapers.**

| arXiv | Titel / Inhalt | Jahr | Relevanz |
|---|---|---|---|
| [2512.04805](https://arxiv.org/abs/2512.04805) | Xiao, Liu, Fan, Deng, 2D LR-Ising bis L = 8192, σ\* = 2 | 2025 | Hauptgegenstand des Sak-Streits |
| [2608.20750](https://arxiv.org/abs/2608.20750) | Liu, Xiao, Fan, Deng, 2D LR-Perkolation bis L = 16384 | 2026 | Präzise ρ_c-Benchmarks |
| [2606.22407](https://arxiv.org/abs/2606.22407) | „Perturbative Renormalization and Universality Diagram for Long-Range Quantum Criticality" | 2026 | Quanten-Gegenstück σ\* = 2 |
| [2508.18808](https://arxiv.org/abs/2508.18808) | Hutchcroft, „Critical LRP II" | 2025 | Offene Frage α_c(d), nur Aussagen bis auf Konstanten |
| [2604.07130](https://arxiv.org/abs/2604.07130) | Okuyama & Ohzeki, 1D LR-Spinglas auf der Nishimori-Linie | 2026 | Offenes Fenster 3/2 ≤ α ≤ 2 |
| [2403.00421](https://arxiv.org/abs/2403.00421) | „Monte Carlo based techniques for quantum magnets with long-range interactions" | 2024 | Die Pipeline, die GZL im Quantenfall ersetzt |

**Bewertung.**

| Kriterium | Punkte | Begründung |
|---|---|---|
| Neuheit | 4 | Keine aktuellen klassischen LR-Reihen gefunden. Alte Arbeiten sind ungeprüft. |
| Passung | 3 | `evaluate_graph` passt direkt. Korpus, Vertexgewichte und Hardcore-Abbildung sind Eigenbau. |
| Hebel | 4 | Deterministische, nicht-universelle Größen über ganze σ-Bereiche, dazu die ganze Brillouin-Zone |
| Machbarkeit | 2 | Graphenenumeration und Reihenanalyse sind Monate Arbeit. |
| Wow | 3 | β_c(σ)-Kurven gegen MC-Punkte, aber abstrakt |
| Anschluss | 5 | Mitten in einem aktiven Streit zwischen drei Communities |

### Idee 3: Modulare Graphfunktionen als duale Graph-Zetas, GZL als Numerik-Labor der Stringtheorie

**Pitch.** Modulare Graphfunktionen (MGFs) aus Genus-1-Stringamplituden sind für planare Graphen exakt GZL-Graph-Zetas des dualen Graphen auf dem Impulsgitter Λ_τ = ℤ + τℤ. Die Physik der langreichweitigen Quantenmagnete und die Welt der Stringamplituden teilen damit dasselbe numerische Objekt.

**Das Problem im Zielfeld.** MGFs sind „non-holomorphic modular functions associated with Feynman graphs for a conformal scalar field theory on a two-dimensional torus" **[belegt]** ([arXiv:1512.06779](https://arxiv.org/abs/1512.06779)). Der Stand der Werkzeuge:

- Alle Identitäten bis Gewicht 6 und alle dihedralen Identitäten bei Gewicht 7 sind bewiesen ([arXiv:1608.04393](https://arxiv.org/abs/1608.04393)).
- Ein Mathematica-Paket überführt MGFs in iterierte Eisenstein-Integrale, deckt aber nur Topologien **bis vier Vertices** ab ([arXiv:2502.05531](https://arxiv.org/abs/2502.05531)).
- Elliptische MGFs werden über äquivariante iterierte Integrale gelöst ([arXiv:2511.15883](https://arxiv.org/abs/2511.15883)).
- Objekte mit komplexen Indizes tauchen neu auf ([arXiv:2512.21413](https://arxiv.org/abs/2512.21413)).

Numerische Werte bei beliebigem innerem τ für größere planare Graphen und für nicht-ganzzahlige Exponenten gibt es als Werkzeug nicht allgemein.

**Abbildung auf GZL [abgeleitet, numerisch gestützt].** Ein MGF ist C_Γ(τ) = Σ_{p_e ∈ Λ_τ∖0} Π_v δ(Σ ±p_e) Π_e (τ₂/(π|p_e|²))^{a_e}. Das Impulserhaltungsgitter ist das Zyklengitter von Γ. Nach Whitney lässt es sich genau dann als Differenzen von Vertexvariablen schreiben, wenn Γ planar ist, und diese Vertices sind die Flächen von Γ, also die Knoten des Duals Γ\*. Daraus folgt:

C_Γ(τ) = (τ₂/π)^{Σa_e} · ζ_{Λ_τ, Γ\*}(k = 0) mit ν_e = 2a_e.

Die Zuordnung im Einzelnen:

- Dihedrale („Banana"-)MGFs werden zu Kreis-Zetas (`zeta_circle`).
- Trihedrale werden zu Theta-artigen SP-Graphen.
- Der selbstduale Tetraeder K₄ läuft über das Tensornetz.
- Die Kante mit Impuls, Z_{Λ_τ,2a}(k), ist der Baustein aus Green-Funktion bzw. Kronecker-Eisenstein-Reihe.
- GZLs externes k entspräche dem Torus-Punkt z elliptischer MGFs **[spekulativ, nur im dihedralen Fall plausibilisiert]**.

```python
A = np.array([[1, tau.real], [0, tau.imag]])                 # |A(m,n)| = |m + n*tau|
C_abc = (tau.imag/np.pi)**(a+b+c) * gzl.zeta_circle([2*a, 2*b, 2*c], A)   # dihedral, 2 Vertices
C_K4  = (tau.imag/np.pi)**12 * gzl.evaluate_graph(K4, np.full(6, 4.0), A, n_points=64)  # tetraedrisch
```

**Pilotergebnisse [Pilot]** (`mgf_test.py`, `mgf_test2.py`). Die Nachprüfung ohne GZL bestätigt beide Identitäten auf 10⁻¹⁴ bis 10⁻¹⁵ (Abschnitt 7.4):

- **C₂,₂,₁₊ε/₂ → (2/5)E₅ + ζ(5)/30.** Nach Richardson-Extrapolation in ε liegt der relative Fehler bei **1,1–1,2×10⁻⁸** an allen vier getesteten τ (i, e^{iπ/3}, 0,21 + 1,37i, −0,4 + 0,95i).
- Die **Zagier-Identität C₁,₁,₁ = E₃ + ζ(3)** ist auf etwa 9×10⁻⁷ reproduziert.
- Jede Kreis-Auswertung braucht ~0,1 s.
- Die tetraedrische MGF (alle ν = 4) bei τ = i konvergiert mit n_points = 32/48/64 gegen 2,29641310909 / …926 / …927, in 0,1–0,3 s.
- Die Zuschreibung von C₂,₂,₁ an D'Hoker–Green–Vanhove ist *nicht verifiziert*.

**Was neu herauskäme.**

- Double-Precision-Werte planarer MGFs mit **≥ 5 Vertices** an beliebigem innerem τ, jenseits des 4-Vertex-Pakets.
- **Verallgemeinerte MGFs mit reellen Exponenten** C_{s₁,s₂,s₃}(τ).
- Elliptische MGFs auf einem vollen z-Gitter per FFT (spekulativ).
- Ein neuer Satz über modulare Funktionen: Mit dem Wörterbuch aus Idee 1 wandert **das Minimum von C_{s,s,s}(τ) über der Fundamentaldomäne bei s\* = 1,95918 vom hexagonalen zum quadratischen Punkt und für s > 4,303 zu einem Rechteckpunkt** [nachgeprüft, Abschnitt 7]. Nach den Notizen hat diese Frage noch niemand gestellt.

**Erster Meilenstein (3–5 Wochen).** Die dihedralen, trihedralen und tetraedrischen Identitäten an vielen τ systematisch reproduzieren. Die inhomogenen Laplace-Gleichungen numerisch per finiter Differenzen in τ prüfen.

**Ausbaupfad.**

- **Monat 2–4:** Tabellen von 5- und 6-Vertex-MGFs. Der Vergleich mit [arXiv:2502.05531](https://arxiv.org/abs/2502.05531) ist nur bis 4 Vertices möglich.
- **Monat 4–6:** C_{s,s,s} für reelle s, Minimierungsstudie und Anschluss an Idee 1.
- **Monat 6–12:** elliptische MGFs gegen publizierte Werte prüfen. Optional eine Multiprecision-Portierung (mpmath oder Arb) der Kreis-Zeta für PSLQ.

**Vorwissen und Einstieg.** Modulformen, Eisenstein-Reihen, Grundlagen der Genus-1-Stringamplituden, Planarität und Dualität von Graphen. Einstieg: [arXiv:1512.06779](https://arxiv.org/abs/1512.06779), [arXiv:1608.04393](https://arxiv.org/abs/1608.04393), [arXiv:2502.05531](https://arxiv.org/abs/2502.05531).

**Risiken.**

1. Die MGF-Community hat exakte Methoden. GZL erschließt dort **kein blockiertes Problem**, sondern liefert Kreuzchecks und Erweiterungen.
2. Die meisten stringrelevanten MGFs haben Kanten mit a_e = 1, also ν = d = 2. Das liegt außerhalb der GZL-Domäne und ist nur per ε-Extrapolation mit ~10⁻⁸ erreichbar.
3. Nicht-planare MGFs (K₃,₃-Topologien) sind keine GZL-Objekte.
4. Modulare Graph*formen* mit a ≠ b brauchen anisotrope Kerne.
5. Mit doppelter Genauigkeit ist PSLQ nur für winzige Basen möglich.

**Scoop-Risiko: gering.** Kein Bezug zur Roadmap der Autoren.

**Schlüsselpapers.**

| arXiv | Titel / Inhalt | Jahr | Relevanz |
|---|---|---|---|
| [1512.06779](https://arxiv.org/abs/1512.06779) | Modulare Graphfunktionen und single-valued elliptische Polylogarithmen | 2015 | Definition. Impuls an nur zwei Vertices entspricht GZLs Zweiterminal-Struktur |
| [1608.04393](https://arxiv.org/abs/1608.04393) | D'Hoker & Kaidi, alle Identitäten bis Gewicht 6 | 2016 | Prüfstein für die Numerik |
| [2502.05531](https://arxiv.org/abs/2502.05531) | Claasen & Doroudiani, iterierte Eisenstein-Integrale, Paket bis 4 Vertices | 2025 | Stand der Technik, Grenze bei 4 Vertices |
| [2511.15883](https://arxiv.org/abs/2511.15883) | Schlotterer, Sohnle, Tao, elliptische MGFs über äquivariante iterierte Integrale | 2025 | Ziel für k-aufgelöste GZL-Werte |
| [2603.20550](https://arxiv.org/abs/2603.20550) | Dobrowolski, symmetrisierte Mordell–Tornheim-Zeta | 2026 | 1D-Kreis-Zeta = ζ̄_n, Validierungsbank |
| [2512.21413](https://arxiv.org/abs/2512.21413) | Fedosova & Klinger-Logan, Divisor-Summen mit komplexen Indizes | 2025 | Motivation für nicht-ganzzahlige Exponenten |

**Bewertung.**

| Kriterium | Punkte | Begründung |
|---|---|---|
| Neuheit | 5 | Das Wörterbuch MGF ↔ duale Graph-Zeta ist in den Quellen nirgends ausgesprochen. |
| Passung | 3 | Planar und ν > d funktioniert. Die typischen a = 1-Kanten sind Grenzfälle. |
| Hebel | 2 | Exakte Konkurrenzmethoden existieren. |
| Machbarkeit | 4 | Die Piloten laufen, das Wörterbuch ist einfach. |
| Wow | 4 | „Magnetismus-Software berechnet Stringamplituden" ist eine starke Geschichte. |
| Anschluss | 3 | Interesse als Kreuzcheck. Neue s\*-Aussage für math.NT |

### Idee 4: Simulator-Material-Zwillinge mit vollem Schwanz – Rydberg-vdW-Arrays und Tm-Dreiecksmagnete

**Pitch.** Eine Referenz im thermodynamischen Limes für Rydberg-Quantenprozessoren mit 256 Qubits: volle 1qp-Dispersion, Gap und kritischer Punkt mit dem *kompletten* 1/r⁶-Schwanz. Dazu wird exakt beziffert, wie stark der Schwanz die Abbildung „Simulator ↔ Material" für TmMgGaO₄ und NaTmSe₂ verschiebt.

**Das Problem im Zielfeld [belegt].**

- Pasqal fährt das 2D-TFIM auf 16×16 Atomen mit U(R/r)⁶. Den kritischen Punkt h\*/J ≈ 2,6 für L = 16 übernimmt die Arbeit aus QMC. Tensornetze verlieren bei tJ ≳ 1 die Kontrolle, und Spektren werden nicht gezeigt ([arXiv:2608.07178](https://arxiv.org/abs/2608.07178)).
- Die „one-to-one"-Simulation von TmMgGaO₄ auf einem Dreiecks-Rhombus mit 256 Qubits bildet ein J₁–J₂-Modell mit J₂ ≈ 0,05 J₁ ab. Sie behauptet nur „almost the same neighbour dependence". Die MPS-Vergleichsrechnungen dauerten etwa zwei Wochen gegenüber etwa einem Tag auf dem QPU ([arXiv:2603.20372](https://arxiv.org/abs/2603.20372)).
- NaTmSe₂ realisiert das TFIM mit quantitativ bestimmten Austauschparametern und Neutronendaten ([arXiv:2505.09884](https://arxiv.org/abs/2505.09884)).
- Für KTmSe₂ zeigte die GZL-Physik-Arbeit, dass „J₁ + voller Dipol" die INS-Dispersion am besten beschreibt ([arXiv:2609.18761](https://arxiv.org/html/2609.18761), Daten aus [arXiv:2306.03544](https://arxiv.org/abs/2306.03544)).
- Keines dieser Experimente wird mit einer konvergierten Rechnung im thermodynamischen Limes über die volle Brillouin-Zone mit komplettem Schwanz verglichen.

**Abbildung auf GZL [abgeleitet].**

- **Gitter:** `triangular` (auch `square` und `chain`).
- **Korpus:** `tfim1qp` bzw. `tfim0qp`, AFM-Vorzeichen.
- **Kern:** Interaction = a(x) + C₆|x|^{−6} für Rydberg. Für Tm-Materialien mit Momenten ∥ c ist die Dipol-zz-Kopplung in der Ebene exakt isotrop, also a(x) = J₁, J₂ kompakt plus D|x|^{−3}.
- **Observablen:** 1qp-Dispersion ω(k) über die ganze Brillouin-Zone, Gap am K-Punkt, λ_c und zν per DlogPadé.
- Der vdW-Schwanz liefert auf dem Dreiecksgitter J₂/J₁ = 1/27 ≈ 0,037 und J₃/J₁ = 1/64. Das Material hat J₂ ≈ 0,05 J₁ und kein J₃. Die Differenz lässt sich genau als kompakte Korrektur a(x) einspeisen.

```bash
# Rydberg-Referenz (Z2-symmetrische Detuning-Linie), Pseudokonfiguration
gzl series --corpus tfim1qp --A triangular --nu 6 --n-points 48 --order-max 11
# Materialvariante: Interaction.from_shells(...)  (J1, Delta J2 kompakt) + b*|x|^-3 bzw. |x|^-6
```

**Was neu herauskäme.**

- Die ersten Spektren im thermodynamischen Limes mit vollem vdW-Schwanz für die Gitter der 256-Qubit-Maschinen.
- Eine quantifizierte „Schwanzkorrektur" der TmMgGaO₄-Abbildung: Um wie viel verschieben 1/27 statt 0,05 und der J₃-Schwanz Gap und kritisches Feld?
- Eine Modellwahl J₁–J₂ gegen J₁ + Dipol für NaTmSe₂ nach dem KTmSe₂-Muster, samt Abstand zum Quantenkritischen Punkt.

**Erster Meilenstein (3–6 Wochen).**

1. Das 2D-KTmSe₂-Ergebnis aus [arXiv:2609.18761](https://arxiv.org/html/2609.18761) reproduzieren (Validierung).
2. Reines vdW-Modell gegen J₁–J₂-Modell auf dem Dreiecksgitter vergleichen: Dispersion, Gap und λ_c.

**Ausbaupfad.**

- **Monat 2–4:** Vergleich mit Pasqal (h\*/J ≈ 2,6) nach Klärung der Vorzeichen- und Untergitter-Konvention. Die Notizen markieren diesen Punkt als *nicht verifiziert*.
- **Monat 4–6:** NaTmSe₂ und 1D-Rydberg-Ketten für die Caltech-Messungen des dynamischen Strukturfaktors ([arXiv:2601.16275](https://arxiv.org/abs/2601.16275)).
- **Später:** 3D-Stapelung und Spektralgewichte, sobald die Roadmap-Features vorliegen.

**Vorwissen und Einstieg.** pCUT und Linked-Cluster-Grundlagen, DlogPadé, Rydberg-Hamiltonian, Kristallfeld-Physik nicht-Kramers-Ionen. Einstieg: [arXiv:2609.18761](https://arxiv.org/html/2609.18761), [arXiv:2603.20372](https://arxiv.org/abs/2603.20372).

**Risiken.**

1. Der mitgelieferte Korpus gilt nur auf der **Z₂-symmetrischen Linie** δ = Σ_j U_ij/2. Ein Longitudinalfeld, also das volle (Ω, δ)-Diagramm, braucht einen neuen Korpus.
2. pCUT aus dem Transversalfeld-Limes gilt nur in der paramagnetischen Phase. TmMgGaO₄ ordnet bei tiefer Temperatur.
3. Ordnung 11 begrenzt die Präzision nahe λ_c.
4. Die Stapelung in 3D braucht anisotrope Dipolkerne.

**Scoop-Risiko: hoch.** Die Autoren nennen Rydberg-vdW-LRTFIM und DSF-Vergleiche ausdrücklich als nächstes Ziel und haben KTmSe₂ selbst bearbeitet. Empfehlung: Kooperation suchen oder sich auf das konkrete Paar TmMgGaO₄-Simulator und NaTmSe₂ konzentrieren.

**Schlüsselpapers.**

| arXiv | Titel / Inhalt | Jahr | Relevanz |
|---|---|---|---|
| [2603.20372](https://arxiv.org/abs/2603.20372) | 256-Qubit-„one-to-one"-Simulation von TmMgGaO₄ | 2026 | Zielexperiment für die Schwanzkorrektur |
| [2608.07178](https://arxiv.org/abs/2608.07178) | Pasqal/Browaeys, 2D-TFIM auf 16×16 mit vdW | 2026 | Kritischer Punkt aus QMC bei L = 16 |
| [2505.09884](https://arxiv.org/abs/2505.09884) | NaTmSe₂ als TFIM-Realisierung | 2025 | Nächstes Material nach KTmSe₂ |
| [2306.03544](https://arxiv.org/abs/2306.03544) | KTmSe₂-INS, dispersiver Kristallfeldzweig | 2023 | Validierungsdaten |
| [2609.18761](https://arxiv.org/abs/2609.18761) | Duft et al., „Exact and fast series expansions for quantum models with long-range interactions" | 2026 | Methode und KTmSe₂-Präzedenz |
| [2601.16275](https://arxiv.org/abs/2601.16275) | Caltech, DSF einer Rydberg-Kette am kritischen Punkt | 2026 | 1D-Vergleichsziel |

**Bewertung.**

| Kriterium | Punkte | Begründung |
|---|---|---|
| Neuheit | 3 | Neue Systeme, aber dieselbe Methode wie bei den Autoren |
| Passung | 4 | Läuft heute, auf der Z₂-Linie und im 2D-Dipolmodell |
| Hebel | 4 | Voller Schwanz und volle Brillouin-Zone im thermodynamischen Limes statt DMRG/QMC bei L = 16 |
| Machbarkeit | 4 | Mitgelieferte Korpora, das Kern-API reicht aus. |
| Wow | 5 | Direkter Vergleich mit Experiment und Material |
| Anschluss | 4 | Relevant für Rydberg- und Neutronen-Community. Kooperation ist naheliegend. |

### Idee 5: Quanten-Sak-Test und Log-Korrekturen im LRTFIM für d = 1, 2, 3

**Pitch.** Mit rauschfreien Reihen auf einem sehr dichten σ-Gitter prüfen, ob die Exponenten des LRTFIM ihren kurzreichweitigen Wert bei σ = 2 − η_SR (Sak) oder erst bei σ = 2 erreichen, wie es das RG-Papier von 2026 behauptet. Dazu die multiplikativen Logarithmen an der oberen kritischen Grenze σ_uc = 2d/3 und im 3D-Kubikgitter vermessen.

**Das Problem im Zielfeld [belegt].** [arXiv:2606.22407](https://arxiv.org/abs/2606.22407) schlägt ein (d, σ)-Universalitätsdiagramm für Quanten-O(n) vor, mit σ\* = 2 und Vorhersagen für ν, η_ω und η_k. Klassisch behauptet die Deng-Gruppe dasselbe per MC ([arXiv:2512.04805](https://arxiv.org/abs/2512.04805)). Für die antiferromagnetische LR-Ising-Kette weichen Korrelationsexponent und zentrale Ladung laut NQS bei kleinem Zerfallsexponenten ab, und der Zusammenbruch der konformen Invarianz ist „yet to be determined" ([arXiv:2308.09709](https://arxiv.org/abs/2308.09709)). Frühere pCUT+MC-Reihen waren verrauscht und nur dünn in σ abgetastet ([arXiv:1802.06684](https://arxiv.org/abs/1802.06684), [arXiv:2203.08081](https://arxiv.org/abs/2203.08081)).

**Abbildung auf GZL.** Direkt mit den mitgelieferten Korpora **[belegt]**: ν = d + σ auf dem Gitter `chain`, `square`, `triangular` oder `cubic`, FM- und AFM-Vorzeichen, 1qp-Gap als Observable. Die Reihenanalyse umfasst DlogPadé und Ratio-Methoden mit Log-Termen.

```bash
for sigma in $(seq 1.50 0.005 2.20); do
  gzl series --corpus tfim1qp --A chain --nu $(echo "1+$sigma" | bc) --order-max 11 ...
done   # -> z*nu(sigma); d=1: Sak 7/4 vs 2; d=2: ~1,96 vs 2
```

**Was neu herauskäme.**

- Deterministische Kurven zν(σ) in d = 1, 2, 3, direkt gegen die Vorhersagen aus [arXiv:2606.22407](https://arxiv.org/abs/2606.22407) gestellt.
- Die erste quantitative Extraktion multiplikativer Log-Exponenten aus exakten Koeffizienten in d = 2 und 3. In d = 3 gehört die kurzreichweitige Klasse zu 4D-Ising und trägt Logarithmen im gesamten SR-Regime **[abgeleitet]**.

**Erster Meilenstein (2–4 Wochen).** Zuerst den Crossover-Abschnitt von [arXiv:2609.18761](https://arxiv.org/html/2609.18761) lesen. Die Notizen konnten dessen Zahlen nicht extrahieren, und davon hängt die Neuheit ab. Danach zν(σ) für die FM-Kette auf dem Gitter σ ∈ [1,5; 2,2] mit Schrittweite 0,005 berechnen.

**Ausbaupfad.** Zuerst d = 2 (Quadrat, Dreieck AFM), dann das 3D-Kubikgitter, dann Log-bewusste Extrapolation. Langfristig O(2) und O(3) über XXZ-Korpora, die heute fehlen.

**Vorwissen und Einstieg.** Kritische Phänomene langreichweitiger Modelle, Reihenanalyse. Einstieg: [arXiv:2203.08081](https://arxiv.org/abs/2203.08081), [arXiv:2606.22407](https://arxiv.org/abs/2606.22407).

**Risiken.**

1. In d = 2 liegen die beiden Kandidaten nur **~0,04** auseinander. Bei Ordnung 11 und langsamen Crossovers begrenzt der Extrapolationsfehler, nicht das Rauschen.
2. Die Autoren haben σ bereits dicht gescannt. Die Neuheit muss aus der expliziten Konfrontation mit RG und MC und aus der Log-Analyse kommen.

**Scoop-Risiko: mittel–hoch.** Die GZL-Arbeit charakterisiert „the crossover between universality regimes" bereits.

**Schlüsselpapers.**

| arXiv | Titel / Inhalt | Jahr | Relevanz |
|---|---|---|---|
| [2606.22407](https://arxiv.org/abs/2606.22407) | „Perturbative Renormalization and Universality Diagram for Long-Range Quantum Criticality" | 2026 | Zu prüfende Vorhersage σ\* = 2 |
| [2203.08081](https://arxiv.org/abs/2203.08081) | „Scaling at quantum phase transitions above the upper critical dimension" | 2022 | Anker für Exponenten und Log-Korrekturen |
| [1802.06684](https://arxiv.org/abs/1802.06684) | „Quantum criticality of two-dimensional quantum magnets with long-range interactions" | 2019 | MC-basierte Vorgängerreihen |
| [2512.04805](https://arxiv.org/abs/2512.04805) | Klassisch σ\* = 2 per MC | 2025 | Klassisches Gegenstück |
| [2308.09709](https://arxiv.org/abs/2308.09709) | NQS für die AF-LR-Ising-Kette | 2023 | Dokumentierte Uneinigkeit |
| [2609.18761](https://arxiv.org/abs/2609.18761) | GZL-Physik-Arbeit | 2026 | Ausgangspunkt, Überlappung prüfen |

**Bewertung.**

| Kriterium | Punkte | Begründung |
|---|---|---|
| Neuheit | 3 | Die Überlappung mit dem σ-Scan der Autoren ist wahrscheinlich. |
| Passung | 5 | Nur mitgelieferte Korpora |
| Hebel | 3 | Das Rauschen verschwindet, das Ordnungslimit bleibt. |
| Machbarkeit | 4 | Rechenzeit in Minuten, die Analyse ist die Arbeit. |
| Wow | 3 | Exponentenkurven |
| Anschluss | 4 | Aktiver Streit |

---

## Weitere Ideen: vierzehn Kandidaten mit ehrlichem Erweiterungsbedarf

N = Neuheit, P = Passung, H = Hebel, M = Machbarkeit, W = Wow, A = Anschlussfähigkeit (je 1–5).

| Idee | Spur | Kern der Abbildung | Was fehlt | N | P | H | M | W | A | Scoop |
|---|---|---|---|---|---|---|---|---|---|---|
| ATM-Elastizitätskonstanten, Cauchy-Verletzung C₁₂ − C₄₄, Phonon-Zweitmoment Tr D(q) mit Dreikörperkräften ([2504.07338](https://arxiv.org/abs/2504.07338), [2012.05413](https://arxiv.org/abs/2012.05413)) | B | Summe von Dreiecks-Zetas, Ableitung nach A, Bloch-Phase auf einer Kante | Fortsetzung für ν ≤ d in Dreiecksblöcken, nur Bravais (kein hcp) | 4 | 2 | 4 | 3 | 3 | 4 | mittel |
| Vierkörper-Dispersion (Dipolring) | B | 4-Zyklus bzw. K₄ | Diagonalen mit positiven Potenzen, anisotrope Kontraktion | 4 | 1 | 3 | 2 | 3 | 3 | mittel |
| 1D-Rydberg-Kette am Ising-Punkt für Caltech-CFT-Spektren ([2601.16275](https://arxiv.org/abs/2601.16275)) | A | `chain`, ν = 6, AFM | nichts auf der Z₂-Linie; tricritisch und E8 brauchen neue Korpora | 3 | 4 | 3 | 4 | 4 | 4 | hoch |
| Langreichweitiger Benchmark für NLCE+QA auf Ionen-QPUs ([2605.28599](https://arxiv.org/abs/2605.28599)) | A | LR-Referenzdispersion für α = 1–3 | nichts | 3 | 5 | 3 | 5 | 2 | 3 | mittel |
| Bilayer- und Dimer-Kerne (r² + d²)^{−α/2} als asymptotische Potenzreihe plus a(x) ([2408.13145](https://arxiv.org/abs/2408.13145), [2605.13969](https://arxiv.org/abs/2605.13969)) | A | Kerntrick ohne mehratomige Basis | Dimer-Korpus, Konvergenzprüfung | 4 | 3 | 3 | 3 | 2 | 4 | hoch |
| α_c der LR-XY-Kette (LSW 3,0 gegen ≈ 2,7) und der Spin-1-Kette (QMC 2,49(1)) ([2601.20058](https://arxiv.org/abs/2601.20058), [2604.20831](https://arxiv.org/abs/2604.20831), [2604.12754](https://arxiv.org/abs/2604.12754)) | A | Reihen über α | XXZ- und Spin-1-Korpora, gestaffelte Vorzeichen nicht als Kern darstellbar | 2 | 2 | 4 | 2 | 3 | 4 | hoch |
| Dipolares XY/XXZ auf dem Quadratgitter, AFM-Dämpfungsrätsel ([2311.11726](https://arxiv.org/abs/2311.11726), [2605.07685](https://arxiv.org/abs/2605.07685)) | A | Isotropes 1/r³ bei senkrechter Quantisierungsachse | XXZ-Korpus, gapped Entwicklungslimes | 3 | 2 | 4 | 2 | 4 | 4 | hoch |
| 1/S-Spinwellen-Schleifen als Brillouin-Zonen-Faltungen über die Zeta-Algebra | A/B | `graph_convolve`-Algebra auf J(k) | Methodenentwicklung [spekulativ] | 4 | 2 | 3 | 2 | 3 | 3 | gering |
| Magnon-Bindungszustände im 2D-LRTFIM ([2512.09037](https://arxiv.org/abs/2512.09037)) | A | 2qp-Sektor | zwei Impulse, nicht in v1.0.0 | 3 | 1 | 4 | 1 | 4 | 4 | hoch |
| EA-Hochtemperaturreihe des 1D/2D-LR-Spinglases ([2604.07130](https://arxiv.org/abs/2604.07130)) | B | gerade Graphen mit J² ∝ r^{−2σ} | eigener Korpus [spekulativ] | 4 | 3 | 3 | 2 | 2 | 3 | gering |
| Lévy-Flug-Konstanten: P_n(0) = Kreis-Zeta/Zⁿ, Torus-Plateaus für FSS oberhalb d_c ([2412.08814](https://arxiv.org/abs/2412.08814)) | B | Kreis-Zeta bzw. D̂(k) auf dem Gitter | singuläres BZ-Integral, meist genügt EpsteinLib | 2 | 3 | 2 | 4 | 2 | 2 | gering |
| Multipräzisions-Graph-Zeta für PSLQ | B | Port des BZ-Integrals nach mpmath/Arb | Neuimplementierung | 4 | 1 | 3 | 2 | 3 | 3 | gering |
| 3D-skalare Gitterstörungstheorie, Vakuumgraphen | B | kompakter Kern plus Potenzschwanz wie bei Lüscher–Weisz | anisotrope Subleading-Terme, IR-Regulator (ν = 1 < d) | 3 | 1 | 2 | 1 | 2 | 2 | gering |
| Log-Korrekturen im Dicke-Ising-Kontext ([2605.27484](https://arxiv.org/abs/2605.27484)) | A | nur der LR-Ising-Anteil | Photonmode ist all-to-all und keine Gittersumme | 2 | 2 | 2 | 3 | 2 | 3 | hoch |

---

## Überraschende Zusammenhänge: eine Zahl, drei Communities

Der auffälligste Befund dieser Recherche ist eine **Drei-Welten-Identität** **[abgeleitet, durch Pilot gestützt]**. Dieselbe Zahl, die Dreiecks-Graph-Zeta auf einem 2D-Gitter mit Einheitsvolumen, erscheint an drei Stellen:

- als **Dreikörper-Riesz-Energie** in der Kristallchemie;
- als **zweischleifige dihedrale MGF** π^{3s}·C_{s,s,s}(τ) in der Stringtheorie;
- als **Einbettungssumme der Zyklen dritter Ordnung** in der Linked-Cluster-Reihe eines langreichweitigen Quantenmagneten.

Dahinter steckt eine strukturelle Übereinstimmung: GZLs Darstellung der Kreis-Zeta als Brillouin-Zonen-Integral über Produkte von Epstein-Zetas ist genau die „position-space"-Darstellung von MGFs. Die Torus-Green-Funktion G_a(z|τ) ist eine Epstein-Zeta am Wellenvektor z ([arXiv:2609.18918](https://arxiv.org/abs/2609.18918), [arXiv:1512.06779](https://arxiv.org/abs/1512.06779)). In 1D wird dasselbe Objekt zur symmetrisierten Mordell–Tornheim-Zeta aus der Zahlentheorie ([arXiv:2603.20550](https://arxiv.org/abs/2603.20550)). Diese Brücke ist exakt, aber wenig ergiebig, weil 1D bereits über MZVs verstanden ist.

Der **Sak-Streit wird gleichzeitig in vier Communities ausgetragen**, und keine nutzt Reihen mit exakten Gitterkoeffizienten:

- klassisches MC ([arXiv:2512.04805](https://arxiv.org/abs/2512.04805));
- Perkolation ([arXiv:2608.20750](https://arxiv.org/abs/2608.20750));
- Kontinuums-ε-Entwicklungen ([arXiv:2602.07818](https://arxiv.org/abs/2602.07818), [arXiv:2608.15120](https://arxiv.org/abs/2608.15120));
- Quanten-RG ([arXiv:2606.22407](https://arxiv.org/abs/2606.22407)).

Klassische Hochtemperaturreihe und Quanten-Linked-Cluster-Entwicklung teilen dieselben Gittersummen. GZL ist deshalb der gemeinsame Nenner, der die klassische und die Quantenseite mit *einem* Werkzeug prüfen könnte (Ideen 2 und 5).

Zwei weniger offensichtliche Identitäten gelten exakt **[abgeleitet]**. Erstens ist die Rückkehrwahrscheinlichkeit eines Lévy-Flugs auf dem Gitter nach n Schritten P_n(0) = ζ^{(n)}_Λ(ν,…,ν)/Z_{Λ,ν}(0)ⁿ, also eine Kreis-Zeta. Zweitens ist das phononische Zweitmoment Tr D(q) eines Bravais-Kristalls mit ATM-Kräften eine Summe von Dreiecks-Zetas mit Bloch-Phase auf einer Kante. Das ist genau GZLs Modus mit externem Impuls, und damit wird die Gitterdynamik der Kristallchemie zu einem Anwendungsfall der „ganzen Brillouin-Zone per FFT". Spekulativ ist die Verbindung zu Feynman-Perioden: Das Residuum einer Graph-Zeta am Gesamtpol ν|E| = d(|V|−1) sollte die Ortsraum-Periode sein. Das verlangt allerdings ν < d pro Kante und liegt damit außerhalb von GZL.

Zwei aktuelle Arbeiten liefern schließlich die **Theorie dafür, warum der volle Schwanz zählt**, und damit ein Verkaufsargument für alle Ideen der Spur A. Das Abschneiden eines r^{−p}-Schwanzes bei Reichweite R verschiebt lokale Observablen um O(R^{−(p−d)}), und diese Schranke ist optimal ([arXiv:2608.15576](https://arxiv.org/abs/2608.15576)). Surrogate aus Exponentialsummen in iMPS verzerren kritische Diagnosen ([arXiv:2606.20522](https://arxiv.org/abs/2606.20522)). Die DMRG-Studie dipolarer XY-Modelle schneidet bei R_max = W/2 ab ([arXiv:2605.07685](https://arxiv.org/html/2605.07685)).

---

## Verworfene Ideen: wo die Brücke nur mit Gewalt hält

| Idee | Warum sie nicht trägt |
|---|---|
| Photonik, Atomarrays, kollektive Lamb-Verschiebung ([2512.24596](https://arxiv.org/abs/2512.24596), [2609.17313](https://arxiv.org/abs/2609.17313)) | Der Kern ist die komplexe, oszillierende dyadische Green-Funktion e^{ik₀r}/r. GZLs k ist nur eine Bloch-Phase auf einem reellen geraden Kern. Die Bandstruktur ist ohnehin eine Einfachsumme. |
| Lüscher-Zeta und ihre LR-Modifikation ([2507.18399](https://arxiv.org/abs/2507.18399)) | Polkern 1/(p² − q²) mit winkelabhängigen Vertexgewichten, kein Differenzkern |
| QED_L-Koeffizienten, Madelung-Konstanten ([1702.01296](https://arxiv.org/abs/1702.01296), [2608.10041](https://arxiv.org/abs/2608.10041), [2503.00977](https://arxiv.org/abs/2503.00977)) | Einfachsummen. Mehrschleifen-Punktladungssummen faktorisieren. Das richtige Werkzeug ist EpsteinLib, nicht GZL. |
| Ortsraum-Photonpropagator in Gitter-QCD+QED ([2603.13086](https://arxiv.org/abs/2603.13086)) | Die Vertexfunktionen sind MC-Korrelatordaten, keine analytischen Kerne. Außerdem d = 4. |
| Langreichweitige ML-Potenziale ([2606.06617](https://arxiv.org/abs/2606.06617), [2512.18029](https://arxiv.org/abs/2512.18029)) | Nur Zweikörper-Ewald-Summen. Kein ML-Ansatz braucht nicht-baumartige Graphsummen. |
| Many-Body-Dispersion (MBD, [2308.03140](https://arxiv.org/abs/2308.03140)) | Auf Bravais-Gittern ein 3×3-T(k) pro k: anisotrope *Einfach*summen. Reale Systeme sind mehratomig. |
| LR-Kitaev-Ketten ([2412.01076](https://arxiv.org/abs/2412.01076), [2609.28735](https://arxiv.org/abs/2609.28735)) | Freie Fermionen, nur eindimensionale Polylog-Transformierte |
| Perkolations-Dreiecksdiagramm nahe s ≈ d | Der Kern wäre das unbekannte kritische τ(x) ≍ ‖x‖^{−d+α}, das kein Eingabekern und nicht in ℓ¹ ist ([2404.07276](https://arxiv.org/abs/2404.07276)). |
| k-Punkt-Baumsummen der Hochdimensions-Perkolation ([2607.22387](https://arxiv.org/abs/2607.22387)) | Die Green-Funktions-Kerne sind nicht summierbar. |
| Feynman-Perioden ([2607.25595](https://arxiv.org/abs/2607.25595)) | ν < d pro Kante. Graphical functions und HyperInt sind weit überlegen. |
| Epidemien, Ausbreitungskerne | Nur Normierungen niedriger Ordnung, per FFT oder im Kontinuum billig. 2025–26 wurde kein Engpass gefunden. |
| Ionen mit modenbasierten Kopplungen ([2606.13499](https://arxiv.org/abs/2606.13499)) | J_ij ∝ b_ik b_jk ist kein Potenzgesetz. |
| 1D-Tornheim- und MZV-Summen | Theoretisch bereits exakt bekannt, GZL taugt nur als Prüfbank. |
| PSLQ-Entdeckungen mit GZL, wie es heute ist | Doppelte Genauigkeit reicht nur für winzige Basen. |
| LiHoF₄, Pyrochlor-Spin-Eis, NaGdS₂ | Nicht-Bravais **und** anisotroper Dipol **und** ausdrücklich Ziel der Autoren. Erst nach der Roadmap. |
| Finanzen (Ising-Marktmodelle) | In der Recherche fand sich keine Quelle, die eine Gittersumme mit Potenzgesetz als Engpass zeigt. Nicht weiterverfolgt. |

---

## Empfehlung: morgen mit Idee 1 anfangen

**Starte mit „Optimale Gitter für Mehrkörper-Potenzgesetzenergien".** Die Idee hat die beste Kombination aus Neuheit und sofortiger Machbarkeit:

- Sie läuft mit GZL 1.0.0 ohne jede Erweiterung.
- Die Piloten liegen schon vor.
- Sie verbindet drei Welten (Gitteroptimierung, Kristallchemie, modulare Graphfunktionen).
- Das Ergebnis ist ein zeigbares Phasendiagramm.
- Sie berührt die Roadmap der Autoren nur am Rand.

Idee 2 ist der natürliche zweite, langfristige Strang mit dem größten Hebel. Die Infrastruktur für Graphenenumeration und Reihenanalyse, die man dafür baut, nützt später auch den Ideen 4 und 5.

**Schritt 1: Piloten härten (Woche 1–2). Erledigt, siehe Abschnitt 7.** In `/home/user/GapYas/research/pilots/` `opt_test.py`, `opt2.py` und `opt3.py` erneut laufen lassen (`mgf_test2.py` nur aus diesem Verzeichnis). Dann ν\* ≈ 3,91836 und die Barriere bei ~72° mit variiertem n_points, mit `evaluate_graph` statt `zeta_circle` und mit Brute-Force-Summen bei drei ν-Werten absichern. Die Fundamentaldomäne bis Im τ ≈ 5 und mit dichteren Starts absuchen. Ergebnis: eine belastbare 2D-Karte E₃(τ; ν).

**Schritt 2: globale 3D-Suche (Woche 3–8).** Die 3D-Gitter mit Einheitsvolumen über reduzierte Gram-Matrizen parametrisieren (5 Parameter) und die Dreiecks-Zeta für ν ∈ [3,2; 12] global minimieren. Ziel ist, festzustellen, ob BCC das globale Minimum ist oder ob ein anderes Gitter gewinnt. Parallel C_{s,s,s}(τ) für reelle s gegen die Laplace-Gleichungen der ganzzahligen Fälle prüfen (Brücke zu Idee 3).

**Schritt 3: Neuheit sichern und Kontakt aufnehmen (parallel, ab Woche 2).** Die offenen Fragen von Bétermin und die Arbeiten von Luo und Wei auf Mehrkörperaussagen durchsehen. Danach Buchheit bzw. den Schwerdtfeger-Kreis kontaktieren: GZL kommt von ihnen, und die offene BCC-Stabilitätsfrage aus [arXiv:2504.07338](https://arxiv.org/abs/2504.07338) ist ihr Thema. So wird das Scoop-Risiko zur Kooperation, und die für ATM nötige Fortsetzung für ν ≤ d kann gemeinsam in GZL einfließen.

**Schluss.** Durch die Recherche verschiebt sich der Blick: Das eigentliche Potenzial von GZL liegt weniger in „mehr Quantenmagnetismus" als darin, dass GZL als **erstes allgemeines Werkzeug für Graph-strukturierte Potenzgesetzsummen** Gebiete verbindet, die dieselben Zahlen unter verschiedenen Namen berechnen. In Spur A haben die Autoren fast jede naheliegende Richtung öffentlich beansprucht. Neuheit entsteht dort nur über konkrete Experimentpaare oder Kooperation. In Spur B ist das Feld leer: Seit dem Release am 25.09.2026 gibt es keine externe Nutzung, und klassische LR-Reihen, Mehrkörper-Gitteroptimierung und das MGF-Wörterbuch tauchen in keiner Roadmap auf. Die kritischste Erweiterung für den Transfer ist dabei nicht die Anisotropie, sondern die **meromorphe Fortsetzung für Nicht-Brücken-Blöcke mit ν ≤ d**. Sie würde ATM-Kristallphysik und stringrelevante MGFs mit a = 1-Kanten zugleich erschließen.

---

## Nachprüfung der Pilotrechnungen: bestätigt, präzisiert und um eine zweite Phase ergänzt

Die Pilotrechnungen aus den Ideen 1 und 3 wurden nachträglich unabhängig nachgerechnet. Alle Skripte, Rohdaten und Logs liegen in `research/verification/` und laufen mit `gzl 1.0.0` und `epsteinlib`. Die Kennzeichnung **[nachgeprüft]** heißt: mit zwei unabhängigen Verfahren übereinstimmend berechnet. Das ersetzt keine Begutachtung.

**Methode.** Als Referenz dient eine eigene Implementierung ohne GZL-Code (`lattice_sums.py`). Die Dreieckssumme T_ν(Λ) = Σ'_{x,y} |x|^{−ν}|y|^{−ν}|x−y|^{−ν} wird direkt über alle Gitterpunkte in einem Würfel [−R, R]^d summiert. Die innere Faltung läuft per FFT. Danach wird in R extrapoliert (Richardson mit dem bekannten Abschneidefehler ~R^{d−2ν}). Die FFT-Version stimmt mit einer naiven Doppelschleife auf Maschinengenauigkeit überein. Für 2D-Gitter ist die Referenz sogar schneller als GZL (0,04 s gegenüber 0,1 s pro Gitter).

### GZL ist korrekt (V1)

In 55 Testfällen wurde `gzl.zeta_circle` mit der Referenz verglichen: fünf 2D-Gitter und vier 3D-Gitter (FCC, BCC, SC, triklin), mit ν von 2,5 bis 9 und mit ungleichen Exponenten. **Die größte relative Abweichung beträgt 1,9·10⁻¹⁰, typisch sind 10⁻¹⁵** **[nachgeprüft]**. Die größeren Abweichungen treten nur nahe ν = d auf, wo die Extrapolation der Referenz schwieriger wird. `evaluate_graph` liefert für das Dreieck dasselbe wie `zeta_circle`, bis auf 3·10⁻¹⁶. In 3D braucht GZL 1–160 s pro Gitter, die FFT-Referenz etwa 1 s. Für die globale Suche wurde deshalb die Referenz verwendet und GZL als Kontrolle.

### 2D: drei Phasen statt zwei (V2)

| ν | globales Minimum (Suche über die ganze Fundamentaldomäne) | T hexagonal | T Quadrat |
|---|---|---|---|
| 2,2 | hexagonal | 29,7948755694 | 30,2601286381 |
| 3,0 | hexagonal | 13,3917941103 | 13,6528937154 |
| 3,9 | hexagonal | 7,9670208300 | 7,9723204026 |
| 3,95 | Quadrat | 7,7829092994 | 7,7737904082 |
| 6,0 | Quadrat | 3,7457712638 | 3,2499743124 |
| 8,0 | Quadrat | 2,2280735200 | 1,5524012126 |
| 10 | **Rechteck, b/a = 1,0998** (T = 0,7597084) | 1,4063994404 | 0,7622381339 |
| 20 | **Rechteck, b/a = 1,1941** (T = 0,0212133) | 0,1603699 | 0,0234490 |

- **Übergang hexagonal → Quadrat bei ν\* = 3,9183649026** **[nachgeprüft]**. Referenz und GZL liefern denselben Wert auf zehn Stellen. Das bestätigt den Pilotwert 3,91836.
- **Der Übergang ist erster Ordnung** **[nachgeprüft]**. Die Hesse-Matrix zeigt: Das Quadratgitter ist ab ν = 3,8063 ein lokales Minimum, das Hexagonalgitter nur bis ν = 4,2784. Dazwischen sind beide lokal stabil. Bei ν\* liegt auf dem Randbogen |τ| = 1 eine Barriere von 7,89851 auf 7,90956 bei arg τ ≈ 72°, also 0,14 %.
- **Neu gegenüber dem Pilot: Ein zweiter, kontinuierlicher Übergang bei ν₂ = 8,606** **[nachgeprüft]**. Dort wird das Quadratgitter in Streckrichtung instabil, und ein Rechteckgitter τ = i·y übernimmt. Sein Seitenverhältnis wächst stetig: 1,056 (ν = 9), 1,100 (10), 1,141 (12), 1,170 (15), 1,194 (20), 1,215 (30). Die Pilotaussage „Quadrat für alle ν > 3,918“ ist damit **falsch** für ν > 8,606. Der Pilot hatte nur bis ν = 6 gerechnet.
- **Grenzfall ν → ∞ [abgeleitet, nicht bewiesen]:** Für große ν zählt nur das kleinste Dreieck. Im Rechteckgitter konkurrieren das rechtwinklige Dreieck (Produkt √(y + 1/y)) und das kollineare Dreieck (Produkt 2y^{−3/2}). Gleichsetzen ergibt y⁴ + y² = 4, also **y∞ = √((√17 − 1)/2) ≈ 1,2496**. Die numerische Folge läuft darauf zu. Für ν ≥ 50 ist die FFT-Summe wegen des Dynamikumfangs nicht mehr genau genug; diese Werte sind verworfen.

![Energieunterschied Quadrat bzw. bestes Rechteck zum Hexagonalgitter. Die Nulldurchgänge markieren die beiden Übergänge.](../verification/figures/fig_2d_energy_difference.pdf)

![Kleinster Hesse-Eigenwert von Hexagonal- und Quadratgitter. Grau: Bereich, in dem beide lokal stabil sind (Koexistenz, Übergang erster Ordnung). Rechts: Instabilität des Quadrats bei ν₂ = 8,606.](../verification/figures/fig_2d_stability.pdf)

![Energielandschaft log₁₀(T/T_min − 1) über der Fundamentaldomäne bei vier Exponenten. Kreis: Hexagonalgitter, Quadrat: Quadratgitter. Bei ν = 10 liegt das Minimum auf der imaginären Achse oberhalb von τ = i (Rechteckgitter).](../verification/figures/fig_2d_landscape.pdf)

![Seitenverhältnis des optimalen Gitters für ν > ν₂.](../verification/figures/fig_2d_rect_ratio.pdf)

### 3D: BCC ist das beste gefundene Gitter, FCC wird instabil (V3, V5)

Die globale Suche lief über alle 3D-Bravais-Gitter mit Einheitsvolumen, also einen 5-dimensionalen Formraum. Pro ν gab es 16 Startpunkte (FCC, BCC, SC und 13 zufällige), jeweils mit Nelder-Mead und LLL-Reduktion.

| ν | FCC | BCC | SC | bestes gefundenes Gitter |
|---|---|---|---|---|
| 3,2 | 52,0983519937 | **52,0006824132** | 56,0330699444 | BCC |
| 3,5 | 38,1993132080 | **38,0743775719** | 41,9368217767 | BCC |
| 4,0 | 25,2265705426 | **25,0532095478** | 28,6610705747 | BCC |
| 4,5 | 17,9623763656 | **17,7406073980** | 21,0980682627 | BCC |
| 5,0 | 13,3808991194 | **13,1146038631** | 16,2165790058 | BCC |
| 6,0 | 8,0433635309 | **7,7101880730** | 10,2962433924 | BCC |
| 7,0 | 5,1456778058 | **4,7819587144** | 6,8781813217 | BCC |
| 9,0 | 2,3132634600 | **1,9784839245** | 3,2722630842 | BCC |
| 12,0 | 0,7708588009 | **0,5646023215** | 1,1328205539 | BCC |

- **BCC ist bei allen neun Exponenten das beste gefundene Gitter** **[nachgeprüft, numerisch]**. Alle Startpunkte, die tief genug kamen, liefen nach BCC. Die Werte bei ν = 3,5 bis 9 stimmen mit GZL auf 10⁻¹⁰ bis 10⁻¹⁵ überein (V1).
- **BCC ist ein lokales Minimum** gegenüber *allen* volumenerhaltenden Verzerrungen: Alle fünf Hesse-Eigenwerte sind positiv, geprüft bei ν = 3,5; 4,5; 6; 9. Auf dem Bain-Pfad (tetragonale Verzerrung) liegt das Minimum bei c/a = 1, also genau bei BCC. Das betrifft direkt die offene Stabilitätsfrage aus [arXiv:2504.07338](https://arxiv.org/abs/2504.07338), dort allerdings für Lennard-Jones plus ATM und nicht für die reine Dreieckssumme.
- **Neu: FCC verliert bei ν = 3,752 seine lokale Stabilität** **[nachgeprüft]**. Der kleinste Hesse-Eigenwert ist 2,02 bei ν = 3,5, 0,013 bei 3,75 und −5,97 bei 4,0. Oberhalb davon ist FCC ein Sattelpunkt.
- **Einschränkung:** Eine Suche mit 16 Starts ist kein Beweis des globalen Minimums. Unterhalb von ν = 3,2 wurde nicht gerechnet.

![Abstand von FCC und SC zu BCC in Prozent.](../verification/figures/fig_3d_lattices.pdf)

### Modulare Graphfunktionen: Identitäten auf Maschinengenauigkeit (V4)

Beide Identitäten wurden mit der Referenzsumme im Impulsraum geprüft, ohne GZL und an fünf τ-Werten:

- **C₂,₂,₁ = (2/5)E₅ + ζ(5)/30** gilt bis auf **10⁻¹⁵** relativ **[nachgeprüft]**. Der Pilot kam über die ε-Extrapolation in GZL nur auf 10⁻⁸.
- **C₁,₁,₁ = E₃ + ζ(3)** gilt bis auf **10⁻¹⁴** **[nachgeprüft]**. Das gelingt aber nur, wenn man in der Extrapolation die Terme R^{−2} log R berücksichtigt. Mit reinen Potenzen bleibt ein systematischer Fehler von 10⁻⁶. Der Grund: Alle drei Kanten liegen bei ν = d, dadurch divergiert die Nachbarsumme logarithmisch.
- Für C₂,₂,₂, wo alle ν > d sind, stimmen GZL und Referenz auf 10⁻¹⁵ überein. Das Wörterbuch „planare MGF = GZL-Zeta des dualen Graphen“ ist damit am Dreieck/Theta-Graph bestätigt.
- **Folgerung für Idee 1:** Wegen T_{2s}(Λ_τ) = π^{3s} C_{s,s,s}(τ) übersetzen sich die 2D-Ergebnisse direkt. Das Minimum von C_{s,s,s} über der Fundamentaldomäne springt bei **s\* = 1,95918** von e^{iπ/3} nach i und wird für **s > 4,303** zu einem Rechteckpunkt.

### Was sich am Bericht ändert

- Idee 1 ist bestätigt und wird stärker. Das 2D-Phasendiagramm hat **drei** Phasen (hexagonal → Quadrat → Rechteck), mit Übergängen erster und zweiter Ordnung. Dazu kommt in 3D ein klares Bild: BCC ist das beste Gitter, FCC wird ab ν = 3,752 instabil.
- Die Pilotaussage „Quadrat für ν > 3,918“ ist für ν > 8,606 korrigiert.
- Die Präzisionsangabe für die MGF-Identitäten steigt von 10⁻⁸ auf 10⁻¹⁴ bis 10⁻¹⁵. Diese Präzision stammt aber aus der eigenen Referenzsumme und nicht aus GZL. Die Aussage in Idee 3, dass GZL für a = 1-Kanten eine Erweiterung braucht, bleibt richtig.

### Neuheitsprüfung

NOVELTY_PLACEHOLDER

---

## Quellenliste

**Software und Dokumentation**: [GZL GitHub (README)](https://github.com/graph-zeta/gzl), [DOCUMENTATION.md](https://github.com/graph-zeta/gzl/blob/main/DOCUMENTATION.md), [gzl auf PyPI](https://pypi.org/project/gzl/), [epsteinlib auf PyPI](https://pypi.org/project/epsteinlib/).

**Eigene Pilotrechnungen (nicht begutachtet)**: `/home/user/GapYas/research/pilots/`:

- `mgf_test.py`: C₂,₂,₁- und Zagier-Test mit ε-Verschiebung;
- `mgf_test2.py`: Richardson-Extrapolation und tetraedrische K₄-MGF;
- `opt_test.py`: 2D-Gitterscan und 3D FCC/BCC/SC;
- `opt2.py`: Nelder-Mead und Brute-Force-Check;
- `opt3.py`: Brent-Wurzel ν\* und Barriere entlang des Bogens.

**Nachprüfung (unabhängige Referenz, Abschnitt 7)**: `research/verification/`:

- `lattice_sums.py`: Referenzimplementierung der Dreieckssumme ohne GZL-Code;
- `v1_gzl_vs_reference.py`: GZL gegen Referenz, 55 Fälle;
- `v2_2d_landscape.py`, `v2b_high_nu.py`, `v2c_large_nu.py`: 2D-Phasendiagramm, Stabilität, zweiter Übergang;
- `v3_3d_lattices.py`, `v5_3d_bain_stability.py`, `v5b_fcc_threshold.py`: globale 3D-Suche, Bain-Pfad, Hesse-Matrizen;
- `v4_modular_graph_functions.py`, `v4b_c111_log_fit.py`: MGF-Identitäten;
- `make_figures.py`: Abbildungen; Rohdaten und Logs in `results/`.

| arXiv | Kurzangabe (laut Notizen) | Jahr |
|---|---|---|
| [2609.18918](https://arxiv.org/abs/2609.18918) | Buchheit & Rupp, Methode der Graph-Zetafunktionen | 2026 |
| [2609.18761](https://arxiv.org/abs/2609.18761) | Duft, Adelhardt, Koziol, Buchheit, Schmidt, „Exact and fast series expansions for quantum models with long-range interactions" | 2026 |
| [2609.28282](https://arxiv.org/abs/2609.28282) | Buchheit & Busse, anisotrope singuläre Summen aus Ableitungen der Epstein-Zeta | 2026 |
| [2504.11989](https://arxiv.org/abs/2504.11989) | Buchheit & Busse, „Epstein zeta method for many-body lattice sums" | 2025 |
| [2412.16317](https://arxiv.org/abs/2412.16317) | Buchheit, Busse, Gutendorf, EpsteinLib | 2024 |
| [2504.07338](https://arxiv.org/abs/2504.07338) | Robles-Navarro et al., LJ + ATM auf kuboidalen Gittern | 2025 |
| [2012.05413](https://arxiv.org/abs/2012.05413) | Analytische LJ-Schwingungseffekte, Einstein-Modell | 2020 |
| [2403.00421](https://arxiv.org/abs/2403.00421) | Adelhardt, Koziol, Langheld, Schmidt, MC-Techniken für LR-Quantenmagnete | 2024 |
| [1802.06684](https://arxiv.org/abs/1802.06684) | Fey, Kapfer, Schmidt, 2D-LR-Quantenmagnete | 2018 |
| [2203.08081](https://arxiv.org/abs/2203.08081) | Langheld et al., Skalierung oberhalb der oberen kritischen Dimension | 2022 |
| [2408.13145](https://arxiv.org/abs/2408.13145) | Adelhardt, Duft, Schmidt, LR-XXZ-Bilayer | 2024 |
| [2604.12754](https://arxiv.org/abs/2604.12754) | Adelhardt, Muleady, Schmidt, Gorshkov, LR-Spin-1-Kette | 2026 |
| [2604.20831](https://arxiv.org/abs/2604.20831) | LR-Spin-1-Kette, QMC α_c = 2,49(1) | 2026 |
| [2601.20058](https://arxiv.org/abs/2601.20058) | LR-XY-Kette, α_c ≈ 2,7 | 2026 |
| [2606.22407](https://arxiv.org/abs/2606.22407) | Li, Fan, Chen, Deng, Universalitätsdiagramm LR-Quantenkritikalität | 2026 |
| [2512.04805](https://arxiv.org/abs/2512.04805) | Xiao, Liu, Fan, Deng, 2D LR-Ising σ\* = 2 | 2025 |
| [2512.01956](https://arxiv.org/abs/2512.01956) | Xiao et al., 2D LR-Heisenberg | 2025 |
| [2608.20750](https://arxiv.org/abs/2608.20750) | Liu, Xiao, Fan, Deng, 2D LR-Perkolation | 2026 |
| [2508.18807](https://arxiv.org/abs/2508.18807) | Hutchcroft, Critical LRP I | 2025 |
| [2508.18808](https://arxiv.org/abs/2508.18808) | Hutchcroft, Critical LRP II | 2025 |
| [2404.07276](https://arxiv.org/abs/2404.07276) | Hutchcroft, Dreiecksbedingung LRP | 2024 |
| [2607.22387](https://arxiv.org/abs/2607.22387) | Blanc-Renaudie & Hutchcroft, k-Punkt-Funktion | 2026 |
| [2602.07818](https://arxiv.org/abs/2602.07818) | Zweischleifen-ε-Entwicklung LR-O(n) | 2026 |
| [2608.15120](https://arxiv.org/abs/2608.15120) | 6 − ε für LR-Perkolation | 2026 |
| [2510.02458](https://arxiv.org/abs/2510.02458) | Pagni et al., FRG für die 1D-LR-Ising-Kette | 2025 |
| [2604.07130](https://arxiv.org/abs/2604.07130) | Okuyama & Ohzeki, 1D-LR-Spinglas | 2026 |
| [2412.08814](https://arxiv.org/abs/2412.08814) | Liu, Park, Slade, FSS oberhalb d_c | 2024 |
| [2308.09709](https://arxiv.org/abs/2308.09709) | NQS für die AF-LR-Ising-Kette | 2023 |
| [2609.17356](https://arxiv.org/abs/2609.17356) | Luo & Wei, Sarnak–Strömbergsson | 2026 |
| [2605.07580](https://arxiv.org/abs/2605.07580) | Luo & Wei, Minimierer von Theta- und Epstein-Verhältnissen | 2026 |
| [2411.17199](https://arxiv.org/abs/2411.17199) | Deng & Luo, hexagonale Optimalität | 2024 |
| [2407.20762](https://arxiv.org/abs/2407.20762) | Bétermin & Furlanetto, p-Normen | 2024 |
| [2312.01395](https://arxiv.org/abs/2312.01395) | Bétermin, Šamaj, Travěnec, Strukturübergänge | 2023 |
| [1512.06779](https://arxiv.org/abs/1512.06779) | MGFs und elliptische Polylogarithmen | 2015 |
| [1608.04393](https://arxiv.org/abs/1608.04393) | D'Hoker & Kaidi, MGF-Identitäten | 2016 |
| [2502.05531](https://arxiv.org/abs/2502.05531) | Claasen & Doroudiani, iterierte Eisenstein-Integrale | 2025 |
| [2511.15883](https://arxiv.org/abs/2511.15883) | Schlotterer, Sohnle, Tao, elliptische MGFs | 2025 |
| [2512.21413](https://arxiv.org/abs/2512.21413) | Fedosova & Klinger-Logan, Konvolutionsidentitäten | 2025 |
| [2603.20550](https://arxiv.org/abs/2603.20550) | Dobrowolski, Mordell–Tornheim | 2026 |
| [2607.25595](https://arxiv.org/abs/2607.25595) | Portner, Feynman-Perioden und svMZVs | 2026 |
| [2603.20372](https://arxiv.org/abs/2603.20372) | 256-Qubit-Simulation TmMgGaO₄ | 2026 |
| [2608.07178](https://arxiv.org/abs/2608.07178) | Pasqal/Browaeys, 2D-TFIM mit 256 Qubits | 2026 |
| [2601.16275](https://arxiv.org/abs/2601.16275) | Caltech, DSF einer Rydberg-Kette | 2026 |
| [2505.09884](https://arxiv.org/abs/2505.09884) | NaTmSe₂ als TFIM | 2025 |
| [2306.03544](https://arxiv.org/abs/2306.03544) | KTmSe₂-INS | 2023 |
| [2311.11726](https://arxiv.org/abs/2311.11726) | Chen et al., Quench-Spektroskopie dipolar XY | 2023 |
| [2605.07685](https://arxiv.org/abs/2605.07685) | DMRG dipolares XY auf Archimedischen Gittern | 2026 |
| [2605.28599](https://arxiv.org/abs/2605.28599) | NLCE+QA auf Ionen-QPU | 2026 |
| [2512.09037](https://arxiv.org/abs/2512.09037) | Magnon-Bindungszustände im 2D-LRTFIM | 2025 |
| [2605.13969](https://arxiv.org/abs/2605.13969) | Squeezing-Übergang in Potenzgesetz-Bilayern | 2026 |
| [2605.27484](https://arxiv.org/abs/2605.27484) | Koziol, ferromagnetisches Dicke-Ising-Modell | 2026 |
| [2604.22743](https://arxiv.org/abs/2604.22743) | Fitzner, Lesanovsky, Sbierski, Thermometrie Kagome-Rydberg | 2026 |
| [2608.15576](https://arxiv.org/abs/2608.15576) | Optimale Abschneideschranke O(R^{−(p−d)}) | 2026 |
| [2606.20522](https://arxiv.org/abs/2606.20522) | iMPS ohne Surrogat für algebraische Schwänze | 2026 |
| [2512.24596](https://arxiv.org/abs/2512.24596) | Xie & Schotland, Photonbandstruktur atomarer Gitter | 2025 |
| [2609.17313](https://arxiv.org/abs/2609.17313) | Exzeptionelle Topologie in 2D-Atomarrays | 2026 |
| [2507.18399](https://arxiv.org/abs/2507.18399) | Modifizierte Lüscher-Zeta mit langreichweitiger Kraft | 2025 |
| [1702.01296](https://arxiv.org/abs/1702.01296) | Matzelle & Tiburzi, QED-Finite-Volume-Summen | 2017 |
| [2608.10041](https://arxiv.org/abs/2608.10041) | He & Hu, Madelung-Randterm | 2026 |
| [2503.00977](https://arxiv.org/abs/2503.00977) | Calara & Miller, Madelung per Multipolsummation | 2025 |
| [2603.13086](https://arxiv.org/abs/2603.13086) | Erb, Meyer, Ottnad, Photonpropagator Gitter-QCD | 2026 |
| [2606.06617](https://arxiv.org/abs/2606.06617) | PSWF-LR, langreichweitige ML-Potenziale | 2026 |
| [2512.18029](https://arxiv.org/abs/2512.18029) | Latent Ewald | 2025 |
| [2308.03140](https://arxiv.org/abs/2308.03140) | libMBD | 2023 |
| [2412.01076](https://arxiv.org/abs/2412.01076) | LR-Kitaev-Ketten nach Quench | 2024 |
| [2609.28735](https://arxiv.org/abs/2609.28735) | Erweiterte LR-Kitaev-Modelle | 2026 |
| [2606.13499](https://arxiv.org/abs/2606.13499) | 2D-LR-XY mit 200 Ionen | 2026 |
