# Rechercheauftrag: Neue Anwendungsfelder für die Graph Zeta Library (GZL)

## 1. Rolle und Ziel

Du bist ein Forschungs-Scout mit Hintergrund in mathematischer Physik, statistischer Physik
und numerischer Mathematik. Deine Aufgabe ist **nicht**, GZL zusammenzufassen, sondern
**Brücken zu finden**: Wo in der aktuellen Forschung (arXiv, Schwerpunkt 2025–2026) tauchen
genau die mathematischen Objekte auf, die GZL schnell und präzise berechnet, und wo wurde
diese Verbindung **noch nicht** gesehen?

Ziel ist eine Liste konkreter, **neuartiger** Projektideen. „Neuartig“ heißt: so etwas wie ein
Hackathon-Projekt beim ETHack. Es verbindet zwei Welten, die vorher niemand verbunden hat,
und liefert ein Ergebnis, das es vorher nicht gab: eine neue Zahl, eine neue Kurve, ein neues
Phasendiagramm oder einen Vergleich mit einem Experiment.

## 2. Rahmen des Projekts

- **Langfristiges Forschungsprojekt**, kein Wochenend-Hack. Einarbeitungsaufwand ist kein
  Ausschlusskriterium; entscheidend sind Neuheit, wissenschaftlicher Hebel und Passung zu GZL.
- Ideen dürfen auch eigene Erweiterungen von GZL erfordern (neue Korpora, neue Kerne,
  eigene Graphen-Enumeration), solange der Aufwand ehrlich benannt wird.
- Ziel: ein Ergebnis, das Forschenden auffällt (Paper, Kooperation, Abschlussarbeit, Open-Source-Beitrag).

## 3. Was GZL kann (Fakten, bitte gegen die Quellen prüfen)

**Quellen, die du zuerst lesen musst:**
- Methode (Mathematik): arXiv:2609.18918 (Buchheit, Rupp)
- Anwendung (Physik, Linked-Cluster-Expansionen, LRTFIM, Material KTmSe₂): arXiv:2609.18761
- Code und Doku: https://github.com/graph-zeta/gzl (README, DOCUMENTATION.md, AGENTS.md, CHANGELOG.md)
- PyPI: `pip install gzl` (v1.0.0, AGPL-3.0, Python ≥ 3.11)

**Kernobjekt (vereinfacht, exakte Definition im Paper nachlesen):** Zu einem Graphen G mit
Knoten und Kanten wird eine **Gittersumme** berechnet: Man summiert über alle Platzierungen
der Knoten auf einem Bravais-Gitter Λ = A·ℤ^d das Produkt von Kanten-Kernen K(x_i − x_j),
optional mit einer Phase e^{ik·(…)} für einen externen Impuls k. Das nennt sich
„Graph-Zetafunktion“ ζ_G(k; ν).

**Unterstützte Kerne** (`gzl.Interaction`):
K(x) = a(x) + Σ_j b_j·|x|^(−ν_j),
also Potenzgesetze (auch Summen mehrerer Exponenten) plus ein reeller, gerader, kompakt
getragener Kurzreichweiten-Anteil a(x). Dieser Anteil darf anisotrop sein, der Potenzgesetz-Anteil
ist isotrop. Verschiedene Kanten können verschiedene Kerne haben.

**Numerische Stärken:**
- Faktorisierung des Graphen in Blöcke. Brücken und Zyklen werden analytisch über
  verallgemeinerte (Epstein-)Zetafunktionen behandelt. Serien-parallele Blöcke (Baumweite ≤ 2)
  kosten linear in der Knotenzahl, Blöcke mit höherer Baumweite laufen über
  Tensornetz-Bucket-Elimination (Polynomial in der Gittergröße, Exponent hängt nur von der
  Baumweite ab).
- Präzise auch für ν nahe der Dimension d, wo direkte Summation extrem langsam konvergiert.
- **Das gesamte Impulsgitter (Brillouin-Zone) per FFT zum Preis eines einzelnen Impulspunkts.**
  Monte Carlo bräuchte einen Lauf pro Impuls.
- Selbst getestet: `gzl series --corpus tfim1qp --A cubic --nu 4 --n-points 8 --order-max 9`
  liefert 9 Ordnungen × 512 Impulspunkte in ~17 s auf einem Laptop.

**Öffentliche API (Auswahl):** `evaluate_graph` (beliebiger Graph als Kantenliste oder
NetworkX-MultiGraph), `Interaction`, `zeta_circle`, `graph_zeta_general`, `slab_zeta`,
`hybrid_zeta`, `direct_sum_*`, `graph_convolve`/`graph_multiply`, CLI `gzl series` mit
TOML-Konfiguration. Gitter: `chain`, `square`, `triangular`, `cubic` **oder jede beliebige
Gittermatrix A** (d = 1, 2, 3).

**Grenzen (heute):**
- Nur Bravais-Gitter (keine mehratomige Basis). Das steht auf der Roadmap.
- Potenzgesetz-Teil nur isotrop. Anisotrope Langreichweite (z. B. echte Dipol-Dipol-Winkelabhängigkeit)
  steht auf der Roadmap.
- Nur reelle, gerade Kerne.
- Mitgelieferte physikalische Daten („Korpora“) bisher nur für das Transversalfeld-Ising-Modell
  (0qp bis Ordnung 13, 1qp bis Ordnung 11). Heisenberg u. a. sind angekündigt.
- `evaluate_graph` selbst ist aber **modellunabhängig**: Jeder, der seine eigenen Graphen und
  Vorfaktoren mitbringt, kann es nutzen. Genau hier liegt das Potenzial für neue Felder.

Wenn du feststellst, dass eine dieser Angaben nicht stimmt, korrigiere sie ausdrücklich.

## 4. Die zentrale Frage

> **Wo in Physik, Mathematik oder anderen Feldern müssen heute Summen der Form
> „Summe über Gitterplätze von Produkten paarweiser Potenzgesetz- bzw. Kurzreichweiten-Kerne,
> strukturiert durch einen Graphen“ berechnet werden, und wird das dort bisher mühsam gemacht
> (Monte Carlo, Abschneiden bei endlicher Reichweite, Ewald-Summation, kleine Systeme, grobe
> Impulsgitter), sodass GZL einen echten Sprung ermöglichen würde?**

Arbeite in **zwei Spuren**:

**Spur A: „In der Nähe, aber weiter als die Autoren“.** Aktuelle Physik-Themen mit
langreichweitigen Wechselwirkungen auf Gittern, in denen GZL neue Ergebnisse liefern könnte,
die in arXiv:2609.18761 noch nicht vorkommen (andere Gitter, andere Exponenten, andere Observablen,
Vergleich mit neuen Experimenten).

**Spur B: „Transfer in ein anderes Feld“.** Bereiche, in denen dieselbe Mathematik auftaucht,
aber unter anderem Namen, und wo niemand an Graph-Zetafunktionen denkt.

## 5. Kandidatenfelder (Hypothesen, prüfen und nicht einfach übernehmen)

Prüfe mindestens diese Felder, und ergänze eigene:

**Spur A (Quantenvielteilchen und Quantensimulation):**
1. Rydberg-Atom-Arrays und optische Pinzetten (van-der-Waals ν = 6, Dipol ν = 3): aktuelle
   Experimente zu Phasendiagrammen und Anregungsspektren auf 2D/3D-Gittern.
2. Gefangene Ionen mit einstellbarem Exponent 0 < ν < 3: Spektroskopie von Quasiteilchen,
   Lichtkegel, Vergleich mit Theorie.
3. Polare Moleküle, magnetische Atome, Dipolare Quantenmagnete.
4. Langreichweitige Quantenkritikalität: Wo wechselt die Universalitätsklasse mit ν,
   gibt es offene Streitfragen zu kritischen Exponenten, die Reihenentwicklungen entscheiden könnten?
5. Neue Quantenmaterialien mit langreichweitiger Kopplung (KTmSe₂ und verwandte Seltenerd-
   Verbindungen, Neutronenstreu-Daten, frustrierte Dreiecksgitter).
6. Langreichweitige Kitaev-Ketten und Topologie, Magnonen und Spinwellen in höherer Ordnung.

**Spur B (Transfer):**
7. Klassische statistische Physik: Mayer-/Virial- und Hochtemperatur-Entwicklungen langreichweitiger
   Modelle, klassisches Langreichweiten-Ising, Lévy-Flüge und Zufallsbewegungen auf Gittern,
   Gitter-Green-Funktionen.
8. Gitter-Feldtheorie: Störungstheorie auf dem Gitter, Feynman-Diagramme mit Gitterpropagatoren.
9. Elektrostatik und Materialwissenschaft: Madelung-Konstanten, Kristallenergien, Casimir- und
   Van-der-Waals-Summen, Ewald-Alternativen.
10. Zahlentheorie und Analysis: Epstein-Zeta, mehrfache Gittersummen, modulare Formen,
    Gitter-Summen-Identitäten. Kann GZL als numerisches Experimentierwerkzeug neue Vermutungen stützen?
11. Graphentheorie, Kombinatorik und Netzwerke: Zählprobleme, Einbettungen, Baumweite-Algorithmen.
12. Überraschende Felder: maschinelles Lernen (Kernmethoden auf Gittern, periodische Faltungen),
    Signal- und Bildverarbeitung, Kosmologie oder Gravitation auf Gittern, Akustik und Photonik
    (Metamaterialien, periodische Strukturen), Epidemiologie oder Ökologie mit Potenzgesetz-
    Ausbreitungskernen, Finanzen (Ising-Marktmodelle). Hier bitte besonders ehrlich sein, ob die
    Verbindung trägt oder konstruiert ist.

## 6. Vorgehen

1. **Grundlage verstehen:** Lies Abstract, Einleitung und Methodenteil beider GZL-Papers sowie
   README und AGENTS.md. Formuliere in 5 Sätzen, welche Klasse von Summen GZL berechnet und welche nicht.
2. **Suche:** arXiv, Schwerpunkt Veröffentlichungen **2025 bis heute**, ältere Grundlagenpapers nur,
   wenn sie die Verbindung begründen. Kategorien u. a.: cond-mat.str-el, cond-mat.stat-mech,
   cond-mat.quant-gas, quant-ph, physics.atom-ph, hep-lat, math-ph, math.NT, math.CO, cs.DS.
   Suchbegriffe (Englisch), kombinieren und variieren:
   `long-range interactions lattice`, `power-law interactions series expansion`,
   `linked-cluster expansion`, `perturbative continuous unitary transformation long-range`,
   `lattice sum Epstein zeta`, `Rydberg array excitation spectrum`, `trapped ion quasiparticle dispersion`,
   `long-range transverse-field Ising critical exponents`, `dipolar quantum magnet`,
   `Madelung lattice sum`, `lattice Green function power law`, `Lévy flight lattice`,
   `lattice perturbation theory Feynman integrals`, `multiple lattice sums`, `graph zeta function`.
3. **Pro Kandidat prüfen:**
   - Welche Summe bzw. welches Integral wird dort konkret berechnet? Formel angeben.
   - Lässt sie sich auf das GZL-Kernobjekt abbilden? Welcher Graph, welches Gitter, welcher Kern?
     Falls nicht direkt: was fehlt (mehratomige Basis, Anisotropie, komplexe Kerne, d > 3)?
   - Wie wird es dort bisher gerechnet, und wo sind die Grenzen (Ordnung, Genauigkeit, Rechenzeit,
     Impulsauflösung)?
   - **Neuheitscheck:** Hat schon jemand GZL oder die Methode aus arXiv:2609.18918 darauf angewandt?
     Zitierende Arbeiten prüfen (arXiv-Listing, Google Scholar, Semantic Scholar).
   - Aufwand und nötige Erweiterungen (siehe Abschnitt 2).
4. **Querverbindungen suchen:** Gibt es Paare von Feldern, die dieselbe Summe benutzen, ohne
   voneinander zu wissen? Das sind die wertvollsten Funde.

## 7. Bewertung jeder Idee (je 1–5 Punkte, mit kurzer Begründung)

- **Neuheit:** Hat das so noch niemand gemacht?
- **Passung zu GZL:** Funktioniert es mit der heutigen Version, oder braucht es Features der Roadmap?
- **Hebel:** Wie viel besser als der Status quo (Genauigkeit, Geschwindigkeit, volle Brillouin-Zone)?
- **Machbarkeit:** Aufwand bis zum ersten belastbaren Ergebnis, und welche Erweiterungen nötig sind.
- **Wow-Faktor:** Lässt sich das Ergebnis zeigen (Plot, interaktive Demo, Vergleich mit Experiment)?
- **Anschlussfähigkeit:** Interessiert das Forschende (Paper, Abschlussarbeit, Kooperation)?

## 8. Ausgabeformat

1. **Kurzfazit (max. 10 Zeilen):** die drei stärksten Ideen in je einem Satz.
2. **Was GZL wirklich kann:** 5 Sätze plus Liste der Grenzen, inkl. Korrekturen an Abschnitt 3.
3. **Top 5 Projektideen**, jeweils mit:
   - Titel und Ein-Satz-Pitch
   - Das Problem im Zielfeld und warum es heute schwer ist
   - Die konkrete Abbildung auf GZL (Graph, Gitter, Kern, Observable), gerne als Pseudocode
   - Was neu herauskäme
   - Erster Meilenstein (Minimalergebnis) und Ausbaupfad über Monate
   - Nötiges Vorwissen und die besten Einstiegsquellen
   - Risiken und was schiefgehen kann
   - 3–6 Schlüsselpapers (arXiv-ID, Titel, Jahr, Link, ein Satz Relevanz)
   - Bewertung nach Abschnitt 7
4. **Weitere Ideen** (Kurzform, Tabelle mit Bewertung).
5. **Überraschende Zusammenhänge:** Felder, die dieselbe Mathematik teilen, auch wenn noch
   keine fertige Projektidee daraus folgt.
6. **Verworfene Ideen:** kurz, warum sie nicht tragen. Das ist genauso wichtig.
7. **Empfehlung:** eine Idee, mit der ich morgen anfangen sollte, und die ersten drei Schritte.
8. **Quellenliste** aller zitierten Arbeiten.

## 9. Qualitätsregeln

- **Keine erfundenen Papers.** Jede arXiv-ID muss existieren und inhaltlich zum Satz passen.
  Im Zweifel weglassen oder als „nicht verifiziert“ markieren.
- Trenne klar zwischen **belegt** (steht im Paper), **abgeleitet** (deine Schlussfolgerung) und
  **spekulativ** (Idee, die geprüft werden muss).
- Keine konstruierten Verbindungen schönreden. Wenn eine Brücke nur mit Gewalt hält, sag es.
- Wenn eine Idee Features braucht, die GZL heute nicht hat, markiere das deutlich.
- Aktualität zählt: Neueste Arbeiten (2025–2026) bevorzugen, weil dort die offenen Fragen liegen.
- Schreibe die Ergebnisse auf Deutsch, Fachbegriffe und Suchbegriffe dürfen Englisch bleiben.
