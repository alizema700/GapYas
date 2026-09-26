# Optimal lattices for a three-body power-law energy — code and certificates

Supplementary material for the preprint

> A. Suleman, *Optimal lattices for a three-body power-law energy: steep-decay limits, computer-assisted
> proofs and the minima of the modular graph functions C_{s,s,s}* (2026).

DOI of this archive: 10.5281/zenodo.22975868. The computations were carried out with the assistance of the
AI system Claude (Anthropic); see the Declarations of the paper.

## Requirements
Python ≥ 3.11 with `numpy 2.4`, `scipy 1.17`, `python-flint 0.9` (Arb), `mpmath 1.3`; for the numerics
also `epsteinlib 0.6`, `gzl 1.0`, `matplotlib`. Lean 4 + Mathlib for `lean/` (see `lean/README.md`).

```
pip install numpy scipy python-flint mpmath epsteinlib gzl matplotlib
```

## Rigorous proofs (`rigorous/`)
All certificates use outward-rounded IEEE interval arithmetic (`iv.py`) and Arb ball arithmetic.

| Statement (paper) | Command | Log |
|---|---|---|
| Thm 3: BCC maximises P, B&B (≈4 min) | `python bnb3d_rig.py 0.0069` | `bnb3d_rig.log`, `bnb3d_rig_result.json` |
| Thm 3: local lemma, exact c² = 27/1400 | `python bcc_local_exact.py`, `python local3d_rig.py` | `bcc_local_exact.txt`, `local3d_rig.json` |
| Thm 4: square minimiser, ν = 4,5,6,7,8 (20–40 min each) | `python prove_min_rig.py <nu> square` | `prove_rig_square_nu<nu>.log` |
| Thm 4: hexagonal minimiser, ν = 3.5 | `python prove_min_rig.py 3.5 hex` | `prove_rig_hex_nu3.5.log` |
| Thm 4: crossing ν* ∈ (3.918364, 3.918366) | `python certify_crossing_rig.py` | `certify_crossing_rig.log` |
| Thm 4: ν₂ ∈ (8.604, 8.608) | `python certify_nu2.py` | `certify_nu2.log` |
| Thm 4: rectangular minimiser, ν = 10,12,…,20 | `python prove_rect_rig.py <nu> <y0>` (y0 from the log) | `prove_rig_rect_nu<nu>.log` |

Run the commands from inside `rigorous/`.

## Large-exponent asymptotics (`largenu/`)
Balance equation, 1/ν law and figure data (`asym2d.py`, `balance2d.py`, `asym_data.py`), active sets and hull margins (`active_sets.py`), Hessian at BCC (`hess3d.py`), working notes (`largenu.tex`).

## Formal proof (`lean/`)
Lean 4 formalisation of the planar max–min theorem (Theorem 1); no `sorry`, standard axioms only.

## Numerics (`verification/`, `extended/`, `theorem*/`)
- `verification/`: reference lattice sums (`lattice_sums.py`), comparison with GZL, 2D/3D landscapes, modular graph function identities.
- `extended/`: 3D global search, FCC Bain curvature, unequal exponents, hcp and mixed energies, low ν, error bars, dense λ_c.
- `theorem/`, `theorem3d/`, `theorem_mgf/`: earlier versions of the provers with floating-point safety margins (superseded by `rigorous/`, kept for reference).

## Paper (`paper/`)
LaTeX source, figure scripts (`make_paper_figs.py`) and the compiled PDF.

## License
MIT.
