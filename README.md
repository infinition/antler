# ANTLER

Reproducible quantum-many-body study of digital statistical-phase transport in
a number-conserving correlated-hopping SSH ladder, using exact diagonalization.

## Why ANTLER?

ANTLER stands for **Anyonic Numerical Transport for Logical Exchange and
Robustness**. It names the research program: numerical transport in a
correlated-hopping ladder, logical exchange protocols and quantitative
robustness tests. The name describes the motivation of the model, not a claim
that physical anyons or non-Abelian braiding have been demonstrated.

## Scope

This public release contains the validated Abelian digital phase primitive,
the Phase 0 to 4.7 audit trail used to establish it, and the separately
curated Phase 5 to 8 research archive.
It does not claim non-Abelian braiding, universality, fault tolerance, or a
topological quantum computer.

For the digital protocol, the all-fine-step data at `dt=0.125` and
`D=4,6,8,10,12` fit

```text
odd_slope(D) = -1 + 0.66493 D^(-2.28624)
```

with `R^2 = 0.999983`. All retained points have leakage below `5e-5` and
minimum singular value above `0.999975`.

## Quick start

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -r requirements.txt
python SMOKE_TEST.py
python scripts\verify_public_results.py
python paper\generate_figures.py
python paper\verify_paper.py
```

## Reproduce the Phase 4.7 closure

```powershell
python experiments\phase4_7\run_phase47_publication_closure.py
```

The command reads the final numerical JSON files under `results/phase4_7/`
and regenerates `results/phase4_7/publication_closure.json`.

## Repository layout

- `antler/`: basis and Hamiltonian utilities used by the public protocol.
- `experiments/phase0` to `experiments/phase4_7/`: frozen construction,
  logical-frame, transport and validation scripts for the Abelian result.
- `experiments/phase5` to `experiments/phase7/`: constrained non-Abelian
  searches, reference controls and scoped no-go audits.
- `results/phase4` to `results/phase4_7/`: numerical data supporting the
  Abelian result.
- `results/phase5` to `results/phase8c/`: numerical outputs for the
  separately curated Phase 5 to 8 archive.
- `docs/public/`: result summary, reproducibility notes and claim boundary.
- `paper/`: LaTeX manuscript, figures and bibliography for the Abelian paper.

## Claim boundary

The observed ANTLER gate family is Abelian. The public data support a
controlled logical phase primitive in the specified model only. They do not
establish a physical anyon platform or non-Abelian topological computation.
The Phase 5 to 8 archive records scoped negative results, external controls
and an explicit neutral-link extension contract. Its claim boundary is in
`docs/public/PHASE5_8_RESEARCH_ARCHIVE.md`.

## License and citation

No reuse license has been selected yet. Until one is added, the repository is
publicly visible but no reuse permission is granted. Citation metadata will be
maintained in `CITATION.cff`.
