# Reproducibility

Install the dependencies and run the smoke test.

```powershell
pip install -r requirements.txt
python SMOKE_TEST.py
python scripts\verify_public_results.py
python paper\generate_figures.py
python paper\verify_paper.py
```

The public closure can be rebuilt from the included JSON files.

```powershell
python experiments\phase4_7\run_phase47_publication_closure.py
```

The script checks the five fine-step depth points and writes the fitted
summary to `results/phase4_7/publication_closure.json`.

The Phase 5 to 8 archive is supplied as exact scripts, JSON outputs and
scoped documents. Its public index is `PHASE5_8_RESEARCH_ARCHIVE.md`.
