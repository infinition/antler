# Paper build

From the repository root, regenerate the numerical figures and tables with:

```powershell
python paper\generate_figures.py
```

Compile the manuscript from `paper/` with a LaTeX distribution:

```powershell
pdflatex main.tex
bibtex main
pdflatex main.tex
pdflatex main.tex
```

The figure generator reads only the tracked JSON files under `results/`.
The manuscript reports the Abelian Phase 4.7 result and its explicit claim
boundary. It is not a paper claiming non-Abelian braiding.
