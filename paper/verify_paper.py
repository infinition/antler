"""Check manuscript inputs that do not require a local LaTeX compiler."""
from __future__ import annotations

import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PAPER = ROOT / "paper"


def main() -> None:
    main_tex = (PAPER / "main.tex").read_text(encoding="utf-8")
    generated = (PAPER / "generated_results.tex").read_text(encoding="utf-8")
    bibliography = (PAPER / "references.bib").read_text(encoding="utf-8")
    closure = json.loads((ROOT / "results" / "phase4_7" / "publication_closure.json").read_text())

    assert "Fabien POLLY" in main_tex
    assert "https://github.com/infinition/antler" in main_tex
    assert f"{closure['fit']['a']:.5f}" in main_tex
    assert f"{closure['fit']['p']:.5f}" in main_tex
    assert f"{closure['fit']['r2']:.6f}" in main_tex

    figures = re.findall(r"\\includegraphics(?:\[[^]]*\])?\{([^}]+)\}", main_tex)
    assert figures
    for figure in figures:
        assert (PAPER / "figures" / figure).exists(), figure

    labels = set(re.findall(r"\\label\{([^}]+)\}", main_tex + generated))
    references = set(re.findall(r"\\ref\{([^}]+)\}", main_tex))
    assert references <= labels, sorted(references - labels)

    citations = {
        key.strip()
        for group in re.findall(r"\\cite\{([^}]+)\}", main_tex)
        for key in group.split(",")
    }
    keys = set(re.findall(r"@\w+\{([^,]+),", bibliography))
    assert citations <= keys, sorted(citations - keys)

    for path in PAPER.rglob("*"):
        if path.is_file() and path.suffix in {".tex", ".bib", ".md", ".py"}:
            assert "\u2014" not in path.read_text(encoding="utf-8", errors="ignore"), path

    print(f"PASS: paper inputs verified, figures={len(figures)}, citations={len(citations)}")


if __name__ == "__main__":
    main()
