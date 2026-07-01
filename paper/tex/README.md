# SEDE — LaTeX build

## Current, canonical source: `SEDE.tex`
- **`SEDE.tex`** is regenerated from the canonical Markdown (`../SEDE.md`) by **`../build_tex.sh`**
  (mdmath + pandoc + `../template.tex` + `../header.tex`), so it never silently diverges from the paper.
- Compile with **XeLaTeX** (Unicode Greek/math via STIX Two Text + `header.tex`'s `newunicodechar` maps):
  ```
  cd paper && bash build_tex.sh      # -> tex/SEDE.tex and tex/SEDE.pdf
  ```
- **Figures are vector**: `\includegraphics` uses extension-less basenames under `../../output`, so
  XeLaTeX prefers the `*.pdf` (vector) over `*.png`. Regenerate the vector PDFs with
  `python paper/make_vector_figures.py` (or the per-figure `run_*.py`).
- References are **inline** in `SEDE.md` (no BibTeX needed for this build).
- `SEDE.pdf` is the file to submit (arXiv / JCAP initial submission). The Paper II companion builds the
  same way via `../build_foundations.sh` → `SEDE_foundations.tex` / `SEDE_foundations.pdf`.
