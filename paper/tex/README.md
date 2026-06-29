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
- Current: 37 pp, vector figures, includes the no-SH0ES robustness (§5.4 ii′), Fig 5 (CMB) and
  Fig 9 (E1), and the round-4 closes. This is the file to submit (arXiv / JCAP initial submission).

## Superseded / acceptance-stage WIP (do not treat as current)
- `main.tex` (RevTeX) + `body.tex` + `main_jcap.tex` (jcappub) are an earlier **hand-port** that has
  **diverged** (older structure: still has `fig_chr_experiments` and an "Honest verdict" subsection that
  `SEDE.md` removed; missing this session's additions). They are kept only as the **journal-class
  (jcappub, pdflatex) target for the acceptance stage**, which needs a manual re-port with macro-based
  (not literal-Unicode) math. JCAP accepts the initial submission as a single PDF, so `SEDE.pdf` suffices
  to submit; the jcappub version is for revision/acceptance.
