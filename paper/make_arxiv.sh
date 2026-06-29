#!/usr/bin/env bash
# Build the arXiv source tarball from the CURRENT manuscript (tex/SEDE.tex, regenerated
# from the canonical SEDE.md by build_tex.sh). Flattens figure paths to ./ and bundles the
# 11 VECTOR figure PDFs. References are inline in SEDE.md (no bibtex/.bbl needed).
#
#   Engine: XeLaTeX. The source loads fontspec + STIX Two Text (both ship with TeX Live,
#   so arXiv has them); arXiv's AutoTeX detects fontspec and processes with xelatex.
#
# Run from paper/:  bash make_arxiv.sh
set -e
cd "$(dirname "$0")"

# 0. ensure the LaTeX is in sync with SEDE.md (single source of truth)
bash build_tex.sh >/dev/null

rm -rf arxiv && mkdir -p arxiv

# 1. manuscript: prepend an engine hint, flatten figure paths ../../output[/paper]/figX -> figX
{ echo '% arXiv: process with XeLaTeX (fontspec + STIX Two Text; references inline, no BibTeX).'
  sed -E 's#\{\.\./\.\./output/(paper/)?#{#g' tex/SEDE.tex
} > arxiv/SEDE.tex

# 2. the 11 vector figures the manuscript references, copied as flat basenames
FIGS=(
  paper/fig1_mechanism paper/fig2_inputs paper/fig3_datafit paper/fig4_eos
  paper/fig5_oos paper/fig6_forecast
  fig_cmb_earlyde fig_F5_calibration fig_e1_mechanism fig_delta_orthogonality fig_syk_scrambling
)
for f in "${FIGS[@]}"; do
  cp "../output/$f.pdf" "arxiv/$(basename "$f").pdf"
done

# 3. human-readable submission note (not processed by arXiv)
cat > arxiv/ARXIV_SUBMIT.txt <<'RDM'
SEDE — arXiv source package
  Top file : SEDE.tex   (compile with XeLaTeX; run twice for the table of contents)
  Figures  : 11 vector PDFs (flat, ./), referenced extension-less
  Refs     : inline in the manuscript (no BibTeX, no .bbl)
  Fonts    : STIX Two Text via fontspec (ships with TeX Live / arXiv)
If arXiv asks for a processor, choose XeLaTeX.
RDM

# 4. tarball
( cd arxiv && tar czf ../SEDE_arxiv.tar.gz \
    SEDE.tex ARXIV_SUBMIT.txt \
    fig1_mechanism.pdf fig2_inputs.pdf fig3_datafit.pdf fig4_eos.pdf \
    fig5_oos.pdf fig6_forecast.pdf fig_cmb_earlyde.pdf fig_F5_calibration.pdf \
    fig_e1_mechanism.pdf fig_delta_orthogonality.pdf fig_syk_scrambling.pdf )

echo "wrote paper/SEDE_arxiv.tar.gz ($(tar tzf SEDE_arxiv.tar.gz | wc -l | tr -d ' ') files)"
