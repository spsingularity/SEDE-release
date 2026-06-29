#!/usr/bin/env bash
# Build SEDE.pdf from SEDE.md.
#   - xelatex engine (Unicode text: Greek, dashes, sub/superscripts via STIX Two Text)
#   - custom template.tex = pandoc default with unicode-math swapped for plain fontspec
#     (so the informal inline math notation in prose, e.g. ∫ ∝ ⟹, maps cleanly via
#      header.tex's \newunicodechar instead of being forced into math mode)
#   - header.tex maps the Unicode math operators used in prose to LaTeX math commands
#   - superscript/subscript markdown extensions disabled so a bare '^' is escaped, not parsed
# Figures live in ../output/ (regenerate via make_all_figures.py and
# make_paper_figures.py). Run from the paper/ directory:  bash build_pdf.sh
set -e
cd "$(dirname "$0")"
# Assemble the build markdown: wrap informal inline math (mdmath) + lift title/abstract into
# YAML metadata and add the author (makedoc), so --toc lands after the abstract.
# temp lives in paper/ so the figures' ../output/ relative paths still resolve.
python3 makedoc.py SEDE.md .SEDE.proc.md
trap 'rm -f .SEDE.proc.md' EXIT
pandoc -f markdown-superscript-subscript .SEDE.proc.md -o SEDE.pdf \
  --pdf-engine=xelatex \
  --template=template.tex \
  -H header.tex \
  -V mainfont="STIX Two Text" \
  -V geometry:margin=1in \
  -V fontsize=10pt \
  -V colorlinks=true -V linkcolor=RoyalBlue -V urlcolor=RoyalBlue \
  --toc --toc-depth=2 \
  2> build.log || { echo "BUILD FAILED — see paper/build.log"; grep -A3 '! ' build.log | head; exit 1; }
echo "built paper/SEDE.pdf"
grep -ci "missing character" build.log | awk '{if($1>0) print "note: "$1" missing-glyph warnings (see build.log)"; else print "no missing glyphs"}'
