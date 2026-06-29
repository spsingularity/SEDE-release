#!/usr/bin/env bash
# Re-sync the editable LaTeX source from the canonical SEDE.md (single source of truth).
#   - Same mdmath + pandoc(+template.tex,+header.tex) toolchain as build_pdf.sh, so it
#     compiles under xelatex exactly like the PDF build.
#   - Emits a STANDALONE, current tex/SEDE.tex (front matter + body + inline references).
#   - Rewrites figure paths to extension-less basenames under ../../output so xelatex
#     picks up the VECTOR *.pdf (not the *.png) — JCAP/arXiv-preferred.
# Run from paper/:  bash build_tex.sh
set -e
cd "$(dirname "$0")"
mkdir -p tex

# 1. assemble the build markdown (inline-math wrap + title/abstract -> YAML), as in build_pdf.sh
python3 makedoc.py SEDE.md .SEDE.proc.md
trap 'rm -f .SEDE.proc.md' EXIT

# 2. markdown -> standalone LaTeX (NO pdf-engine; keep the .tex). graphicspath added via header.
pandoc -f markdown-superscript-subscript .SEDE.proc.md -o tex/SEDE.tex \
  --standalone \
  --shift-heading-level-by=-1 \
  --template=template.tex \
  -H header.tex \
  -V mainfont="STIX Two Text" \
  -V geometry:margin=1in \
  -V fontsize=10pt \
  -V colorlinks=true -V linkcolor=RoyalBlue -V urlcolor=RoyalBlue \
  --toc --toc-depth=2

# 3. figure paths: pandoc writes ../output/...png (relative to paper/). From tex/ that is
#    ../../output/...; strip the .png so the vector .pdf is preferred by graphicspath/xelatex.
perl -0pi -e 's#\{\.\./output/([^}]+?)\.png\}#{../../output/\1}#g' tex/SEDE.tex

# 4. compile (xelatex twice: ToC + cross-refs). References are inline in SEDE.md (no bibtex).
( cd tex && xelatex -interaction=nonstopmode SEDE.tex >SEDE.build.log 2>&1 \
            && xelatex -interaction=nonstopmode SEDE.tex >>SEDE.build.log 2>&1 ) \
  || { echo "TEX BUILD FAILED — see tex/SEDE.build.log"; grep -A2 '^!' tex/SEDE.build.log | head -20; exit 1; }

echo "built paper/tex/SEDE.tex  ->  tex/SEDE.pdf"
grep -ci "missing character" tex/SEDE.build.log | awk '{if($1>0) print "note: "$1" missing-glyph warnings"; else print "no missing glyphs"}'
grep -cE '\\includegraphics' tex/SEDE.tex | awk '{print $1" figures (extension-less -> vector PDF)"}'
