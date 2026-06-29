"""Regenerate every paper figure as a vector PDF (alongside the existing PNGs).

matplotlib writes a true vector PDF when handed a ``*.pdf`` path, so this driver
monkeypatches ``Figure.savefig`` to emit a ``*.pdf`` next to every ``*.png`` it
writes, then runs the four figure-producing scripts as ``__main__``. The LaTeX
manuscript (paper/tex/main.tex) references the figures without an extension via
``\\graphicspath``, so pdflatex picks up the PDFs automatically.

Run from anywhere:
    python paper/make_vector_figures.py
"""
import os
import sys
import runpy

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)                       # for `import sede`
sys.path.insert(0, os.path.join(ROOT, "scripts"))  # for flat sibling imports between drivers
os.chdir(ROOT)  # the run-scripts write to relative output/ paths

import matplotlib
matplotlib.use("Agg")
import matplotlib.figure as mfig

_orig_savefig = mfig.Figure.savefig


def _savefig(self, fname, *args, **kwargs):
    _orig_savefig(self, fname, *args, **kwargs)
    if isinstance(fname, (str, os.PathLike)):
        s = os.fspath(fname)
        if s.lower().endswith(".png"):
            pdf = s[:-4] + ".pdf"
            # dpi is irrelevant for vector output; keep bbox_inches etc.
            kw = {k: v for k, v in kwargs.items() if k != "dpi"}
            _orig_savefig(self, pdf, *args, **kw)
            print(f"  + vector  {pdf}")


mfig.Figure.savefig = _savefig

SCRIPTS = [
    "paper/figures.py",                    # fig1_mechanism .. fig6_forecast
    "scripts/make_paper_figures.py",       # fig_F5_calibration (+ fig_F1, fig_F3, appendix)
    "scripts/run_cmb_earlyde.py",          # fig_cmb_earlyde        (Fig 5)
    "scripts/run_e1_figure.py",            # fig_e1_mechanism       (Fig 9)
    "scripts/run_delta_orthogonality.py",  # fig_delta_orthogonality (Fig 10)
    "scripts/run_chr_experiments.py",      # fig_chr_experiments
    "scripts/run_syk_scrambling.py",       # fig_syk_scrambling     (Fig G1)
]

for sc in SCRIPTS:
    print(f"=== {sc} ===")
    try:
        runpy.run_path(os.path.join(ROOT, sc), run_name="__main__")
    except SystemExit:
        pass
    except Exception as e:  # one bad script should not block the rest
        print(f"  {sc} FAILED: {e}")

print("\nDone. Vector PDFs written next to the PNGs in output/ and output/paper/.")
