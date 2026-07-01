"""
make_foundations_figures.py — figures for Paper II (SEDE_foundations.md).

Produces a coherent 5-figure set into output/:
  foundations_fig1_reduction.png   — the reduction-chain schematic (the spine)
  (fig2 = criticality_soc.png, fig3 = deposition_drive.png — reused, already built)
  foundations_fig4_socsignature.png — mean-field SOC signature (avalanche + slowing)
  foundations_fig5_inferences.png   — ensemble inequivalence (C<0) + GREA viscosity

Run from repo root:  python paper/make_foundations_figures.py
"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
OUT = os.path.join(ROOT, "output")
os.makedirs(OUT, exist_ok=True)

GREEN, BLUE, ORANGE, GREY = "#2e7d32", "#1565c0", "#e65100", "#555555"


# ---------------------------------------------------------------------------
# Fig 1 — the reduction chain (schematic)
# ---------------------------------------------------------------------------
def fig_reduction():
    fig, ax = plt.subplots(figsize=(9.2, 6.2))
    ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")

    def box(x, y, w, h, text, color, fc=None):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.04,rounding_size=0.12",
                     linewidth=1.6, edgecolor=color, facecolor=fc or "white"))
        ax.text(x + w/2, y + h/2, text, ha="center", va="center", fontsize=8.6, color="black")

    def arrow(x1, y1, x2, y2, color=GREY):
        ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
                     mutation_scale=13, linewidth=1.3, color=color))

    box(3.3, 9.0, 3.4, 0.8, "Volume-law postulate\n(horizon counts its volume, Δ=1)", GREY, "#f2f2f2")
    # four sub-claims
    subs = [("state\nDERIVED (SYK)", GREEN, 0.3),
            ("form\nREDUCED (thermal)", BLUE, 2.7),
            ("scale\nDERIVED (CKN)", GREEN, 5.1),
            ("count\nOPEN — the residue", ORANGE, 7.5)]
    for txt, col, x in subs:
        box(x, 7.4, 2.2, 0.85, txt, col)
        arrow(5.0, 9.0, x + 1.1, 8.25, col)
    # chain from count
    arrow(8.6, 7.4, 8.6, 6.7, ORANGE)
    box(5.6, 5.9, 4.0, 0.8, "count ⟺ connectivity ⟸ range:\ngravity 1/r (α≤d) ⇒ non-additive ⇒ volume", BLUE)
    arrow(7.6, 5.9, 6.2, 5.3, BLUE)
    box(3.0, 4.5, 5.0, 0.8, "driven-NESS conjecture  (equilibrium→area, driven→volume)\nmaps: J≥J_c (jcoup) + deposition→drive (deposdrv)", BLUE)
    arrow(5.5, 4.5, 5.5, 3.9, BLUE)
    box(2.8, 3.1, 5.4, 0.8, "ratchet through spinodal at z*; GSL-locked at volume;\nrobust over a broad untuned basin (socatt)", BLUE)
    arrow(5.5, 3.1, 5.5, 2.5, GREY)
    box(3.4, 1.7, 4.2, 0.8, "RESIDUE: is the microcanonical dS\nstate count volume?  (dS holography)", ORANGE)
    arrow(5.5, 1.7, 5.5, 1.1, GREEN)
    box(2.9, 0.3, 5.2, 0.7, "DECIDED EMPIRICALLY: Δ via DESI DR3 + Euclid (~0.09)", GREEN, "#eaf5ea")

    # legend
    for i, (c, t) in enumerate([(GREEN, "derived / decided"), (BLUE, "reduced"), (ORANGE, "open")]):
        ax.add_patch(FancyBboxPatch((0.2, 0.3 + i*0.5), 0.3, 0.3, boxstyle="round,pad=0.02",
                     edgecolor=c, facecolor="white", linewidth=1.6))
        ax.text(0.6, 0.45 + i*0.5, t, fontsize=7.5, va="center")
    ax.set_title("The reduction of SEDE's one postulate", fontsize=11)
    fig.tight_layout()
    p = os.path.join(OUT, "foundations_fig1_reduction.png")
    fig.savefig(p, dpi=140); plt.close(fig); print("wrote", p)


# ---------------------------------------------------------------------------
# Fig 4 — mean-field SOC signature
# ---------------------------------------------------------------------------
def fig_socsignature():
    rng = np.random.default_rng(0)
    def avalanche(rng, cap=200000):
        active, size = 1, 0
        while active > 0 and size < cap:
            b = rng.poisson(active); size += active; active = b
        return size
    sizes = np.array([avalanche(rng) for _ in range(80000)])
    sizes = sizes[(sizes > 0) & (sizes < 200000)]
    edges = np.unique(np.round(np.logspace(0, np.log10(sizes.max()), 24)).astype(int))
    cnt, _ = np.histogram(sizes, bins=edges)
    ctr = np.sqrt(edges[:-1]*edges[1:]); pdf = cnt/np.diff(edges)/cnt.sum()
    ok = pdf > 0

    from sede import chr_mechanism as chr
    from scipy.integrate import solve_ivp
    J, TH = 6.0, 3.0; sig = lambda x: 1/(1+np.exp(-x))
    z_star = chr.transition_redshift()
    ng = np.linspace(-4, 1.5, 700); zg = np.exp(-ng)-1
    fz = np.clip(chr.f_eq(chr.control_variance(np.clip(zg, 0, None))), 1e-6, 1)
    def rhs(n, m):
        m = float(np.clip(m, 0, 1)); return 25*(-m+sig(J*m+float(np.interp(n, ng, fz))-TH))
    m = solve_ivp(rhs, (ng[0], ng[-1]), [sig(-TH)], t_eval=ng, rtol=1e-8, atol=1e-10).y[0]
    arg = J*m+fz-TH; lam = 1-J*sig(arg)*(1-sig(arg)); chi = 1/np.maximum(np.abs(lam), 1e-3)
    m_sp = 0.5 - np.sqrt(0.25 - 1.0/J)               # area-branch spinodal
    chi_area = np.where(m <= m_sp, chi, np.nan)       # show only the metastable area branch

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    ax[0].loglog(ctr[ok], pdf[ok], "o", ms=4, color=BLUE, label="critical branching")
    s = np.array([2, ctr[ok].max()/15]); ax[0].loglog(s, 0.5*s**-1.5, "k--", label=r"$P(s)\propto s^{-3/2}$")
    ax[0].set_xlabel("avalanche size s"); ax[0].set_ylabel("P(s)")
    ax[0].set_title("(a) mean-field SOC universality"); ax[0].legend(fontsize=8)
    sel = (zg > 0) & (zg < 3.5)
    ax[1].plot(zg[sel], chi_area[sel], color=ORANGE, lw=2)
    ax[1].axvline(z_star, ls=":", c="0.6", label=f"z*={z_star:.2f}")
    ax[1].set_xlim(3.5, 0); ax[1].set_xlabel("redshift z")
    ax[1].set_ylabel(r"susceptibility / $\tau_{ac}$  ($1/\lambda_{stab}$)")
    ax[1].set_title("(b) area-phase critical slowing-down at z*"); ax[1].legend(fontsize=8)
    ax[1].annotate("area phase\nloses stability", xy=(z_star, 0.85*np.nanmax(chi_area[sel])),
                   xytext=(2.3, 0.6*np.nanmax(chi_area[sel])), fontsize=7.5, color=GREY,
                   arrowprops=dict(arrowstyle="->", color="0.6"))
    fig.tight_layout()
    p = os.path.join(OUT, "foundations_fig4_socsignature.png")
    fig.savefig(p, dpi=140); plt.close(fig); print("wrote", p)


# ---------------------------------------------------------------------------
# Fig 5 — literature inferences: ensemble (C<0) + GREA viscosity
# ---------------------------------------------------------------------------
def fig_inferences():
    from sede.friedmann import E_SEDE_lambda, compute_growth_factor
    Om, Or, g, lam = 0.30, 9e-5, 1.4964, 0.5; ODE0 = 1-Om-Or
    z = np.linspace(0, 3, 300); H = E_SEDE_lambda(z, Om, g, lam, Or)
    S = 1/H**2; T = H/(2*np.pi); C = T*np.gradient(S, T)
    a = 1/(1+z); D = compute_growth_factor(z, Om, Or)
    f = np.clip((1-np.exp(-g*D**2))/(1-np.exp(-g)), 1e-9, 1)
    lna = np.log(a); eps = -np.gradient(np.log(H), lna); dlnf = np.gradient(np.log(f), lna)
    rho = ODE0*f*H**(2*lam); zeta = rho*dlnf/(9*H**2); w = -1+(1/3)*(2*lam*eps-dlnf)

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    ax[0].plot(z, C/S, color=BLUE, lw=2)
    ax[0].axhline(-2, ls="--", c="0.6", label="C = −2S (analytic)")
    ax[0].set_xlabel("redshift z"); ax[0].set_ylabel("C / S")
    ax[0].set_title("(a) horizon negative specific heat\n⇒ microcanonical (volume) selected")
    ax[0].legend(fontsize=8); ax[0].set_ylim(-2.6, 0)
    ax[0].axhspan(-2.6, 0, color="#fdecea", alpha=0.5, zorder=0)
    ax[0].text(1.5, -0.45, "C < 0  ⇒  canonical ensemble ill-defined\n(ensemble inequivalence, long-range)",
               ha="center", fontsize=7.8, color=GREY)
    ax2 = ax[1]; ln1 = ax2.plot(z, zeta, color=ORANGE, lw=2, label=r"$\zeta(z)$ (bulk viscosity)")
    ax2.set_xlabel("redshift z"); ax2.set_ylabel(r"effective bulk viscosity $\zeta$", color=ORANGE)
    ax2.tick_params(axis="y", labelcolor=ORANGE); ax2.set_xlim(3, 0)
    axb = ax2.twinx(); ln2 = axb.plot(z, w, color=GREY, lw=1.5, ls="--", label="w(z)")
    axb.axhline(-1, ls=":", c="0.7"); axb.set_ylabel("w(z)", color=GREY)
    ax2.set_title("(b) SEDE as structure-gated GREA\n(ζ≥0, peaks at gate-activation z≲0.5)")
    lns = ln1+ln2; ax2.legend(lns, [l.get_label() for l in lns], fontsize=8, loc="center right")
    fig.tight_layout()
    p = os.path.join(OUT, "foundations_fig5_inferences.png")
    fig.savefig(p, dpi=140); plt.close(fig); print("wrote", p)


if __name__ == "__main__":
    fig_reduction()
    fig_socsignature()
    fig_inferences()
    print("Paper II figures done (fig2=criticality_soc.png, fig3=deposition_drive.png reused).")
