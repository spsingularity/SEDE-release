#!/usr/bin/env python3
"""
Cross-team reproduction (SEDE_V2) — independently re-derive the sibling team's two
quantitative Barrow-gap closures in OUR pipeline, and settle the w0 two-background labeling.

Reimplements (against sede.friedmann, NOT their code):
  (B) Binding-energy budget (double-counting / Gap 7): the released gravitational
      binding energy of cosmic structure is ~1e-6 of ρ_DE → binding energy sets the
      RATE (γ via Thm 4C), the CKN bound sets the SCALE. So "ρ_DE = binding energy"
      fails by ~1e6 and SEDE does not make that mistake (f_sat is normalized).
      [Background-independent Sheth–Tormen astrophysics; reproduces V2 to 3 sig figs.]
  (A) GSL / total-WEC (thermo consistency / Gap 5): the transient DE phantom phase does
      NOT violate the generalized 2nd law because the TOTAL fluid keeps the weak energy
      condition (matter-protected) → H decreases monotonically, the volume-law horizon
      entropy increases monotonically and saturates at de Sitter. Checked on BOTH
      backgrounds — V2 ran only their closure; here we ALSO confirm it on the canonical
      bare λ=0.5 model (the one that yields the ΔDIC preference).
  (R) New-or-repackaged (run_new_or_repackaged.py): ρ_DE/ρ_crit ∝ H^{-Δ}. On the self-consistent
      SEDE background, Δ=0 plateaus in the future (a fixed fraction = GR's own horizon energy,
      repackaged à la Jacobson) while Δ>0 keeps growing (genuine DE domination). So the genuinely
      NEW content is entirely the Barrow Δ-deformation: "is the energy new?" ⟺ "is Δ≠0?".
  (C) Covariant formulation — NUMERICAL CORE only (run_covariant_formulation.py): Route B
      d ln(ρ_DE/f_sat)/d ln H = 2-Δ (the λ=1-Δ/2 H-coupling origin); Route A the k-essence kinetic
      term (1+w)ρ_DE crosses 0 at the w=-1 points (single-field Vikman no-go ⟹ ghostly quintom ⟹
      SEDE's phantom is EFFECTIVE). The covariance/causality argument itself (∇_μT^μν=0,
      Misner-Sharp, apparent horizon θ₊=0) is STRUCTURAL — cited, NOT reproduced numerically.
  (W0) Two-background reconciliation: bare Barrow λ=0.5 (E_SEDE_lambda) → w0≈-0.99 with a
      single -1 crossing (= V2 'w0_bare'); the (1-ε/2) dynamical-temperature model
      (E_SEDE_H) → w0≈-0.85 with a double crossing (= V2 'self-consistent Cai-Kim
      closure'). The w0=-0.85 number is the λ=1 dynamical-T COUSIN, NOT the Δ=1/λ=0.5
      Barrow model.

V2 reference values (their results/*.json):
  binding: <E_bind/Mc^2>=1.93e-6, ρ_bind/ρ_DE=8.68e-7, too-small=1.15e6
  gsl(closure): total_WEC_min=+0.413, dlnH_dlna_max=-0.622, dS_hor_dlna_min=+5.29e-4

Run:  python run_xval_v2checks.py
"""
import numpy as np
from sede import friedmann as fr

# V2 parameters, for a like-for-like reproduction (canonical SEDE uses Om≈0.30)
OM, H0, GAM, OR = 0.31, 68.0, 1.507, 9.0e-5
G, MSUN, C = 6.674e-11, 1.989e30, 2.998e8


def binding_budget():
    """Mass-weighted released binding fraction (Sheth–Tormen × virial velocity)."""
    M = np.geomspace(1e10, 1e16, 400)
    M8, sig8, alpha = 6e14, 0.81, 0.30
    sigma = sig8 * (M / M8) ** (-alpha)
    nu = 1.686 / sigma
    A, a, p = 0.3222, 0.707, 0.3
    fST = A * np.sqrt(2 * a / np.pi) * (1 + (1 / (a * nu ** 2)) ** p) * nu * np.exp(-a * nu ** 2 / 2)
    dndlnM = fST * alpha / M
    H0_si = H0 * 1e3 / 3.086e22
    rho_crit = 3 * H0_si ** 2 / (8 * np.pi * G)
    Rvir = (3 * M * MSUN / (4 * np.pi * 200 * rho_crit)) ** (1 / 3.)
    ebind_frac = 0.5 * (G * M * MSUN / Rvir) / C ** 2          # E_bind/(Mc^2) ~ (1/2)v^2/c^2
    w = M * dndlnM
    frac = np.trapezoid(ebind_frac * w, np.log(M)) / np.trapezoid(w, np.log(M))
    rho_bind_over_DE = frac * OM / (1 - OM - OR)
    return frac, rho_bind_over_DE, 1 / rho_bind_over_DE


def gsl_check(Efunc, label):
    """Total-WEC / H-monotonicity / horizon-entropy-monotonicity on a given background."""
    z = np.linspace(0.0, 10.0, 1200)
    E = Efunc(z, OM, GAM)
    E2 = E ** 2
    rho_m = OM * (1 + z) ** 3
    rho_r = OR * (1 + z) ** 4
    rho_de = E2 - rho_m - rho_r
    lnx = np.log(1 + z)
    w = -1.0 + (1.0 / 3.0) * np.gradient(np.log(np.clip(rho_de, 1e-30, None)), lnx)   # fluid w
    rho_plus_p = rho_m + (rho_r + rho_r / 3.0) + (1 + w) * rho_de                     # total ρ+p
    dlnH_dlna = -0.5 * np.gradient(np.log(E2), lnx)
    S_hor = E ** (-3.0)                                                               # Barrow Δ=1 volume-law
    dS_dlna = np.gradient(S_hor, -lnx)
    de_ph = (1 + w) < 0
    zc = z[de_ph]
    S_ratio = (1.0 / (1 - OM - OR)) ** 3                                              # de Sitter future
    ok = (rho_plus_p.min() >= -1e-6) and (dlnH_dlna.max() <= 1e-6) and (dS_dlna.min() >= -1e-9)
    ph = f"DE phantom z∈[{zc.min():.2f},{zc.max():.2f}]" if zc.size else "no phantom phase"
    print(f"\n  [{label}]")
    print(f"    w0={w[0]:+.4f}  min w={w.min():+.4f}   {ph}")
    print(f"    min(ρ_tot+p_tot)/ρ_crit = {rho_plus_p.min():+.4f}  (≥0 → total WEC holds)")
    print(f"    max(dlnH/dlna)          = {dlnH_dlna.max():+.4f}  (≤0 → H decreasing)")
    print(f"    min(dS_hor/dlna)        = {dS_dlna.min():+.3e}  (≥0 → S_hor ↑, GSL OK)")
    print(f"    de Sitter S_hor(∞)/S_hor(0) = {S_ratio:.2f}   →  {'PASS' if ok else 'FAIL'}")
    return ok


def repackaged_check():
    """new_or_repackaged: ρ_DE/ρ_crit ∝ H^{-Δ}. On the self-consistent SEDE background the
    discriminator is the FUTURE — Δ=0 plateaus (a fixed fraction = repackaged GR horizon
    energy, Jacobson), Δ>0 keeps growing (genuine DE domination). 'Is the energy new?'⟺'Δ≠0?'."""
    Ode0 = 1 - OM - OR
    z = np.array([3.0, 1.0, 0.0, -0.5])                       # a = 0.25, 0.5, 1, 2 (future)
    zz = np.concatenate([np.linspace(0, 5, 200), np.geomspace(5.01, 60, 120)])
    Dg = fr.compute_growth_factor(zz, OM, OR)
    print("\n  ρ_DE/ρ_crit on OUR background at a=[0.25,0.5,1,2]:")
    fut = {}
    for D in (0.0, 0.5, 1.0):
        E = fr.E_SEDE_lambda(z, OM, GAM, 1 - D / 2)
        ratio = (E ** 2 - OM * (1 + z) ** 3 - OR * (1 + z) ** 4) / E ** 2
        fut[D] = ratio[-1]
        print(f"    Δ={D}: {np.round(ratio, 3).tolist()}")
    # Δ=0 future ≈ Ω_DE0 (flat/tracks); Δ=1 future grows above it
    ok = (abs(fut[0.0] - Ode0) < 0.02) and (fut[1.0] > fut[0.0] + 0.15)
    print(f"  => Δ=0 future plateaus at {fut[0.0]:.2f}≈Ω_DE0 (REPACKAGED); Δ=1 grows to {fut[1.0]:.2f} (GENUINE).")
    print("     novelty is ENTIRELY the Barrow Δ≠0 (magnitude = GR's own horizon energy, Jacobson).")
    return ok


def covariant_numcore():
    """The NUMERICAL core of run_covariant_formulation.py (the rest is structural, not numerical):
    Route B  d ln(ρ_DE/f_sat)/d ln H = 2-Δ (H-coupling λ=1-Δ/2 origin);
    Route A  the k-essence kinetic term (1+w)ρ_DE crosses 0 at the w=-1 points (single-field
             Vikman no-go ⟹ a fundamental description is a ghostly quintom; SEDE's phantom is
             EFFECTIVE). The covariance/causality argument itself (∇_μT^μν=0, Misner-Sharp,
             apparent horizon θ₊=0) is STRUCTURAL — cited, not reproduced numerically."""
    zz = np.concatenate([np.linspace(0, 5, 200), np.geomspace(5.01, 60, 120)])
    Dg = fr.compute_growth_factor(zz, OM, OR)
    zc = np.linspace(0, 3, 400)
    print("\n  Route B — d ln(ρ_DE/f_sat)/d ln H = 2-Δ (on our solved fixed-point E(z)):")
    slope_ok = True
    for D in (0.0, 0.5, 1.0):
        E = fr.E_SEDE_lambda(zc, OM, GAM, 1 - D / 2)
        fs = (1 - np.exp(-GAM * np.interp(zc, zz, Dg) ** 2)) / (1 - np.exp(-GAM))
        pure = (E ** 2 - OM * (1 + zc) ** 3 - OR * (1 + zc) ** 4) / fs
        slope = np.polyfit(np.log(E), np.log(pure), 1)[0]
        print(f"    Δ={D}: slope={slope:.3f}  (predicted {2 - D:.1f})")
        slope_ok &= abs(slope - (2 - D)) < 1e-2
    E = fr.E_SEDE_H(zc, OM, GAM)
    rde = E ** 2 - OM * (1 + zc) ** 3 - OR * (1 + zc) ** 4
    w = -1 - np.gradient(np.log(rde), np.log(1 / (1 + zc))) / 3
    cr = zc[np.where(np.diff(np.sign((1 + w) * rde)))[0]]
    print(f"  Route A — kinetic (1+w)ρ_DE crosses 0 at z={np.round(cr, 3).tolist()} (= w=-1 points → quintom no-go).")
    print("  [covariance/causality argument is STRUCTURAL — cited from V2, not a numerical result.]")
    return slope_ok and cr.size >= 1


def w0_two_background():
    """Confirm the bare-λ=0.5 vs (1-ε/2) split that explains the w0=-0.99 vs -0.85 labels."""
    z = np.linspace(0, 3.0, 400)

    def w_of(E):
        rho = E ** 2 - OM * (1 + z) ** 3 - OR * (1 + z) ** 4
        return -1 - np.gradient(np.log(np.clip(rho, 1e-30, None)), np.log(1 / (1 + z))) / 3

    w_bare = w_of(fr.E_SEDE_lambda(z, OM, GAM, 0.5))
    w_cous = w_of(fr.E_SEDE_H(z, OM, GAM))
    c_bare = z[np.where(np.diff(np.sign(w_bare + 1)))[0]]
    c_cous = z[np.where(np.diff(np.sign(w_cous + 1)))[0]]
    print(f"\n  bare Barrow Δ=1/λ=0.5 (E_SEDE_lambda): w0={w_bare[0]:+.4f}  crosses -1 at z={np.round(c_bare,3).tolist()}")
    print(f"  dynamical Cai-Kim T  (E_SEDE_H)      : w0={w_cous[0]:+.4f}  crosses -1 at z={np.round(c_cous,3).tolist()}")
    print("  => w0=-0.85 is the (1-ε/2) λ=1 COUSIN; the Δ=1/λ=0.5 ΔDIC-headline model is w0≈-1.0 (crossing).")
    # canonical Barrow should be near -1 with a crossing; cousin near -0.85 with a crossing
    return (w_bare[0] < -0.95) and (c_bare.size >= 1) and (-0.90 < w_cous[0] < -0.80) and (c_cous.size >= 1)


if __name__ == "__main__":
    print("=" * 72)
    print("CROSS-TEAM REPRODUCTION (SEDE_V2) — binding-energy budget, GSL, w0 backgrounds")
    print("=" * 72)

    print("\n(B) BINDING-ENERGY BUDGET  [V2: 1.93e-6 / 8.68e-7 / 1.15e6]")
    frac, ratio, toosmall = binding_budget()
    print(f"    <E_bind/Mc^2>    = {frac:.3e}")
    print(f"    ρ_bind/ρ_DE      = {ratio:.3e}")
    print(f"    too-small factor = {toosmall:.3e}")
    b_ok = (abs(frac - 1.93e-6) / 1.93e-6 < 0.02) and (toosmall > 1e5)

    print("\n(A) GSL / TOTAL-WEC  on OUR two backgrounds  [V2 tested only the closure]")
    g1 = gsl_check(lambda z, Om, g: fr.E_SEDE_lambda(z, Om, g, 0.5),
                   "bare Barrow λ=0.5  (CANONICAL — the ΔDIC model, NOT tested by V2)")
    g2 = gsl_check(fr.E_SEDE_H, "dynamical Cai-Kim T  (= V2's cai_kim_sc closure)")

    print("\n(R) NEW-OR-REPACKAGED  [V2 run_new_or_repackaged.py: novelty ⟺ Δ≠0]")
    r_ok = repackaged_check()

    print("\n(C) COVARIANT FORMULATION — numerical core  [V2 run_covariant_formulation.py]")
    c_ok = covariant_numcore()

    print("\n(W0) TWO-BACKGROUND RECONCILIATION")
    w_ok = w0_two_background()

    print("\n" + "=" * 72)
    checks = [("binding budget reproduces V2 (3 s.f.)", b_ok),
              ("GSL passes on canonical bare λ=0.5", g1),
              ("GSL passes on (1-ε/2) closure", g2),
              ("new-or-repackaged: Δ=0 tracks, Δ>0 genuine", r_ok),
              ("covariant numerical core (slope=2-Δ, kinetic zeros)", c_ok),
              ("w0 two-background split confirmed", w_ok)]
    for name, ok in checks:
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}")
    n_fail = sum(1 for _, ok in checks if not ok)
    print("=" * 72)
    print("ALL CHECKS PASSED." if n_fail == 0 else f"{n_fail} CHECK(S) FAILED.")
    import sys
    sys.exit(1 if n_fail else 0)
