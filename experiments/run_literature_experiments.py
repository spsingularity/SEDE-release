"""
Literature-inspired falsifiability experiments for SEDE, computed and compared to ΛCDM and to
Gough's Information Dark Energy. These CHARACTERISE how distinguishable SEDE is on each probe; they
do not by themselves confirm/refute it (the decisive tests remain Δ via DESI DR3/Euclid and the E1
mechanism). Sources: Gough 2022 (H(z) dip falsifiability); DESI 2024 physics-focused DE taxonomy;
Sandage-Loeb redshift drift (ELT-HIRES/SKA).

  E1: H(z) deviation vs ΛCDM   — SEDE vs Gough's published −2.2% @ z≈1.7 falsifiable signature
  E2: redshift drift Δv(z)     — z≈2–5 "redshift desert" (ELT/SKA)
  E3: emergent/thawing/mirage  — DE-fraction falloff + CPL (w0,wa) class

Output: console, results/literature_experiments.json, output/fig_literature_experiments.png
"""
import json, os
import numpy as np
from sede import friedmann as fr

OMR = 9.0e-5
C = 299792.458


def _smd(Om):
    zh = np.linspace(0, 30, 8000)
    El = np.sqrt(Om * (1 + zh)**3 + OMR * (1 + zh)**4 + (1 - Om - OMR))
    sfrd = 0.015 * (1 + zh)**2.7 / (1 + ((1 + zh) / 2.9)**5.6)
    integ = sfrd / (El * (1 + zh))
    return zh, np.cumsum((integ[::-1] * np.gradient(zh)[::-1]))[::-1]


def E_lcdm(z, Om): return np.sqrt(Om * (1 + z)**3 + OMR * (1 + z)**4 + (1 - Om - OMR))
def E_sede(z, Om): return fr.E_SEDE_barrow(np.atleast_1d(z), Om, 1.4964, Delta=1.0)
def E_gough(z, Om):
    zh, cum = _smd(Om); g = np.sqrt(np.clip(np.interp(z, zh, cum), 0, None) / cum[0])
    return np.sqrt(Om * (1 + z)**3 + OMR * (1 + z)**4 + (1 - Om - OMR) * g)


def rho_de(Efun, z, Om):           # ρ_DE/ρ_crit0 = E² − matter − rad
    return np.atleast_1d(Efun(z, Om))**2 - Om * (1 + z)**3 - OMR * (1 + z)**4


def main():
    Om = 0.30
    z = np.linspace(0, 5, 400)
    El, Es, Eg = E_lcdm(z, Om), E_sede(z, Om), E_gough(z, Om)
    out = {"Om": Om}

    # ── E1: H(z) deviation vs ΛCDM ──────────────────────────────────────────────
    dHs = 100 * (Es / El - 1); dHg = 100 * (Eg / El - 1)
    print("=" * 74)
    print("E1  H(z) deviation vs ΛCDM (matched Ω_m)  [Gough's falsifiable signature]")
    print("=" * 74)
    for zi in (0.5, 1.0, 1.7, 2.3, 3.0):
        i = abs(z - zi).argmin()
        print(f"   z={zi:.1f}:  SEDE {dHs[i]:+5.2f}%   Gough {dHg[i]:+5.2f}%")
    print(f"   extremum: SEDE {dHs.min():+.2f}% @ z={z[dHs.argmin()]:.2f}  |  "
          f"Gough {dHg.min():+.2f}% @ z={z[dHg.argmin()]:.2f}  (Gough 2022 quotes −2.2% @ z≈1.7)")
    print("   => ΛCDM(0) < SEDE(~−0.5%) < Gough(~−2.2%): SEDE ~5× subtler in H(z).")
    out["E1"] = dict(z=z.tolist(), dH_sede=dHs.tolist(), dH_gough=dHg.tolist(),
                     sede_min=float(dHs.min()), gough_min=float(dHg.min()))

    # ── E2: redshift drift Δv(z) over 20 yr ─────────────────────────────────────
    H0, dt = 70.0, 20.0
    fac = C * 1e5 * dt / 9.778e11   # → cm/s
    dv = lambda E: fac * (H0 - H0 * E / (1 + z))
    dvl, dvs, dvg = dv(El), dv(Es), dv(Eg)
    print("\n" + "=" * 74)
    print("E2  redshift drift Δv over 20 yr (cm/s)  [ELT-HIRES/SKA, z≈2–5]")
    print("=" * 74)
    for zi in (2.0, 3.0, 4.0, 5.0):
        i = abs(z - zi).argmin()
        print(f"   z={zi:.1f}:  ΛCDM {dvl[i]:+6.2f}  SEDE {dvs[i]:+6.2f}  Gough {dvg[i]:+6.2f}  "
              f"| SEDE−ΛCDM {dvs[i]-dvl[i]:+.3f}  SEDE−Gough {dvs[i]-dvg[i]:+.3f}")
    print("   => ELT precision ~1–3 cm/s: catches Gough−ΛCDM (~0.5–0.7), SEDE−ΛCDM (~0.1–0.2) below it.")
    out["E2"] = dict(z=z.tolist(), dv_lcdm=dvl.tolist(), dv_sede=dvs.tolist(), dv_gough=dvg.tolist())

    # ── E3: emergent/thawing/mirage classification ──────────────────────────────
    rs = rho_de(E_sede, z, Om); rs0 = rs[0]
    rg = rho_de(E_gough, z, Om); rg0 = rg[0]
    OmDE_s = rs / (Es**2); OmDE_g = rg / (Eg**2)
    # CPL w0,wa from a fit of w(z) over z<1
    zc = z[z < 1.0]
    ws = -1 - (1 / 3.) * (1 + zc) * np.gradient(np.log(np.maximum(rho_de(E_sede, zc, Om), 1e-6)), zc)
    a = 1 / (1 + zc)
    A = np.vstack([np.ones_like(a), 1 - a]).T
    w0, wa = np.linalg.lstsq(A, ws, rcond=None)[0]
    mirage = -3.66 * (1 + w0)
    print("\n" + "=" * 74)
    print("E3  emergent/thawing/mirage class  [DESI 2024 physics-focused DE]")
    print("=" * 74)
    print(f"   SEDE CPL: w0={w0:+.3f}, wa={wa:+.3f}  (mirage line wa≈{mirage:+.2f}; thawing wa<0,w0>−1)")
    print(f"   DE fraction left at z=2 / z=3:  SEDE ρ_DE/ρ_DE0 = {rs[abs(z-2).argmin()]/rs0:.3f} / "
          f"{rs[abs(z-3).argmin()]/rs0:.3f}   Gough = {rg[abs(z-2).argmin()]/rg0:.3f} / "
          f"{rg[abs(z-3).argmin()]/rg0:.3f}")
    cls = "EMERGENT" if (rs[abs(z-2).argmin()]/rs0 < 0.85 or w0 > -0.9) else "thawing/near-Λ"
    print(f"   => SEDE sits in the {cls} quadrant (DE largely intact at z~3 ⇒ closer to near-Λ/mirage,")
    print(f"      Gough falls off faster ⇒ more 'emergent').")
    out["E3"] = dict(w0=float(w0), wa=float(wa), mirage_wa=float(mirage),
                     sede_rhoDE_z3=float(rs[abs(z-3).argmin()]/rs0),
                     gough_rhoDE_z3=float(rg[abs(z-3).argmin()]/rg0))

    os.makedirs("results", exist_ok=True)
    json.dump(out, open("results/literature_experiments.json", "w"), indent=2)
    print("\nwrote results/literature_experiments.json")
    _fig(z, dHs, dHg, dvl, dvs, dvg, rs/rs0, rg/rg0)
    return out


def _fig(z, dHs, dHg, dvl, dvs, dvg, fs, fg):
    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    except Exception as e:
        print(f"(figure skipped: {e})"); return
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.3))
    ax[0].plot(z, dHs, "C3", label="SEDE"); ax[0].plot(z, dHg, "C0", label="Gough IDE")
    ax[0].axhline(0, color="k", lw=1, ls=":"); ax[0].plot(1.7, -2.2, "k*", ms=10, label="Gough 2022 quoted")
    ax[0].set_xlabel("z"); ax[0].set_ylabel("ΔH/H vs ΛCDM [%]"); ax[0].set_title("E1: H(z) deviation"); ax[0].legend(fontsize=8)
    ax[1].plot(z, dvs - dvl, "C3", label="SEDE − ΛCDM"); ax[1].plot(z, dvg - dvl, "C0", label="Gough − ΛCDM")
    ax[1].axhspan(-2, 2, color="gray", alpha=0.2, label="ELT ~±2 cm/s")
    ax[1].set_xlim(1.5, 5); ax[1].set_xlabel("z"); ax[1].set_ylabel("Δv difference [cm/s, 20 yr]")
    ax[1].set_title("E2: redshift drift"); ax[1].legend(fontsize=8)
    ax[2].plot(z, fs, "C3", label="SEDE"); ax[2].plot(z, fg, "C0", label="Gough"); ax[2].axhline(1, color="k", lw=1, ls=":")
    ax[2].set_xlabel("z"); ax[2].set_ylabel("ρ_DE(z)/ρ_DE(0)"); ax[2].set_title("E3: emergent character (DE falloff)"); ax[2].legend(fontsize=8)
    fig.suptitle("Literature-inspired falsifiability experiments — SEDE vs Gough vs ΛCDM", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95]); os.makedirs("output", exist_ok=True)
    fig.savefig("output/fig_literature_experiments.png", dpi=130); plt.close(fig)
    print("wrote output/fig_literature_experiments.png")


if __name__ == "__main__":
    main()
