"""
Observational test of the STRUCTURE → DARK ENERGY link, on real data:
  ΛCDM (no structure)  vs  SEDE (DE ∝ gravitational growth D)  vs  Gough IDE (DE ∝ baryonic SMD^0.5).

Two regimes:
  (A) SIMPLE  : DESI DR2 BAO + Gold-2018 fσ8, with a FREE sound horizon (α=c/H0r_d marginalised).
  (B) JOINT   : DESI BAO (PHYSICAL r_d from the sound-horizon integral) + Pantheon+ SN
                + compressed CMB (shift R and acoustic scale ℓ_A) + fσ8. This removes the
                free-r_d freedom and is much closer to the repo's CMB-pinned joint conditions.
Nuisances (H0·r_d in A; SN offset M_B; growth σ8) are marginalised analytically. H0 fixed at 67.4
in (B) — the compressed BAO/CMB observables are ~H0-independent (ratios), so this is a fair Ω_m fit.
"""
import numpy as np
from sede import data_loader as dl
from sede import friedmann as fr
from sede import cmb as cmbmod

C = 299792.458
OMR = 9.0e-5
Z_DRAG, Z_STAR = 1059.94, 1089.92
H0_FIX = 67.4
R_PL, R_PL_E = 1.7502, 0.0046          # Planck 2018 shift parameter (load_planck)
LA_PL, LA_PL_E = 301.71, 0.09          # Planck 2018 acoustic scale

z_b, types, mean, cov = dl.load_desi_dr2()
Cinv = np.linalg.inv(cov)
z_f, fs8_obs, fs8_err = (np.asarray(x, float) for x in dl.load_fss8())
z_sn, m_sn, Csn = dl.load_pantheon_plus()
Lsn = np.linalg.cholesky(Csn)
u1 = np.linalg.solve(Lsn, np.ones(len(z_sn)))      # for analytic M marginalisation
den1 = u1 @ u1


def sfrd(zp):
    return 0.015 * (1 + zp)**2.7 / (1 + ((1 + zp) / 2.9)**5.6)


def smd_curve(Om):
    zh = np.linspace(0, 30, 8000)
    El = np.sqrt(Om * (1 + zh)**3 + OMR * (1 + zh)**4 + (1 - Om - OMR))
    integ = sfrd(zh) / (El * (1 + zh))
    cum = np.cumsum((integ[::-1] * np.gradient(zh)[::-1]))[::-1]
    return zh, cum


def E_lcdm(z, Om):
    z = np.atleast_1d(z)
    return np.sqrt(Om * (1 + z)**3 + OMR * (1 + z)**4 + (1 - Om - OMR))


def E_sede(z, Om):
    return fr.E_SEDE_barrow(np.atleast_1d(z), Om, 1.4964, Delta=1.0)


def E_gough(z, Om):
    z = np.atleast_1d(z)
    zh, cum = smd_curve(Om)
    g = np.sqrt(np.clip(np.interp(z, zh, cum), 0, None) / cum[0])
    return np.sqrt(Om * (1 + z)**3 + OMR * (1 + z)**4 + (1 - Om - OMR) * g)


def cum_invE(Efun, Om, zmax, n=4000):
    zg = np.linspace(0, zmax, n)
    invE = 1.0 / Efun(zg, Om)
    Ig = np.concatenate([[0], np.cumsum(0.5 * (invE[1:] + invE[:-1]) * np.diff(zg))])
    return zg, Ig


# ───────────────────────── (A) simple: BAO+fσ8, free r_d ─────────────────────────
def bao_chi2_freerd(Efun, Om):
    zg, Ig = cum_invE(Efun, Om, 2.6)
    s = np.zeros(len(z_b))
    for k, (z, t) in enumerate(zip(z_b, types)):
        I = np.interp(z, zg, Ig); E = float(np.atleast_1d(Efun(np.array([z]), Om))[0])
        s[k] = I if t == 'DM/rd' else (1.0 / E if t == 'DH/rd' else (z * I**2 / E)**(1 / 3))
    sCs, sCm, mCm = s @ Cinv @ s, s @ Cinv @ mean, mean @ Cinv @ mean
    return mCm - sCm**2 / sCs


def fs8_chi2(Efun, Om):
    D, f = fr.compute_growth_model(z_f, Om, lambda z: Efun(z, Om))
    shape = f * D; w = 1.0 / fs8_err**2
    s8 = np.sum(w * shape * fs8_obs) / np.sum(w * shape**2)
    return float(np.sum(w * (fs8_obs - s8 * shape)**2))


# ───────────────────────── (B) joint: physical r_d + SN + CMB ────────────────────
def bao_chi2_phys(Efun, Om, H0):
    rd = cmbmod.sound_horizon_rs(Om, H0, 0.02237, Z_DRAG)
    zg, Ig = cum_invE(Efun, Om, 2.6)
    pred = np.zeros(len(z_b))
    for k, (z, t) in enumerate(zip(z_b, types)):
        I = np.interp(z, zg, Ig); E = float(np.atleast_1d(Efun(np.array([z]), Om))[0])
        DM = (C / H0) * I; DH = C / (H0 * E)
        val = DM if t == 'DM/rd' else (DH if t == 'DH/rd' else (z * DM**2 * DH)**(1 / 3))
        pred[k] = val / rd
    d = mean - pred
    return float(d @ Cinv @ d)


def sn_chi2(Efun, Om):
    zg, Ig = cum_invE(Efun, Om, 2.6)
    I = np.interp(z_sn, zg, Ig)
    shape = 5.0 * np.log10((1 + z_sn) * np.maximum(I, 1e-8))      # M_B + 5log(c/H0) marginalised
    d = m_sn - shape
    ud = np.linalg.solve(Lsn, d)
    return float(ud @ ud - (ud @ u1)**2 / den1)


def cmb_chi2(Efun, Om, H0):
    zg, Ig = cum_invE(Efun, Om, Z_STAR, n=6000)
    I_star = np.interp(Z_STAR, zg, Ig)
    DM_star = (C / H0) * I_star
    R = np.sqrt(Om) * H0 * DM_star / C
    rs_star = cmbmod.sound_horizon_rs(Om, H0, 0.02237, Z_STAR)
    lA = np.pi * DM_star / rs_star
    return ((R - R_PL) / R_PL_E)**2 + ((lA - LA_PL) / LA_PL_E)**2


def scan(Efun, joint):
    best = (1e18, None, {})
    for Om in np.linspace(0.26, 0.36, 41):
        if joint:
            cb = bao_chi2_phys(Efun, Om, H0_FIX); cs = sn_chi2(Efun, Om)
            cc = cmb_chi2(Efun, Om, H0_FIX); cg = fs8_chi2(Efun, Om)
            tot = cb + cs + cc + cg; parts = dict(BAO=cb, SN=cs, CMB=cc, fs8=cg)
        else:
            cb = bao_chi2_freerd(Efun, Om); cg = fs8_chi2(Efun, Om)
            tot = cb + cg; parts = dict(BAO=cb, fs8=cg)
        if tot < best[0]:
            best = (tot, Om, parts)
    return best


def main():
    models = [("ΛCDM ", E_lcdm), ("SEDE ", E_sede), ("Gough", E_gough)]
    for joint, label, npts, npar in [(False, "(A) SIMPLE — DESI BAO + fσ8, FREE r_d", 29, 3),
                                     (True, "(B) JOINT — BAO(phys r_d) + Pantheon+ SN + CMB(R,ℓ_A) + fσ8",
                                      13 + len(z_sn) + 2 + 16, 3)]:
        print("=" * 80); print(label, f"   [{npts} data pts]"); print("=" * 80)
        base = None
        for name, Efun in models:
            tot, Om, parts = scan(Efun, joint)
            if name.strip() == "ΛCDM":
                base = tot
            pj = "  ".join(f"{k}={v:.1f}" for k, v in parts.items())
            print(f"  {name}: total={tot:8.1f}  Ω_m={Om:.3f}   [{pj}]   Δχ²(vs ΛCDM)={tot-base:+.2f}")
        print()


if __name__ == "__main__":
    main()
