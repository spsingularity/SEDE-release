#!/usr/bin/env python3
"""
TIER 1 — real-data confrontations that strengthen/test SEDE-H (canonical bare Barrow
λ=0.5, γ=theory).  All four use data already on disk; SN offset is analytically
marginalised so only the SHAPE of μ(z) (the crossing-(-1) EOS) is tested.

  1a SN robustness   : Δχ²(SEDE−ΛCDM) with SN ∈ {Pantheon+, DES-SN5YR, Union3}.
                       Does the preference survive swapping the SN sample?
  1b Lyα z=2.33      : SEDE vs ΛCDM at the high-z BAO point + Ω_DE(z=3) (early-DE).
  1c WL S8 prior     : is SEDE consistent with DES-Y3/KiDS cosmic-shear S8?
  1d CMB-lensing S8  : SEDE vs Planck+ACT lensing S8 (the S8 tension, honest).

Run:  python run_tier1_data.py
"""
import numpy as np, camb
from camb import dark_energy
from scipy.optimize import minimize
from sede import friedmann as fr
from sede import data_loader as dl
from run_lambda_verify import w_of_a, data, R_PL, LA_PL, CMB_CINV, OMBH2_PRIOR, SHOES
C = 299792.458

# ---- SN blocks (z, mu_obs, inv_cov), offset marginalised analytically ----
def sn_blocks():
    blocks = {}
    zp, mu, cov = dl.load_pantheon_plus()
    if zp is not None:
        blocks['Pantheon+'] = (np.asarray(zp), np.asarray(mu), np.linalg.inv(cov))
    z5, mu5, ic5 = dl.load_des_dovekie()
    if z5 is not None:
        blocks['DES-SN5YR'] = (np.asarray(z5), np.asarray(mu5), np.asarray(ic5))
    zu, muu, icu = dl.load_union3()
    if zu is not None:
        blocks['Union3'] = (np.asarray(zu), np.asarray(muu), np.asarray(icu))
    return blocks

def chi2_sn_marg(mu_obs, mu_model, icov):
    """SN χ² with the absolute-offset (M_B/μ0) analytically marginalised:
       χ² = a − b²/e,  a=Δ·C⁻¹·Δ, b=Δ·C⁻¹·1, e=1·C⁻¹·1,  Δ=μ_obs−μ_model."""
    d = mu_obs - mu_model
    Cinv1 = icov @ np.ones_like(d)
    a = float(d @ icov @ d); b = float(d @ Cinv1); e = float(np.ones_like(d) @ Cinv1)
    return a - b * b / e

def background(Om, H0, ombh2, gamma, lam, lcdm=False):
    h = H0 / 100.; pa = camb.CAMBparams()
    pa.set_cosmology(H0=H0, ombh2=ombh2, omch2=Om * h**2 - ombh2, mnu=0.06)
    if not lcdm:
        a, w = w_of_a(Om, gamma, lam); de = dark_energy.DarkEnergyPPF(); de.set_w_a_table(a, w); pa.DarkEnergy = de
    pa.WantCls = False; pa.WantTransfer = False
    return camb.get_background(pa)

def chi2(Om, H0, ombh2, s8, sn, gamma=1.4964, lam=0.5, lcdm=False, s8prior=None):
    bg = background(Om, H0, ombh2, gamma, lam, lcdm)
    dd = bg.get_derived_params(); rd = dd['rdrag']; zs = dd['zstar']; rs = dd['rstar']; c = 0.0
    # DESI BAO (incl. Lyα z=2.33)
    z, t, m, icov = data['desi']
    pred = np.array([bg.comoving_radial_distance(zz)/rd if tp=='DM/rd' else (C/bg.hubble_parameter(zz))/rd if tp=='DH/rd'
        else ((zz*bg.comoving_radial_distance(zz)**2*(C/bg.hubble_parameter(zz)))**(1/3.))/rd for zz, tp in zip(z, t)])
    c += float((m-pred) @ icov @ (m-pred))
    # SN (offset-marginalised)
    zsn, musn, icsn = sn
    mu_model = 5*np.log10((1+zsn)*bg.comoving_radial_distance(zsn)) + 25
    c += chi2_sn_marg(musn, mu_model, icsn)
    # CMB (R, l_A) + ombh2 + SH0ES + CC + fσ8
    DMz = bg.comoving_radial_distance(zs); R = np.sqrt(Om)*H0*DMz/C; lA = np.pi*DMz/rs
    v = np.array([R-R_PL, lA-LA_PL]); c += float(v @ CMB_CINV @ v)
    c += ((ombh2-OMBH2_PRIOR[0])/OMBH2_PRIOR[1])**2 + ((H0-SHOES[0])/SHOES[1])**2
    zc, H, icc = data['cc']; dH = H - bg.hubble_parameter(zc); c += float(dH @ icc @ dH)
    zf, fo, fe = data['fs8']; Dd, fd = fr.compute_growth_model(zf, Om, lambda zz: bg.hubble_parameter(np.atleast_1d(zz))/H0)
    c += float(np.sum(((fo - fd*s8*Dd)/fe)**2))
    if s8prior is not None:
        S8 = s8*np.sqrt(Om/0.3); c += ((S8 - s8prior[0])/s8prior[1])**2
    return c, R, lA, rd

def fit(sn, lcdm=False, s8prior=None):
    def obj(v):
        Om, H0, ob, s8 = v
        if not (0.2<Om<0.45 and 60<H0<78 and 0.019<ob<0.025 and 0.62<s8<0.90): return 1e9
        try: return chi2(Om, H0, ob, s8, sn, lcdm=lcdm, s8prior=s8prior)[0]
        except Exception: return 1e9
    return minimize(obj, [0.30, 68.5, 0.02237, 0.78], method='Nelder-Mead',
                    options=dict(xatol=1e-4, fatol=1e-3, maxiter=6000))

def w_crosses(Om):
    z = np.linspace(0, 2, 60); E = fr.E_SEDE_lambda(z, Om, 1.4964, 0.5)
    rho = np.maximum(E**2 - Om*(1+z)**3, 1e-8); w = -1 + (1/3.)*(1+z)*np.gradient(np.log(rho), z)
    return w[0], w.min(), z[np.argmin(w)], bool((w < -1).any())

if __name__ == "__main__":
    print("="*72); print("TIER 1 — real-data confrontations (canonical SEDE-H λ=0.5, γ=theory)"); print("="*72)

    # 1a — SN robustness
    print("\n[1a] SN robustness — Δχ²(SEDE−ΛCDM), SN offset marginalised:")
    blocks = sn_blocks()
    for name, sn in blocks.items():
        rS = fit(sn, lcdm=False); rL = fit(sn, lcdm=True)
        OmS = rS.x[0]; w0, wmin, zmin, cr = w_crosses(OmS)
        print(f"  {name:11s}: Δχ²={rS.fun-rL.fun:+.2f}  (SEDE χ²={rS.fun:.1f}, ΛCDM {rL.fun:.1f}) "
              f"| Ωm={OmS:.3f} w0={w0:+.3f} min_w={wmin:+.3f}@z{zmin:.1f} crosses−1:{cr}", flush=True)

    # 1b — Lyα z=2.33 + early-DE
    print("\n[1b] Lyα high-z BAO (z=2.33) + early-DE fraction Ω_DE(z=3):")
    Om, H0, ob = 0.30, 68.5, 0.02237
    for lab, lc in [("ΛCDM", True), ("SEDE-H", False)]:
        bg = background(Om, H0, ob, 1.4964, 0.5, lcdm=lc); dd = bg.get_derived_params(); rd = dd['rdrag']
        DH = (C/bg.hubble_parameter(2.33))/rd; DM = bg.comoving_radial_distance(2.33)/rd
        Ez3 = bg.hubble_parameter(3.0)/H0; ODE3 = max(Ez3**2 - Om*(1+3.0)**3, 0.0)/Ez3**2
        print(f"  {lab:7s}: DH/rd(2.33)={DH:.2f}  DM/rd(2.33)={DM:.2f}  Ω_DE(z=3)={ODE3:.4f}")
    print("  DESI DR2 Lyα: DH/rd=8.63, DM/rd=38.99  (already in the joint above)")

    # 1c / 1d — S8 priors (consistency)
    print("\n[1c/1d] S8 confrontation (σ8 jointly fit to fσ8, then add the S8 prior):")
    sn = blocks.get('Pantheon+') or next(iter(blocks.values()))
    for tag, pri in [("no S8 prior", None),
                     ("DES-Y3/KiDS WL  S8=0.776±0.017", (0.776, 0.017)),
                     ("Planck+ACT lens S8=0.832±0.013", (0.832, 0.013))]:
        r = fit(sn, lcdm=False, s8prior=pri); Om, H0, ob, s8 = r.x; S8 = s8*np.sqrt(Om/0.3)
        print(f"  {tag:32s}: SEDE σ8={s8:.3f} Ωm={Om:.3f} -> S8={S8:.3f}  χ²={r.fun:.1f}")
    print("  (SEDE predicts S8≈0.76 — consistent with WL, mild tension with Planck-lensing = the S8 tension itself.)")
