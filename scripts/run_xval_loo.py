#!/usr/bin/env python3
"""
#9 LEAVE-ONE-PROBE-OUT cross-validation (ML-style robustness, rarely done for a DE model).

Re-fit the canonical parameter-free SEDE-H (λ=0.5, γ=theory) AND ΛCDM with each probe
dropped in turn, in the CAMB-pinned joint. If the ΔDIC≈−2.9 preference is real (not
driven by one dataset), Δχ²(SEDE−ΛCDM) should stay negative and the parameters stable
across every holdout.  Probes: DESI BAO, SN (Pantheon+), CMB (R,l_A), CC, fσ8.

Run:  python run_xval_loo.py     (~10 min; 6 holdouts × 2 models, CAMB in the loop)
"""
import numpy as np, camb
from camb import dark_energy
from scipy.optimize import minimize
from scipy.linalg import cho_solve
from sede import friedmann as fr
from run_lambda_verify import w_of_a, data, R_PL, LA_PL, CMB_CINV, OMBH2_PRIOR, SHOES
C = 299792.458
PROBES = ["DESI", "SN", "CMB", "CC", "fs8"]

def chi2(Om, H0, ombh2, MB, s8, drop=set(), lcdm=False):
    h = H0/100.; pa = camb.CAMBparams(); pa.set_cosmology(H0=H0, ombh2=ombh2, omch2=Om*h**2-ombh2, mnu=0.06)
    if not lcdm:
        a, w = w_of_a(Om, 1.4964, 0.5); de = dark_energy.DarkEnergyPPF(); de.set_w_a_table(a, w); pa.DarkEnergy = de
    pa.WantCls = False; pa.WantTransfer = False; bg = camb.get_background(pa)
    dd = bg.get_derived_params(); rd = dd['rdrag']; zs = dd['zstar']; rs = dd['rstar']; c = 0.0
    if "DESI" not in drop:
        z, t, m, icov = data['desi']
        pred = np.array([bg.comoving_radial_distance(zz)/rd if tp=='DM/rd' else (C/bg.hubble_parameter(zz))/rd if tp=='DH/rd'
            else ((zz*bg.comoving_radial_distance(zz)**2*(C/bg.hubble_parameter(zz)))**(1/3.))/rd for zz, tp in zip(z, t)])
        c += float((m-pred) @ icov @ (m-pred))
    if "SN" not in drop and data['pan'] is not None:
        zp, mu, chol = data['pan']; dmu = mu - (5*np.log10((1+zp)*bg.comoving_radial_distance(zp))+25+MB)
        c += float(dmu @ cho_solve(chol, dmu))
    if "CMB" not in drop:
        DMz = bg.comoving_radial_distance(zs); R = np.sqrt(Om)*H0*DMz/C; lA = np.pi*DMz/rs
        v = np.array([R-R_PL, lA-LA_PL]); c += float(v @ CMB_CINV @ v)
    c += ((ombh2-OMBH2_PRIOR[0])/OMBH2_PRIOR[1])**2 + ((H0-SHOES[0])/SHOES[1])**2   # priors always on
    if "CC" not in drop:
        zc, H, icc = data['cc']; dH = H - bg.hubble_parameter(zc); c += float(dH @ icc @ dH)
    if "fs8" not in drop:
        zf, fo, fe = data['fs8']; Dd, fd = fr.compute_growth_model(zf, Om, lambda zz: bg.hubble_parameter(np.atleast_1d(zz))/H0)
        c += float(np.sum(((fo - fd*s8*Dd)/fe)**2))
    return c

def fit(drop=set(), lcdm=False):
    def obj(v):
        Om, H0, ob, MB, s8 = v
        if not (0.2<Om<0.45 and 60<H0<78 and 0.019<ob<0.025 and -20.5<MB<-18.5 and 0.62<s8<0.90): return 1e9
        try: return chi2(Om, H0, ob, MB, s8, drop, lcdm)
        except Exception: return 1e9
    return minimize(obj, [0.30, 68.5, 0.02237, -19.40, 0.78], method='Nelder-Mead',
                    options=dict(xatol=1e-4, fatol=1e-3, maxiter=6000))

if __name__ == "__main__":
    print("="*72); print("#9 LEAVE-ONE-PROBE-OUT — is the SEDE preference robust to dropping any probe?"); print("="*72)
    print(f"\n  {'holdout':>12s} {'Δχ²(SEDE−ΛCDM)':>16s} {'Ωm':>7s} {'H0':>7s} {'σ8':>7s}  (SEDE best-fit)")
    rows = []
    for drop in [set()] + [{p} for p in PROBES]:
        rS = fit(drop, lcdm=False); rL = fit(drop, lcdm=True)
        Om, H0, ob, MB, s8 = rS.x; d = rS.fun - rL.fun
        tag = "none (all)" if not drop else "− " + list(drop)[0]
        rows.append(d)
        print(f"  {tag:>12s} {d:>16.2f} {Om:>7.3f} {H0:>7.2f} {s8:>7.3f}", flush=True)
    rows = np.array(rows)
    print(f"\n  Δχ² range over all holdouts: [{rows.min():+.2f}, {rows.max():+.2f}]  mean={rows.mean():+.2f}")
    allneg = bool(np.all(rows < 0))
    print(f"  SEDE preferred (Δχ²<0) in EVERY holdout: {allneg}  ->",
          "robust, not driven by one dataset." if allneg else "NOT robust — one probe drives it.")
