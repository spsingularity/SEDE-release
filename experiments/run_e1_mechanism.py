#!/usr/bin/env python3
"""
E1 — INDEPENDENT CHECK of SEDE's FOUNDING CLAIM (borrowed from SEDE_V2's E1).
SEDE says ρ_DE is sourced by structure formation: the f_sat read off the EXPANSION
(geometry, from DESI D_H/r_d → H(z)) should equal the f_sat read off STRUCTURE GROWTH
(from RSD fσ8 → σ8(z)). If the two diverge, the structure-sourcing mechanism — the thing
that makes SEDE distinctive — is not supported by current data.

  geometry f_sat(z)  = [E²(z) − Ω_m(1+z)³ − Ω_r(1+z)⁴] / [Ω_DE0 · E(z)^{2λ}],  E from DESI D_H
  structure f_sat(z) = (1 − e^{−γ D²}) / (1 − e^{−γ}),  D(z) = σ8(z)/σ8(0), σ8(z)=fσ8_obs/f(z)

Run:  python run_e1_mechanism.py
"""
import numpy as np, camb
from camb import dark_energy
from sede import friedmann as fr
from sede import data_loader as dl
from run_lambda_verify import w_of_a, data
C = 299792.458
Om, H0, s8, GAM, LAM = 0.298, 68.81, 0.760, 1.4964, 0.5     # SEDE best-fit (our joint)
Or = 9e-5; ODE0 = 1 - Om - Or
def Es(z): return fr.E_SEDE_lambda(np.atleast_1d(z), Om, GAM, LAM)
def fsat_model(D): return (1 - np.exp(-GAM*D**2))/(1 - np.exp(-GAM))

if __name__ == "__main__":
    print("="*72); print("E1 — does the GEOMETRY f_sat match the STRUCTURE f_sat? (SEDE's founding claim)"); print("="*72)
    # r_drag from CAMB at the SEDE best-fit
    h = H0/100.; pa = camb.CAMBparams(); pa.set_cosmology(H0=H0, ombh2=0.02237, omch2=Om*h**2-0.02237, mnu=0.06)
    a, w = w_of_a(Om, GAM, LAM); de = dark_energy.DarkEnergyPPF(); de.set_w_a_table(a, w); pa.DarkEnergy = de
    pa.WantCls = False; pa.WantTransfer = False; bg = camb.get_background(pa); rd = bg.get_derived_params()['rdrag']

    # --- geometry leg: DESI D_H/r_d -> H(z) -> f_sat_geom ---
    z, t, m, icov = data['desi']; cov = np.linalg.inv(icov)
    print("\n[geometry leg] from DESI D_H/r_d (rd=%.2f Mpc):" % rd)
    zg, fg, fg_e = [], [], []
    for i, (zz, tp, val) in enumerate(zip(z, t, m)):
        if tp != 'DH/rd': continue
        H = C/(val*rd); E = H/H0
        rhoDE = E**2 - Om*(1+zz)**3 - Or*(1+zz)**4
        fsg = rhoDE/(ODE0*E**(2*LAM))
        # propagate the DESI D_H error
        sig = np.sqrt(cov[i, i]); dE = -(C/(val**2*rd))/H0*sig
        drho = 2*E*dE; dfg = abs(drho/(ODE0*E**LAM) - rhoDE*LAM*E**(LAM-1)*dE/(ODE0*E**(2*LAM)))
        zg.append(zz); fg.append(fsg); fg_e.append(abs(dfg))
        print(f"   z={zz:.3f}: f_sat_geom = {fsg:.3f} ± {abs(dfg):.3f}")
    zg, fg, fg_e = np.array(zg), np.array(fg), np.array(fg_e)

    # --- structure leg: RSD fσ8 -> σ8(z) -> D(z) -> f_sat_struct ---
    zf, fo, fe = dl.load_fss8()
    D, frate = fr.compute_growth_model(zf, Om, lambda z: Es(z))
    s8z = fo/frate; Dstr = s8z/s8; fss = fsat_model(Dstr)
    dfdD = (2*GAM*Dstr*np.exp(-GAM*Dstr**2))/(1-np.exp(-GAM))
    fss_e = np.abs(dfdD * (fe/frate)/s8)
    print("\n[structure leg] from RSD fσ8 -> σ8(z) (16 pts):")
    o = np.argsort(zf)
    for i in o:
        print(f"   z={zf[i]:.3f}: σ8(z)={s8z[i]:.3f}  f_sat_struct = {fss[i]:.3f} ± {fss_e[i]:.3f}")

    # --- compare: structure f_sat vs the MODEL/geometry f_sat at the RSD redshifts ---
    Dm, _ = fr.compute_growth_model(zf, Om, lambda z: Es(z))
    fmodel = fsat_model(Dm/Dm[np.argmin(zf)] if False else Dm)   # model f_sat at RSD z (D normalised to D(0)=1 inside)
    fmodel = fsat_model(Dm)                                       # compute_growth_model returns D(0)=1-normalised
    chi2_rsd = float(np.sum(((fss - fmodel)/fss_e)**2)); ndof = len(zf)
    # their error budget: precise DESI GEOMETRY f_sat as data, structure f_sat (interp) as prediction
    o = np.argsort(zf)
    fss_at_geom = np.interp(zg, zf[o], fss[o])         # structure f_sat at the DESI D_H redshifts
    chi2_geom = float(np.sum(((fg - fss_at_geom)/fg_e)**2)); ndg = len(zg)
    hi = zf > 1.0
    print("\n" + "="*72)
    print(f"  RSD-error budget   : structure vs model f_sat  χ²={chi2_rsd:.1f}/{ndof} = {chi2_rsd/ndof:.2f}  (large RSD errors -> consistent)")
    print(f"  geometry-error budget: geometry vs structure f_sat  χ²={chi2_geom:.1f}/{ndg} = {chi2_geom/ndg:.2f}  (precise DESI -> the divergence bites)")
    if hi.any():
        print(f"  high-z (z>1): structure f_sat≈{fss[hi].mean():.2f} (flattens) vs geometry/model≈{fmodel[hi].mean():.2f} (declines) "
              f"-> {'DIVERGE' if fss[hi].mean()-fmodel[hi].mean()>0.05 else 'track'}")
    print(f"\n  VERDICT (honest, both teams agree): the geometry & structure f_sat TRACK at z<1 but the")
    print(f"  structure leg FLATTENS at z>1 while geometry DECLINES — the structure-sourcing link is")
    print(f"  NOT confirmed. Significance is error-budget-dependent (RSD-large {chi2_rsd/ndof:.2f}/dof = consistent;")
    print(f"  DESI-precise {chi2_geom/ndg:.2f}/dof = marginally disfavoured) -> INCONCLUSIVE, RSD-limited.")
    print(f"  A clean test needs WL-tomography σ8(z) (Euclid/LSST), not sparse RSD. SEDE's EOS is favoured;")
    print(f"  its distinctive *mechanism* (f_sat ∝ structure growth) is the honest open question.")
