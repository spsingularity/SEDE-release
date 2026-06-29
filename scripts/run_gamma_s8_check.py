#!/usr/bin/env python3
"""(b) Re-check the sibling team's two pushbacks under OUR CAMB-pinned full joint
(DESI+Pantheon+ + CMB(R,l_A) + CC + fσ8 + ombh2 prior + SH0ES), with σ8 fit JOINTLY
to the 16-pt fσ8 RSD data (not hand-set):
  1. γ_data for SEDE-H λ=0.5  (they claim ~1.7 matching p=5/3; we got ~0.74-0.78)
  2. S8 easing: S8(SEDE-H) vs S8(ΛCDM) with σ8 jointly fit
     (they claim no easing -> both ~0.78; our memory says SEDE eases to 0.74)
S8 = σ8 * sqrt(Ωm/0.3).
"""
import numpy as np, camb
from camb import dark_energy
from scipy.optimize import minimize
from scipy.linalg import cho_solve
from sede import friedmann as fr
from run_camb_joint import load, R_PL, LA_PL, CMB_CINV, OMBH2_PRIOR, SHOES
from run_lambda_verify import chi2, fit, data
C=299792.458

def S8(Om, s8): return s8*np.sqrt(Om/0.3)

# ---- LCDM full-joint fit, returning the parameter vector ----
def lcdm_fit():
    def obj(v):
        Om,H0,ob,MB,s8=v
        if not(0.2<Om<0.45 and 60<H0<78 and 0.019<ob<0.025 and -20.5<MB<-18.5 and 0.62<s8<0.90): return 1e9
        h=H0/100.; pa=camb.CAMBparams(); pa.set_cosmology(H0=H0,ombh2=ob,omch2=Om*h**2-ob,mnu=0.06)
        pa.WantCls=False; pa.WantTransfer=False
        bg=camb.get_background(pa); d=bg.get_derived_params(); rd=d['rdrag']; zs=d['zstar']; rs=d['rstar']; c=0.0
        z,t,m,icov=data['desi']
        pred=np.array([bg.comoving_radial_distance(zz)/rd if tp=='DM/rd' else (C/bg.hubble_parameter(zz))/rd if tp=='DH/rd'
            else ((zz*bg.comoving_radial_distance(zz)**2*(C/bg.hubble_parameter(zz)))**(1/3.))/rd for zz,tp in zip(z,t)])
        dd=m-pred; c+=float(dd@icov@dd)
        zp,mu,chol=data['pan']; dmu=mu-(5*np.log10((1+zp)*bg.comoving_radial_distance(zp))+25+MB); c+=float(dmu@cho_solve(chol,dmu))
        DMz=bg.comoving_radial_distance(zs); R=np.sqrt(Om)*H0*DMz/C; lA=np.pi*DMz/rs; v=np.array([R-R_PL,lA-LA_PL]); c+=float(v@CMB_CINV@v)
        c+=((ob-OMBH2_PRIOR[0])/OMBH2_PRIOR[1])**2 + ((H0-SHOES[0])/SHOES[1])**2
        zc,H,icc=data['cc']; dH=H-bg.hubble_parameter(zc); c+=float(dH@icc@dH)
        zf,fo,fe=data['fs8']; Dd,fd=fr.compute_growth_model(zf,Om,lambda zz:bg.hubble_parameter(np.atleast_1d(zz))/H0); c+=float(np.sum(((fo-fd*s8*Dd)/fe)**2))
        return c
    return minimize(obj,[0.30,68.5,0.02237,-19.40,0.78],method='Nelder-Mead',options=dict(xatol=1e-4,fatol=1e-3,maxiter=6000))

if __name__=="__main__":
    print("Fitting LCDM (full joint, sigma8 jointly fit to fσ8)...", flush=True)
    rL=lcdm_fit(); Om,H0,ob,MB,s8=rL.x
    print(f"  ΛCDM:  Ωm={Om:.3f} H0={H0:.2f} σ8={s8:.3f} -> S8={S8(Om,s8):.3f}  χ²={rL.fun:.2f}", flush=True)

    print("Fitting SEDE-H λ=0.5, γ=theory(1.4964) [parameter-free]...", flush=True)
    rP=fit(1.4964,0.5,gfree=False,use_shoes=True); Om,H0,ob,MB,s8=rP.x
    print(f"  SEDE-H γ=th: Ωm={Om:.3f} H0={H0:.2f} σ8={s8:.3f} -> S8={S8(Om,s8):.3f}  Δχ²={rP.fun-rL.fun:+.2f}", flush=True)

    print("Fitting SEDE-H λ=0.5, γ FREE [γ_data]...", flush=True)
    rF=fit(None,0.5,gfree=True,use_shoes=True); Om,H0,ob,MB,s8,g=rF.x
    print(f"  SEDE-H γ-free: γ_data={g:.2f}  Ωm={Om:.3f} H0={H0:.2f} σ8={s8:.3f} -> S8={S8(Om,s8):.3f}  Δχ²={rF.fun-rL.fun:+.2f}", flush=True)

    print("\n=== VERDICT ===")
    print(f"  γ_data (ours, CAMB-pinned, λ=0.5) = {g:.2f}   [their claim: ~1.7 (p=5/3); their cai_kim γ-free=0.84]")
    print(f"  S8: ΛCDM={S8(*[rL.x[0],rL.x[4]]):.3f}  vs  SEDE-H(γ-free)={S8(rF.x[0],rF.x[4]):.3f}  (σ8 jointly fit)")
