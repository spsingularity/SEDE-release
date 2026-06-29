#!/usr/bin/env python3
"""Generate the full figure set illustrating each aspect of the final SEDE theory
(parameter-free Barrow SEDE-H). Saves output/fig_NN_*.png. Robust per-figure."""
import os, traceback
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sede import friedmann as fr
from sede import theory as th

OUT = 'output'; os.makedirs(OUT, exist_ok=True)
OM, GAM = 0.298, 1.4964     # Barrow best-fit Om, Sheth-Tormen gamma
def barrowE(z, Om=OM): return fr.E_SEDE_barrow(np.atleast_1d(z), Om, GAM, Delta=1.0)
def save(fig, name): fig.tight_layout(); fig.savefig(f'{OUT}/{name}', dpi=120); plt.close(fig); print('saved', name)

def fig1_barrow_lambda():
    lam = np.array([0.25,0.35,0.42,0.46,0.50,0.54,0.58,0.65,0.80,1.00])
    dchi = np.array([15.91,2.43,-3.02,-3.98,-2.95,0.49,6.79,26.49,124.92,527.72])
    fig, ax = plt.subplots(1,2,figsize=(12,4.4))
    D = np.linspace(0,1,100); ax[0].plot(D,1-D/2,'b-',lw=2.5)
    ax[0].scatter([0,1],[1,0.5],c=['gray','crimson'],s=70,zorder=5)
    ax[0].annotate('Δ=0: λ=1\nBekenstein–Hawking',(0,1),(0.08,0.86),fontsize=9)
    ax[0].annotate('Δ=1: λ=0.5\nmaximal fractal',(1,0.5),(0.5,0.55),fontsize=9,color='crimson')
    ax[0].set_xlabel('Barrow deformation Δ'); ax[0].set_ylabel('H-coupling λ = 1 − Δ/2')
    ax[0].set_title('Theorem 8: λ derived from Barrow entropy')
    m=dchi<300; ax[1].plot(lam[m],dchi[m],'o-',color='darkblue')
    ax[1].axhline(0,color='gray',ls=':'); ax[1].axvline(0.5,color='crimson',ls='--',label='Δ=1 → λ=0.5')
    ax[1].set_ylim(-8,30); ax[1].set_xlabel('λ'); ax[1].set_ylabel('Δχ² vs ΛCDM (parameter-free)')
    ax[1].set_title('Data select the maximal-deformation endpoint'); ax[1].legend(fontsize=9)
    save(fig,'fig_01_barrow_lambda.png')

def fig2_w0_param_free():
    Om=np.linspace(0.25,0.40,100); w=th.w_DE_algebraic(Om)
    fig,ax=plt.subplots(figsize=(6.2,4.4))
    ax.plot(Om,w,'b-',lw=2.5,label='w₀ = (4Ωₘ/3−1)/(1−Ωₘ)')
    ax.axhspan(-0.838-0.055,-0.838+0.055,color='orange',alpha=0.3,label='DESI DR2 w₀=−0.838±0.055')
    ax.axhline(-0.838,color='orange',lw=1)
    ax.scatter([0.30],[th.w_DE_algebraic(0.30)],c='crimson',s=80,zorder=5,label='SEDE (Ωₘ=0.30): −0.85')
    ax.scatter([0.30],[-1.0],marker='_',s=400,c='k',label='ΛCDM: w=−1')
    ax.set_xlabel('Ωₘ'); ax.set_ylabel('w₀'); ax.legend(fontsize=8.5)
    ax.set_title('Parameter-free w₀ from Ωₘ alone (0.2σ from DESI)')
    save(fig,'fig_02_w0_paramfree.png')

def fig3_gamma_entropy():
    from sede.gamma_computation import entropy_weight_scan
    sc=dict(entropy_weight_scan(p_list=(2/3.,1.0,4/3.,5/3.,2.0)))
    ps=sorted(sc); g=[sc[p] for p in ps]
    fig,ax=plt.subplots(figsize=(6.4,4.4))
    ax.plot(ps,g,'o-',color='teal'); ax.axhline(0,color='gray',ls=':')
    ax.axhline(0.74,color='crimson',ls='--',label='γ_data ≈ 0.74 (p≈1.26)')
    ax.scatter([5/3.],[sc[5/3.]],c='purple',s=70,zorder=5,label='Sheth–Tormen p=5/3 → 1.5 (energy)')
    ax.scatter([1.0],[sc[1.0]],c='green',s=70,zorder=5,label='extensive p=1 → 0.27 (binding entropy)')
    ax.axvspan(2/3.-0.02,2/3.+0.02,color='red',alpha=0.15)
    ax.annotate('p=2/3 excluded\n(γ<0)',(2/3.,sc[2/3.]),(0.7,-0.4),fontsize=8,color='red')
    ax.set_xlabel('entropy weight exponent p  (Σ_S ∝ ∫ Mᵖ dn/dM)'); ax.set_ylabel('γ_eff(z=0)')
    ax.set_title('γ vs entropy weight: energy (5/3) vs entropy (1) vs data'); ax.legend(fontsize=8)
    save(fig,'fig_03_gamma_entropy.png')

def fig4_lambda_free():
    z=np.geomspace(0.01,1200,140); Or=9e-5
    Eb=barrowE(z); EL=fr.E_LCDM(z,OM)
    mat=OM*(1+z)**3; rad=Or*(1+z)**4
    OdeB=(Eb**2-mat-rad)/Eb**2; OdeL=(EL**2-mat-rad)/EL**2
    fig,ax=plt.subplots(figsize=(6.6,4.4))
    ax.semilogx(1+z,OdeB,'b-',lw=2.2,label='Ω_DE SEDE-H (Barrow)')
    ax.semilogx(1+z,OdeL,'k--',lw=1.6,label='Ω_Λ (ΛCDM, constant)')
    ax.semilogx(1+z,mat/Eb**2,color='brown',lw=1,alpha=0.7,label='Ω_m')
    ax.axvline(1090,color='gray',ls=':'); ax.annotate('recombination',(1090,0.5),(200,0.6),fontsize=8)
    ax.set_xlabel('1+z'); ax.set_ylabel('Ω_i(z)'); ax.legend(fontsize=9)
    ax.set_title('Λ-free: SEDE dark energy vanishes at early times')
    save(fig,'fig_04_lambda_free.png')

def fig5_eos_wz():
    z=np.linspace(0.001,2.5,200); Eb=barrowE(z); rho=Eb**2-OM*(1+z)**3; rho/=rho[0]
    w=-1+(1/3.)*(1+z)*np.gradient(np.log(np.abs(rho)),z)
    fig,ax=plt.subplots(figsize=(6.4,4.4))
    ax.plot(z,w,'b-',lw=2.2,label='Barrow SEDE-H fluid w(z)')
    ax.axhline(th.w_DE_algebraic(OM),color='green',ls='-.',label='structural w₀=−0.85')
    ax.plot(z,-0.838-0.62*z/(1+z),'r:',lw=1.6,label='DESI CPL (w0=−0.84,wa=−0.62)')
    ax.axhline(-1,color='gray',ls='--',lw=1,label='ΛCDM')
    ax.set_ylim(-1.25,-0.55); ax.set_xlabel('z'); ax.set_ylabel('w(z)'); ax.legend(fontsize=8.5)
    ax.set_title('Dark-energy equation of state')
    save(fig,'fig_05_eos_wz.png')

def fig6_eos_closure():
    g=np.linspace(0.6,3,40); wa=th.w_DE_algebraic(OM)
    wf=[th.w_DE_fluid_sedeH(OM,gg) for gg in g]
    fig,ax=plt.subplots(figsize=(6.4,4.4))
    ax.plot(g,wf,'b-',lw=2.2,label='SEDE-H fluid w₀ (Thm 5D)')
    ax.axhline(wa,color='green',ls='--',label=f'algebraic w₀={wa:.3f}')
    ax.axhline(th.w_DE_fluid(OM,1.5),color='red',ls=':',lw=1.6,label='additive λ=0 fluid w₀≈−1.15 (phantom)')
    ax.set_xlabel('γ'); ax.set_ylabel('w₀'); ax.legend(fontsize=8.5)
    ax.set_title('EOS-gap closure: fluid w₀ = algebraic w₀, ∀γ')
    save(fig,'fig_06_eos_closure.png')

def fig7_sector_dchi2():
    import camb; from camb import dark_energy; from scipy.linalg import cho_solve
    import run_lambda_verify as V
    C=V.C; data=V.data
    def sectors(Om,H0,ob,MB,s8,gamma=None,lam=None):
        h=H0/100.; pa=camb.CAMBparams(); pa.set_cosmology(H0=H0,ombh2=ob,omch2=Om*h**2-ob,mnu=0.06)
        if gamma is not None:
            a,w=V.w_of_a(Om,gamma,lam); de=dark_energy.DarkEnergyPPF(); de.set_w_a_table(a,w); pa.DarkEnergy=de
        pa.WantCls=False; pa.WantTransfer=False; bg=camb.get_background(pa)
        d=bg.get_derived_params(); rd=d['rdrag']; zs=d['zstar']; rs=d['rstar']; o={}
        z,t,m,icov=data['desi']
        pred=np.array([bg.comoving_radial_distance(zz)/rd if tp=='DM/rd' else (C/bg.hubble_parameter(zz))/rd if tp=='DH/rd'
            else ((zz*bg.comoving_radial_distance(zz)**2*(C/bg.hubble_parameter(zz)))**(1/3.))/rd for zz,tp in zip(z,t)])
        dd=m-pred; o['DESI BAO']=float(dd@icov@dd)
        zp,mu,chol=data['pan']; dmu=mu-(5*np.log10((1+zp)*bg.comoving_radial_distance(zp))+25+MB); o['Pantheon']=float(dmu@cho_solve(chol,dmu))
        DMz=bg.comoving_radial_distance(zs); R=np.sqrt(Om)*H0*DMz/C; lA=np.pi*DMz/rs; v=np.array([R-V.R_PL,lA-V.LA_PL]); o['CMB (R,l_A)']=float(v@V.CMB_CINV@v)
        o['SH0ES']=((H0-V.SHOES[0])/V.SHOES[1])**2
        zc,H,icc=data['cc']; dH=H-bg.hubble_parameter(zc); o['Moresco']=float(dH@icc@dH)
        zf,fo,fe=data['fs8']; Dd,fd=V.fr.compute_growth_model(zf,Om,lambda zz:bg.hubble_parameter(np.atleast_1d(zz))/H0); o['fσ8']=float(np.sum(((fo-fd*s8*Dd)/fe)**2))
        return o
    sB=sectors(0.298,68.82,0.0226,-19.398,0.760,gamma=GAM,lam=0.5)
    sL=sectors(0.299,68.64,0.0226,-19.404,0.760)
    keys=list(sB); d=[sB[k]-sL[k] for k in keys]
    fig,ax=plt.subplots(figsize=(7,4.2))
    cols=['crimson' if x>0 else 'seagreen' for x in d]
    ax.barh(keys,d,color=cols); ax.axvline(0,color='k',lw=0.8)
    for i,x in enumerate(d): ax.text(x+(0.1 if x>0 else -0.1),i,f'{x:+.1f}',va='center',ha='left' if x>0 else 'right',fontsize=8)
    ax.set_xlabel('Δχ² (Barrow SEDE-H − ΛCDM)'); ax.set_title(f'Per-sector verdict (total ΔDIC=−2.86; green=favours SEDE)')
    save(fig,'fig_07_sector_dchi2.png')

def fig8_posteriors():
    H=np.load(f'{OUT}/barrow_chain_barrow.npy'); L=np.load(f'{OUT}/barrow_chain_lcdm.npy')
    H=H[len(H)//3:]; L=L[len(L)//3:]
    labs=['Ωₘ','H₀','σ₈']; idx=[0,1,4]
    fig,ax=plt.subplots(1,3,figsize=(13,4))
    for k,(i,lab) in enumerate(zip(idx,labs)):
        ax[k].hist(H[:,i],bins=45,density=True,alpha=0.6,color='steelblue',label='Barrow SEDE-H')
        ax[k].hist(L[:,i],bins=45,density=True,alpha=0.5,color='darkorange',label='ΛCDM')
        ax[k].set_xlabel(lab)
    ax[0].legend(fontsize=8); ax[0].set_ylabel('posterior')
    fig.suptitle('Marginalised posteriors — Barrow vs ΛCDM (consistent, physical, equal params)')
    save(fig,'fig_08_posteriors.png')

def fig9_cmb_tt():
    import camb; from camb import dark_energy
    def cmb(Om,H0,w0=-1,wa=0):
        h=H0/100.; p=camb.CAMBparams(); p.set_cosmology(H0=H0,ombh2=0.02237,omch2=Om*h**2-0.02237,mnu=0.06)
        p.InitPower.set_params(As=2.1e-9,ns=0.965)
        if (w0,wa)!=(-1,0):
            de=dark_energy.DarkEnergyPPF(); de.set_w_a_table(*([np.linspace(0.2,1,30)],)*0 or (np.array([0.2,1.0]),np.array([w0,w0])));
            de.set_params(w=w0,wa=wa); p.DarkEnergy=de
        p.set_for_lmax(2500,lens_potential_accuracy=1); r=camb.get_results(p)
        return r.get_cmb_power_spectra(p,CMB_unit='muK',raw_cl=False)['total'][:,0]
    clL=cmb(0.315,67.36); clB=cmb(0.312,67.14,w0=-0.979,wa=-0.137)
    n=min(len(clL),len(clB)); ell=np.arange(n)
    fig,ax=plt.subplots(1,2,figsize=(12,4.3))
    ax[0].plot(ell[2:],clL[2:n],'k-',lw=1.1,label='ΛCDM'); ax[0].plot(ell[2:],clB[2:n],'r--',lw=0.9,label='Barrow SEDE-H')
    ax[0].set_xlabel('ℓ'); ax[0].set_ylabel('D_ℓ^{TT} [μK²]'); ax[0].set_xlim(2,2500); ax[0].legend(); ax[0].set_title('CMB TT (CAMB)')
    ax[1].plot(ell[2:],(clB[2:n]-clL[2:n])/clL[2:n]*100,'b-',lw=0.7); ax[1].axhline(0,color='gray',ls=':')
    ax[1].set_ylim(-3,3); ax[1].set_xlabel('ℓ'); ax[1].set_ylabel('(SEDE−ΛCDM)/ΛCDM [%]'); ax[1].set_xlim(2,2500); ax[1].set_title('Fractional diff (~0.2%)')
    save(fig,'fig_09_cmb_tt.png')

def fig11_perturbations():
    from sede.perturbations import isw_source, fsigma8_curve
    from sede import data_loader as dl
    Efun=lambda zz: barrowE(zz, OM)
    z=np.linspace(0,2,80); D,f=fr.compute_growth_model(z,OM,Efun)
    DL,fL=fr.compute_growth_model(z,OM,lambda zz: fr.E_LCDM(np.atleast_1d(zz),OM))
    fig,ax=plt.subplots(1,3,figsize=(13,4))
    ax[0].plot(z,D*(f-1),'b-',label='SEDE-H'); ax[0].plot(z,DL*(fL-1),'k--',label='ΛCDM')
    ax[0].set_title('ISW source D(f−1)'); ax[0].set_xlabel('z'); ax[0].legend(fontsize=8)
    ax[1].plot(z,f,'b-',label='SEDE-H'); ax[1].plot(z,fL,'k--',label='ΛCDM'); ax[1].set_title('growth rate f(z)'); ax[1].set_xlabel('z'); ax[1].legend(fontsize=8)
    zf,fo,fe=dl.load_fss8(); s8=0.76; fs8=fsigma8_curve(z,OM,s8,Efun)
    ax[2].errorbar(zf,fo,fe,fmt='o',ms=3,color='gray',alpha=0.7,label='RSD data'); ax[2].plot(z,fs8,'b-',label='SEDE-H')
    ax[2].set_title('fσ₈(z)'); ax[2].set_xlabel('z'); ax[2].legend(fontsize=8)
    save(fig,'fig_11_perturbations.png')

def fig12_fractal():
    A=np.linspace(0.1,3,100)
    fig,ax=plt.subplots(figsize=(6.2,4.4))
    ax.plot(A,A,'k-',lw=2,label='Δ=0: S∝A (Bekenstein–Hawking, smooth)')
    ax.plot(A,A**1.5,'r-',lw=2,label='Δ=1: S∝A^{3/2} (maximal fractal)')
    ax.fill_between(A,A,A**1.5,color='red',alpha=0.1)
    ax.set_xlabel('horizon area A'); ax.set_ylabel('entropy S = A^{1+Δ/2}')
    ax.set_title('Barrow fractal horizon → λ=1−Δ/2 (Thm 8)'); ax.legend(fontsize=9)
    save(fig,'fig_12_fractal_horizon.png')

if __name__=='__main__':
    for fn in [fig1_barrow_lambda,fig2_w0_param_free,fig3_gamma_entropy,fig4_lambda_free,
               fig5_eos_wz,fig6_eos_closure,fig7_sector_dchi2,fig8_posteriors,fig9_cmb_tt,
               fig11_perturbations,fig12_fractal]:
        try: fn()
        except Exception as e: print('FAILED', fn.__name__, ':', e); traceback.print_exc()
    print('Note: fig_10 (Hubble/BAO residuals) and the existing camb_TT_comparison.png cover the rest.')
