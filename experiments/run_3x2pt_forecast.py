"""
JOINT 3×2pt + CMB-lensing Fisher with the FULL Stage-IV nuisance vector — the publication-grade version
of the E1 mechanism forecast (no deflation hack).

Fields: N_s shear + N_l clustering + CMB-lensing in one per-ℓ covariance matrix (all cross-spectra
κκ/δκ/δδ and φδ/φκ/φφ + covariance exact). Marginalised parameters:
  cosmology  : Ω_m, σ8, Δ, h, n_s, Ω_b            (Planck priors on h,n_s,Ω_b)
  astrophysics: galaxy bias b_i (per lens bin, free); intrinsic alignment A_IA, η_IA (NLA)
  systematics : photo-z shifts Δz_i (per source & lens bin, σ=0.002); shear mult. bias m_i (σ=0.005)
All three nuisances act on the field WEIGHTS: photo-z shifts n_i(z); IA adds an alignment kernel to the
shear weight (so GI/II + IA×{clustering,CMB-lens} all appear automatically); m_i scales the shear field.
SEDE's Δ enters only via growth D(z;Δ) and distances. Output: console + results/3x2pt_forecast.json.

Honest scope: linear P(k) to ℓ_max=1500; representative CMB-lensing reconstruction noise. This is the
full-nuisance forecast the earlier 'deflation' estimated.
"""
import json, os
import numpy as np
from scipy.special import erf
from sede import friedmann as fr

C = 299792.458; H0H = 2997.92458; GAM = 1.4964
ZSTAR = 1089.9; C1RHO = 0.0134; ZPIV = 0.62
NS, NL = 5, 5
FSKY, NGAL_S, NGAL_L, SIG_E = 0.4, 27.0, 18.0, 0.26
N_KCMB = 2.0e-8
Z0 = 0.70; LMIN, LMAX, NL_ELL = 20, 1500, 16
ARC2SR = (180 * 60 / np.pi)**2


def lam(D): return 1.0 - D / 2.0
def Efun(z, Om, D): return np.asarray(fr.E_SEDE_lambda(np.atleast_1d(z), Om, GAM, lam(D)))


def _T(k, Om, h, Ob):
    G = Om * h * np.exp(-Ob - np.sqrt(2 * h) * Ob / Om); q = k / G
    return (np.log(1 + 2.34*q)/(2.34*q) * (1 + 3.89*q + (16.1*q)**2 + (5.46*q)**3 + (6.71*q)**4)**-0.25)


def Plin(k, Om, s8, h, Ob, ns):
    x = k*8.0; Wt = 3*(np.sin(x)-x*np.cos(x))/x**3
    s2 = np.trapezoid(k**2*(k**ns*_T(k,Om,h,Ob)**2)*Wt**2, k)/(2*np.pi**2)
    return (s8**2/s2)*k**ns*_T(k, Om, h, Ob)**2


def base_bins(zg, nbin):
    nz = zg**2*np.exp(-(zg/Z0)**1.5)
    cdf = np.concatenate([[0], np.cumsum(0.5*(nz[1:]+nz[:-1])*np.diff(zg))]); cdf/=cdf[-1]
    e = np.interp(np.linspace(0,1,nbin+1), cdf, zg); sigz=0.05*(1+zg); out=[]
    for a,b in zip(e[:-1],e[1:]):
        w=0.5*(erf((b-zg)/(np.sqrt(2)*sigz))-erf((a-zg)/(np.sqrt(2)*sigz)))
        ni=nz*w; ni/=np.trapezoid(ni,zg); out.append(ni)
    return out


def shift(n, dz, zg):
    ns = np.interp(zg - dz, zg, n, left=0.0, right=0.0); s = np.trapezoid(ns, zg)
    return ns/s if s > 0 else n


# parameter layout
NCOS = 6        # Om s8 D h ns Ob
IDX = dict(Om=0, s8=1, D=2, h=3, ns=4, Ob=5, A_IA=6, eta=7)
B0 = 8; DZS0 = B0+NL; DZL0 = DZS0+NS; MS0 = DZL0+NL; NPAR = MS0+NS


def cl_all(th, ells, zg, bs, bl):
    Om,s8,D,h,ns,Ob,A_IA,eta = th[:8]
    bvec = th[B0:B0+NL]; dzs = th[DZS0:DZS0+NS]; dzl = th[DZL0:DZL0+NL]; mvec = th[MS0:MS0+NS]
    E = Efun(zg, Om, D)
    chi = H0H*np.concatenate([[0], np.cumsum(0.5*(1/E[1:]+1/E[:-1])*np.diff(zg))])
    Dg,_ = fr.compute_growth_model(zg, Om, lambda z: Efun(z, Om, D))
    cz = np.linspace(zg[-1], ZSTAR, 300); chistar = chi[-1]+H0H*np.trapezoid(1/Efun(cz,Om,D), cz)
    pref = 1.5*Om/H0H**2*(1+zg)*chi
    sb = [shift(bs[i], dzs[i], zg) for i in range(NS)]
    lb = [shift(bl[i], dzl[i], zg) for i in range(NL)]
    Wg = []
    for i in range(NS):
        g = np.array([np.trapezoid((sb[i]*(chi-chi[j])/np.where(chi>0,chi,1))[j:], zg[j:]) for j in range(len(zg))])
        Wk = pref*np.clip(g,0,None)                                            # lensing
        WI = -A_IA*((1+zg)/(1+ZPIV))**eta * C1RHO*Om/np.clip(Dg,1e-3,None) * sb[i]*(E/H0H)  # IA (NLA)
        Wg.append((1+mvec[i])*(Wk+WI))                                         # ×(1+m) shear calibration
    Wd = [bvec[i]*lb[i]*(E/H0H) for i in range(NL)]
    Wphi = pref*(chistar-chi)/chistar
    Wall = Wg + Wd + [Wphi]; nf=len(Wall); m=chi>1e-3
    kg=np.logspace(-4,2,1800); P0=Plin(kg,Om,s8,h,Ob,ns); meas=(H0H/E)[m]
    Cl=np.zeros((len(ells),nf,nf))
    for li,l in enumerate(ells):
        k=(l+0.5)/chi[m]; Pk=np.interp(k,kg,P0)*Dg[m]**2; w=meas/chi[m]**2*Pk
        Wm=[Wx[m] for Wx in Wall]
        for a in range(nf):
            for b in range(a,nf):
                Cl[li,a,b]=Cl[li,b,a]=np.trapezoid(w*Wm[a]*Wm[b], zg[m])
    return Cl


def noise_matrix():
    nf=NS+NL+1; N=np.zeros((nf,nf))
    for i in range(NS): N[i,i]=SIG_E**2/((NGAL_S/NS)*ARC2SR)
    for i in range(NL): N[NS+i,NS+i]=1.0/((NGAL_L/NL)*ARC2SR)
    N[-1,-1]=N_KCMB; return N


def main():
    zg=np.linspace(1e-3,3.2,170); bs=base_bins(zg,NS); bl=base_bins(zg,NL)
    ells=np.unique(np.logspace(np.log10(LMIN),np.log10(LMAX),NL_ELL).astype(int)).astype(float)
    fid=np.zeros(NPAR)
    fid[:8]=[0.30,0.80,1.0,0.67,0.965,0.049, 0.5, 0.0]
    fid[B0:B0+NL]=[1.2+0.4*i/NL for i in range(NL)]
    names=['Om','s8','Delta','h','ns','Ob','A_IA','eta_IA']+[f'b{i}' for i in range(NL)]+\
          [f'dzs{i}' for i in range(NS)]+[f'dzl{i}' for i in range(NL)]+[f'm{i}' for i in range(NS)]
    step=np.array([0.004,0.01,0.03,0.01,0.01,0.001,0.05,0.1]+[0.03]*NL+[0.002]*NS+[0.002]*NL+[0.004]*NS)
    prior={'h':0.005,'ns':0.004,'Ob':0.0005,'eta_IA':3.0}
    for i in range(NS): prior[f'dzs{i}']=0.002
    for i in range(NL): prior[f'dzl{i}']=0.002
    for i in range(NS): prior[f'm{i}']=0.005
    N=noise_matrix(); Cf=cl_all(fid,ells,zg,bs,bl)
    dC=[]
    for i in range(NPAR):
        tp=fid.copy(); tm=fid.copy(); tp[i]+=step[i]; tm[i]-=step[i]
        dC.append((cl_all(tp,ells,zg,bs,bl)-cl_all(tm,ells,zg,bs,bl))/(2*step[i]))
    dl=np.gradient(ells); F=np.zeros((NPAR,NPAR))
    for li in range(len(ells)):
        Cinv=np.linalg.inv(Cf[li]+N); pf=(2*ells[li]+1)/2*FSKY*dl[li]
        for a in range(NPAR):
            for b in range(a,NPAR):
                F[a,b]=F[b,a]=F[a,b]+pf*np.trace(Cinv@dC[a][li]@Cinv@dC[b][li])
    for nm,s in prior.items(): F[names.index(nm),names.index(nm)]+=1/s**2
    sig=np.sqrt(np.diag(np.linalg.inv(F))); sD=sig[IDX['D']]; sG=0.192
    sE1=np.sqrt(sG**2+sD**2)
    print("="*80); print("JOINT 3×2pt + CMB-LENSING, FULL NUISANCE VECTOR (photo-z + IA + shear bias)"); print("="*80)
    print(f"  {NPAR} params marginalised: cosmology(6)+IA(2)+bias({NL})+photo-z({NS+NL})+shear-m({NS})")
    print(f"  σ(Ω_m)={sig[0]:.4f}  σ(σ8)={sig[1]:.4f}  σ(A_IA)={sig[6]:.3f}  σ(Δ)_growth={sD:.3f}")
    print("-"*80)
    print(f"  σ(Δ)_growth: WL 0.25 | joint(no nuisance) 0.05 | JOINT+FULL NUISANCES {sD:.3f}")
    print(f"  σ(Δ)_geom (DESI DR3) {sG} ⟹ σ(Δ_geom−Δ_growth)={sE1:.3f} ⟹ mechanism test {1/sE1:.1f}σ")
    print("  (σ(σ8) is now realistic Stage-IV ⟹ this REPLACES the earlier deflation estimate.)")
    out=dict(npar=NPAR, sigma_Om=float(sig[0]), sigma_s8=float(sig[1]), sigma_A_IA=float(sig[6]),
             sigma_Delta_growth=float(sD), sigma_Delta_geom=sG, mechanism_sigma=float(1/sE1))
    os.makedirs("results",exist_ok=True); json.dump(out,open("results/3x2pt_forecast.json","w"),indent=2)
    print("\nwrote results/3x2pt_forecast.json"); return out


if __name__=="__main__":
    main()
