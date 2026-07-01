"""
REALISTIC weak-lensing tomographic Fisher for σ(Δ)_growth — the "further step" beyond run_wl_fisher.py:
  + nonlinear P(k,z) (Halofit-Takahashi), ℓ_max=3000   [statistical GAIN]
  + full cosmology marginalised (Ω_m,σ8,Δ,h,n_s,Ω_b) with Planck priors on h,n_s,Ω_b   [degradation]
  + intrinsic alignments (NLA, A_IA) and baryon feedback (A_b) nuisances marginalised   [degradation]

Question: does the nonlinear GAIN survive realistic systematics, and what is the net σ(Δ)_growth (hence
the E1 structure-sourcing test significance)? SEDE's Δ enters only via growth D(z;Δ) and distances.
Output: console + results/wl_fisher_realistic.json.  (Stage-IV: f_sky=0.4, n_gal=30, 10 bins.)
"""
import json, os
import numpy as np
from scipy.special import erf
from sede import friedmann as fr
from sede import halofit as hf

C = 299792.458; GAM = 1.4964; H0H = 2997.92458
FSKY, NGAL, SIG_E, NBIN, Z0 = 0.4, 30.0, 0.26, 10, 0.70
LMIN, LMAX, NL = 20, 3000, 20
C1RHO = 0.0134                       # NLA normalisation C1·ρ_crit
kgrid = np.logspace(-4, 2, 2500)


def lam(D): return 1.0 - D / 2.0
def Efun(z, Om, D): return np.asarray(fr.E_SEDE_lambda(np.atleast_1d(z), Om, GAM, lam(D)))


def T_bbks(k, Om, h, Ob):
    G = Om * h * np.exp(-Ob - np.sqrt(2 * h) * Ob / Om); q = k / G
    return (np.log(1 + 2.34*q)/(2.34*q) * (1 + 3.89*q + (16.1*q)**2 + (5.46*q)**3 + (6.71*q)**4)**-0.25)


def Plin0(k, Om, s8, h, Ob, ns):
    Psh = k**ns * T_bbks(k, Om, h, Ob)
    x = k * 8.0; W = 3*(np.sin(x)-x*np.cos(x))/x**3
    s2 = np.trapezoid(k**2 * (k**ns*T_bbks(k,Om,h,Ob)**2) * W**2, k)/(2*np.pi**2)
    return (s8**2 / s2) * k**ns * T_bbks(k, Om, h, Ob)**2


def nz_bins(zg):
    nz = zg**2*np.exp(-(zg/Z0)**1.5)
    cdf = np.concatenate([[0], np.cumsum(0.5*(nz[1:]+nz[:-1])*np.diff(zg))]); cdf/=cdf[-1]
    edges = np.interp(np.linspace(0,1,NBIN+1), cdf, zg); sigz=0.05*(1+zg)
    out=[]
    for a,b in zip(edges[:-1],edges[1:]):
        w=0.5*(erf((b-zg)/(np.sqrt(2)*sigz))-erf((a-zg)/(np.sqrt(2)*sigz)))
        ni=nz*w; ni/=np.trapezoid(ni,zg); out.append(ni)
    return out


def cl_matrix(th, ells, zg, bins):
    Om, s8, D, h, ns, Ob, A_IA, A_b = th
    E = Efun(zg, Om, D)
    chi = H0H*np.concatenate([[0], np.cumsum(0.5*(1/E[1:]+1/E[:-1])*np.diff(zg))])
    Dg,_ = fr.compute_growth_model(zg, Om, lambda z: Efun(z, Om, D))
    dchidz = H0H/E
    OmZ = Om*(1+zg)**3 / E**2
    P0k = Plin0(kgrid, Om, s8, h, Ob, ns)
    bary = 1 - A_b*( (kgrid/5.0)**2/(1+(kgrid/5.0)**2) )      # baryon high-k suppression
    # lensing efficiency q_i and IA kernel wI_i
    q=[]; wI=[]
    for ni in bins:
        g=np.array([np.trapezoid((ni*(chi-chi[j])/np.where(chi>0,chi,1))[j:], zg[j:]) for j in range(len(zg))])
        q.append(1.5*Om/H0H**2*(1+zg)*np.clip(g,0,None)*chi)
        wI.append(-A_IA*C1RHO*Om/np.clip(Dg,1e-3,None) * ni / dchidz)   # NLA, per dχ
    q=np.array(q); wI=np.array(wI)
    nb=len(bins); Cl=np.zeros((len(ells),nb,nb)); m=chi>1e-3
    for li,l in enumerate(ells):
        k=(l+0.5)/chi[m]
        D2lin = k**3 * (np.interp(k,kgrid,P0k)*Dg[m]**2) /(2*np.pi**2)
        # nonlinear via halofit per-l is costly; approximate boost from z=0 table scaled by growth is
        # inaccurate, so apply halofit on the (k,z) along the ray using local Ω_m(z):
        Pnl = np.interp(k,kgrid,P0k*bary)*Dg[m]**2
        # nonlinear boost: evaluate halofit on a coarse representative — use linear*boost(k,zeff)
        D2 = k**3*Pnl/(2*np.pi**2)
        boost = hf.nonlinear_D2(k, D2, np.mean(OmZ[m]))/np.where(D2>0,D2,1)
        Pnl = Pnl*np.clip(boost,1,None)
        w = dchidz[m]/chi[m]**2 * Pnl
        for i in range(nb):
            for j in range(i,nb):
                GG=q[i,m]*q[j,m]; GI=q[i,m]*wI[j,m]+wI[i,m]*q[j,m]; II=wI[i,m]*wI[j,m]
                Cl[li,i,j]=Cl[li,j,i]=np.trapezoid(w*(GG+GI+II), zg[m])
    return Cl


def main():
    zg=np.linspace(1e-3,3.2,200); bins=nz_bins(zg)
    ells=np.unique(np.logspace(np.log10(LMIN),np.log10(LMAX),NL).astype(int)).astype(float)
    n_sr=(NGAL/NBIN)*(180*60/np.pi)**2; N=np.eye(NBIN)*SIG_E**2/n_sr
    fid=[0.30,0.80,1.0, 0.67,0.965,0.049, 1.0, 0.0]
    step=[0.004,0.01,0.03, 0.01,0.01,0.001, 0.1,0.1]
    names=['Om','s8','Delta','h','ns','Ob','A_IA','A_b']
    prior=dict(h=0.005, ns=0.004, Ob=0.0005, A_b=0.3)        # Planck-ish + baryon theory; A_IA free
    Cf=cl_matrix(fid,ells,zg,bins)
    dC=[]
    for i in range(len(fid)):
        tp=fid[:]; tm=fid[:]; tp[i]+=step[i]; tm[i]-=step[i]
        dC.append((cl_matrix(tp,ells,zg,bins)-cl_matrix(tm,ells,zg,bins))/(2*step[i]))
    dl=np.gradient(ells); n=len(fid); F=np.zeros((n,n))
    for li in range(len(ells)):
        Cinv=np.linalg.inv(Cf[li]+N); pref=(2*ells[li]+1)/2*FSKY*dl[li]
        for a in range(n):
            for b in range(n):
                F[a,b]+=pref*np.trace(Cinv@dC[a][li]@Cinv@dC[b][li])
    for nm,s in prior.items():
        F[names.index(nm),names.index(nm)]+=1/s**2
    sig=np.sqrt(np.diag(np.linalg.inv(F)))
    sD=sig[names.index('Delta')]
    print("="*76); print("REALISTIC WL FISHER (Halofit ℓ≤3000 + full cosmology + IA + baryons, marginalised)"); print("="*76)
    print("  marginalised 1σ: "+"  ".join(f"{nm}={sig[i]:.4f}" for i,nm in enumerate(names)))
    print(f"\n  σ(Δ)_growth: linear/3-param {0.241:.2f}  →  REALISTIC {sD:.2f}")
    import run_e1_forecast as E1
    e1=E1.main(); sG=e1['sigma_Delta_geom']; sE1=np.sqrt(sG**2+sD**2)
    print(f"  σ(Δ)_geom {sG:.2f} | σ(Δ_geom−Δ_growth) realistic = {sE1:.2f} ⟹ decoupling test {1/sE1:.1f}σ")
    out=dict(sigma=dict(zip(names,sig.tolist())), sigma_Delta_growth_realistic=float(sD),
             sigma_Delta_geom=sG, decoupling_sigma=float(1/sE1))
    os.makedirs("results",exist_ok=True); json.dump(out,open("results/wl_fisher_realistic.json","w"),indent=2)
    print("\nwrote results/wl_fisher_realistic.json")
    return out


if __name__=="__main__":
    main()
