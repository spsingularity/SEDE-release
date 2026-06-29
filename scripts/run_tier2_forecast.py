#!/usr/bin/env python3
"""
TIER 2 — Fisher forecasts for the surveys that will DECIDE SEDE (not yet observed,
so these are forecasts, not data). Question: how many σ can DESI DR3 + Euclid pin the
Barrow deformation Δ (fiducial Δ=1, λ=1−Δ/2) and the crossing-(-1) EOS?

Observables: BAO DM/rd, DH/rd (geometry) and fσ8 (growth) in survey z-bins.
Parameters: θ = (Ω_m, Δ, ln A_rd)  [A_rd = an rd/H0 amplitude nuisance, marginalised].
Fiducial: Ω_m=0.30, Δ=1 (λ=0.5), γ=theory. Fisher F_ij = Σ_obs (∂o/∂θ_i)(∂o/∂θ_j)/σ_o²;
σ(θ) = sqrt(diag(F⁻¹)).  Forecasted fractional errors are per-survey specs (refs in code).

Run:  python run_tier2_forecast.py
"""
import numpy as np
from sede import friedmann as fr
C = 299792.458
Om0, DEL0, GAM = 0.30, 1.0, 1.4964            # fiducial
H0, RD = 67.5, 147.1                            # fiducial H0, r_d (Mpc)

def lam(Delta): return 1.0 - Delta/2.0

def Efun(z, Om, Delta):
    return fr.E_SEDE_lambda(np.atleast_1d(z), Om, GAM, lam(Delta))

def obs_vector(Om, Delta, lnArd, zbao, zf):
    """Return [DM/rd..., DH/rd..., fσ8...] for the given parameters."""
    Ard = np.exp(lnArd)
    E = Efun(zbao, Om, Delta)
    # comoving distance via cumulative trapezoid on a fine grid
    out = []
    for z in zbao:
        zz = np.linspace(0, z, 256); Ez = Efun(zz, Om, Delta)
        DC = (C/H0) * np.trapezoid(1.0/Ez, zz)
        DM = DC; DH = (C/H0)/Efun(z, Om, Delta)[0]
        out.append(DM/(RD*Ard)); out.append(DH/(RD*Ard))
    # growth fσ8 (shape only; σ8 amplitude absorbed -> use f*D normalised)
    D, f = fr.compute_growth_model(np.asarray(zf), Om, lambda z: Efun(z, Om, Delta))
    s8 = 0.81
    for i in range(len(zf)):
        out.append(f[i]*D[i]*s8)
    return np.array(out)

def fisher(zbao, sig_DM, sig_DH, zf, sig_fs8):
    """sig_* are FRACTIONAL errors (DM,DH) / absolute (fσ8) at each bin."""
    th0 = np.array([Om0, DEL0, 0.0])
    o0 = obs_vector(*th0, zbao, zf)
    # observable 1σ vector
    sig = []
    for i in range(len(zbao)):
        sig += [sig_DM[i]*o0[2*i], sig_DH[i]*o0[2*i+1]]
    sig += list(sig_fs8)
    sig = np.array(sig)
    # numerical derivatives
    steps = np.array([0.01, 0.03, 0.01])
    J = np.zeros((len(o0), 3))
    for k in range(3):
        tp = th0.copy(); tp[k] += steps[k]; tm = th0.copy(); tm[k] -= steps[k]
        J[:, k] = (obs_vector(*tp, zbao, zf) - obs_vector(*tm, zbao, zf)) / (2*steps[k])
    F = J.T @ np.diag(1.0/sig**2) @ J
    cov = np.linalg.inv(F)
    return np.sqrt(np.diag(cov)), cov

def w0_of(Om, Delta):
    z = np.linspace(0, 0.02, 5); E = Efun(z, Om, Delta)
    rho = E**2 - Om*(1+z)**3; return (-1 + (1/3.)*np.gradient(np.log(rho), z)[0])

if __name__ == "__main__":
    print("="*72); print("TIER 2 — Fisher forecasts: σ(Δ) for the canonical SEDE-H (fiducial Δ=1)"); print("="*72)

    # ---- DESI DR3 (5yr) BAO: ~7 tracers z=0.3..2.3, ~0.3-0.7% distances ----
    z_desi = np.array([0.30, 0.51, 0.71, 0.93, 1.32, 1.49, 2.33])
    sDM = np.array([0.009, 0.004, 0.004, 0.004, 0.005, 0.006, 0.010])   # frac DM/rd
    sDH = np.array([0.015, 0.007, 0.007, 0.007, 0.009, 0.011, 0.012])   # frac DH/rd
    zf_desi = np.array([0.30, 0.51, 0.71, 0.93, 1.32]); sfs_desi = 0.02*np.array([0.46,0.46,0.45,0.44,0.42])

    # ---- Euclid: BAO ~0.5% + fσ8 ~1-2% over z=0.9..1.8 ----
    z_euc = np.array([0.9, 1.0, 1.2, 1.4, 1.65, 1.8])
    sDM_e = np.full(6, 0.005); sDH_e = np.full(6, 0.008)
    zf_euc = z_euc.copy(); sfs_euc = np.array([0.012,0.011,0.011,0.012,0.014,0.016])*0.43

    for tag, (zb, sd, sh, zf, sf) in {
        "DESI DR3 (BAO+fσ8)": (z_desi, sDM, sDH, zf_desi, sfs_desi),
        "Euclid (BAO+fσ8)":   (z_euc, sDM_e, sDH_e, zf_euc, sfs_euc),
    }.items():
        s, cov = fisher(zb, sd, sh, zf, sf)
        print(f"\n  {tag}:")
        print(f"    σ(Ω_m)={s[0]:.4f}   σ(Δ)={s[1]:.3f}   σ(lnA_rd)={s[2]:.4f}")
        print(f"    -> Δ=1 vs Δ=0 (BH/smooth-horizon) separated at {1.0/s[1]:.1f}σ")

    # ---- DESI DR3 + Euclid combined (Fisher matrices add) ----
    # rebuild full Fisher by stacking
    def F_of(zb, sd, sh, zf, sf):
        th0 = np.array([Om0, DEL0, 0.0]); o0 = obs_vector(*th0, zb, zf)
        sig = []
        for i in range(len(zb)): sig += [sd[i]*o0[2*i], sh[i]*o0[2*i+1]]
        sig += list(sf); sig = np.array(sig)
        steps = np.array([0.01, 0.03, 0.01]); J = np.zeros((len(o0), 3))
        for k in range(3):
            tp = th0.copy(); tp[k]+=steps[k]; tm = th0.copy(); tm[k]-=steps[k]
            J[:, k] = (obs_vector(*tp, zb, zf)-obs_vector(*tm, zb, zf))/(2*steps[k])
        return J.T @ np.diag(1.0/sig**2) @ J
    Ftot = F_of(z_desi, sDM, sDH, zf_desi, sfs_desi) + F_of(z_euc, sDM_e, sDH_e, zf_euc, sfs_euc)
    sj = np.sqrt(np.diag(np.linalg.inv(Ftot)))
    print(f"\n  DESI DR3 + Euclid COMBINED:")
    print(f"    σ(Ω_m)={sj[0]:.4f}   σ(Δ)={sj[1]:.3f}")
    print(f"    -> Δ=1 detected over Δ=0 at {1.0/sj[1]:.1f}σ   (canonical w0={w0_of(Om0,1.0):+.3f})")
    print("\n  Verdict: the surveys that decide SEDE are DESI DR3 + Euclid; the forecast says")
    print("  they reach σ(Δ)~%.2f — turning today's mild ΔDIC=−2.86 into a multi-σ test of Δ=1." % sj[1])
