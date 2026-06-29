"""
Full weak-lensing tomographic Fisher forecast (C_ℓ^κκ, cross-bin correlations) for σ(Δ)_growth —
the proper version of the E1 mechanism forecast (refines the fσ8/σ8 proxy in run_e1_forecast.py).

Stage-IV cosmic shear (Euclid+LSST-class): N tomographic bins, Limber C_ℓ^{ij} from the 3D matter
power, Gaussian covariance with shape noise, trace-Fisher over θ=(Ω_m, σ8, Δ). SEDE's Δ enters ONLY
through the growth factor D(z;Δ) and the distances χ(z;Δ),H(z;Δ) (dark energy is smooth, c_s²=1, so
the transfer-function shape is Δ-independent). Then recompute the E1 consistency σ(Δ_geom−Δ_growth)
with the geometry σ(Δ)_geom from run_e1_forecast.

Honest scope: linear P(k) (BBKS+Sugiyama Γ) to ℓ_max=1500 (quasi-linear; conservative — nonlinear
scales add growth info but also baryonic/systematic uncertainty); θ=(Ω_m,σ8,Δ) with h,n_s,Ω_b fixed
(matched to the proxy's 3-parameter space). Output: console + results/wl_fisher.json.
"""
import json, os
import numpy as np
from scipy.special import erf
from sede import friedmann as fr

C = 299792.458
GAM = 1.4964
H0H = 2997.92458   # c/H0 in Mpc/h  (distances in Mpc/h are H0-independent)

# fixed (marginalised-away) cosmology
H, OB, NS = 0.67, 0.049, 0.965
FID = dict(Om=0.30, s8=0.80, D=1.0)

# survey (Stage-IV combined)
FSKY, NGAL, SIG_E, NBIN, Z0 = 0.4, 30.0, 0.26, 10, 0.70
LMIN, LMAX, NL = 20, 1500, 18


def lam(D): return 1.0 - D / 2.0
def Efun(z, Om, D): return np.asarray(fr.E_SEDE_lambda(np.atleast_1d(z), Om, GAM, lam(D)))


def T_bbks(k, Om):                       # k in h/Mpc; Sugiyama-corrected BBKS transfer
    G = Om * H * np.exp(-OB - np.sqrt(2 * H) * OB / Om)
    q = k / G
    return (np.log(1 + 2.34 * q) / (2.34 * q) *
            (1 + 3.89 * q + (16.1 * q)**2 + (5.46 * q)**3 + (6.71 * q)**4)**-0.25)


def Pk_norm(Om, s8):                     # amplitude A so σ8 matches, for P=A k^ns T²
    k = np.logspace(-4, 2, 3000)
    Pshape = k**NS * T_bbks(k, Om)**2
    x = k * 8.0
    W = 3 * (np.sin(x) - x * np.cos(x)) / x**3
    s2 = np.trapezoid(k**2 * Pshape * W**2, k) / (2 * np.pi**2)
    return s8**2 / s2


def Plin(k, Om, s8):
    return Pk_norm(Om, s8) * k**NS * T_bbks(k, Om)**2


def nz_bins(zg):                         # equal-number tomographic bins + photo-z smoothing
    nz = zg**2 * np.exp(-(zg / Z0)**1.5)
    cdf = np.concatenate([[0], np.cumsum(0.5 * (nz[1:] + nz[:-1]) * np.diff(zg))])
    cdf /= cdf[-1]
    edges = np.interp(np.linspace(0, 1, NBIN + 1), cdf, zg)
    sigz = 0.05 * (1 + zg)
    bins = []
    for a, b in zip(edges[:-1], edges[1:]):
        wsel = 0.5 * (erf((b - zg) / (np.sqrt(2) * sigz)) - erf((a - zg) / (np.sqrt(2) * sigz)))
        ni = nz * wsel
        ni /= np.trapezoid(ni, zg)       # ∫ n_i dz = 1
        bins.append(ni)
    return bins


def cl_matrix(Om, s8, D, ells, zg, bins):
    E = Efun(zg, Om, D)
    chi = H0H * np.concatenate([[0], np.cumsum(0.5 * (1 / E[1:] + 1 / E[:-1]) * np.diff(zg))])  # Mpc/h
    Dg, _ = fr.compute_growth_model(zg, Om, lambda z: Efun(z, Om, D))   # D(z), D(0)=1
    # convert n_i(z) → n_i(χ): n_i(χ)dχ = n_i(z)dz, dχ/dz = H0H/E
    dchidz = H0H / E
    # lensing efficiency q_i(χ) = (3/2)Ω_m (1/H0H)² (1+z) χ ∫_z^∞ n_i(z')(χ'-χ)/χ' dz'
    q = []
    for ni in bins:
        g = np.array([np.trapezoid((ni * (chi - chi[j]) / np.where(chi > 0, chi, 1))[j:], zg[j:])
                      for j in range(len(zg))])
        g = np.clip(g, 0, None)
        qi = 1.5 * Om / H0H**2 * (1 + zg) * chi * g
        q.append(qi)
    q = np.array(q)
    # Limber: C_ℓ^{ij} = ∫ dz (dχ/dz) q_i q_j / χ² P(k=(ℓ+0.5)/χ, z)
    nb = len(bins); Cl = np.zeros((len(ells), nb, nb))
    kgrid = np.logspace(-4, 2, 3000); P0 = Plin(kgrid, Om, s8)
    m = chi > 1e-3
    base = dchidz[m] / chi[m]**2                                # dχ/dz · 1/χ²
    for li, l in enumerate(ells):
        k = (l + 0.5) / chi[m]
        Pk = np.interp(k, kgrid, P0) * Dg[m]**2                 # P(k,z) along the line of sight
        w = base * Pk                                          # common weight
        for i in range(nb):
            for j in range(i, nb):
                Cl[li, i, j] = Cl[li, j, i] = np.trapezoid(w * q[i, m] * q[j, m], zg[m])
    return Cl


def fisher():
    zg = np.linspace(1e-3, 3.2, 240)
    bins = nz_bins(zg)
    ells = np.unique(np.logspace(np.log10(LMIN), np.log10(LMAX), NL).astype(int)).astype(float)
    # shape noise N^{ii} = σ_ε² / n_i  (n_i in sr^-1)
    n_sr = (NGAL / NBIN) * (180 * 60 / np.pi)**2
    N = np.eye(NBIN) * SIG_E**2 / n_sr
    th = [FID['Om'], FID['s8'], FID['D']]
    step = [0.004, 0.01, 0.03]
    Cfid = cl_matrix(*th, ells, zg, bins)
    dC = []
    for i in range(3):
        tp = th[:]; tm = th[:]; tp[i] += step[i]; tm[i] -= step[i]
        dC.append((cl_matrix(*tp, ells, zg, bins) - cl_matrix(*tm, ells, zg, bins)) / (2 * step[i]))
    # ℓ-bin widths
    dl = np.gradient(ells)
    F = np.zeros((3, 3))
    for li in range(len(ells)):
        Ctot = Cfid[li] + N
        Cinv = np.linalg.inv(Ctot)
        pref = (2 * ells[li] + 1) / 2 * FSKY * dl[li]
        for a in range(3):
            for b in range(3):
                F[a, b] += pref * np.trace(Cinv @ dC[a][li] @ Cinv @ dC[b][li])
    return F


def main():
    print("=" * 76)
    print("WEAK-LENSING TOMOGRAPHIC FISHER (C_ℓ^κκ) — σ(Δ)_growth, the proper E1 forecast")
    print(f"  Stage-IV: f_sky={FSKY}, n_gal={NGAL}/arcmin², {NBIN} bins, ℓ={LMIN}-{LMAX}, θ=(Ω_m,σ8,Δ)")
    print("=" * 76)
    F = fisher()
    cov = np.linalg.inv(F)
    sOm, sS8, sD_wl = np.sqrt(np.diag(cov))
    print(f"  σ(Ω_m)={sOm:.4f}   σ(σ8)={sS8:.4f}   σ(Δ)_growth[WL]={sD_wl:.3f}")
    # geometry σ(Δ) from the proxy forecast
    import run_e1_forecast as E1
    e1 = E1.main()
    sD_geom = e1["sigma_Delta_geom"]; sD_proxy = e1["sigma_Delta_growth"]
    sD_E1 = np.sqrt(sD_geom**2 + sD_wl**2)
    print("-" * 76)
    print(f"  σ(Δ)_geom  [DESI DR3]               = {sD_geom:.3f}")
    print(f"  σ(Δ)_growth: fσ8 proxy {sD_proxy:.3f}  →  full WL C_ℓ {sD_wl:.3f}  ({'tighter' if sD_wl<sD_proxy else 'looser'})")
    print(f"  σ(Δ_geom − Δ_growth)  (E1 consistency, full WL) = {sD_E1:.3f}")
    print(f"  ⇒ geometry/growth decoupling Δ_growth−Δ_geom=1 detectable at {1/sD_E1:.1f}σ "
          f"(was {1/np.sqrt(sD_geom**2+sD_proxy**2):.1f}σ with the proxy)")
    print("  HONEST: linear P(k)≤ℓ1500, h/n_s/Ω_b fixed → indicative; nonlinear scales would tighten σ(Δ)_growth.")
    out = dict(sigma_Om=sOm, sigma_s8=sS8, sigma_Delta_growth_WL=float(sD_wl),
               sigma_Delta_geom=sD_geom, sigma_Delta_E1=float(sD_E1),
               detect_decoupling_1=float(1 / sD_E1))
    os.makedirs("results", exist_ok=True)
    json.dump(out, open("results/wl_fisher.json", "w"), indent=2)
    print("\nwrote results/wl_fisher.json")
    return out


if __name__ == "__main__":
    main()
