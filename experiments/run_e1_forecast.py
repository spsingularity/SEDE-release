"""
E1 MECHANISM FORECAST — can Euclid/LSST test SEDE's distinctive claim (dark energy SOURCED by
structure)?  SEDE predicts the deformation inferred from GEOMETRY (BAO/SN → H(z)) equals the one
inferred from GROWTH (RSD fσ8 + weak-lensing σ8(z)): Δ_geom = Δ_growth. A model where DE is NOT
structure-sourced can have Δ_growth ≠ Δ_geom. We forecast σ(Δ)_geom and σ(Δ)_growth SEPARATELY for
DESI DR3 (geometry) and Euclid+LSST (growth), and the consistency precision σ(Δ_geom−Δ_growth) — i.e.
how large a geometry/growth DECOUPLING the surveys could detect. This targets what makes SEDE *SEDE*,
not just another w(z).

Honest scope: simplified Fisher — the weak-lensing-tomography information is represented through its
delivered growth observables (fσ8 + σ8(z) in z-bins at forecast precision), not a full C_ℓ^κκ shear
likelihood. Captures the geometry-vs-growth split; absolute numbers are indicative.

Output: console + results/e1_forecast.json.
"""
import json, os
import numpy as np
from sede import friedmann as fr

C = 299792.458
GAM = 1.4964


def lam(D): return 1.0 - D / 2.0
def Efun(z, Om, D): return fr.E_SEDE_lambda(np.atleast_1d(z), Om, GAM, lam(D))


def bao_vec(Om, D, lnA, zb):
    A = np.exp(lnA)
    out = []
    for z in zb:
        zz = np.linspace(0, z, 200)
        I = np.trapezoid(1.0 / Efun(zz, Om, D), zz)
        out.append(A * I)                                  # DM/rd ∝ A·∫dz/E
    for z in zb:
        out.append(A / float(Efun(np.array([z]), Om, D)[0]))   # DH/rd ∝ A/E
    return np.array(out)


def fs8_vec(Om, D, s8, zf):
    Dg, f = fr.compute_growth_model(np.asarray(zf), Om, lambda z: Efun(z, Om, D))
    return s8 * f * Dg                                     # fσ8 = σ8·f·D


def deriv(func, theta, i, dx):
    tp = list(theta); tm = list(theta); tp[i] += dx; tm[i] -= dx
    return (func(*tp) - func(*tm)) / (2 * dx)


def fisher(func, theta, sig, steps):
    n = len(theta)
    d = [deriv(func, theta, i, steps[i]) for i in range(n)]
    w = 1.0 / np.asarray(sig)**2
    return np.array([[np.sum(d[i] * w * d[j]) for j in range(n)] for i in range(n)])


def sigma_param(F, idx):
    return float(np.sqrt(np.linalg.inv(F)[idx, idx]))


def main():
    Om0, D0, lnA0, s80 = 0.30, 1.0, 0.0, 0.78

    # ── GEOMETRY: DESI DR3 BAO (θ = Ω_m, Δ, lnA) ──────────────────────────────
    zb = np.array([0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.8, 2.3])
    fid_b = bao_vec(Om0, D0, lnA0, zb)
    sig_b = np.concatenate([0.005 * fid_b[:len(zb)], 0.009 * fid_b[len(zb):]])   # 0.5% DM, 0.9% DH
    Fg = fisher(lambda Om, D, lnA: bao_vec(Om, D, lnA, zb), (Om0, D0, lnA0), sig_b,
                (0.003, 0.02, 0.01))
    sD_geom = sigma_param(Fg, 1)

    # ── GROWTH: Euclid+LSST fσ8 / σ8(z) tomography (θ = Ω_m, Δ, σ8) ────────────
    zf = np.array([0.3, 0.5, 0.7, 0.9, 1.0, 1.2, 1.4, 1.65, 1.8])               # LSST(z<1)+Euclid(z>0.9)
    frac = np.array([0.015, 0.015, 0.015, 0.013, 0.020, 0.020, 0.020, 0.022, 0.022])  # ~1.5% LSST, ~2% Euclid
    fid_f = fs8_vec(Om0, D0, s80, zf)
    sig_f = frac * fid_f
    Fr = fisher(lambda Om, D, s8: fs8_vec(Om, D, s8, zf), (Om0, D0, s80), sig_f,
                (0.003, 0.02, 0.01))
    sD_growth = sigma_param(Fr, 1)

    # ── E1 consistency: σ(Δ_geom − Δ_growth), and detectability of a decoupling ──
    sD_E1 = np.sqrt(sD_geom**2 + sD_growth**2)
    # combined Δ measurement IF structure-sourcing assumed (Δ shared): add Δ-blocks
    sD_comb = 1.0 / np.sqrt(1.0 / sD_geom**2 + 1.0 / sD_growth**2)

    print("=" * 76)
    print("E1 MECHANISM FORECAST — geometry-Δ vs growth-Δ (structure-sourcing test)")
    print("=" * 76)
    print(f"  σ(Δ)_geom    [DESI DR3 BAO]            = {sD_geom:.3f}")
    print(f"  σ(Δ)_growth  [Euclid+LSST fσ8/σ8(z)]   = {sD_growth:.3f}")
    print(f"  σ(Δ_geom − Δ_growth)  (E1 consistency) = {sD_E1:.3f}")
    print(f"  σ(Δ) combined (structure-sourcing assumed) = {sD_comb:.3f}")
    print("-" * 76)
    print(f"  ⇒ a geometry/growth DECOUPLING of Δ_growth−Δ_geom = 1 (a maximally non-structure-sourced")
    print(f"    rival) is detectable at {1.0/sD_E1:.1f}σ; a decoupling of 0.5 at {0.5/sD_E1:.1f}σ.")
    print(f"  ⇒ Δ=0 (area-law) vs Δ=1 separated by geometry alone at {1.0/sD_geom:.0f}σ, by growth alone at {1.0/sD_growth:.0f}σ.")
    print("  HONEST: this is the test of SEDE's *mechanism* (DE tied to growth), the one probe that")
    print("  distinguishes SEDE from a generic w(z). Simplified Fisher (fσ8/σ8 proxy for WL tomography).")

    out = dict(sigma_Delta_geom=sD_geom, sigma_Delta_growth=sD_growth,
               sigma_Delta_E1=sD_E1, sigma_Delta_combined=sD_comb,
               detect_decoupling_1=1.0 / sD_E1, detect_decoupling_0p5=0.5 / sD_E1)
    os.makedirs("results", exist_ok=True)
    json.dump(out, open("results/e1_forecast.json", "w"), indent=2)
    print("\nwrote results/e1_forecast.json")
    return out


if __name__ == "__main__":
    main()
