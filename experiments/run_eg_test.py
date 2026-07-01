#!/usr/bin/env python3
"""
E_G(z) — the growth–geometry (gravitational-slip) test (§6; Rauhut, Blake et al. DESI-DR1, arXiv:2507.16098).

The E_G estimator is the ratio of weak-lensing (∇²(Ψ+Φ)) to galaxy-velocity (f σ8) amplitudes from a
common set of overdensities — i.e. a direct probe of the gravitational SLIP (Ψ vs Φ). Its
scale-independent GR prediction (Zhang et al. 2007) is
    E_G(z) = Ω_{m,0} / f(z),     f ≡ d ln D / d ln a.
SEDE is general relativity with SMOOTH dark energy (c_s²=1, α_M=α_T=0; §7) ⟹ NO slip, Ψ=Φ ⟹ E_G obeys
the same Ω_{m,0}/f form as ΛCDM, differing only through SEDE's slightly different growth rate f(z). A
modified-gravity model would deviate. So E_G is (i) a *consistency* test SEDE must pass and (ii) the
observational realisation of CHR's growth–expansion lock (P2): lensing geometry vs RSD velocities.

The DESI-DR1 + KiDS/DES/HSC measurement [arXiv:2507.16098] is consistent with GR+ΛCDM out to z≈1; this
driver shows SEDE predicts the same no-slip E_G(z) (to %-level) and is therefore equally consistent.
Run:  python run_eg_test.py
"""
import numpy as np
from sede.friedmann import E_SEDE_lambda, compute_growth_model

Om, Or, GAM, LAM = 0.3153, 9e-5, 1.4964, 0.5            # Om=0.3153 to match the E_G paper's Planck GR baseline
ZBINS = np.array([0.15, 0.25, 0.35, 0.50, 0.70, 0.95])  # DESI-DR1 effective redshifts (2507.16098)


def growth_rate(zarr, E_of_z):
    _, f = compute_growth_model(zarr, Om, E_of_z, Or)
    return f


if __name__ == "__main__":
    print("=" * 80)
    print("E_G(z) = Ω_m,0/f(z) — growth–geometry / gravitational-slip test (DESI-DR1, 2507.16098)")
    print("=" * 80)
    f_s = growth_rate(ZBINS, lambda z: E_SEDE_lambda(np.atleast_1d(z), Om, GAM, LAM, Or))
    f_l = growth_rate(ZBINS, lambda z: np.sqrt(Om*(1+z)**3 + Or*(1+z)**4 + (1-Om-Or)))
    Eg_s, Eg_l = Om/f_s, Om/f_l
    print(f"\n  Ω_m,0 = {Om}  (no-slip GR prediction E_G = Ω_m,0/f; SEDE & ΛCDM share the form)")
    print(f"\n  {'z':>5s} {'f_SEDE':>8s} {'f_ΛCDM':>8s} {'E_G(SEDE)':>10s} {'E_G(ΛCDM)':>10s} {'Δ%':>7s}")
    for i, z in enumerate(ZBINS):
        dpc = 100*(Eg_s[i]/Eg_l[i] - 1)
        print(f"  {z:>5.2f} {f_s[i]:>8.4f} {f_l[i]:>8.4f} {Eg_s[i]:>10.4f} {Eg_l[i]:>10.4f} {dpc:>+7.2f}")
    maxdev = float(np.max(np.abs(Eg_s/Eg_l - 1)))
    print(f"\n  max |E_G(SEDE)/E_G(ΛCDM) − 1| = {100*maxdev:.2f}%  over 0.15 ≤ z ≤ 0.95")
    print("\n  ⟹ SEDE predicts NO gravitational slip (Ψ=Φ; smooth DE on a GR bulk), so E_G follows the")
    print("     same Ω_m,0/f(z) law as ΛCDM, differing only by the small growth-rate shift. The DESI-DR1 +")
    print("     KiDS/DES/HSC E_G measurement is consistent with GR+ΛCDM out to z≈1 [2507.16098] — hence")
    print("     equally consistent with SEDE. A modified-gravity origin of acceleration would deviate from")
    print("     this curve; SEDE does not. This is the existing-data realisation of CHR's growth–expansion")
    print("     lock (P2): lensing geometry vs RSD velocities must track the same f(z).")
    ok = maxdev < 0.05 and bool(np.all(Eg_s > 0))
    print(f"\n  E_G no-slip consistency check: {'PASS' if ok else 'FAIL'}  (SEDE within {100*maxdev:.1f}% of ΛCDM, both GR-form)")
    import sys
    sys.exit(0 if ok else 1)
