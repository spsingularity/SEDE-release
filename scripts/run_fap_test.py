#!/usr/bin/env python3
"""
F_AP(z) — the r_d-INDEPENDENT (and H0-independent) BAO-shape test (Turyshev 2026, arXiv:2602.05368).

Turyshev's review proposes F_AP(z) ≡ D_M(z)/D_H(z) as a calibration-free diagnostic built directly
from the published DESI (D_M/r_d, D_H/r_d): the sound horizon r_d cancels in the ratio. It is also
H0-independent — F_AP(z) = E(z)·∫₀^z dz'/E(z') depends ONLY on the expansion *shape* E(z)=H/H0. So it
tests SEDE's w(z) shape with NO reliance on r_d, H0, the CMB sound-horizon, or SN calibration — the
exact systematics the post-DESI robustness debate (Wang & Mota 2025; Liu et al 2024) centres on.

We compute F_AP from the DESI DR2 tracer bins that report both D_M/r_d and D_H/r_d (LRG1, LRG2,
LRG+ELG, ELG, QSO, Lyα), propagating each bin's 2×2 covariance, and compare SEDE vs ΛCDM (matched Ω_m).
Run:  python run_fap_test.py
"""
import numpy as np
from scipy.integrate import quad
from sede import friedmann as fr
from sede.data_loader import load_desi_dr2

Om, Or, GAM, LAM = 0.30, 9e-5, 1.4964, 0.5
TRACER = {0.510: 'LRG1', 0.706: 'LRG2', 0.934: 'LRG+ELG', 1.321: 'ELG', 1.484: 'QSO', 2.330: 'Lya'}


def E_sede(z): return float(fr.E_SEDE_lambda(np.atleast_1d(z), Om, GAM, LAM, Or)[0])
def E_lcdm(z): return float(np.sqrt(Om*(1+z)**3 + Or*(1+z)**4 + (1-Om-Or)))


def fap_model(z, E):
    """F_AP(z) = D_M/D_H = E(z)·∫₀^z dz'/E(z')  (r_d- and H0-independent; pure shape)."""
    integ = quad(lambda zp: 1.0/E(zp), 0.0, z)[0]
    return E(z) * integ


def fap_data():
    """Build F_AP = (D_M/r_d)/(D_H/r_d) per tracer bin, with error from the 2×2 sub-covariance."""
    z, t, m, cov = load_desi_dr2()
    out = []
    for z0, name in TRACER.items():
        idx = [i for i in range(len(z)) if abs(z[i]-z0) < 2e-3]
        if len(idx) != 2:
            continue
        # identify which row is DM/rd vs DH/rd
        iDM = [i for i in idx if t[i] == 'DM/rd'][0]; iDH = [i for i in idx if t[i] == 'DH/rd'][0]
        a, b = m[iDM], m[iDH]                          # a = DM/rd, b = DH/rd
        F = a/b                                        # r_d cancels
        ca, cb = cov[iDM, iDM], cov[iDH, iDH]; cab = cov[iDM, iDH]
        dFda, dFdb = 1.0/b, -a/b**2                    # error propagation F=a/b
        varF = dFda**2*ca + dFdb**2*cb + 2*dFda*dFdb*cab
        out.append((z0, name, F, np.sqrt(varF)))
    return out


if __name__ == "__main__":
    print("=" * 82)
    print("F_AP(z) = D_M/D_H — r_d- & H0-INDEPENDENT BAO-shape test (Turyshev 2026)")
    print("=" * 82)
    data = fap_data()
    print(f"\n  {'z':>5s} {'tracer':>8s} {'F_AP(data)':>12s} {'±':>7s} {'SEDE':>8s} {'ΛCDM':>8s} "
          f"{'(SEDE−d)/σ':>11s} {'(Λ−d)/σ':>9s}")
    chi2_s = chi2_l = 0.0
    for z0, name, F, sF in data:
        fs = fap_model(z0, E_sede); fl = fap_model(z0, E_lcdm)
        ps, pl = (fs-F)/sF, (fl-F)/sF; chi2_s += ps**2; chi2_l += pl**2
        print(f"  {z0:>5.2f} {name:>8s} {F:>12.3f} {sF:>7.3f} {fs:>8.3f} {fl:>8.3f} "
              f"{ps:>11.2f} {pl:>9.2f}")
    n = len(data)
    print(f"\n  χ²(SEDE)  = {chi2_s:.2f}  (n={n} r_d-independent points, χ²/n={chi2_s/n:.2f})")
    print(f"  χ²(ΛCDM)  = {chi2_l:.2f}  (χ²/n={chi2_l/n:.2f})")
    print(f"  Δχ²(SEDE−ΛCDM) = {chi2_s-chi2_l:+.2f}  (negative ⟹ SEDE's *shape* fits the calibration-free")
    print(f"                   F_AP better; both Ω_m matched, so this isolates the w(z) shape)")
    print(f"\n  ⟹ F_AP removes r_d, H0, the CMB sound horizon, and SN calibration — the systematics the")
    print(f"     post-DESI robustness debate (Wang & Mota; Liu et al) centres on. SEDE's evolving-w shape")
    print(f"     is {'consistent with' if chi2_s/n < 2 else 'in tension with'} the DESI DR2 BAO shape on this clean diagnostic (per-bin 2×2 cov;")
    print(f"     cross-bin correlations neglected — small for independent tracers).")
    import sys
    sys.exit(0)
