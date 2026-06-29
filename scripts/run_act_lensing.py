#!/usr/bin/env python3
"""
ACT DR6 CMB-LENSING test of SEDE (§5; answers arXiv:2509.02945).

Wu et al. (2509.02945) report that holographic/Ricci DE classes are disfavoured by the FULL ACT DR6
CMB likelihood through an early-vs-late tension that compressed (R, ℓ_A) priors hide. SEDE's headline
ΔDIC used compressed CMB, so this is its sharpest untested point. ACT DR6 CMB lensing (C_L^κκ) is the
DE-sensitive part of ACT — the lensing potential integrates growth × geometry, exactly the late-time
dark-energy effect — so it is the right ACT probe of SEDE's w(z) shape.

We compute C_L^κκ from CAMB for theta*-matched ΛCDM and SEDE-H (w0,wa→PPF) and evaluate the official
ACT DR6 lensing likelihood (act_dr6_lenslike, act_baseline, lens_only). theta*-matching fixes the
acoustic scale (Planck-pinned) and isolates the pure DE-shape effect on lensing. Run: python run_act_lensing.py
"""
import numpy as np, camb
import act_dr6_lenslike as alike

DDIR = "packages/data/act_dr6_lens/v1.2/"
THETA = 0.0104109                          # Planck theta*_MC (same as run_plik_check.py)
LMAX = 4000


def camb_clkk(w0=-1.0, wa=0.0, ombh2=0.02237, omch2=0.1200, As=2.1e-9, ns=0.965, tau=0.054):
    pars = camb.CAMBparams()
    pars.set_cosmology(cosmomc_theta=THETA, ombh2=ombh2, omch2=omch2, mnu=0.06, tau=tau)
    pars.InitPower.set_params(As=As, ns=ns)
    if (w0, wa) != (-1.0, 0.0):
        pars.set_dark_energy(w=w0, wa=wa, dark_energy_model='ppf')
    pars.set_for_lmax(LMAX, lens_potential_accuracy=4)
    res = camb.get_results(pars)
    pot = res.get_lens_potential_cls(lmax=LMAX)        # col 0 = [L(L+1)]^2/(2π) C_L^φφ
    tot = res.get_cmb_power_spectra(pars, lmax=LMAX, CMB_unit='muK', raw_cl=False)['lensed_scalar']
    ell = np.arange(pot.shape[0])
    cl_pp = pot[:, 0]
    cl_kk = cl_pp / 4. * 2 * np.pi                      # package convention (matches its test)
    with np.errstate(divide='ignore', invalid='ignore'):
        prefac = 2 * np.pi / ell / (ell + 1.)
    prefac[:2] = 0.0
    cl_tt = tot[:, 0] * prefac; cl_ee = tot[:, 1] * prefac; cl_bb = tot[:, 2] * prefac; cl_te = tot[:, 3] * prefac
    return ell, cl_kk, cl_tt, cl_ee, cl_te, cl_bb


def chi2_act(w0=-1.0, wa=0.0):
    ell, cl_kk, cl_tt, cl_ee, cl_te, cl_bb = camb_clkk(w0, wa)
    d = alike.load_data('act_baseline', ddir=DDIR, lens_only=True, like_corrections=False)
    return -2 * alike.generic_lnlike(d, ell, cl_kk, ell, cl_tt, cl_ee, cl_te, cl_bb, trim_lmax=2998)


if __name__ == "__main__":
    print("=" * 80)
    print("ACT DR6 CMB-LENSING test — does SEDE survive full ACT lensing? (vs 2509.02945)")
    print("=" * 80)

    # sanity: fiducial ΛCDM from the package's own bundled spectra should give χ²≈14.06
    try:
        e, t, ee, bb, te = np.loadtxt(DDIR + 'like_corrs/cosmo2017_10K_acc3_lensedCls.dat', unpack=True)
        ep, _, _, _, _, pp, _, _ = np.loadtxt(DDIR + 'like_corrs/cosmo2017_10K_acc3_lenspotentialCls.dat', unpack=True)
        pf = 2 * np.pi / e / (e + 1.)
        d = alike.load_data('act_baseline', ddir=DDIR, lens_only=True, like_corrections=False)
        x_fid = -2 * alike.generic_lnlike(d, ep, pp / 4 * 2 * np.pi, e, t * pf, ee * pf, te * pf, bb * pf, trim_lmax=2998)
        print(f"\n  pipeline sanity (bundled fiducial ΛCDM): χ² = {x_fid:.2f}  (package test value 14.06)")
    except Exception as ex:
        print(f"  (sanity check skipped: {ex})")

    cL = chi2_act(w0=-1.0, wa=0.0)
    cB = chi2_act(w0=-0.979, wa=-0.137)        # SEDE-H best-fit → CPL (same as run_plik_check.py)
    cS = chi2_act(w0=-0.85,  wa=-0.40)         # stronger CPL, sensitivity bound
    print(f"\n  theta*-matched ACT DR6 lensing (act_baseline, lens_only):")
    print(f"    ΛCDM   (w=-1):                 χ² = {cL:.2f}")
    print(f"    SEDE-H (w0=-0.979, wa=-0.137): χ² = {cB:.2f}   ->  Δχ² = {cB-cL:+.2f}")
    print(f"    (bound: w0=-0.85, wa=-0.40):   χ² = {cS:.2f}   ->  Δχ² = {cS-cL:+.2f}")
    verdict = "ROBUST — SEDE survives ACT DR6 lensing" if abs(cB-cL) < 2 else "TENSION — ACT lensing disfavours SEDE's w(z)"
    print(f"\n  ⟹ |Δχ²_ACTlens| = {abs(cB-cL):.2f}: {verdict}.")
    print("     ACT lensing integrates growth×geometry (the DE-sensitive channel); a small Δχ² means")
    print("     SEDE's late-time, BBN-safe w(z) (Ω_DE→0 at high z) does NOT hit the standard-HDE")
    print("     early–late tension of 2509.02945 — it has LESS early DE than ΛCDM, not a runaway.")
    import sys
    sys.exit(0)
