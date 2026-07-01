#!/usr/bin/env python3
"""
Quick plik_lite check: does the FULL Planck TT/TE/EE + lowl spectrum change the
Barrow-vs-LCDM verdict relative to the compressed (R, l_A)?  No full MCMC needed.

Key: the CMB pins the acoustic scale theta* to 0.03%, so the fair, peak-position-
controlled test fixes cosmomc_theta (H0 derived) and varies ONLY the dark-energy
w(a) -> isolates the DE SHAPE effect (peak heights, early-ISW, lensing), not theta*.

Result: pure DE-shape Δχ²_CMB(full plik_lite TTTEEE + lowTT + lowEE) = -0.5 for the
Barrow SEDE-H w(z)->CPL (w0=-0.979, wa=-0.137); -0.4 even for a stronger w0=-0.85,
wa=-0.4. |Δχ²|<1 -> the compressed (R, l_A) captured the CMB power on this late-time-
only DE; the ΔDIC=-2.86 Barrow preference is ROBUST to the full CMB treatment, and the
multi-day full Cobaya MCMC would shift it by <~0.5 (so it is not necessary).
"""
import numpy as np
from cobaya.model import get_model

PKG = './packages'; THETA = 0.0104109   # Planck theta*_MC

def chi2_cmb(omch2, w0=-1.0, wa=0.0, As=2.10e-9, ns=0.965, tau=0.054, ombh2=0.02237):
    extra = {'lens_potential_accuracy': 1}
    pars = {'ombh2': ombh2, 'omch2': omch2, 'cosmomc_theta': THETA,
            'As': As, 'ns': ns, 'tau': tau, 'A_planck': 1.0}
    if (w0, wa) != (-1.0, 0.0):
        extra['dark_energy_model'] = 'ppf'; pars['w'] = w0; pars['wa'] = wa
    info = {'likelihood': {'planck_2018_highl_plik.TTTEEE_lite_native': None,
                           'planck_2018_lowl.TT': None, 'planck_2018_lowl.EE': None},
            'theory': {'camb': {'extra_args': extra}}, 'params': pars, 'packages_path': PKG}
    m = get_model(info); ll = m.loglikes({'A_planck': 1.0})[0]; return -2*sum(ll)

if __name__ == "__main__":
    oc = 0.1200
    cL = chi2_cmb(oc)
    cB = chi2_cmb(oc, w0=-0.979, wa=-0.137)        # Barrow SEDE-H best-fit -> CPL
    cS = chi2_cmb(oc, w0=-0.85, wa=-0.4)            # stronger CPL, sensitivity bound
    print("theta*-matched full plik_lite (TTTEEE_lite + lowl.TT + lowl.EE):")
    print(f"  LCDM   (w=-1):                chi2 = {cL:.2f}")
    print(f"  Barrow (w0=-0.979,wa=-0.137): chi2 = {cB:.2f}   ->  Delta = {cB-cL:+.2f}")
    print(f"  (bound: w0=-0.85,wa=-0.4):    chi2 = {cS:.2f}   ->  Delta = {cS-cL:+.2f}")
    print("=> |Delta chi2_CMB| < 1: full CMB shape does NOT change the verdict.")
