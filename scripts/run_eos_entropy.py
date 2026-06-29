#!/usr/bin/env python3
"""
Theorem 13 — the dark-energy EOS as a readout of horizon entropy production.
Derivation: ρ_DE ∝ H^{2λ} f_sat ⟹  1+w_DE = (1/3)[ 2λ·ε  −  d ln f_sat/d ln a ],  ε ≡ −d ln H/d ln a.
So w_DE−(−1) is the COMPETITION between cosmic expansion (2λε, pushes quintessence, w>−1) and
horizon ENTROPY PRODUCTION (d ln f_sat/d ln a, pushes phantom, w<−1). The −1 crossing is exactly
where they balance; both de Sitter brackets (ε→0, df_sat→0) give w→−1. This driver verifies the
identity against the direct fluid EOS and shows the decomposition. Run: python run_eos_entropy.py
"""
import numpy as np
from sede import friedmann as fr

Om, Or, gamma, lam = 0.30, 9e-5, 1.4964, 0.5
def E(z): return fr.E_SEDE_lambda(np.atleast_1d(z), Om, gamma, lam)
def fsat(z):
    D = fr.compute_growth_factor(np.atleast_1d(z), Om, Or)
    return np.clip((1-np.exp(-gamma*D**2))/(1-np.exp(-gamma)), 1e-30, 1)

if __name__ == "__main__":
    print("="*78); print("THEOREM 13 — dark-energy EOS = expansion vs horizon entropy production"); print("="*78)
    z = np.linspace(0, 3, 400)
    Ez = E(z); lnE = np.log(Ez); eps = (1+z)*np.gradient(lnE, z)          # ε = -dlnH/dlna = (1+z) dlnE/dz
    lnf = np.log(fsat(z)); dlnf_dlna = -(1+z)*np.gradient(lnf, z)         # dlnf/dlna = -(1+z) dlnf/dz
    expansion = (2*lam*eps)/3.0                                          # pushes w > -1 (quintessence)
    entropy   = dlnf_dlna/3.0                                            # pushes w < -1 (phantom)
    w_formula = -1 + expansion - entropy                                 # Theorem 13
    rho = np.maximum(Ez**2 - Om*(1+z)**3 - Or*(1+z)**4, 1e-12)
    w_direct = -1 + (1/3.)*(1+z)*np.gradient(np.log(rho), z)             # direct fluid EOS
    resid = np.max(np.abs(w_formula - w_direct)[5:-5])

    print(f"\n  identity check: max|w_formula − w_direct| over z∈[0,3] = {resid:.2e}  "
          f"{'✓ verified' if resid<1e-2 else 'FAIL'}")
    print(f"\n  {'z':>4s} {'ε':>6s} {'expansion(+)':>12s} {'entropyProd(−)':>14s} {'1+w':>8s} {'w':>8s}")
    for zz in (0.0, 0.3, 0.5, 1.0, 2.0, 3.0):
        i = np.argmin(np.abs(z-zz))
        print(f"  {zz:>4.1f} {eps[i]:>6.2f} {expansion[i]:>12.3f} {entropy[i]:>14.3f} "
              f"{(w_formula[i]+1):>8.3f} {w_formula[i]:>8.3f}")
    zc = z[np.argmin(np.abs(w_formula+1))]
    print(f"\n  w crosses −1 at z≈{zc:.2f}: where 2λε = dlnf_sat/dlna (expansion = entropy production).")
    print("  ⟹ PHYSICAL READING: the dark-energy EOS is a direct meter of the universe's entropy")
    print("  production. w<−1 (phantom) = horizon entropy growing FASTER than expansion dilutes it")
    print("  (structure-formation era); w>−1 = expansion winning; w=−1 = the two balance (the de Sitter")
    print("  brackets, ε→0 & df_sat→0). Measuring w(z) measures dS_horizon/dt — the GSL made observable.")
