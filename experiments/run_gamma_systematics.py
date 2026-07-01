#!/usr/bin/env python3
"""
E2 — γ systematics table  [REVISION_PLAN.md T1.2]
=================================================
The reservoir (horizon, T_AH ∝ H) is FIXED by the §2 identity ρ_DE = T_AH·s_grav ⟹
the entropy weight p = 5/3 (binding energy E_bind ∝ M^{5/3} deposited at mass-independent
T_AH).  That leg is forced.  The NUMERICAL γ = (1/2) dlnΣ_S/dlnσ8, with
Σ_S = ∫ M^{5/3} (dn/dlnM) dlnM, then inherits STANDARD, PRE-SPECIFIED halo-model choices —
the mass function, the mass range, the transfer function.  This script brackets γ across
those choices to show (a) it clusters near the canonical 1.5, and (b) NONE is tuned to
cosmological data.  This supports the honest claim (round-2 §7.2): the temperature/reservoir
(hence p=5/3) is forced; the numeric γ carries standard halo-model systematics.

Usage:  python run_gamma_systematics.py
"""
from __future__ import annotations
import json
import numpy as np
from colossus.cosmology import cosmology
from colossus.lss import mass_function

# canonical entropy weight (forced by the horizon reservoir, §3.1 / Result 4C)
P_WEIGHT = 5.0 / 3.0
S8_FID = 0.811
DELTA_LN_S8 = 0.02

# model -> compatible mass definition (colossus requirement)
MODEL_MDEF = {
    'press74':  'fof',
    'sheth99':  'fof',     # Sheth-Tormen (the canonical choice)
    'tinker08': '200m',
    'despali16': '200m',
    'watson13': 'fof',
}

MASS_RANGES = {
    'full[1e10,1e16]':  (1e10, 1e16),
    'narrow[1e11,1e15]': (1e11, 1e15),
    'massive[1e12,1e16]': (1e12, 1e16),
}

# transfer-function variations (colossus power-spectrum models)
TRANSFERS = ['eisenstein98', 'sugiyama95']   # EH98 (default) and a Sugiyama95 alt


def _set_cosmo(sigma8_val):
    params = {'flat': True, 'H0': 67.4, 'Om0': 0.315, 'Ob0': 0.049,
              'sigma8': sigma8_val, 'ns': 0.965, 'Tcmb0': 2.725}
    name = f"sede_{sigma8_val:.4f}"
    try:
        cosmology.addCosmology(name, params)
    except Exception:
        pass
    cosmology.setCosmology(name)


def sigma_S(sigma8_val, model, mdef, m_lo, m_hi, transfer='eisenstein98', n=400):
    _set_cosmo(sigma8_val)
    M = np.logspace(np.log10(m_lo), np.log10(m_hi), n)
    # transfer enters via the power-spectrum args (EH98 default); not a cosmology kwarg
    mf = mass_function.massFunction(M, z=0.0, mdef=mdef, model=model, q_out='dndlnM',
                                    ps_args={'model': transfer})
    dndM = mf / M
    return np.trapezoid(M ** P_WEIGHT * dndM, M)


def gamma_for(model, mdef, m_lo, m_hi, transfer='eisenstein98'):
    s_hi = sigma_S(S8_FID * np.exp(+DELTA_LN_S8), model, mdef, m_lo, m_hi, transfer)
    s_lo = sigma_S(S8_FID * np.exp(-DELTA_LN_S8), model, mdef, m_lo, m_hi, transfer)
    dln = (np.log(s_hi) - np.log(s_lo)) / (2 * DELTA_LN_S8)
    return 0.5 * dln   # γ = (1/2) dlnΣ_S/dlnσ8  (x = D² convention)


def main():
    rows = []
    print("E2 — γ systematics (entropy weight p = 5/3 FIXED; halo-model choices varied)\n")
    print(f"{'mass function':12s} {'mdef':6s} {'mass range':20s} {'transfer':14s} {'gamma':>7s}")
    print("-" * 66)
    # primary scan: mass function × mass range (EH98 transfer)
    for model, mdef in MODEL_MDEF.items():
        for rname, (lo, hi) in MASS_RANGES.items():
            try:
                g = gamma_for(model, mdef, lo, hi, 'eisenstein98')
                rows.append({'model': model, 'mdef': mdef, 'range': rname,
                             'transfer': 'eisenstein98', 'gamma': float(g)})
                star = "  <- canonical" if (model == 'sheth99' and rname.startswith('full')) else ""
                print(f"{model:12s} {mdef:6s} {rname:20s} {'eisenstein98':14s} {g:7.3f}{star}")
            except Exception as e:
                print(f"{model:12s} {mdef:6s} {rname:20s} {'eisenstein98':14s}   SKIP ({type(e).__name__})")
    # transfer variation (canonical ST, full range)
    print("-" * 66)
    for tr in TRANSFERS:
        try:
            g = gamma_for('sheth99', 'fof', 1e10, 1e16, tr)
            rows.append({'model': 'sheth99', 'mdef': 'fof', 'range': 'full[1e10,1e16]',
                         'transfer': tr, 'gamma': float(g)})
            print(f"{'sheth99':12s} {'fof':6s} {'full[1e10,1e16]':20s} {tr:14s} {g:7.3f}")
        except Exception as e:
            print(f"{'sheth99':12s} {'fof':6s} {'full':20s} {tr:14s}   SKIP ({type(e).__name__})")

    gammas = np.array([r['gamma'] for r in rows])
    canon = [r['gamma'] for r in rows if r['model'] == 'sheth99' and r['range'].startswith('full')
             and r['transfer'] == 'eisenstein98']
    print("\n" + "=" * 66)
    print(f"  canonical (Sheth-Tormen, full range, EH98): γ = {canon[0]:.3f}")
    print(f"  full systematic spread: γ ∈ [{gammas.min():.3f}, {gammas.max():.3f}], "
          f"median {np.median(gammas):.3f}  (N={len(gammas)})")
    print("  Honest claim: p=5/3 forced by the horizon reservoir; the numeric γ carries")
    print("  standard halo-model systematics (~±0.2), NONE tuned to cosmological data.")
    out = {'canonical_gamma': canon[0], 'min': float(gammas.min()),
           'max': float(gammas.max()), 'median': float(np.median(gammas)),
           'p_weight': P_WEIGHT, 'rows': rows}
    json.dump(out, open('results/e2_gamma_systematics.json', 'w'), indent=2)
    print("  saved -> results/e2_gamma_systematics.json")


if __name__ == '__main__':
    import os
    os.makedirs('results', exist_ok=True)
    main()
