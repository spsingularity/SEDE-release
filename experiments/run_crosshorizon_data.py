#!/usr/bin/env python3
"""
TIER 3 — cross-horizon DATA confrontation. Does the cosmologically-fitted Δ=1 survive
contact with REAL black-hole observations?  We test SEDE's universal-Δ prediction against:
  (1) the GWTC binary-black-hole catalog — Hawking area theorem + the Barrow generalised
      second law, event by event (real masses/spins);
  (2) the primordial-black-hole evaporation channel — the ONE place Δ is observable.

Honest framing built in: the area/entropy channel turns out to be CONSISTENT-but-not-
CONSTRAINING (S_Barrow is monotonic in area, so the 2nd law holds for every Δ≥0); the
discriminating observable is PBH evaporation.  Run: python run_crosshorizon_data.py
"""
import numpy as np
from sede import barrow_bh as bh

# Real GWTC binary-black-hole mergers: (name, m1, m2, M_final, a*_final) [M_sun], LVK/GWTC.
GWTC = [
    ("GW150914", 35.6, 30.6, 63.1, 0.69),
    ("GW151226", 13.7,  7.7, 20.5, 0.74),
    ("GW170104", 31.0, 20.1, 49.1, 0.66),
    ("GW170814", 30.7, 25.3, 53.4, 0.72),
    ("GW170729", 50.2, 34.0, 79.5, 0.81),
    ("GW190521", 85.0, 66.0, 142.0, 0.72),
    ("GW190412", 30.1,  8.3, 38.0, 0.67),
]
MPL_G = 2.176e-5          # Planck mass in grams
MSUN_G = 1.989e33         # solar mass in grams

if __name__ == "__main__":
    print("="*74); print("TIER 3 (DATA) — is the cosmological Δ=1 consistent with REAL black-hole horizons?"); print("="*74)

    # (1) GWTC catalog: area theorem + Barrow GSL, event by event
    print("\n[1] GWTC binary-black-hole catalog (area theorem A_f≥A_1+A_2, Barrow GSL ΔS_B≥0):")
    print(f"    {'event':>10s} {'ΔA/(A1+A2)':>11s} {'area thm':>9s}  GSL holds for Δ=0/0.5/1")
    area_ok = gsl_all = 0
    for name, m1, m2, Mf, af in GWTC:
        A1, A2, Af, dA, frac = bh.area_theorem_margin(m1, m2, Mf, af)
        athm = dA > 0; area_ok += athm
        gsl = []
        for D in (0.0, 0.5, 1.0):
            Si = bh.bh_entropy(m1, 0.0, D) + bh.bh_entropy(m2, 0.0, D)
            Sf = bh.bh_entropy(Mf, af, D); gsl.append(Sf > Si)
        gsl_all += all(gsl)
        print(f"    {name:>10s} {100*frac:>10.0f}% {'✓' if athm else '✗':>9s}  {'/'.join('✓' if g else '✗' for g in gsl)}")
    N = len(GWTC)
    print(f"\n    Area theorem holds: {area_ok}/{N};  Barrow GSL holds (Δ=0,0.5,1): {gsl_all}/{N} for EVERY Δ.")
    print("    => the catalog is CONSISTENT with Δ=1 — but ALSO with Δ=0. The entropy/area channel")
    print("       does NOT constrain Δ_BH: S_Barrow∝A^{1+Δ/2} is monotonic in A, so the 2nd law is")
    print("       satisfied for any Δ≥0 whenever the (Δ-independent) area theorem holds. Consistent,")
    print("       not constraining — exactly why ringdown areas alone can't measure horizon entropy.")

    # (2) the discriminating channel: PBH evaporation
    print("\n[2] Primordial-black-hole evaporation — the one place Δ IS observable:")
    M_star_g = 5.1e14                      # Hawking: PBH mass evaporating now (~10^15 g)
    M_star_sun = M_star_g / MSUN_G
    Tr = bh.bh_temperature_ratio(M_star_sun, 1.0)        # T_Barrow/T_Hawking at Δ=1
    tr = bh.evaporation_time_ratio_log10(M_star_sun, 1.0)
    print(f"    Standard (Δ=0): M* ≈ {M_star_g:.1e} g evaporates within a Hubble time (γ-ray signature).")
    print(f"    Barrow (Δ=1):   T_B/T_H ≈ {Tr:.1e}  ⟹  evaporation time ×10^{tr:.0f}  -> NEVER evaporates.")
    print("    => Δ=1 predicts NO PBH evaporation in the observable mass window. The current ABSENCE of")
    print("       a PBH evaporation γ-ray signal is CONSISTENT with Δ=1; a future DETECTION (standard")
    print("       Hawking spectrum from a ~10^15 g PBH) would favour Δ_BH≈0 and FALSIFY universal Δ=1.")

    # verdict
    print("\n" + "="*74)
    print("VERDICT (cross-horizon, real data):")
    print(f"  • GWTC catalog: {gsl_all}/{N} events satisfy the Barrow GSL for Δ=1 — CONSISTENT, but the")
    print("    area/entropy channel cannot DISTINGUISH Δ=1 from Δ=0 (monotonicity). Not yet a constraint.")
    print("  • PBH evaporation: Δ=1 forbids it; the observed absence is consistent; a detection falsifies.")
    print("  • Decisive future probes: PBH/Hawking-radiation searches, and ringdown 'area-quantization'")
    print("    spectroscopy with ET/CE. Until then the cosmological Δ=1 is UNREFUTED but UNCONFIRMED on")
    print("    black-hole horizons — the cross-horizon test remains the open, genuine decider.")
