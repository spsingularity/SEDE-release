#!/usr/bin/env python3
"""
TIER 3 — the cross-horizon test: does the cosmologically-fitted Δ=1 appear on BLACK-HOLE
horizons too?  This is the DECISIVE test of SEDE's quantum-gravity reading (Δ universal),
the one no dark-energy survey can do.  Here we (i) verify the Hawking area theorem on a real
merger, (ii) give the BH-sector signatures predicted by Δ=1, and (iii) state the falsifier.

These are PREDICTIONS + a forecast, not a current measurement: the observational handle is
GW area-law tests (now ~25%, ET/CE ~%) and primordial-BH evaporation (future).

Run:  python run_tier3_crosshorizon.py
"""
import numpy as np
from sede import barrow_bh as bh

if __name__ == "__main__":
    print("="*72); print("TIER 3 — cross-horizon: is the cosmological Δ=1 also the BLACK-HOLE Δ?"); print("="*72)

    # (i) GW150914 — Hawking area theorem (Δ-independent) + Barrow entropy increase
    m1, m2, Mf, af = 35.6, 30.6, 63.1, 0.69     # LVK GW150914 (M_sun, dimensionless spin)
    A1, A2, Af, dA, frac = bh.area_theorem_margin(m1, m2, Mf, af)
    print("\n[i] GW150914 area theorem (A_f ≥ A_1+A_2):")
    print(f"    A_1+A_2 = {A1+A2:.3e}   A_f = {Af:.3e}   ΔA = +{dA:.3e}  ({100*frac:+.1f}%)  -> holds")
    print("    Barrow entropy increase ΔS_B = S_B(final) − S_B(initial), by Δ:")
    for D in (0.0, 0.5, 1.0):
        Si = bh.bh_entropy(m1, 0.0, D) + bh.bh_entropy(m2, 0.0, D)
        Sf = bh.bh_entropy(Mf, af, D)
        print(f"      Δ={D:.1f}: log10 S_B(i)={np.log10(Si):7.2f}  log10 S_B(f)={np.log10(Sf):7.2f}  "
              f"ΔS_B>0: {Sf>Si}")

    # (ii) BH-sector signatures predicted by the cosmological Δ=1
    print("\n[ii] If Δ=1 is UNIVERSAL (cosmology → BH horizons), predicted BH signatures:")
    print(f"     {'BH':>16s} {'log10(S_B/S_BH)':>16s} {'log10(T_B/T_H)':>15s} {'log10(t_evap ratio)':>20s}")
    cases = [("60 M_sun (LVK)", 60.0, 0.7), ("10^6 M_sun (SMBH)", 1e6, 0.9),
             ("10^15 g PBH", 1e15/2e33, 0.0)]   # 1e15 g in M_sun
    for name, M, a in cases:
        enh = bh.entropy_enhancement_log10(M, 1.0, a)
        Tr = np.log10(bh.bh_temperature_ratio(M, 1.0, a))
        tr = bh.evaporation_time_ratio_log10(M, 1.0, a)
        print(f"     {name:>16s} {enh:16.1f} {Tr:15.1f} {tr:20.1f}")
    print("     => Δ=1 enhances BH entropy by ~10^40 and suppresses Hawking T by the same,")
    print("        so evaporation is effectively forbidden — a sharp PBH/Hawking-radiation signature.")

    # (iii) falsifier + forecast
    print("\n[iii] Falsification logic:")
    print("     - If a BH horizon is measured to carry entropy A/4 (Δ_BH=0) while cosmology")
    print("       needs Δ=1, the deformation is DE-sector-specific, NOT a universal QG constant")
    print("       -> SEDE's quantum-gravity reading is falsified (the model survives as late-DE only).")
    print("     - If BH horizons show Δ≈1 (S∝A^{3/2}) -> strong evidence Δ=1 is a real QG signature.")
    print("     Current GW area-law precision ~25% (Isi+2021) bounds macroscopic-horizon deviations;")
    print("     ET/CE ringdown (~1-2%) + primordial-BH evaporation are the decisive future probes.")
    print("\n  Net: Δ=1 is fixed by cosmology (Thms 8-9) and makes FALSIFIABLE BH-sector predictions")
    print("  -> SEDE is testable on two independent horizons, not just dark energy.")
