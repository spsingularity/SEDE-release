#!/usr/bin/env python3
"""
Verlinde cross-check (§7): does the SAME volume-law normalisation that fixes SEDE's dark ENERGY
(ρ_DE ~ ρ_crit, via CKN/flatness, §8.2) also yield the dark-MATTER acceleration scale a₀ ~ cH₀
that Verlinde (2016) derives from the same volume-law de Sitter entropy?

This is an ORDER-OF-MAGNITUDE consistency check, NOT a derivation of a₀ and NOT an endorsement of
Verlinde's dark-matter proposal (which is contested at cluster/lensing scales and lacks a CMB/
covariant completion; §7). The point: both ρ_DE ~ ρ_crit and a₀ ~ cH₀ descend from the single
cosmological-horizon scale, so a volume-law dark sector should give both with consistent O(1) factors.

Distinctive SEDE input: for the canonical λ = 1−Δ/2 = 0.5 background, the late-time de Sitter
attractor has H_∞ = Ω_DE0·H₀ (E² = Ω_DE0·(E²)^{1/2} ⟹ E_∞ = Ω_DE0), versus ΛCDM's H_∞ = √Ω_Λ·H₀.
So SEDE predicts a slightly different de Sitter acceleration than ΛCDM.

Run:  python run_verlinde_crosscheck.py
"""
import numpy as np

# --- constants (SI) ---
c = 2.99792458e8                       # m/s
Mpc = 3.0856775814913673e22            # m
H0_kms = 67.4                          # km/s/Mpc (fiducial; result is insensitive at O(1))
H0 = H0_kms * 1e3 / Mpc                # s^-1
Om = 0.30
Or = 9e-5
Ode = 1.0 - Om - Or                    # Ω_DE0 ≈ 0.70

a0_obs = 1.2e-10                       # observed MOND/RAR scale (McGaugh–Lelli–Schombert 2016), m/s²
TWO_PI = 2.0 * np.pi                   # Verlinde's geometric factor a₀ ≈ cH/(2π)

# --- candidate de Sitter / horizon acceleration scales ---
a_H0    = c * H0                                   # Hubble acceleration today
a_DE    = c * np.sqrt(Ode) * H0                    # from current ρ_DE: H_DE = √Ω_DE0·H₀
a_inf   = c * Ode * H0                             # SEDE λ=0.5 attractor: H_∞ = Ω_DE0·H₀ (distinctive)


def show(label, a_raw):
    a_v = a_raw / TWO_PI
    print(f"  {label:<38s}  raw = {a_raw:.2e}  | ×1/2π = {a_v:.2e}  | (×1/2π)/a₀_obs = {a_v/a0_obs:.2f}")


if __name__ == "__main__":
    print("=" * 86)
    print("VERLINDE CROSS-CHECK (§7) — one volume-law scale for dark energy AND apparent dark matter?")
    print("=" * 86)
    print(f"\n  H₀ = {H0_kms} km/s/Mpc = {H0:.3e} s⁻¹ ;  Ω_DE0 = {Ode:.3f}")
    print(f"  cH₀ = {a_H0:.3e} m/s²    observed a₀ (RAR/SPARC) = {a0_obs:.2e} m/s²")
    print(f"  cH₀ / a₀_obs = {a_H0/a0_obs:.2f}   (the raw Hubble-acceleration coincidence)")

    print("\n  Candidate de Sitter acceleration scales (Verlinde: a₀ ≈ (de Sitter accel)/2π):")
    print(f"  {'scale':<38s}  {'value':>10s}     {'/2π':>10s}     ratio to a₀_obs")
    show("a = cH₀            (Hubble today)", a_H0)
    show("a = c√Ω_DE0·H₀     (from ρ_DE)", a_DE)
    show("a = cΩ_DE0·H₀      (SEDE λ=0.5 H_∞)", a_inf)

    best = a_DE / TWO_PI / a0_obs
    print("\n  VERDICT (honest, order-of-magnitude):")
    print(f"  • The horizon/de Sitter acceleration scale set by the SAME normalisation that gives")
    print(f"    ρ_DE ~ ρ_crit lands within an O(1) factor of the observed MOND scale a₀: with Verlinde's")
    print(f"    geometric 1/2π, the ρ_DE-based estimate is {best:.2f}×a₀_obs (i.e. ~{abs(1-best)*100:.0f}% off).")
    print(f"  • This is a CONSISTENCY check — both ρ_DE~ρ_crit and a₀~cH₀ descend from one horizon scale —")
    print(f"    NOT a derivation of a₀'s precise value, and NOT an adoption of Verlinde's dark-matter")
    print(f"    mechanism (its cluster/lensing failures and missing CMB/covariant completion are untouched).")
    print(f"  • Distinctive: SEDE's λ=0.5 attractor H_∞=Ω_DE0·H₀ predicts a de Sitter accel ~{a_inf/a_DE:.2f}×")
    print(f"    the ρ_DE-based one (ΛCDM uses √Ω_Λ instead of Ω_DE0) — a future discriminator if the")
    print(f"    volume-law dark-matter link is ever made quantitative.")

    ok = 0.3 < best < 3.0
    print(f"\n  O(1)-consistency check: {'PASS' if ok else 'FAIL'}  (0.3 < (1/2π)c√Ω_DE0 H₀ / a₀ < 3)")
    import sys
    sys.exit(0 if ok else 1)
