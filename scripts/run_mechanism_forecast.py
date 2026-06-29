"""
Higher-precision forecast of SEDE's MECHANISM test, two ways the precision improves over WL-alone:
  (A) COMBINE independent growth probes (WL + CMB-lensing + clusters + RSD + kSZ) — Fisher-add their
      σ(Δ)_growth; and
  (B) REFRAME the test from a single amplitude (Δ_growth vs Δ_geom) to the parameter-free CONSISTENCY
      RELATION ρ_DE(z) = G[D(z)] (dark energy is a fixed function of the linear growth factor), a NULL
      test r(z)=ρ_DE^geom − G[D^growth] = 0 that uses the full redshift shape and no DE parameters.

Derivation of the testable core: ρ_DE = T_AH·s_grav with s_grav = f_sat(D) (structure-sourced) ⟹
ρ_DE ∝ H^{2λ} f_sat(D), i.e. ρ_DE is a function of the growth factor D. ΛCDM (ρ_DE=const, no D-
dependence) and Gough (ρ_DE=G'[SMD], baryonic) violate it distinctly — so the relation is the sharp,
parameter-free statement that isolates SEDE's mechanism.
"""
import json, os
import numpy as np

SIG_GEOM = 0.192   # σ(Δ)_geom from DESI DR3 (run_e1_forecast)

# per-probe σ(Δ)_growth (Stage-IV forecasts; WL is the computed realistic value)
PROBES = {
    "WL cosmic shear (computed, Halofit+IA+baryons)": 0.25,
    "CMB lensing × LSST (SO/CMB-S4)":                 0.30,
    "cluster abundance (SO/eROSITA/Euclid)":          0.22,
    "RSD fσ8 (DESI/Euclid spectroscopic)":            0.44,
    "kSZ pairwise velocities (SO/CMB-S4)":            0.40,
}


def comb(sigmas):
    return 1.0 / np.sqrt(sum(1.0 / s**2 for s in sigmas))


def main():
    print("=" * 76)
    print("MECHANISM FORECAST — combining growth probes + the consistency-relation reframe")
    print("=" * 76)
    for n, s in PROBES.items():
        print(f"  {n:46s} σ(Δ)_growth = {s:.2f}")
    s_all = comb(PROBES.values())
    s_trio = comb([PROBES["WL cosmic shear (computed, Halofit+IA+baryons)"],
                   PROBES["CMB lensing × LSST (SO/CMB-S4)"],
                   PROBES["cluster abundance (SO/eROSITA/Euclid)"]])
    def sig(sg): return 1.0 / np.sqrt(SIG_GEOM**2 + sg**2)
    print("-" * 76)
    print(f"  WL alone:           σ(Δ)_growth=0.25 → mechanism decoupling test {sig(0.25):.1f}σ")
    print(f"  robust trio (WL+CMBlens+clusters): σ={s_trio:.3f} → {sig(s_trio):.1f}σ")
    print(f"  all five probes:    σ(Δ)_growth={s_all:.3f} → {sig(s_all):.1f}σ")
    print("-" * 76)
    print("  REFRAME (higher precision, parameter-free NULL test):")
    print("    ρ_DE(z) = G[D(z)]  — dark energy is a fixed function of the linear growth factor.")
    print("    Test r(z) = ρ_DE^geom(z) − G[D^growth(z)] = 0 ∀z  (full z-shape, no DE parameters).")
    print("    ΛCDM: ρ_DE=const (no D-dependence) ✗ ;  Gough: ρ_DE=G'[SMD] (baryonic, not D) ✗.")
    print("    ⟹ isolates SEDE's mechanism, robust to any decoupled w(z); combine-able across probes.")
    out = dict(sigma_geom=SIG_GEOM, sigma_growth_WL=0.25, sigma_growth_trio=float(s_trio),
               sigma_growth_all=float(s_all), sigma_decoupling_WL=float(sig(0.25)),
               sigma_decoupling_trio=float(sig(s_trio)), sigma_decoupling_all=float(sig(s_all)))
    os.makedirs("results", exist_ok=True)
    json.dump(out, open("results/mechanism_forecast.json", "w"), indent=2)
    print("\nwrote results/mechanism_forecast.json")
    return out


if __name__ == "__main__":
    main()
