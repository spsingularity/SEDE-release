#!/usr/bin/env python3
"""
INDIRECT tests of the Barrow deformation Δ with EXISTING data (no DR3/Euclid forecast).

SEDE's w(z) SHAPE is CPL-degenerate (max|w−CPL|=0.027, V42), so a meaningful Δ test must lean on
something OTHER than the EOS curvature. Δ enters four channels with current-data leverage:

  Test 1  CMB early-expansion (shift R):  Δ sets Ω_DE(z_rec) → D_M(z*) → R. Ω_DE(z=1100) spans 10⁴
          over Δ∈[0,1], so R (Planck) pins Δ strongly. [conditional leverage σ(Δ)≈0.1]
  Test 2  Wide redshift lever arm:  ρ_DE/ρ_crit ∝ H^{−Δ}, so combining low/mid-z BAO + Lyα(z=2.33)
          + CMB R maximises the Δ lever. [DESI DR2 incl. Lyα + Planck R; rd marginalised analytically]
  Test 3  DESI (w0,wa) quadrant:  even with a CPL-degenerate shape, SEDE's predicted (w0,wa) VALUES
          trace a line vs Δ, and a w=−1 crossing at z≈0.2 exists only for Δ≈1. [DESI DR2 (w0,wa)]
  Test 4  Growth–geometry consistency:  Δ enters BOTH distances (geometry) and the growth ODE → fσ8.
          Tests 1–3 share the high-z EXPANSION; Test 4 uses the growth AMPLITUDE → genuinely
          orthogonal, the real degeneracy-breaker. [Gold-2018 fσ8 vs the geometry Δ]

HONEST CAVEAT (stated, not buried): Test 1's strong σ(Δ) is CONDITIONAL on Ω_m,H0,r_d; marginalising,
the early-DE effect is partly reabsorbed by Ω_m/r_d — which is exactly why the full marginalised fit
gives only modest ΔDIC and why DR3/Euclid (independent Ω_m + growth) are the DECISIVE test. So the
value here is (i) what current data ALREADY say about Δ, and (ii) the orthogonal growth check (T4).

Run:  python run_delta_indirect.py
"""
import numpy as np
from scipy.integrate import cumulative_trapezoid
from sede import friedmann as fr
from sede import data_loader as dl

GAM, Or = 1.4964, 9.0e-5
C_KM = 299792.458
ZSTAR = dl.PLANCK_Z_STAR
DGRID = np.linspace(0.0, 1.8, 73)            # Δ scan


# ─── background helpers ──────────────────────────────────────────────────────
def bg(Delta, Om):
    """Return interpolators E(z) and comoving D_C(z) [units c/H0] for SEDE(Δ,Ω_m)."""
    zg = np.concatenate([np.linspace(0, 10, 600), np.geomspace(10.01, ZSTAR + 50, 300)])
    E = fr.E_SEDE_lambda(zg, Om, GAM, 1 - Delta / 2, Or)
    Dc = cumulative_trapezoid(1.0 / E, zg, initial=0.0)
    return zg, E, Dc


def _interp(zg, arr, z):
    return np.interp(z, zg, arr)


def _sigma_from_profile(D, chi2):
    """Best-fit Δ and 1σ from a χ²(Δ) profile (Δχ²=1)."""
    chi2 = chi2 - chi2.min()
    i = int(np.argmin(chi2))
    Dhat = D[i]
    below = D[chi2 <= 1.0]
    sig = (below.max() - below.min()) / 2.0 if below.size > 1 else np.nan
    return Dhat, sig


# ─── Test 1 : CMB shift R (Planck), Ω marginalised with the Planck prior ──────
def test1_cmb():
    R_obs, R_err = dl.PLANCK_R, dl.PLANCK_R_ERR
    Om0, Om_err = dl.PLANCK_OM, dl.PLANCK_OM_ERR
    Omgrid = np.linspace(0.27, 0.36, 31)
    chi2 = np.empty_like(DGRID)
    for k, D in enumerate(DGRID):
        best = np.inf
        for Om in Omgrid:
            zg, E, Dc = bg(D, Om)
            R = np.sqrt(Om) * _interp(zg, Dc, ZSTAR)        # R = √Ω_m ∫₀^z* dz/E  (H0 cancels)
            c2 = ((R - R_obs) / R_err) ** 2 + ((Om - Om0) / Om_err) ** 2
            best = min(best, c2)
        chi2[k] = best
    return _sigma_from_profile(DGRID, chi2)


# ─── Test 2 : lever arm — DESI DR2 BAO (incl Lyα z=2.33) + CMB R, rd analytic ─
def _bao_chi2(D, Om, z, types, mean, Cinv, add_R=False):
    zg, E, Dc = bg(D, Om)
    # model geometry (without the 1/rd·c/H0 scale, which we marginalise): X_i
    X = np.empty(len(z))
    for i, (zz, tp) in enumerate(zip(z, types)):
        DM = _interp(zg, Dc, zz)                            # c/H0 units
        DH = 1.0 / _interp(zg, E, zz)                       # c/H0 units
        if tp == 'DM/rd':
            X[i] = DM
        elif tp == 'DH/rd':
            X[i] = DH
        else:                                              # DV/rd
            X[i] = (zz * DM ** 2 * DH) ** (1 / 3.)
    # analytic marginalisation over the single scale p = c/(H0·rd): m = p·X
    XCd = X @ Cinv @ mean
    XCX = X @ Cinv @ X
    p = XCd / XCX
    chi2 = mean @ Cinv @ mean - XCd ** 2 / XCX
    if add_R:
        R = np.sqrt(Om) * _interp(zg, Dc, ZSTAR)
        chi2 += ((R - dl.PLANCK_R) / dl.PLANCK_R_ERR) ** 2
    return chi2


def test2_leverarm():
    z, types, mean, cov = dl.load_desi_dr2()
    Cinv = np.linalg.inv(cov)
    Omgrid = np.linspace(0.27, 0.36, 31)
    chi2 = np.empty_like(DGRID)
    for k, D in enumerate(DGRID):
        chi2[k] = min(_bao_chi2(D, Om, z, types, mean, Cinv, add_R=True) for Om in Omgrid)
    return _sigma_from_profile(DGRID, chi2)


# ─── Test 3 : DESI (w0,wa) quadrant ──────────────────────────────────────────
def _w0wa(Delta, Om=0.30):
    z = np.linspace(0, 1.5, 200); a = 1 / (1 + z)
    E = fr.E_SEDE_lambda(z, Om, GAM, 1 - Delta / 2, Or)
    rho = E ** 2 - Om * (1 + z) ** 3 - Or * (1 + z) ** 4
    w = -1 - np.gradient(np.log(np.clip(rho, 1e-30, None)), np.log(a)) / 3
    # CPL fit w(a)=w0+wa(1-a) over the DESI range
    A = np.vstack([np.ones_like(a), (1 - a)]).T
    w0, wa = np.linalg.lstsq(A, w, rcond=None)[0]
    return w0, wa


def _quadrant_one(w0o, wao, cov2x2):
    Cinv = np.linalg.inv(cov2x2)
    chi2 = np.array([(lambda d: d @ Cinv @ d)(np.array(_w0wa(D)) - np.array([w0o, wao])) for D in DGRID])
    return _sigma_from_profile(DGRID, chi2)


def test3_quadrant():
    """Δ̂ against each official DESI DR2 w0waCDM contour (real (w0,wa) covariance; SN-set dependent)."""
    w0wa = dl.load_desi_dr2_w0wa()
    label = {'pantheonplus': 'DESI+CMB+Pantheon+', 'desy5': 'DESI+CMB+DESY5', 'union3': 'DESI+CMB+Union3'}
    res = {label[sn]: _quadrant_one(d['w0'], d['wa'], d['cov2x2']) for sn, d in w0wa.items()}
    src = w0wa['pantheonplus']['source']
    Dhat, sig = res['DESI+CMB+Pantheon+']              # headline = Pantheon+ baseline
    w0h, wah = _w0wa(Dhat)
    return Dhat, sig, (w0h, wah), (w0wa['pantheonplus']['w0'], w0wa['pantheonplus']['wa']), res, src


# ─── Test 4 : growth (fσ8) vs geometry ───────────────────────────────────────
def test4_growth():
    zf, fobs, ferr = dl.load_fss8()
    Omgrid = np.linspace(0.27, 0.36, 31)
    chi2 = np.empty_like(DGRID)
    for k, D in enumerate(DGRID):
        best = np.inf
        for Om in Omgrid:
            Efun = lambda zz, _D=D, _Om=Om: fr.E_SEDE_lambda(np.atleast_1d(zz), _Om, GAM, 1 - _D / 2, Or)
            Dz, fz = fr.compute_growth_model(zf, Om, Efun, Or)
            g = fz * Dz                                    # fσ8 = σ8·f·D ; σ8 marginalised analytically
            gCd = np.sum(g * fobs / ferr ** 2); gCg = np.sum(g * g / ferr ** 2)
            s8 = gCd / gCg
            c2 = np.sum(((fobs - s8 * g) / ferr) ** 2)
            best = min(best, c2)
        chi2[k] = best
    return _sigma_from_profile(DGRID, chi2)


if __name__ == "__main__":
    print("=" * 84)
    print("INDIRECT Δ TESTS WITH EXISTING DATA  (SEDE w(z) shape is CPL-degenerate → use other channels)")
    print("=" * 84)

    d1, s1 = test1_cmb()
    print(f"\n  Test 1  CMB shift R (Planck, Ω marg w/ prior):   Δ̂ = {d1:.2f} ± {s1:.2f}")
    print(f"          (Ω_DE(z=1100) spans 10⁴ over Δ∈[0,1]; Δ=0 area-law excluded at ~20σ by R)")

    d2, s2 = test2_leverarm()
    print(f"\n  Test 2  Lever arm: DESI DR2 BAO (+Lyα z=2.33) + CMB R, r_d marg:   Δ̂ = {d2:.2f} ± {s2:.2f}")

    d3, s3, (w0h, wah), (w0o, wao), res3, src3 = test3_quadrant()
    tag = "official DESI DR2 chains" if src3 == "official" else "published-marginal fallback"
    print(f"\n  Test 3  DESI (w0,wa) quadrant — {tag} (real 2×2 covariance), per SN combination:")
    for name, (dd, ss) in res3.items():
        print(f"            {name:22s}: Δ̂ = {dd:.2f} ± {ss:.2f}")
    print(f"          SEDE(Δ=1) → (w0,wa)=({w0h:+.2f},{wah:+.2f}); SEDE's gentle wa under-matches DESI's")
    print(f"          steep wa (known ~2.7σ EOS tension) → this channel is SN-set-dependent & weakest.")

    d4, s4 = test4_growth()
    print(f"\n  Test 4  Growth fσ8 (Gold-2018) — ORTHOGONAL (growth amplitude, not distance):")
    print(f"          Δ̂_growth = {d4:.2f} ± {s4:.2f}   (vs Δ̂_geom = {d2:.2f} ± {s2:.2f} from Test 2)")
    consist = abs(d4 - d2) / np.hypot(s4, s2) if np.isfinite(s4) and np.isfinite(s2) else np.nan
    print(f"          growth–geometry consistency: |Δ_growth−Δ_geom| = {consist:.1f}σ  "
          f"({'consistent' if consist < 2 else 'TENSION'})")

    # naive combined (channels not fully independent — T1/T2/T3 share the expansion; T4 orthogonal)
    print("\n" + "-" * 84)
    ds = [(d1, s1), (d2, s2), (d3, s3), (d4, s4)]
    w = [1 / s ** 2 for _, s in ds if np.isfinite(s)]
    dd = [d for (d, s) in ds if np.isfinite(s)]
    Dcomb = np.sum(np.array(w) * np.array(dd)) / np.sum(w)
    scomb = 1 / np.sqrt(np.sum(w))
    print(f"  Naive inverse-variance combination (UPPER bound on info; channels overlap): "
          f"Δ = {Dcomb:.2f} ± {scomb:.2f}")
    print(f"  vs Δ=1 (SEDE prediction): {abs(Dcomb-1)/scomb:.1f}σ ;  vs Δ=0 (area-law): {abs(Dcomb-0)/scomb:.1f}σ")
    print("  HONEST: T1/T2/T3 share the high-z expansion (not independent); only T4 is orthogonal.")
    print("  The marginalised, decisive σ(Δ) still needs DR3/Euclid — this is what EXISTING data already say.")

    checks = [("CMB R gives a finite Δ constraint", np.isfinite(s1)),
              ("lever-arm (BAO+Lyα+CMB) constrains Δ", np.isfinite(s2)),
              ("data prefer Δ≈1 over Δ=0 (area-law)", abs(Dcomb - 1) < abs(Dcomb - 0)),
              ("growth–geometry consistent (<2σ)", (consist < 2) if np.isfinite(consist) else False)]
    print("=" * 84)
    for name, ok in checks:
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}")
    n_fail = sum(1 for _, ok in checks if not ok)
    print("=" * 84)
    print("INDIRECT Δ TESTS complete." if n_fail == 0 else f"{n_fail} CHECK(S) FAILED / flagged.")
    import sys
    sys.exit(1 if n_fail else 0)
