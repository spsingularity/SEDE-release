#!/usr/bin/env python3
"""
APPLES-TO-APPLES marginalised CAMB-in-the-loop MCMC: ΛCDM vs SEDE (Barrow Δ=1) vs GOUGH IDE.

Gough's Information Dark Energy (ρ_DE ∝ SMD(a)^0.5, baryonic stellar-mass-density driver) is run
through the SAME pipeline as SEDE — its background EOS w(a) is derived from ρ_DE ∝ SMD^0.5 and fed to
CAMB via DarkEnergyPPF.set_w_a_table, exactly as SEDE's w(a) is. Same data (DESI DR2 BAO + Pantheon+
+ cosmic chronometers + compressed CMB R,ℓ_A + fσ8 + ω_b/SH0ES priors), same 5 params
(Ω_m, H0, ω_b, M_B, σ8), same DIC. This makes the SEDE-vs-Gough comparison fully fair.

Run with miniconda python (has camb):  /opt/miniconda3/bin/python3 run_gough_mcmc.py [--steps N --burn N]
"""
import argparse, time
import multiprocessing as mp
import numpy as np
import camb
from camb import dark_energy
import run_lambda_verify as V        # SEDE chi2 + data + consts (R_PL, LA_PL, CMB_CINV, priors)
import run_barrow_mcmc as B          # B.lcdm_chi2

fr = V.fr; C = V.C; data = V.data
GAMMA_TH, LAM = 1.4964, 0.5


# ── Gough: ρ_DE ∝ SMD(z)^0.5, SMD = ∫_z^∞ SFRD/((1+z')H) dz' (Madau-Dickinson, fiducial bg) ──
def _smd_curve():
    zh = np.linspace(0, 30, 6000); Om = 0.31
    El = np.sqrt(Om * (1 + zh)**3 + 9e-5 * (1 + zh)**4 + (1 - Om))
    sfrd = 0.015 * (1 + zh)**2.7 / (1 + ((1 + zh) / 2.9)**5.6)
    integ = sfrd / (El * (1 + zh))
    cum = np.cumsum((integ[::-1] * np.gradient(zh)[::-1]))[::-1]      # ∫_z^∞
    return zh, cum

_ZH, _CUM = _smd_curve()


def w_of_a_gough(n=250):
    z = np.linspace(0, 8, n)
    smd = np.interp(z, _ZH, _CUM)
    rho = np.maximum((smd / _CUM[0])**0.5, 1e-8)                      # ρ_DE ∝ SMD^0.5, ρ(0)=1
    w = -1 + (1 / 3.) * (1 + z) * np.gradient(np.log(rho), z)
    a = 1 / (1 + z); i = np.argsort(a)
    return a[i], w[i]


def gough_chi2(Om, H0, ombh2, MB, s8):
    h = H0 / 100.
    pa = camb.CAMBparams(); pa.set_cosmology(H0=H0, ombh2=ombh2, omch2=Om * h**2 - ombh2, mnu=0.06)
    a, w = w_of_a_gough(); de = dark_energy.DarkEnergyPPF(); de.set_w_a_table(a, w); pa.DarkEnergy = de
    pa.WantCls = False; pa.WantTransfer = False; bg = camb.get_background(pa)
    d = bg.get_derived_params(); rd = d['rdrag']; zs = d['zstar']; rs = d['rstar']; c = 0.0
    z, t, m, icov = data['desi']
    pred = np.array([bg.comoving_radial_distance(zz) / rd if tp == 'DM/rd'
                     else (C / bg.hubble_parameter(zz)) / rd if tp == 'DH/rd'
                     else ((zz * bg.comoving_radial_distance(zz)**2 * (C / bg.hubble_parameter(zz)))**(1 / 3.)) / rd
                     for zz, tp in zip(z, t)])
    dd = m - pred; c += float(dd @ icov @ dd)
    zp, mu, chol = data['pan']
    from scipy.linalg import cho_solve
    dmu = mu - (5 * np.log10((1 + zp) * bg.comoving_radial_distance(zp)) + 25 + MB)
    c += float(dmu @ cho_solve(chol, dmu))
    DMz = bg.comoving_radial_distance(zs); R = np.sqrt(Om) * H0 * DMz / C; lA = np.pi * DMz / rs
    v = np.array([R - V.R_PL, lA - V.LA_PL]); c += float(v @ V.CMB_CINV @ v)
    c += ((ombh2 - V.OMBH2_PRIOR[0]) / V.OMBH2_PRIOR[1])**2 + ((H0 - V.SHOES[0]) / V.SHOES[1])**2
    zc, H, icc = data['cc']; dH = H - bg.hubble_parameter(zc); c += float(dH @ icc @ dH)
    zf, fo, fe = data['fs8']
    Dd, fd = fr.compute_growth_model(zf, Om, lambda zz: bg.hubble_parameter(np.atleast_1d(zz)) / H0)
    c += float(np.sum(((fo - fd * s8 * Dd) / fe)**2))
    return c


def _inbox(th):
    Om, H0, ob, MB, s8 = th
    return (0.20 < Om < 0.45 and 60 < H0 < 78 and 0.0195 < ob < 0.025
            and -20.5 < MB < -18.5 and 0.62 < s8 < 0.90)


def logp_lcdm(t):
    if not _inbox(t): return -np.inf
    try: return -0.5 * B.lcdm_chi2(*t)
    except Exception: return -np.inf

def logp_sede(t):
    if not _inbox(t): return -np.inf
    try: return -0.5 * V.chi2(t[0], t[1], t[2], t[3], t[4], GAMMA_TH, LAM)[0]
    except Exception: return -np.inf

def logp_gough(t):
    if not _inbox(t): return -np.inf
    try: return -0.5 * gough_chi2(*t)
    except Exception: return -np.inf


def run(name, lp, args):
    import emcee
    theta0 = np.array([0.305, 68.4, 0.02237, -19.40, 0.78])
    scales = np.array([0.008, 0.4, 0.0003, 0.02, 0.015])
    rng = np.random.default_rng(11)
    p0 = theta0 + scales * rng.standard_normal((args.walkers, 5))
    pool = mp.Pool(args.workers) if args.workers > 1 else None
    s = emcee.EnsembleSampler(args.walkers, 5, lp, pool=pool)
    t0 = time.time(); print(f"[{name}] burn {args.burn}...", flush=True)
    st = s.run_mcmc(p0, args.burn, progress=False); s.reset()
    print(f"[{name}] prod {args.steps}...", flush=True)
    s.run_mcmc(st, args.steps, progress=False)
    if pool: pool.close()
    chi2s = -2 * s.get_log_prob(flat=True)
    cmin, cmean = float(chi2s.min()), float(chi2s.mean()); pD = cmean - cmin; DIC = cmean + pD
    med = np.median(s.get_chain(flat=True), axis=0)
    print(f"[{name}] {time.time()-t0:.0f}s acc={np.mean(s.acceptance_fraction):.2f} "
          f"Om={med[0]:.3f} H0={med[1]:.2f} s8={med[4]:.3f} | χ²min={cmin:.2f} <χ²>={cmean:.2f} p_D={pD:.1f} DIC={DIC:.2f}",
          flush=True)
    return dict(cmin=cmin, cmean=cmean, DIC=DIC)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--steps', type=int, default=250); ap.add_argument('--burn', type=int, default=120)
    ap.add_argument('--walkers', type=int, default=24); ap.add_argument('--workers', type=int, default=8)
    args = ap.parse_args()
    try: mp.set_start_method('fork')
    except RuntimeError: pass
    print("=" * 72); print("ΛCDM vs SEDE(Barrow Δ=1) vs GOUGH IDE — apples-to-apples marginalised CAMB MCMC"); print("=" * 72)
    rL = run('LCDM ', logp_lcdm, args)
    rS = run('SEDE ', logp_sede, args)
    rG = run('GOUGH', logp_gough, args)
    print("\n" + "=" * 72 + "\nMARGINALISED COMPARISON (ΔDIC vs ΛCDM; negative = preferred over ΛCDM)\n" + "=" * 72)
    for nm, r in [('SEDE ', rS), ('GOUGH', rG)]:
        print(f"  {nm}: χ²min Δ={r['cmin']-rL['cmin']:+.2f}  <χ²> Δ={r['cmean']-rL['cmean']:+.2f}  "
              f"ΔDIC={r['DIC']-rL['DIC']:+.2f}")
    print(f"  (ΛCDM: χ²min={rL['cmin']:.2f}, DIC={rL['DIC']:.2f})")


if __name__ == "__main__":
    main()
