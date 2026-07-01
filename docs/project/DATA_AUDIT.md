# Data audit (2026-06-23)

Full sweep of every downloaded / hardcoded dataset for discrepancies, duplication,
and anomalies. Automated checks: NaN/inf, duplicate rows, covariance symmetry +
positive-definiteness + condition number, and value-range sanity.

## Summary

| dataset | N | status | notes |
|---|---:|---|---|
| DESI DR2 BAO | 13 | ✅ clean | cov PD (cond 1.2e2); z=2.33 is DH,DM-ordered but correctly labelled |
| Pantheon+ | 1580 | ✅ clean | 94 duplicate-CID = standard multi-survey light curves (handled by cov) |
| Moresco CC | 15 | ⚠️ note | values standard; covariance is a 3% systematic APPROXIMATION |
| **fσ8 RSD** | 16 | ✅ **FIXED** | was 18-pt with bad over-tight entries; now standard Gold-2018 |
| Planck (compressed) | — | ⚠️ note | old (H0,Ω_m,R) path; superseded by CAMB (R,l_A) in the joint |
| SH0ES | — | ✅ **FIXED** | was 73.0±1.0; unified to Riess+2022 73.04±1.04 |
| DES-Dovekie | 1820 | ✅ clean | cross-validation only (not combined with Pantheon+ → no double-count) |
| Union3 | 22 | ✅ clean | cross-validation only |

## Findings in detail

### Fixed this session
1. **fσ8 — erroneous over-tight entries (the big one).** The old 18-point array had
   z=0.85 with σ=0.035 (real eBOSS-ELG error ≈0.095), and non-standard z=1.05 and
   z=2.33 points. Three points dominated χ² (17.5 of 30.9) and added a SPURIOUS
   ~+2.4 to SEDE-H's disfavouring (the fσ8 Δ(SEDE-H−ΛCDM) flips from +2.06 to −0.31
   with clean data). Replaced with the standard Gold-2018 16-point compilation.
   Effect: best-fit joint Δχ² +6.46 → +3.60.
2. **SH0ES prior inconsistency.** `data_loader` had (73.0, 1.0) while the CAMB
   drivers used (73.04, 1.04). Unified to Riess+2022 (73.04 ± 1.04).

### Benign (no action — verified correct)
- **Pantheon+ "duplicates":** 1580 rows / 1466 unique CIDs. 94 CIDs repeat — these
  are the SAME SN observed by multiple surveys (e.g. 1999ac twice, m_b_corr 13.665 /
  13.714), the standard Pantheon+ structure that the full STAT+SYS covariance is
  built to handle. The 181 same-z pairs are 109 same-CID (multi-survey) + 72 distinct
  SNe sharing a Hubble-frame z. Not double-counting.
- **DESI z=2.33 ordering:** the file lists DH/rd then DM/rd at the Lyα point (reverse
  of the other tracers), but each row's type label matches its value and covariance
  row, so it is parsed correctly.
- **DES-Dovekie & Union3:** load cleanly; used only as cross-validation alternatives
  to Pantheon+, never combined with it, so there is no SN double-counting in the joint.

### Notes / methodology (not errors)
- **Moresco covariance** is the diagonal stat error + a 3% correlated-systematic
  floor — an APPROXIMATION to the full Moresco+2020 covariance (which has
  scale-dependent SPS/IMF correlations). This slightly under-weights correlations;
  impact on the joint is ≲1 χ². The 15-point compilation also mixes Moresco CC with
  older Simon+2005/Stern+2010 points (standard but heterogeneous).
- **Two CMB treatments coexist:** the legacy compressed (H0, Ω_m, R) Gaussian priors
  (`data_loader.load_planck`, used by the old `mcmc_pipeline`) and the CAMB-exact
  (R, l_A) used by `run_camb_joint`/`run_camb_mcmc`. The CAMB path is canonical for
  all headline (Phase 1/2) results; the compressed path remains for the fast
  tabulated pipeline.
- Planck z_star = 1089.9 (legacy) vs CAMB z_star = 1090.0 — negligible.

## Verdict
After fixing the fσ8 compilation and the SH0ES value, no remaining data errors. The
fσ8 fix is the only one that moves a result (best-fit Δχ² +6.46 → +3.60; marginalised
re-run in progress). All covariances are symmetric and positive-definite.

## eBOSS DR16 full-shape (pre-DESI BAO+RSD) — added for the blind-DESI out-of-sample test
Genuine pre-DESI data: SDSS BOSS DR12 + eBOSS DR16 "BAOplus" consensus full-shape
(D_M/r_d, D_H/r_d, fσ8) for LRG (z=0.38, 0.51, 0.698) and QSO (z=1.48). Loader:
`data_loader.load_eboss_dr16_fs()` (12 measurements; block-diagonal LRG 9×9 ⊕ QSO 3×3 cov).
Used by `run_xval_oos.py` (blind-DESI: train pre-DESI → predict DESI DR2, Δχ²=−8.15 for SEDE).

Public source (CobayaSampler/bao_data, Alam et al. 2021). Re-fetch with:
```
BASE=https://raw.githubusercontent.com/CobayaSampler/bao_data/master
for f in sdss_DR16_BAOplus_LRG_FSBAO_DMDHfs8.dat sdss_DR16_BAOplus_LRG_FSBAO_DMDHfs8_covtot.txt \
         sdss_DR16_BAOplus_QSO_FSBAO_DMDHfs8.dat sdss_DR16_BAOplus_QSO_FSBAO_DMDHfs8_covtot.txt; do
  curl -s "$BASE/$f" -o "data/$f"; done
```
(data/ is gitignored — files stay local; this records the exact public provenance.)

## DESI DR2 w0waCDM (w0,wa) covariance — official chains

Test 3 of `run_delta_indirect.py` uses the **official DESI DR2 cosmology chains** (not the
published-marginal approximation). We store only the small getdist outputs (`chain.covmat`,
`chain.margestats`) per SN combination, parsed by `data_loader.load_desi_dr2_w0wa()`:

```
BASE=https://data.desi.lbl.gov/public/papers/y3/bao-cosmo-params/cobaya/base_w_wa
# dirs: desi-bao-all_{pantheonplus,desy5sn,union3}_planck2018-lowl-TT-clik_planck2018-lowl-EE-clik\
#       _planck-NPIPE-highl-CamSpec-TTTEEE_planck-act-dr6-lensing
for sn in pantheonplus desy5 union3; do  # (desy5 dir is 'desy5sn')
  curl -sL "$BASE/<dir>/chain.covmat"     -o data/desi_dr2_w0wa_${sn}.covmat
  curl -sL "$BASE/<dir>/chain.margestats" -o data/desi_dr2_w0wa_${sn}.margestats
done
```
Extracted (w0, wa, ρ): Pantheon+ (-0.838, -0.620, -0.895); DESY5 (-0.752, -0.860, -0.905);
Union3 (-0.667, -1.090, -0.932). The full Zenodo archive (DOI 10.5281/zenodo.16644576) is 1.3 GB;
we pulled only the two getdist files per chain (~5 KB each). Loader falls back to published marginals
if the files are absent. (data/ is gitignored — this records the exact public provenance.)
