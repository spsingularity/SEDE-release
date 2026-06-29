# Data

The observational datasets used by the analysis are **not redistributed** in
this repository. Fetch each one from its original source and place it in `data/`
with the filename shown below; `sede/data_loader.py` reads local files only (no
auto-download). Cite the original sources in any derived work.

| File(s) | Probe | Source |
|---|---|---|
| `desi_2024_bao_*`, `desi_dr2_desi_gaussian_bao_*` | DESI BAO (DR2) | DESI Collaboration, DESI DR2 BAO |
| `desi_dr2_w0wa_{desy5,pantheonplus,union3}.*` | DESI DR2 (w0,wa) compressed posteriors | DESI Collaboration |
| `Pantheon+SH0ES.dat`, `Pantheon+SH0ES_STAT+SYS.cov` | Type Ia SNe (Pantheon+ & SH0ES) | Scolnic et al. 2022; Brout et al. 2022 |
| `DES5YR_HD.csv`, `DES_Dovekie_HD.csv`, `DES_STAT+SYS.npz` | DES 5-yr SNe Hubble diagram | DES Collaboration 2024 |
| `union3.npz`, `union3_mu.fits` | Union3 SNe | Rubin et al. 2023 |
| `sdss_DR16_BAOplus_{LRG,QSO}_FSBAO_DMDHfs8*` | eBOSS DR16 full-shape BAO + fσ8 | eBOSS Collaboration 2021 |
| `moresco_Hz*.dat`, `moresco_cov_*` | Cosmic chronometers H(z) | Moresco et al. (compilation) |
| `coma_sdss.csv` | Coma cluster galaxy kinematics (SDSS) | SDSS (cross-field Tsallis test) |

## Where to get them

Re-fetch the originals from the references above and drop them into `data/`
with the exact filenames in the table. The eBOSS DR16 files mirror those in the
public `CobayaSampler/bao_data` repository; DESI and Pantheon+/DES/Union3
products are on their respective collaboration data pages. Each dataset is
distributed under its own terms — consult the originating collaboration.
