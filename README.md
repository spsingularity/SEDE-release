# SEDE — Structural Entropy Dark Energy

A **Λ-free** model of dark energy derived from horizon thermodynamics: dark energy is
the thermodynamic conjugate of gravitational entropy at the cosmic apparent horizon,
modulated by structure formation.

```
ρ_DE = T_AH · s_grav ,   s_grav = s_AH · f_sat(z) ,   f_sat ≤ 1 (Bousso bound).
```

The working model (**SEDE-H**) keeps the dynamical-horizon structure and the Barrow
entropy of the horizon, giving

```
ρ_DE ∝ H^{2λ} f_sat(z) ,   λ = 1 − Δ/2   (Barrow deformation Δ ∈ [0,1], Theorem 8).
```

DE vanishes at early times (Ω_DE(z≈1090) ≈ 0 — no cosmological constant), and the
magnitude ρ_DE ~ M_P²H₀² is small because the universe is large (no 10¹²⁰ fine-tuning).

---

## Headline result (current)

**Parameter-free Barrow SEDE-H is preferred over ΛCDM** in the rigorous, marginalised,
CAMB-in-the-loop analysis — at the **same 5 parameters as ΛCDM** (Ω_m, H₀, ω_b, M_B, σ8),
with **no free dark-energy parameter**:

| metric | Barrow SEDE-H | ΛCDM | Δ |
|---|---:|---:|---:|
| χ²_min | 1437.88 | 1440.86 | −2.98 |
| **DIC** | 1447.80 | 1450.66 | **−2.86 (Barrow preferred)** |
| effective params p_D | 5.0 | 4.9 | ≈ equal |

The **canonical SEDE-H is the bare Barrow holographic model** ρ_DE ∝ H^{2−Δ} f_sat with
**λ = 1 − Δ/2 = 0.5** — the *unique parameter-free* background carrying the derived Δ=1 (a
no-go theorem, §W, shows the alternative — adding the Cai–Kim temperature factor (1−ε/2) —
is not parameter-free at λ=0.5). Every input is derived:
- **λ = 0.5** from the Barrow horizon deformation **Δ = 1** (λ=1−Δ/2, Theorem 8) — and Δ=1
  is itself **derived** (Theorem 9): the unique geometric (fractal dimension d_H=2+Δ ≤ 3,
  space-filling) and 2nd-law (entropy maximal) extremum; the data confirm Δ ≈ 1.0–1.1,
- **γ ≈ 1.50** from the Sheth–Tormen halo mass function,
- the **EOS signature is w(z) crossing −1 with w_a < 0** — w0 ≈ −1.0 today, dipping to
  ≈ −1.08 by z≈2 — i.e. DESI DR2's evolving-DE (thawing/crossing) quadrant, CPL-compatible.

> **Note on w₀ = −0.85.** The closed-form structural value w₀ = (4Ω_m/3 − 1)/(1 − Ω_m) = −0.85
> is a property of the *dynamical-horizon cousin* `E_SEDE_H` (the (1−ε/2) Cai–Kim factor,
> ε0 = 2Ω_m, Theorem 5D), which is a λ=1 background — **not** the canonical λ=0.5 model (whose
> ε0 ≈ 1.5Ω_m gives w0 ≈ −1 instead). The two backgrounds cannot be merged for free (§W). So
> −0.85 is a *structural estimate / cousin value*, while the canonical model's actual EOS is
> the crossing-(−1), w_a<0 form above. Both land in DESI's preferred region.

The preference is **robust to the full Planck likelihood** (plik_lite TT/TE/EE shape
changes it by <0.5, `run_plik_check.py`) and is **independently reproduced** by a sibling
implementation (ΔDIC ≈ −3.6).

**Is the preference real, or just model flexibility?** Three independent cross-validations
(§Z) say real: (i) a **false-preference calibration** — under ΛCDM-truth mocks SEDE does *not*
win (null mean +0.43, 500 mocks), and the real Δχ²=−2.96 is beaten by only 1/500 null draws
(**p≈0.002, ~3σ**, 95% UL 0.005), i.e. not a flexibility artifact; (ii) **out-of-sample
prediction** — SEDE trained on *pre-DESI*
eBOSS DR16 data predicts the held-out DESI BAO better than ΛCDM (**Δχ²=−8.15**), with low-z→high-z
extrapolation −9.90; (iii) **code-independence** — SEDE's geometry agrees between CAMB and CLASS
to 5×10⁻⁵.

**Honest caveat:** the preference is statistically *real* but *modest in size* — ΔDIC=−2.86 is
mild on the Jeffreys scale, and both models sit at χ²/dof≈3 against DESI's internal BAO tension
(SEDE is "much less wrong," not "right"). SEDE-H is **preferred and out-of-sample-confirmed**,
but not yet *established* (strong>5 / decisive>10). The decisive remaining test is empirical
(DESI DR3 / Euclid → σ(Δ)≈0.09), not computational.

**The biggest open question (§AA).** The evidence above favours SEDE *as a w(z) model*. Its
distinctive *mechanism* — f_sat literally sourced by structure growth — is a separate, stronger
claim, and it is **not yet confirmed**: the f_sat reconstructed from the expansion (geometry) and
from RSD growth (structure) **track at z<1 but diverge at z>1** (structure flattens, geometry
declines; marginally disfavoured at χ²≈8/6 with precise DESI errors, consistent within the large
RSD errors). Inconclusive and RSD-limited — it needs Euclid/LSST weak-lensing tomography. So the
*phenomenology* is supported; the *founding physical mechanism* is the honest frontier.

---

## Code

### `sede/` package
| module | purpose |
|---|---|
| `theory.py` | Theorems 1–13 (incl. **5D** EOS-gap closure, **6** smooth-DE perturbations, **7** H₀ scope, **8** Barrow λ, **9** Δ=1 derived, **10** QG-from-first-principles, **11** seam=BBN-Δ, **12** inflation↔DE brackets, **13** EOS=entropy production); f_sat, EOS, reheating |
| `friedmann.py` | backgrounds: `E_SEDE`, `E_SEDE_H` (dynamical), `E_SEDE_lambda` / `E_SEDE_barrow` (λ-family), growth ODEs, table builders |
| `gamma_computation.py` | γ from Sheth–Tormen / COLOSSUS; entropy-weight scan; running γ_eff |
| `halo_entropy.py` | coarse-grained NFW phase-space + binding entropy → the entropy weight p |
| `thermo.py` | dimensional (SI) horizon thermodynamics; derives f_sat(0)=1 |
| `cmb.py` | compressed CMB: sound horizon r_s, shift R, acoustic scale l_A |
| `camb_background.py` | CAMB-exact background (SEDE-H via tabulated w(a) PPF) |
| `perturbations.py` | smooth-DE growth, ISW source, S8 |
| `barrow_bh.py` | cross-horizon Barrow black-hole thermodynamics (Tier 3): S=(A/4)^{1+Δ/2}, T_B, area theorem |
| `lambda_family.py` | λ-family landscape diagnostic |
| `likelihood.py`, `data_loader.py`, `mcmc_pipeline.py` | legacy tabulated joint likelihood + emcee + data |
| `verification.py` | 51 analytical tests (51/51 PASS) |

### Drivers
**Canonical (CAMB-in-the-loop):**
- `run_barrow_mcmc.py` — **the headline run**: marginalised parameter-free Barrow SEDE-H vs ΛCDM
- `run_camb_joint.py` / `run_camb_mcmc.py` — CAMB-exact joint best-fit / marginalised MCMC
- `run_lambda_verify.py` — λ-family verification (Barrow λ=0.5 etc.)
- `run_plik_check.py` — full Planck plik_lite robustness check
- `run_camb_check.py` — CAMB TT-spectrum / θ* comparison
- `run_gamma_s8_check.py` — EOS-background table + γ_data/S8 (CAMB-pinned)

**Experiments to strengthen & decide SEDE (§X):**
- `run_tier1_data.py` — Tier 1 (real data): SN robustness across Pantheon+/DES-SN5YR/Union3, Lyα z=2.33, S8
- `run_tier2_forecast.py` — Tier 2 (forecast): Fisher σ(Δ) for DESI DR3 + Euclid
- `run_tier3_crosshorizon.py` — Tier 3 (cross-horizon): Barrow black-hole predictions + falsifier
- `run_crosshorizon_data.py` — Tier 3 (data): Δ=1 vs the GWTC BBH catalog + PBH-evaporation discriminator
- `run_e1_mechanism.py` — founding-claim test: geometry f_sat vs structure f_sat (the open frontier, §AA)

**QG / perturbation checks:**
- `run_combined_qg_tests.py` — combined-QG-theory junction tests (C-QG1 the scale–form seam, C-QG2 Δ=holographic-vs-volume-law, C-QG3 coincidence, C-QG4 T·S=E)
- `run_class_perturbations.py` — W9: full CLASS Boltzmann perturbations (fσ8/σ8/S8/lensing/ISW), validates smooth-DE

**Cross-validations (§Y):**
- `run_xval_consistency.py` — self-consistency: GSL-along-history, ρ_DE=T·s closed loop, BBN null, age
- `run_xval_loo.py` — leave-one-probe-out robustness (ΔΧ²<0 in every holdout)

---

## Repository layout

```
reproduce_all.py     one entry point — runs every stage, prints a PASS/FAIL summary
sede/                the model package (Friedmann/thermo/likelihood/perturbations/…)
scripts/             the analysis drivers (run_*.py, make_*.py); run from the repo root
docs/                PREREGISTRATION.md, DATA_AUDIT.md, DATA.md
paper/               manuscript source + LaTeX/arXiv build chain (paper/SEDE.pdf)
results/             reference outputs (some are inter-stage inputs)
predictions.json     sha256-locked pre-registration artifact (see docs/PREREGISTRATION.md)
data/                observational datasets — you supply these (gitignored; see docs/DATA.md)
output/, cache/      generated figures / tables (gitignored, regenerated on demand)
```

## Quick start

> **Environment:** verified on Python 3.12 with **numpy 2.2.6, scipy 1.13.1, camb 1.6.6**.
> `requirements.txt` pins `numpy>=2.0,<2.3` for a reason: the code uses `np.trapezoid`
> (needs numpy ≥ 2.0), while numpy ≥ 2.3 makes `float(size-1 array)` a hard error that
> breaks camb 1.6.6's `set_cosmology`. Use a fresh venv to keep these pins isolated.

```bash
pip install -r requirements.txt        # core deps (honors the numpy pin above)
pip install -e .                        # make `sede` importable from anywhere (optional)

# ONE entry point — reproduces every result and figure, in order, with a PASS/FAIL summary:
python reproduce_all.py            # tests, theorems, EOS/γ/S8, headline MCMC, figures (~25 min)
python reproduce_all.py --quick    # same, but short MCMC chains (smoke test, ~5 min)
python reproduce_all.py --fast     # skip the slow MCMC stage (everything else, ~5 min)
python reproduce_all.py --plik     # also run the full Planck plik_lite check (needs cobaya, below)
```

Or run any stage on its own (from the repo root, so relative `data/`/`output/` paths resolve):

```bash
python -c "from sede.verification import run_all_tests; run_all_tests()"   # 68/68 tests
python -c "from sede import theory; theory.print_theorems()"               # Theorems 1–13
python scripts/run_barrow_mcmc.py        # headline: parameter-free Barrow SEDE-H vs ΛCDM (~20 min)
python scripts/run_gamma_s8_check.py     # EOS-background table + γ_data/S8 (CAMB-pinned, ~5 min)
python scripts/run_tier1_data.py         # Tier 1: SN robustness (Pantheon+/DES5YR/Union3) + Lyα + S8
python scripts/run_tier2_forecast.py     # Tier 2: Fisher σ(Δ) forecast for DESI DR3 + Euclid
python scripts/run_tier3_crosshorizon.py # Tier 3: cross-horizon Barrow black-hole predictions
python scripts/make_all_figures.py       # all figures -> output/
python scripts/run_plik_check.py         # full Planck plik_lite robustness check (~1 min)
```

`reproduce_all.py --plik` / `scripts/run_plik_check.py` need the native Planck likelihoods installed once:
`cobaya-install planck_2018_highl_plik.TTTEEE_lite_native planck_2018_lowl.TT planck_2018_lowl.EE -p ./packages`

---

## Data & docs

- **`paper/SEDE.pdf`** — the paper itself (the authoritative description of the model and results).
- **`docs/DATA_AUDIT.md`** — audit of every dataset (DESI DR2, Pantheon+, Moresco, fσ8, etc.).
- **`docs/PREREGISTRATION.md`** — dated, committed falsifiable predictions (Δ=1, w₀≈−0.98, etc.) for DESI DR3 / Euclid / LVK; sha256-locks `predictions.json` (regenerate with `python scripts/gen_predictions.py`, verify with `shasum -a 256 predictions.json` from the repo root).
- **`docs/DATA.md`** — provenance and fetch instructions for the datasets (not redistributed here).
- `data/`, `output/`, and `cache/` are gitignored (you supply `data/`; the rest is regenerated).

## Status of the physics
The canonical SEDE-H (bare Barrow holographic, λ=0.5) reproduces DESI DR2's evolving-DE hint
with an EOS that **crosses −1 (w_a<0)** — w0≈−1 today, ≈−1.08 by z≈2 — sitting in DESI's
preferred thawing/crossing quadrant. (The closed-form w₀=−0.85 is the λ=1 cousin's structural
value, §W, not the canonical model's — both land in DESI's region.) It is perturbatively safe
(smooth DE, c_s²=1), preserves BBN/recombination, and at zero extra parameters is
mildly-to-moderately **preferred** over ΛCDM under the gold-standard CMB. It is H₀-agnostic by
construction (the Hubble tension is out of scope).

**On the Barrow-Δ / BBN question (adopted resolution, §S; Theorem 8):** SEDE-H is Barrow
*holographic dark energy* — ρ_DE ∝ H^{2−Δ} f_sat added to standard GR — **not** Barrow
*modified gravity*. The Saridakis BBN bound (Δ≲1.4×10⁻⁴) is for the modified-gravity case;
in SEDE the Barrow exponent multiplies only the f_sat-gated ρ_DE, so H(BBN) is the standard
matter+radiation value (V22) and the bound does not apply. The **constant Δ=1 is BBN-safe by
construction** — no running Δ needed. This commits SEDE's scope (the horizon entropy sources
only the DE sector), a stated, falsifiable stance. A z_eq-tied running Δ (SEDE_V2) reproduces
the identical observable model as a robustness cross-check.

It is a real, testable, Λ-free model — not yet established — awaiting DR3/Euclid.

## License

This repository is **dual-licensed** by component:

- **Code** (the `sede/` package, `scripts/`, `reproduce_all.py`, and all other source) — **MIT License** (see `LICENSE`).
- **Manuscript text and figures** (everything under `paper/`, and the figures in `output/`) — **Creative Commons Attribution 4.0 International (CC BY 4.0)**, https://creativecommons.org/licenses/by/4.0/ (see `paper/LICENSE`).

If you reuse the code, MIT applies; if you reuse text or figures, please attribute under CC BY 4.0.
