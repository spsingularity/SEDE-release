# E1 + E2 results (Tier-1 experiments) — for the manuscript edits

## E1 — Full (uncompressed) primary-CMB likelihood, MARGINALISED  [T1.1, the gate]

**Driver:** `run_full_cmb_mcmc.py --minimize`  ·  **engine:** cobaya 3.6.2 + camb 1.6.6
**Likelihood:** `planck_2018_highl_plik.TTTEEE_lite_native` + `planck_2018_lowl.TT` + `planck_2018_lowl.EE`
**Marginalised params (both models):** ombh2, omch2, cosmomc_theta, logA(→As), ns, tau, A_planck
**SEDE w(z):** canonical E_SEDE_volume (Barrow Δ=1, λ=0.5), γ=1.5 → CPL (w0,wa)=(−0.976,−0.141),
CPL-residual 0.020 (validates the paper's max|w−CPL|≈0.027 claim).

| quantity | ΛCDM | SEDE (γ=1.5) |
|---|---|---|
| χ²_min (TTTEEE+lowT+lowE) | 1003.55 | 1003.11 |
| best-fit Ω_m | 0.315 | 0.295 |
| best-fit H0 | 67.3 | 69.6 |

**Δχ²(SEDE−ΛCDM) = −0.43**, per channel: TTTEEE_lite **−0.37**, lowl.TT (low-ℓ ISW channel) −0.03,
lowl.EE −0.04. (With γ=1.0 instead: −0.81; both negative, both |Δχ²|<1.)

**ACT DR6 lensing at the marginalised best-fit** (not in the joint; evaluated at the primary best-fit
via `act_dr6_lenslike`): ΛCDM 14.10 → SEDE 13.79, **Δχ²_lens = −0.31** — favourable, consistent with the
standalone fiducial −0.32. The low-ℓ ISW-sensitive channel (lowl.TT) is flat (−0.03).

**Verdict: PASS.** SEDE does NOT inherit the standard-HDE early–late tension at the full
primary-CMB (marginalised) level. The (small) preference sits in the high-ℓ peak structure —
exactly the channel flagged as untested — and is *negative*, not the large *positive* Δχ² the HDE
tension would produce. Both minimizer restarts converged to identical minima (robust, not a local min).

**Robustness:** SEDE (w0,wa) recomputed at the best-fit Ω_m=0.295 → (−0.978,−0.139), shift Δw0=−0.0015
vs the fiducial used in the fit ⟹ fixing (w0,wa) during the marginalisation is sound.

**Honest scope (state in the paper):**
- plik_lite has foregrounds PRE-marginalised; the full non-lite plik (sampled foregrounds) and ACT DR6
  PRIMARY multifrequency are NOT installed → this is plan step (a). ACT DR6 *lensing* run separately
  (−0.32). The genuinely-new number is the MARGINALISED Δχ² on the primary peaks: −0.43.
- This is a profiled (minimize) Δχ²; a full MCMC ΔDIC (`--mcmc`) is the optional next step but p_D is
  shared (same 5 cosmo params), so ΔDIC ≈ Δχ² is expected.

**Cross-CMB consistency:** compressed (R,ℓ_A) headline ≈ −0.5 · plik θ*-proxy −0.50 · ACT lensing
−0.32 · full-primary marginalised −0.43 — all |Δχ²_CMB|<1, all favouring SEDE slightly. The wall isn't there.

---

## E2 — γ systematics  [T1.2]

**Driver:** `run_gamma_systematics.py`  ·  entropy weight p=5/3 FIXED (horizon reservoir, Result 4C);
halo-model choices varied.

| mass function (full range, EH98) | γ |
|---|---|
| Press–Schechter | 1.566 |
| **Sheth–Tormen (canonical)** | **1.496** |
| Tinker08 | 1.504 |
| Despali16 | 1.496 |
| Watson13 | 1.484 |

- **Robust to mass-function choice:** all full-range MFs give γ = 1.48–1.57.
- **Robust to transfer:** EH98 1.496 vs Sugiyama95 1.515.
- **Dominant systematic = mass-range low-cutoff:** narrow [1e11,1e15] pulls γ→1.14–1.22; massive
  [1e12,1e16] keeps γ≈1.5.
- **Full spread: γ ∈ [1.14, 1.60], median 1.496.**

**Claim for §3.1 (round-2 §7.2):** the reservoir (hence p=5/3) is forced by ρ_DE=T_AH·s_grav; the
numerical γ≈1.5 inherits standard, pre-specified halo-model systematics (~±0.2, dominated by the mass
range), NONE tuned to cosmological data.

JSON: `results/e1_minimize.json`, `results/e2_gamma_systematics.json`.

---

## Joint with full primary CMB folded in (T1.1 step 1, `run_joint_fullcmb.py`)

The genuine closure of CL's #1: the full primary CMB (plik_lite TTTEEE + lowℓ) folded **into the joint**
(BAO+SN+CC+fσ8+SH0ES + full CMB), replacing the compressed (R, ℓ_A) priors. Canonical λ=0.5 background
(E_SEDE_lambda, consistent with E1 — NOT the E_SEDE_H cousin); σ8 derived; Nelder–Mead profiled.

| | χ²_min | Om | H0 |
|---|---|---|---|
| ΛCDM | 2443.99 | 0.296 | 68.71 |
| SEDE (Barrow λ=0.5, γ=1.4964) | 2440.83 | 0.295 | 68.91 |

**Joint Δχ²(SEDE−ΛCDM) = −3.17** with full CMB in the loop — vs ≈ −2.9 with compressed CMB. The
preference **survives and marginally strengthens**; the headline is NOT a compressed-CMB artifact.
Written into §5.1. JSON: `results/joint_fullcmb.json`.
