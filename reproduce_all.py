#!/usr/bin/env python3
"""
reproduce_all.py — one entry point that reproduces SEDE's results and figures.

Runs each stage as a subprocess (so a failure in one is isolated and reported),
in dependency order, then prints a PASS/FAIL summary with timings.

Stages
  1. tests    verification suite (26/26)                           ~seconds   [always]
  2. theorems print Theorems 1–9                                   ~seconds   [always]
  3. eos      EOS-background table + γ/S8 check (run_gamma_s8_check) ~minutes  [CAMB]
  4. tier1    real-data: SN robustness + Lyα + S8 (run_tier1_data)  ~2 min     [CAMB]
  5. tier2    Fisher forecast σ(Δ): DESI DR3+Euclid (run_tier2_…)   ~seconds
  6. tier3    cross-horizon Barrow-BH predictions (run_tier3_…)     ~seconds
  7. xval     self-consistency: GSL, closed-loop, BBN, age (run_xval_consistency) ~seconds
  8. loo      leave-one-probe-out robustness (run_xval_loo)         ~10 min    [CAMB]
  9. mcmc     headline marginalised ΔDIC (run_barrow_mcmc)          ~20 min    [CAMB+emcee]
 10. plik     full Planck plik_lite robustness (run_plik_check)     ~1 min     [cobaya; opt-in]
 11. figures  all figures -> output/ (make_all_figures)             ~minutes   [matplotlib]

Usage
  python reproduce_all.py                 # tests, theorems, eos, mcmc, figures (no plik)
  python reproduce_all.py --quick         # mcmc with short chains (smoke test, ~3 min)
  python reproduce_all.py --fast          # skip mcmc entirely (everything else)
  python reproduce_all.py --plik          # also run the plik_lite check (needs cobaya install)
  python reproduce_all.py --only tests,figures
  python reproduce_all.py --skip mcmc
"""
import argparse, os, subprocess, sys, time

ROOT = os.path.dirname(os.path.abspath(__file__))
PY = sys.executable
ALL = ["tests", "theorems", "eos", "tier1", "tier2", "tier3", "xval", "diag", "classx", "oos", "e1", "ch", "chr", "cqg", "fhist", "eosS", "pert", "v2x", "qg1", "qg2", "qg3", "qg4", "syk", "cset", "ebnd", "tnet", "structde", "goughmcmc", "litexp", "e1fore", "wlfish", "wlreal", "mech", "x2pt", "dind", "vol", "loo", "calib",
       "verl", "eft", "fap", "eg", "biref", "tloo", "mcmc", "actlens", "plik",
       "gammasys", "probedecomp", "infocrit", "dyntx", "refcheck", "eftlag", "bayesev", "deltaorth", "noshoes", "cmbfig", "e1fig", "fullcmb", "jointcmb", "figures"]


def run(cmd, env=None):
    """Stream a subprocess live; return (ok, seconds)."""
    t0 = time.time()
    e = dict(os.environ)
    e["PYTHONPATH"] = ROOT + (os.pathsep + e["PYTHONPATH"] if e.get("PYTHONPATH") else "")
    if env:
        e.update(env)
    r = subprocess.run(cmd, cwd=ROOT, env=e)
    return r.returncode == 0, time.time() - t0


def stage_cmd(name, args):
    if name == "tests":
        return [PY, "-c", "from sede.verification import run_all_tests; run_all_tests()"]
    if name == "theorems":
        return [PY, "-c", "from sede import theory; theory.print_theorems()"]
    if name == "eos":
        return [PY, "run_gamma_s8_check.py"]
    if name == "tier1":
        return [PY, "run_tier1_data.py"]        # real-data: SN robustness, Lyα, S8
    if name == "tier2":
        return [PY, "run_tier2_forecast.py"]    # Fisher: σ(Δ) for DESI DR3 + Euclid
    if name == "tier3":
        return [PY, "run_tier3_crosshorizon.py"]  # cross-horizon Barrow-BH predictions
    if name == "xval":
        return [PY, "run_xval_consistency.py"]    # GSL, closed-loop, BBN, age (self-consistency)
    if name == "diag":
        return [PY, "run_xval_diagnostics.py"]    # (w0,wa), statefinder, growth-index, ISW, first law
    if name == "classx":
        return [PY, "run_xval_class.py"]          # CLASS-vs-CAMB background (skips if classy absent)
    if name == "oos":
        return [PY, "run_xval_oos.py"]            # out-of-sample prediction (CV1 blind / CV2 split)
    if name == "e1":
        return [PY, "run_e1_mechanism.py"]        # founding-claim test: geometry vs structure f_sat
    if name == "ch":
        return [PY, "run_crosshorizon_data.py"]   # Tier 3 data: Δ=1 vs GWTC catalog + PBH evaporation
    if name == "chr":
        return [PY, "run_chr_experiments.py"]     # Theorem 14: Critical Horizon Response — A+C+E unified
    if name == "verl":
        return [PY, "run_verlinde_crosscheck.py"] # §7: Verlinde volume-law a₀~cH₀ vs ρ_DE~ρ_crit (O(1) check)
    if name == "eft":
        return [PY, "run_eft_stability.py"]       # §7: EFT-of-DE corner (α_T=α_M=0, c_s²=1) + coupling stability
    if name == "fap":
        return [PY, "run_fap_test.py"]            # §5.7: r_d-independent F_AP=D_M/D_H shape test (Turyshev 2026)
    if name == "eg":
        return [PY, "run_eg_test.py"]             # §6: E_G growth–geometry/no-slip test (DESI-DR1 2507.16098)
    if name == "biref":
        return [PY, "run_birefringence.py"]       # §7: birefringence discriminators (SEDE β=0 vs axion-DE rival)
    if name == "cqg":
        return [PY, "run_combined_qg_tests.py"]   # combined-theory junction tests (CKN scale × volume-law)
    if name == "fhist":
        return [PY, "run_fsat_history.py"]        # cosmic f_sat U-shape: inflation↔DE de Sitter brackets
    if name == "eosS":
        return [PY, "run_eos_entropy.py"]         # EOS = expansion vs horizon entropy production (Thm 13)
    if name == "pert":
        return [PY, "run_class_perturbations.py"] # W9: full CLASS perturbations (validates smooth-DE)
    if name == "v2x":
        return [PY, "run_xval_v2checks.py"]       # cross-team: reproduce SEDE_V2 binding/GSL + w0 split
    if name == "qg1":
        return [PY, "run_qg_route1_ckn.py"]       # QG postulate route 1: CKN fixes scale not form (reduction)
    if name == "qg2":
        return [PY, "run_qg_route2_driving.py"]   # QG postulate route 2: volume-law as driven steady state
    if name == "qg3":
        return [PY, "run_qg_route3_roughening.py"] # QG postulate route 3: Δ as roughening universality class (borrowed S1/KPZ)
    if name == "qg4":
        return [PY, "run_qg_route4_ckn_saturation.py"] # QG postulate route 4: CKN-saturation=Δ=0, ruled out (Δ>0 data-required)
    if name == "syk":
        return [PY, "run_syk_scrambling.py"]      # QG postulate route 6: SYK scrambler — maximal chaos settles the STATE half, RMT+Page+OTOC
    if name == "cset":
        return [PY, "run_causal_set.py"]          # QG postulate route 7: causal-set dof COUNTING — bulk∝volume vs links∝area (Minkowski + de Sitter)
    if name == "structde":
        return [PY, "run_structure_de_test.py"]   # structure→DE test on real data: ΛCDM vs SEDE(growth) vs Gough(SMD) on DESI BAO+fσ8 (compressed)
    if name == "goughmcmc":
        return [PY, "run_gough_mcmc.py"]          # apples-to-apples marginalised CAMB MCMC: ΛCDM vs SEDE vs Gough IDE (needs camb)
    if name == "litexp":
        return [PY, "run_literature_experiments.py"]  # falsifiability exps: H(z) dip / redshift drift / emergent-class vs Gough & ΛCDM
    if name == "e1fore":
        return [PY, "run_e1_forecast.py"]         # E1 mechanism forecast: σ(Δ)geom vs σ(Δ)growth → structure-sourcing testability (Euclid/LSST)
    if name == "wlfish":
        return [PY, "run_wl_fisher.py"]           # full WL C_ℓ tomographic Fisher: σ(Δ)_growth (refines E1 mechanism forecast)
    if name == "wlreal":
        return [PY, "run_wl_fisher_realistic.py"]  # realistic WL Fisher: Halofit+IA+baryons+full cosmology → net σ(Δ)_growth
    if name == "mech":
        return [PY, "run_mechanism_forecast.py"]  # combine growth probes + consistency-relation reframe → mechanism test ~4σ
    if name == "x2pt":
        return [PY, "run_3x2pt_forecast.py"]      # joint 3×2pt + CMB-lensing Fisher: properly-correlated σ(Δ)_growth (mechanism ~4σ)
    if name == "ebnd":
        return [PY, "run_entropy_bounds.py"]      # route C: is Δ=1 allowed? volume vs smooth-area bounds + fractal escape
    if name == "tnet":
        return [PY, "run_tensor_network.py"]      # route C: is Δ=1 realisable? RT min-cut local(area) vs nonlocal(volume)
    if name == "dind":
        return [PY, "run_delta_indirect.py"]      # indirect Δ tests with existing data (CMB R, lever arm, w0wa, growth)
    if name == "vol":
        return [PY, "run_volume_equiv.py"]        # volume formulation (S∝V_AH, no Δ) == Barrow Δ=1 guard
    if name == "loo":
        return [PY, "run_xval_loo.py"]            # leave-one-probe-out robustness (~10 min)
    if name == "tloo":
        return [PY, "run_tracer_loo.py"]          # §5.4/5.7: tracer-level LOO (drop DESI LRG1/LRG2; ~10 min)
    if name == "calib":
        return [PY, "run_xval_calibration.py"]    # false-preference mock calibration + Bayes (~15 min)
    if name == "mcmc":
        c = [PY, "run_barrow_mcmc.py"]
        if args.quick:
            c += ["--steps", "300", "--burn", "100"]
        return c
    if name == "actlens":
        return [PY, "run_act_lensing.py"]         # §5.7: ACT DR6 CMB-lensing test (answers 2509.02945; needs 344MB data)
    if name == "plik":
        return [PY, "run_plik_check.py"]
    # ── revision additions (referee response) ──
    if name == "fullcmb":
        return [PY, "run_full_cmb_mcmc.py", "--minimize"]   # E1: full primary CMB marginalised (needs cobaya/plik)
    if name == "jointcmb":
        return [PY, "run_joint_fullcmb.py", "--minimize"]   # full CMB folded into the joint (~20 min)
    if name == "gammasys":
        return [PY, "run_gamma_systematics.py"]             # γ systematics across halo-model choices
    if name == "probedecomp":
        return [PY, "run_probe_decomposition.py"]           # per-probe Δχ² decomposition at joint best-fit
    if name == "infocrit":
        return [PY, "run_info_criteria.py"]                 # AIC/BIC alongside DIC
    if name == "dyntx":
        return [PY, "run_dynT_crossing.py"]                 # crossing-redshift resolution (dynamical-T test)
    if name == "refcheck":
        return [PY, "verify_refs.py"]                       # arXiv reference verification (needs network)
    if name == "eftlag":
        return [PY, "run_eft_lagrangian.py"]                # field-theory rep + stability through w=−1 crossing
    if name == "bayesev":
        return [PY, "run_bayesian_evidence.py"]             # nested-sampling Bayesian evidence (needs dynesty)
    if name == "deltaorth":
        return [PY, "run_delta_orthogonality.py"]           # Fig 8: SEDE-vs-Barrow-HDE Δ orthogonality (the gate)
    if name == "noshoes":
        return [PY, "run_no_shoes_robustness.py"]           # no-SH0ES robustness: Δχ²/lnB/p with SH0ES dropped
    if name == "cmbfig":
        return [PY, "run_cmb_earlyde.py"]                   # Fig 5: CMB early-DE (why SEDE survives the CMB)
    if name == "e1fig":
        return [PY, "run_e1_figure.py"]                     # Fig 9: E1 mechanism on current data
    if name == "figures":
        return [PY, "make_all_figures.py"]
    raise ValueError(name)


def main():
    ap = argparse.ArgumentParser(description="Reproduce SEDE results + figures.")
    ap.add_argument("--quick", action="store_true", help="short MCMC chains (smoke test)")
    ap.add_argument("--fast", action="store_true", help="skip the slow MCMC stage")
    ap.add_argument("--plik", action="store_true", help="include the plik_lite check (needs cobaya)")
    ap.add_argument("--only", help="comma list of stages to run (subset of %s)" % ",".join(ALL))
    ap.add_argument("--skip", help="comma list of stages to skip")
    args = ap.parse_args()

    stages = list(ALL)
    if not args.plik:
        stages.remove("plik")          # opt-in (needs the native Planck likelihoods)
    if args.fast and "mcmc" in stages:
        stages.remove("mcmc")
    if args.only:
        want = [s.strip() for s in args.only.split(",")]
        stages = [s for s in ALL if s in want]
    if args.skip:
        skip = [s.strip() for s in args.skip.split(",")]
        stages = [s for s in stages if s not in skip]

    print("=" * 70)
    print("SEDE — reproduce_all.py   stages:", " -> ".join(stages))
    if args.quick:
        print("  (quick mode: MCMC uses short chains; DIC numbers are indicative only)")
    print("=" * 70, flush=True)

    results = []
    for s in stages:
        print(f"\n{'#'*70}\n# STAGE: {s}\n{'#'*70}", flush=True)
        cmd = stage_cmd(s, args)
        # drivers live in scripts/; cwd stays ROOT so their relative data/results/
        # output/ paths and flat sibling imports (scripts/ on sys.path[0]) resolve.
        if len(cmd) > 1 and isinstance(cmd[1], str) and cmd[1].endswith(".py"):
            cmd[1] = os.path.join("scripts", cmd[1])
        ok, dt = run(cmd)
        results.append((s, ok, dt))
        print(f"\n--> stage '{s}': {'PASS' if ok else 'FAIL'}  ({dt:.0f}s)", flush=True)

    print("\n" + "=" * 70 + "\nSUMMARY\n" + "=" * 70)
    for s, ok, dt in results:
        print(f"  {s:9s} {'PASS' if ok else 'FAIL':4s}  {dt:7.0f}s")
    n_fail = sum(1 for _, ok, _ in results if not ok)
    print("=" * 70)
    print("Figures (if run) are in:  output/")
    if n_fail == 0:
        print("ALL STAGES PASSED.")
    else:
        print(f"{n_fail} stage(s) FAILED — see logs above.")
    sys.exit(1 if n_fail else 0)


if __name__ == "__main__":
    main()
