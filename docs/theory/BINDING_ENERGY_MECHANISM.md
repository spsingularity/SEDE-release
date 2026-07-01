# The binding-energy → horizon-entropy registration mechanism

**Status:** theory-development scratch / not yet folded into `paper/SEDE.md`.
**Scope:** fills the one piece §3.1 / A.3 currently leaves as a *prescription* — *how* local
halo virialisation/binding energy registers on the **cosmic apparent-horizon** entropy ledger
without literally transferring enough energy to *be* dark energy.

---

## 1. The gap, as the paper already states it

The paper names this gap honestly and bounds it hard:

- **A.3 (Prescription 4C)** calls the binding-energy → horizon-entropy step a *"fixed prescription,
  not a theorem,"* and the construction concedes *"we do not claim a demonstrated transport
  mechanism."*
- **§2 / App. F** set the boundary the mechanism must respect: the released binding energy is
  `E_bind/Mc² ~ σ²/c² ~ 10⁻⁶` of ρ_DE — **six orders of magnitude too small to *be* dark energy** —
  and a literal matter→DE flux of O(ρ_crit) over a Hubble time would spoil `ρ_m ∝ a⁻³` (BBN, CMB).

So the missing content is **entropic bookkeeping, never energy**.

### 1.1 Four bars any acceptable mechanism must clear

1. **No energy transfer.** Entropy / state-counting only.
2. **Fix the *shape*, not the magnitude.** It must motivate the weight `p = 5/3` (→ γ); the magnitude
   is CKN + flatness (§8.2), full stop.
3. **Preserve f_sat as a functional of the *background* growth `D²(z)`** (c_s² = 1, Thm 6B / §3.4).
   Anything that responds to *local* δ must average back to a homogeneous background functional.
4. **Stay holographic, not modified-gravity.** §2 / §8.3: routing the deformation through the
   Einstein/Friedmann side is the `(M_P/H)^Δ ~ 10⁶¹` BBN catastrophe. Registration must live in the
   DE fluid.

### 1.2 What a mechanism can and cannot buy

A mechanism here upgrades the **prose** — from *"we don't claim a transport mechanism"* to
*"registration is entropic bookkeeping via X"* — **without** upgrading the **epistemic status**
(it stays a *motivated prescription*). **None of these derives `p = 5/3` from scratch**; that still
comes from `E_bind ∝ M^{5/3}` plus the choice of `T_AH` as the reservoir temperature. Keep that
honest.

---

## 2. Candidate mechanisms — ranking

Ranked by (defensibility + fit to the paper's existing commitments + whether it touches `p = 5/3` +
referee safety).

| Rank | Mechanism | Fit to paper | Verdict |
|---|---|---|---|
| **1** | **Response-field / order parameter** (`φ = f_sat`, source `J ~ D²`, susceptibility `φ = 1−e^{−γJ}`) | Maps onto the *exact* `f_sat = (1−e^{−γD²})/(1−e^{−γ})` form already in A.2 — that kernel **is** a saturating susceptibility to source `D²` with coupling γ. "Structure is the control field, not the fuel" is §2 verbatim. | **Primary framing.** Cleanest defeat of the energy-transfer objection (magnetism analogy: H controls M, doesn't become it). Least retrofitting. |
| **2** | **Separate-universe patchwise horizons** (`δ_L → H_local → R_AH,local → S_AH,local`; global f_sat = ensemble average) | Only candidate with a real *computational* handle (peak–background split, standard & controlled). Directly reconciles bar 3: DE smooth locally, yet the patch-averaged horizon ledger responds to long-wavelength structure → naturally a functional of `D²(z)`. Stays GR + fluid. | **Computational backbone.** Adopt alongside #1. |
| **3** | **Gravitational-memory ledger** (horizon = causal boundary whose state-count responds to interior stress-energy history) | Matches the paper's existing "ledger"/bookkeeping language and the no-transport stance perfectly. | **Interpretive spine.** Safest, but framing not content — pair with #2's machinery. |
| **4** | **Information-erasure / phase-space irreversibility** (shell-crossing + virialisation → lost phase-space info → ΔS_grav) | Best "why an *entropy*, and why irreversible" content; supports the Clausius `ΔS = E_bind/T_AH` step. Cleanest arrow-of-time story. | **Supporting language.** Strong companion to #3; doesn't alone justify the *cosmic-horizon* T. |
| **5** | **Entanglement-wedge re-encoding** (collapse → bulk entanglement reorganisation → changed horizon entropy capacity) | Deepest fit to SEDE's *actual* postulate (§8.4: s_grav = coarse-grained **entanglement** entropy of the horizon volume) and the only candidate that *naturally* yields `M^{5/3}`. | **One-line holographic reading only.** Must be labelled interpretation, not derivation — but touches the model's real ontology. |
| **6** | **Coarse-graining / Buchert backreaction** (inhomogeneous GR → smoothed FRW → effective horizon pressure → ρ_DE) | Intuitive, referee-recognisable, BUT routes through "effective pressure → ρ_DE," i.e. the **averaging-modifies-Friedmann** framing the paper rejects (and Green–Wald-style critiques target). | **Risky for *this* paper.** Pulls toward the modified-gravity quadrant §8.3 disowns. Use the *concept* (info lost under smoothing), not the Buchert *machinery*. |
| **7** | **Soft-hair / BMS** (horizon carries soft gravitational memory written by structure) | Most speculative; BMS/soft hair is well-defined at null infinity / BH horizons, **not** at the cosmological apparent horizon. | **Avoid as load-bearing.** Footnote at most. |

---

## 3. Recommended stack

The memory-ledger is the right *spine*, but the optimal construction is a three-layer stack:

- **Framing:** #1 (response/control-field) — "structure is the control field, not the fuel"; it is
  already the math (`f_sat` = saturating response to `D²`).
- **Mechanism with teeth:** #2 (separate-universe) — makes "horizon responds to structure while DE
  stays smooth" *calculable* and protects c_s² = 1 by construction (patch average → background
  functional).
- **Language:** #3 + #4 — memory ledger + irreversibility, which is what the paper already gestures at.
- **Depth note:** #5 — a single sentence connecting to the entanglement-volume postulate (§8.4).

### 3.1 Two cautions when this is eventually written into the paper

- The phrase *"removes phase-space information from the smooth FRW description"* is the
  coarse-graining idea (#4/#6). Fine — but the **next clause must keep the registration in the DE
  fluid**, or a referee reads it as Buchert backreaction and invokes §8.3 against the model.
- **Hold the line on status.** Replace *"we do not claim a transport mechanism"* with *"the
  registration is entropic horizon bookkeeping (a patch-averaged horizon response /
  gravitational-memory ledger), not heat transport."* This strengthens the prose while keeping A.3
  honestly a *prescription*.

---

## 4. Distilled statement

> Structure formation writes gravitational memory into the coarse-grained horizon state. Binding
> energy measures the irreversibility of that write operation. The horizon does not *receive* the
> energy; it *updates its activation state*. The released binding energy `E_bind ∝ M^{5/3}` fixes the
> entropy **weight** (`p = 5/3 → γ`, the redshift shape of f_sat); the **magnitude** of ρ_DE remains
> the horizon scale (CKN + flatness, §8.2), not energy borrowed from structure.

In schematic form (bookkeeping, **not** transported heat):

```
virialisation → E_bind ∝ M^{5/3} → ΔS_ledger ~ E_bind / T_AH → γ → f_sat(D²) → ρ_DE ∝ H f_sat
                                       (ΔS_ledger ≠ ΔS_literal-transported)
```

---

## 4a. How does halo collapse "write gravitational information"?

This is the load-bearing question. It has a real answer for the **first half** (collapse genuinely
produces gravitational entropy — established physics) and an honest **postulate** for the second half
(that the cosmic horizon is what counts it). Keep them separate.

### (i) Established physics: collapse increases gravitational entropy

"Writing gravitational information" is not a metaphor for heat — it is the well-understood fact that
gravitational collapse is an **entropy-producing, irreversible** process at the level of the *field
configuration*, not the gas.

- **Violent relaxation / phase mixing (Lynden-Bell 1967).** A halo forms from cold, nearly-uniform
  initial conditions (low velocity dispersion, smooth flow, shallow potential). As it collapses the
  potential fluctuates on the collapse timescale and scatters orbits. The *fine-grained* phase-space
  distribution is conserved (Liouville — no information is truly destroyed), but the **coarse-grained**
  distribution any cosmological observer actually uses spreads and mixes → a strict entropy increase.
  Shell-crossing / caustic formation is the same thing geometrically: the smooth single-valued velocity
  field becomes multi-valued and the laminar pre-collapse flow becomes unrecoverable from the smoothed
  description.
- **It is gravitational, not thermal (Penrose Weyl-curvature hypothesis).** The early universe starts
  in a very *low* gravitational-entropy state (Weyl curvature ≈ 0, smooth, no lumps); structure
  formation drives it toward *high* gravitational entropy (large Weyl curvature, deep wells, clumpy
  field). The "information written" is precisely this growth of structure **in the gravitational field
  itself** — the metric goes from featureless to richly structured, irreversibly.

So *halo collapse writes gravitational information* means: **collapse converts a smooth, low-entropy
gravitational configuration into a structured, high-entropy one, and the coarse-grained observer cannot
undo it.**

### (ii) Why binding energy is the right *measure* of how much got written

The amount of irreversible reorganisation scales with how deep the wells got. By the virial theorem
(`2K + U = 0`) a virialised halo's binding energy is `|E_bind| ~ K ~ Mσ²` — the kinetic energy it had
to develop to settle. With NFW-like `R_vir ∝ M^{1/3}` this gives `E_bind ∝ M^{5/3}`. So `E_bind` is a
clean proxy for "strength of the write operation" — and the `σ²/c² ~ 10⁻⁶` smallness (App. F) is
exactly why this is an *entropy bookkeeping* quantity and never an energy supply.

### (iii) SEDE's postulate: why the *cosmic horizon* is the page

Here the standard physics ends and the model's postulate begins (be honest — the paper is). Steps (i)
and (ii) are textbook. The leap is *why* this locally-produced gravitational entropy registers on the
cosmic apparent horizon, at temperature `T_AH`. SEDE's answer is the volume-law ansatz, and it
dissolves the "how does it travel there?" problem by making transport unnecessary:

> The horizon entropy in SEDE is a **coarse-grained count of the bulk gravitational degrees of freedom
> of the whole causal patch** (the volume-law / entanglement-entropy postulate, §8.4) — *not* a surface
> integral of heat near a 2-sphere.

If that is what horizon entropy *is*, collapse changes it **directly, with no signal propagating
outward**, because collapse reorganises the very bulk gravitational dof the horizon entropy counts. The
"writing" and the "reading" are the *same* degrees of freedom — nothing is deposited *onto* a distant
surface; the global coarse-grained gravitational state the horizon entropy tallies is simply different
after a halo forms than before. This is why §2 can say "entropy bookkeeping, not energy transfer" and
mean it literally. The separate-universe picture (§2 mechanism #2) is the concrete handle: a collapsing
region slightly shifts its patch's local expansion and hence its local horizon state; averaging over
patches yields the homogeneous `f_sat`.

### (iv) The honest boundary

- **Solid:** collapse is irreversible and increases coarse-grained gravitational entropy (violent
  relaxation, Weyl growth, Penrose). `E_bind ∝ M^{5/3}` measures it. Textbook.
- **Postulate:** that the *cosmic-horizon* entropy is the bulk count which this increment updates, at
  temperature `T_AH`. The model's one open foundational item (§8.4) — it makes the registration need
  *no transport*, but it is asserted, not derived, and is what DESI DR3 + Euclid test through Δ.

**Net:** collapse irreversibly restructures the gravitational field (real, measurable as binding
energy); and because SEDE postulates the horizon entropy to be the bulk gravitational count of the
causal patch, that same restructuring *is* a change in the horizon ledger — not a delivery to it. The
first clause is physics; the second is the bet the model is built on.

---

## 5. Draft insert (for §3.1 / A.3 — NOT yet applied)

> The binding-energy-to-horizon-entropy step should not be read as a literal transport of virial heat
> to the apparent horizon. The mechanism is a **gravitational-memory ledger**: halo virialisation
> irreversibly reorganises the gravitational field and removes phase-space information from the smooth
> FRW description. Since SEDE defines the dark-energy sector as a horizon thermodynamic sector
> (a smooth `c_s²=1` fluid, §3.4 — *not* a modified Friedmann equation), this irreversible
> coarse-graining is registered as a change in the horizon **activation coordinate** `f_sat` rather
> than as a local energy flux. Equivalently (separate-universe picture), long-wavelength structure
> shifts each FRW patch's local horizon state `S_AH,local`, and `f_sat` is the homogeneous
> patch-average — smooth on small scales, yet responsive to structure, consistent with c_s²=1. The
> binding energy `E_bind ∝ M^{5/3}` measures the strength of this irreversible reorganisation; dividing
> by `T_AH` assigns the corresponding entropy weight (`p = 5/3 → γ`). Binding energy thus fixes the
> response kernel γ and the redshift shape of `f_sat`, while the magnitude of ρ_DE remains the horizon
> scale (CKN + flatness, §8.2). This is a **fixed prescription** for the registration, not a derived
> transport mechanism.

> **Note:** §6 (route C8) is a stronger replacement for this paragraph — it *derives* the temperature
> and the Clausius form instead of asserting them. Prefer C8 when this is folded into the paper.

---

## 6. Route C8 — modular-energy / entanglement-equilibrium (rigorous version of §4a-iii)

**What it is.** Not a new mechanism — it is the entanglement-wedge route of §4a(iii) (mechanism #5,
absorbing #1's no-transport and #4's irreversibility) promoted from postulate to derivation sketch
using the **entanglement first law** and the **CHM modular Hamiltonian**.

**The move.** The horizon entropy is the entanglement entropy of the field state reduced to the
horizon ball; it changes when the *global* state changes, so structure formation changes it with
nothing travelling to `R_AH`.

- *Step 1 — first law of entanglement.* For any region, `δS_EE = δ⟨K⟩` to first order (`K` the modular
  Hamiltonian, `ρ = e^{−K}/Z`). Exact, general, no dynamics — relates the entropy of the region to the
  modular energy of the *same* region.
- *Step 2 — modular Hamiltonian of the horizon ball (Casini–Huerta–Myers).* For a ball of radius `R`
  in a maximally symmetric vacuum, `K = (2π/ℏ)∫ (R²−r²)/(2R) T₀₀ d³x + (area term)`, and the conformal
  Killing flow fixes the modular temperature to `T_mod = ℏH/2π = T_AH`. **The temperature is derived,
  not inserted** — a second, independent route to `T_AH` beyond the §2 identity.
- *Step 3 — modular energy of a localised halo.* A halo at `r ≪ R_AH` sees a near-constant modular
  weight `≈ R/2 = 1/(2H)`, so `δ⟨K⟩ = δE/T_AH × O(1)` (the profile O(1) is the kind already absorbed
  into `f_sat(0)=1`). With Step 1: `δS_horizon = δE/T_AH` — the Clausius form, **derived** from the
  first law, not asserted via heat flux.

**Identifying δE with binding energy.** Use smooth FRW as reference, virialised as perturbation, at
fixed rest mass. The rest-mass `T₀₀` cancels (why it is *not* `p=1`); the surviving configuration piece
is the dissipated heat `δE ≃ |E_bind| ∝ GM²/R_vir ∝ M^{5/3}` (virial theorem makes the exponent robust
to which energy you name). **Hence `p = 5/3`, sign correct, no transport.**

**Ledger — derived vs assumed (honest).**

| piece | status |
|---|---|
| temperature `T_AH` | **derived** (modular temperature, CHM) |
| Clausius form `δS = δE/T_AH` | **derived** (entanglement first law) |
| `p = 5/3` scaling | **derived** (virial theorem, given δE = dissipated heat) |
| sign `δS > 0` | **assumed** — requires δE = *dissipated* (coarse-grained) heat, not the conservative `δ⟨K⟩` of the actual states (which is more-bound → lower energy → would give δS<0) |
| CHM exactness | **approximate** — exact only for a CFT vacuum + spherical region; the clumpy cosmic diamond is an approximation |

**The soft spot C8 relocates but does not dissolve.** The first law `δS_EE = δ⟨K⟩` is about the
**fine-grained** (von Neumann) entropy, which is *conserved* under collapse (Liouville). The
binding/virialisation entropy is **coarse-grained** (Lynden-Bell). C8's sign fix — switching to the
*dissipated* heat — is exactly the point where the coarse-graining / second-law input re-enters; the
reversible first law does not supply it. This is the *same* fine-vs-coarse subtlety as §4a(i),
relocated, not removed. So C8 honestly delivers **{temperature, Clausius form}** as derived and the
**{sign, fine→coarse identification}** as still-assumed (better-localised) inputs. Do not let a boxed
`δS = δE/T_AH` imply the sign is derived.

**de Sitter modular-temperature check — RESOLVED (2026-06-27, analytic).** Settled by the structure of
the CHM / Kodama–Hayward construction (surface gravity is geometric — no paper-fetch needed):

- *Exact dS (ε=0).* As the CHM ball → the horizon, the conformal Killing vector generating the modular
  flow becomes the dS static-patch boost; `T_mod = H/2π` exactly (Gibbons–Hawking). Since H is constant,
  `T_AH = (H/2π)(1−ε/2) = H/2π`. **No factor — the modular temperature *is* the horizon temperature.**
  This is C8's second independent route to `T_AH` (the first being the §2 identity).
- *Quasi-dS FRW (ε≠0).* Matter isn't conformal, so CHM is approximate, but the near-edge modular
  temperature is the apparent-horizon **surface gravity** → the Kodama–Hayward value `(H/2π)(1−ε/2)` —
  exactly the Cai–Kim factor SEDE already absorbs into `f_sat(0)=1`. C8 introduces nothing new.
- *Why p=5/3 is untouched.* The `(1−ε/2)` factor is a property of the **epoch**, identical for every
  halo at fixed z; `p` comes only from the **mass** dependence of `ΔS_AH = E_bind/T_AH`, and `T_AH`
  (factor and all) is mass-independent. So ε rescales the amplitude of `f_sat`, never the M-scaling.

Caveat to keep: the only ε-dependence is the surface-gravity one (the modular flow is generated by the
horizon Kodama/Killing vector whose normalisation *is* the surface gravity) — confident, not just
"fairly certain." The genuine residual caveat is the **CFT/conformal-matter** assumption in CHM, not a
hidden ε-factor. Note ε ≈ 0.47 today ⇒ `(1−ε/2) ≈ 0.76`: this O(1) factor is real, not negligible, so
the `f_sat(0)=1` normalisation is doing genuine work — but it cancels from `p`.

---

## 7. Route C9 — self-similar reduction of γ  (NUMERICALLY VERIFIED, 2026-06-27)

C9 propagates the per-halo weight `M^{5/3}` (C8) through the mass function:
`Σ_S(σ8) = ∫ M^p (dn/dlnM) dlnM`, `γ ≡ dlnΣ_S/dlnσ8`. With peak height `ν = δ_c/σ(M)`,
`α ≡ −dlnσ/dlnM = (n_eff+3)/6`, the exact reduction is

> **`γ = (p−1)·⟨1/α⟩_Σ = (p−1)·⟨6/(n_eff+3)⟩_Σ`**,  weight `∝ M^p dn/dlnM`,

with scale-free limit `γ_SF = 6(p−1)/(n_eff+3)`.

### 7.1 Verification against the repo's own COLOSSUS computation

Ran `gamma_computation`-style integrals directly (`scratchpad/c9_check.py`, Sheth–Tormen, EH98,
M ∈ [10¹⁰,10¹⁶], Planck-18 fiducial). Results:

| quantity | value |
|---|---|
| `dlnΣ/dlnσ8`   (C9 convention)              | **2.993** |
| `dlnΣ/dlnσ8²` (code / `f_sat` convention)   | **1.496** ← the paper's "1.496" |
| C9 closed form `(p−1)⟨1/α⟩_Σ` (pop. average) | **2.991**  → **closed/exact = 0.999** |
| C9 single-scale `(p−1)/α(M_S)`               | 2.765      → closed/exact = 0.924 |
| deposition peak `M_S`, `ν_S`                 | **10¹⁴·⁵ M⊙/h, ν_S = 2.36** |
| `n_eff(M_S)`                                 | −1.55 |
| `p=1` value `dlnΣ/dlnσ8`                      | 0.531  (= 0.27 in σ8² convention) |

### 7.2 The headline correction: the "30–50 % gap" was a CONVENTION ARTIFACT

C9's hand-estimate (γ ≈ 2.0–2.3) was compared against the code's **1.496** and called "high by
30–50 %." But those are **different conventions** — they differ by exactly ×2:

- the code's `1.496 = dlnΣ/dlnσ8²` (variance, because `f_sat` uses `x = D² = σ8²/σ8²(0)`);
- C9's closed form is `dlnΣ/dlnσ8` (amplitude).

In the **same** convention, C9's `(p−1)⟨1/α⟩_Σ = 2.991` reproduces the exact `dlnΣ/dlnσ8 = 2.993` to
**0.1 %**. **C9 is not "off by 30–50 %" — it is exact.** The reduction is correct; the apparent gap was
self-inflicted by mixing σ8 and σ8² derivatives.

### 7.3 Two corrections to C9's *narrative* (the formula stands; the story was wrong twice)

1. **There is no gap to explain.** C9's reconciliation prose ("closed form high, low-ν shoulder pulls
   toward shallower slope") is moot — once conventions match the closed form *is* the exact value.
2. **Population-averaging pulls γ UP, not down.** Single-scale at the peak `M_S` gives 2.765; the
   population average gives 2.991 (*higher*). The low-ν shoulder is *lower* mass → *more* negative
   `n_eff` → *larger* `1/α` → *larger* γ. (This confirms my spectral-slope critique; C9's stated
   direction was backwards. The effect is small, ~8 %, swamped by the ×2 convention error.)

### 7.4 What C9 gets RIGHT (and it's strong)

- **The reduction is exact** (0.1 %): `γ = (p−1)⟨1/α⟩_Σ` is a genuine analytic skeleton for γ.
- **Deposition scale confirmed:** entropy is deposited by `ν_S ≈ 2.4` (~2σ) clusters at
  `M_S ~ 10¹⁴·⁵ M⊙/h` — exactly C9's "massive clusters above M*" claim. → independent falsifiable
  handle via cluster-abundance evolution (eROSITA / SPT / Euclid clusters).
- **p=1 → ≈0 holds (with a caveat):** `p=1` gives 0.27 (σ8² conv), vs 1.496 at p=5/3 — a factor ~5.5,
  i.e. "weak gate vs real gate," not a fine distinction. But it is **not exactly 0**: the `(p−1)`
  prefactor predicts 0, the true residual ≈0.27 is a finite-mass-range boundary term (ST normalisation
  `∫f dν = 1` only over the *full* range). State it as "≈0.27, ~5× below p=5/3," not "=0." (Matches
  `gamma_computation.entropy_weight_scan`: p=1 → 0.27.)

### 7.5 ⚠ Paper-text issue found (flagged, NOT fixed — paper untouched)

`§3.1` and `A.3` write **"γ = dlnΣ_S/dlnσ8 = 1.496"**. The *number* 1.496 is correct and consistent
with the `f_sat` that the code actually uses, but the *denominator label is wrong*: numerically
`dlnΣ/dlnσ8 = 2.993`, and it is `dlnΣ/dlnσ8²` (≡ `(1/2)dlnΣ/dlnσ8`) that equals 1.496. The formula
should read `dln σ8²` (or `(1/2) dln σ8`) to match the `x = D² = σ8²/σ8²(0)` convention of A.2. Minor
but real — fix when the paper is next edited.

### 7.6 Closed-form status — which form is exact (numerically tested 2026-06-27)

Two closed forms are on the table; **neither is uniformly exact**, so label them honestly:

| form | value (σ8² conv) | vs exact 1.496 | regime |
|---|---|---|---|
| method-A `(p−1)⟨1/α⟩_Σ` | 1.496 | **0.1 %** | exact in the *cutoff-dominated* regime (= p=5/3, our case). **Fails at p=1** (predicts 0, true 0.27 — drops the measure term). |
| single-scale `½(a·ν_S²−1)`, a=0.707, ν_S=2.36 | 1.468 | 1.9 % | good *heuristic*; carries the physics (γ = cluster-cutoff steepness). The ~2 % is partly fortuitous. |
| pop-averaged `½⟨aν²−1+ST⟩` | 1.244 | 17 % | my rigorous-ification of the heuristic — *worse*, so the cutoff form is **not** exact. |

**For the paper:** quote the full-integral **1.496** as *the* value. Use `½(a·ν_S²−1)` for intuition
(label "leading single-scale approximation"); optionally cite `(p−1)⟨1/α⟩` as the exact reduction *in
the cutoff regime*. Do **not** present any closed form as exact in general.

### 7.7 Consistency caution when fixing A.3 (do not break theory↔data)

When the A.3 definition is changed `d ln σ8 → d ln σ8²`, the **data-side** `γ_data ≈ 1.41` must stay in
the same convention. It already is — `γ_data` multiplies `x = D²` in the likelihood's `f_sat`, hence is
intrinsically σ8². Fixing the theory definition makes both sides consistent and the
"`γ_theory = 1.496 ≈ γ_data`" vindication survives. **Halve only the written definition's apparent
value, never one side of the comparison.**

### 7.8 Status

C9 is **paper-ready** once (a) it is stated in the σ8² convention to match `f_sat`/A.3 (fix the
`d ln σ8 → d ln σ8²` typo, §7.5), (b) the "30–50 % gap" reconciliation prose is deleted (no gap
exists), (c) the p=1 claim is softened to "≈0.27, ~5× below p=5/3," and (d) the closed form is labelled
a leading approximation (§7.6). The exact value 1.496 is independently confirmed by two transfer
functions (COLOSSUS + no-wiggle EH98, ~6 % apart).

*(Verification script: `scratchpad/c9_check.py`. To make permanent, port into `run_gamma_systematics.py`
and add the `γ(M)` / deposition-kernel figure for the paper.)*
