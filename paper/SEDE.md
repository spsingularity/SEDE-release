# Structural Entropy Dark Energy: a fixed-parameter, growth-gated holographic dark-energy model without a cosmological constant

*(JCAP/PRD target — prepared for submission. §§1–10 + appendices A–G in full prose; App. A documents the
load-bearing results and assumptions (one defining ansatz, fixed prescriptions, supporting derivations);
Figures 1–7 in the main text, Fig. G1 in App. G (sources in `output/`, regenerable via
`make_all_figures.py`). All 45 references resolve on the live arXiv (verified, `verify_refs.py`).
Numbers are the verified values from the repository; reproduce via `reproduce_all.py`, checked by the
verification suite (App. B). Pending at proof: a final boldface/copyedit pass and journal volume/page
confirmation.)*

## Abstract

We present Structural Entropy Dark Energy (SEDE), a model without a cosmological constant in which dark
energy is identified with the thermodynamic conjugate of the cosmic apparent-horizon entropy,
ρ_DE = T_AH s_grav f_sat(z),
gated by a structure-growth factor f_sat whose redshift dependence follows a halo-binding-entropy
prescription. A single Friedmann fixed-point equation replaces Λ by a structure-gated, horizon-coupled
term that adds **no continuously-sampled dark-energy parameter** relative to ΛCDM: the dark sector is
fixed instead by a few stated, discrete modelling choices — the H-coupling λ = 1/2 from a volume-law
horizon entropy, the equation of state from spatial flatness, the coupling γ ≈ 1.5 from a halo-binding
prescription, the sound speed c_s² = 1 — one of which (the entropy weight p) was revised once from a
halo to a horizon reservoir. In a marginalised, CAMB-in-the-loop joint analysis (DESI DR2 BAO,
Pantheon+/DES-5YR/Union3 supernovae, cosmic chronometers, fσ8, and compressed Planck data), SEDE is
**mildly preferred** over ΛCDM (ΔDIC ≈ −2.9). The preference survives folding the full primary CMB
likelihood into the joint (Δχ² = −3.17) and the full non-lite Planck plik with all foregrounds sampled
(Δχ² = −0.33) — so it is not an artifact of CMB compression; a nested-sampling Bayesian evidence gives
ln B = +1.89 (positive-but-weak, favouring SEDE). A per-probe decomposition shows this preference is
carried by the CMB-distance and H₀ (SH0ES) channels, **not** by the DESI BAO data — which SEDE fits
comparably to ΛCDM — so "mildly preferred" should not be read as a DESI-BAO preference. Dropping the
contested SH0ES prior entirely, the preference weakens but survives and shifts onto the harder channel
(Δχ² = −1.75, ln B = +0.95, ~1.9σ; the CMB-distance carries it with H₀ no longer pulled high). A
false-preference calibration indicates the preference is **not obviously attributable to model
flexibility** (~2–3σ with SH0ES, ~1.9σ without, depending on the CMB-distance treatment), and a pre-DESI → DESI split shows a SEDE realisation trained on
pre-DESI data predicts the DESI DR2 BAO trend better than ΛCDM. The equation
of state crosses −1 with w_a < 0, in DESI's evolving-dark-energy quadrant, and the vacuum-energy-scale
problem is recast via the Cohen–Kaplan–Nelson bound rather than inherited. The construction rests on
**one explicit quantum-gravitational postulate** — a volume-law (thermalised) horizon entropy — testable
by a measurement of the deformation Δ; under the stated Fisher assumptions, DESI DR3 and Euclid reach
σ(Δ) ≈ 0.09, sharply distinguishing the volume-law and area-law branches. SEDE is **suggestive but not
established**; the next generation of surveys will decide.

---

## 1. Introduction

The standard cosmological model, ΛCDM, fits an extraordinary range of data with a cosmological
constant Λ whose value confronts two long-standing puzzles. The **cosmological-constant problem**:
if Λ is the energy of the quantum vacuum, naïve effective field theory overshoots the observed value
by ~10¹²⁰. The **coincidence problem**: the dark-energy and matter densities are comparable today,
ρ_DE ~ ρ_m, despite scaling differently with expansion — so we appear to live at a special epoch.
Recent DESI DR2 results, moreover, give a mild preference for *evolving* dark energy (w crossing −1)
over a constant Λ, reopening the question of whether Λ is fundamental at all.

We approach this through an accounting question. As cosmic structure forms, collapsing halos release
gravitational binding energy, and dark energy activates over the *same* epoch that structure
formation approaches saturation. Is that a coincidence, or is the dark sector tied to structure formation by an energy
ledger? It is tempting to guess that the binding energy released by structure *is* the dark energy —
but a one-line budget settles it: a virialised system has |E_bind|/Mc² ~ σ²/c², and even clusters
have σ²/c² ~ 10⁻⁵ (galaxies ~10⁻⁷), so the binding-energy density is ~10⁶ times too small to *be*
the observed ρ_DE. The link, if there is one, cannot be an energy supply.

The resolution we propose is thermodynamic. We identify dark energy with the **conjugate thermal
energy density of the cosmic horizon** — ρ_DE = T_AH·s_grav, the conjugate
energy of the apparent-horizon entropy (motivated by
the horizon thermodynamics of Jacobson 1995; Cai–Kim 2005) — whose *scale* is the horizon's own,
T_AH s ~ M_P²H² ~ ρ_crit, rather than the QFT vacuum (so the vacuum-energy-scale problem is recast,
not inherited) and not structure. What structure supplies is not
energy but **entropy bookkeeping**: by the horizon first law δQ = T dS, the binding-entropy ledger
associated with halo virialisation fixes the *redshift dependence* of the effective activation
fraction f_sat(z) ∈ [0,1], while its amplitude and the f_sat(0) = 1 normalisation are fixed
separately by the horizon-scale normalisation. Energy is conserved throughout — total covariant conservation is automatic,
and the dark-energy fluid is self-conserved — with the binding energy playing a purely entropic,
shape-setting role (§2). The coincidence is then a consequence of the mechanism, not a tuning: dark
energy turns on with a gate whose redshift dependence tracks structure formation as it approaches
saturation, i.e. near "now."

The idea that dark energy is tied to the entropy/information produced as structure forms is not new —
it originates with Gough's information dark energy [Gough 2008, 2011, 2025] — and we credit that
priority explicitly. What is specific to SEDE is the *mechanism* (a gravitational horizon-entropy
realisation, ρ_DE = T_AH·s_grav, driven by the gravitational growth factor rather than the baryonic
star-formation history) and a parameter-free, first-principles dark sector; the relationship is
quantified in §7.6.

The model is, to a useful approximation, **ΛCDM with a cosmological constant that fades in as
structure forms**: schematically Λ_eff(z) ∝ f_sat(z) (the exact running term is
3Ω_DE0 H H₀ f_sat, §2), with f_sat running from ≈0 in the early universe to 1
today. By **Λ-free** we mean specifically that the late-time acceleration is *not* a constant
stress-energy term (p = −ρ = const.) but a self-conserved horizon fluid whose present amplitude is
fixed by flatness (Ω_DE0 = 1 − Ω_m − Ω_r) — exactly as flat ΛCDM fixes Ω_Λ0 — and whose redshift
dependence is set by the horizon/growth prescription rather than held constant; it does *not* mean the
absence of a dark-energy normalisation. Its distinctive feature is that the dark sector adds **no
fitted parameter** — every input is fixed by a stated theoretical choice (§3), so SEDE has the same
cosmological parameter count as the ΛCDM baseline used in this analysis (the five of the compressed
late-time pipeline of §4). The price
is a single postulate about the quantum-gravitational structure of the horizon (§8), which we make
explicit, plus the accompanying fixed choices (the halo-binding-entropy derivation of γ (Thm 4C), the smooth-DE
closure c_s²=1), which we flag as such rather than as derivations.

This paper is organised around the energy ledger. §2 states the model and establishes energy
conservation as its backbone, fixing the (entropic, not energetic) role of binding energy; §3 fixes
the dark-sector inputs — the structure coupling γ from the halo-binding-entropy derivation (Thm 4C), then w(z), λ, and the sound
speed c_s² — separating the horizon-set *magnitude* from the structure-set *shape*. §4 describes the
data and the CAMB-in-the-loop methodology; §5 presents the results, including a false-preference
calibration and out-of-sample tests; §6 gives the falsifiable predictions; §7 the relation to dark
matter and other sectors. §8 then supplies the quantum-gravity underpinning — the magnitude scale and
the maximal deformation — and isolates the one irreducible postulate; §9 discusses the open problems
and §10 concludes. Throughout we maintain that SEDE is a *falsifiable proposal*, preferred at a modest
level, not an established result.

---

## 2. The model and energy conservation

SEDE is general relativity with a single extra fluid: a horizon-fluid dark-energy density identified
with the conjugate thermal energy density of the cosmic apparent horizon. The whole construction is governed by one bookkeeping rule — that energy is conserved when
this fluid's redshift dependence is gated by structure growth. We first state the model (the conjugate horizon-fluid ansatz, the structure
gate, the Friedmann equation, and how Λ is replaced), and then make energy conservation the backbone:
it is what fixes the role of structure's gravitational binding energy as *entropic*, not energetic, and
what shows the w = −1 crossing to be sourcing rather than a phantom field.

**The conjugate horizon-fluid ansatz.** At the cosmological apparent horizon R_AH = 1/H (flat FRW),
the dynamical Cai–Kim temperature is T_AH = (H/2π)(1 − ε/2), ε ≡ −Ḣ/H². For the dark-energy sector,
which tends to a de Sitter attractor (ε → 0), we adopt the **Gibbons–Hawking temperature T_AH = H/2π**
— the ε → 0 limit, i.e. the appropriate *leading* temperature for a slowly-evolving horizon fluid. This
is an approximation, made deliberately: the dynamical (1 − ε/2) correction is largest in the matter era,
exactly where f_sat → 0 and dark energy is negligible, and including it in full on the volume-law fluid
is in fact *disfavoured* — it yields an unstable, phantom-everywhere background with no w = −1 crossing
(§6, `run_dynT_crossing.py`). So T_AH = H/2π is the physically-appropriate choice for this sector and
it is what enters the calibrated E(z). We *identify* the dark-energy density with the conjugate
horizon-fluid energy density,
$$\rho_{\rm DE} = T_{AH}\, s_{\rm grav}\, f_{\rm sat}(z),$$
where s_grav is the gravitational entropy *density* and f_sat is an effective horizon-entropy
activation **gate**, normalised to f_sat(0) = 1 — i.e. a value *relative to today*, not a literal
fraction bounded by 1 (it reaches ≈ 1.23 in the future, below); its redshift dependence is set by the
growth-weighted halo prescription of §3.1, not derived from the Bousso bound itself. *Range and future behaviour (defining the
convention explicitly):* for z ≥ 0 — the regime of every data fit here — the linear growth obeys
0 ≤ D(z) ≤ 1, so 0 ≤ f_sat ≤ 1. f_sat is **not clipped**; in the future (z < 0) linear growth *saturates*
(D freezes at D_∞ ≈ 1.44 as dark energy suppresses growth — it does not run away), so f_sat saturates to a
**finite** de Sitter value f_sat(D_∞) ≈ 1.23 (modestly above 1, never reaching the formal x → ∞ ceiling
1/(1−e^{−γ}) ≈ 1.29). The background tends to a de Sitter attractor (H → H_∞, ρ_DE → const, w → −1; A.7,
Result 12) — bounded and well-defined for all z, with the ≤ 1 statement specific to z ≥ 0. (The gate is
normalised to its present value; it is *not* the literal fraction of horizon entropy deposited by
collapsed structures, which is far smaller — ∼10⁻⁷; App. F.) Clausius horizon thermodynamics motivates the relation and
supplies the conservation logic, but the absolute identification is the defining SEDE **ansatz**, not
an identity (App. A.1). We take s_grav to be the coarse-grained **entanglement entropy in the
apparent-horizon volume**, S_grav ∝ V_AH ∝ R_AH³ — equivalently a *constant* entropy density
s_grav = s_0 (in contrast to the Bekenstein–Hawking area-law density s_AH ∝ H). With T_AH ∝ H,
$$\rho_{\rm DE} = T_{AH}\,s_0\,f_{\rm sat} \;\propto\; H\,f_{\rm sat}, \qquad
s_0 = \frac{\Omega_{\rm DE0}\,\rho_{\rm crit,0}}{T_{AH,0}}\ \text{(fixed by flatness)},$$
giving **ρ_DE ∝ H f_sat** with no free parameter. This volume law is exactly the maximal Barrow
deformation (S ∝ A^{1+Δ/2} = A^{3/2} = R³ at Δ=1, λ = 1 − Δ/2 = 1/2); we keep Δ only as a *test*
parameter quantifying deviations from the volume law (§8, and the data measurement of Δ), not as a
model input.

**The structure gate.** Structure formation supplies the *redshift dependence* of the activation gate
at a rate set by the collapsed-halo abundance (Result 3, App. A), yielding a closed form for the
normalised activation fraction,
$$f_{\rm sat}(z) = \frac{1 - e^{-\gamma D^2(z)}}{1 - e^{-\gamma}}, \qquad D(z)\ \text{the linear growth factor},$$
where D(z) is the model's **self-consistent** growth factor (computed with SEDE's own E(z), which in
turn depends on f_sat — the background is solved as the joint fixed point; §4.4, App. B).
with f_sat(0) = 1 imposed as the present-day normalisation convention and f_sat → 0 at high z. The
single shape parameter γ is fixed by a halo-binding-entropy derivation (§3.1, Thm 4C), not fitted.

**The Friedmann equation.** Inserting ρ_DE into the (unmodified) Friedmann constraint gives a single
**fixed-point equation** for the dimensionless expansion rate E = H/H₀:
$$\boxed{\;E^2(z) = \Omega_m(1+z)^3 + \Omega_r(1+z)^4 + \Omega_{\rm DE0}\, f_{\rm sat}(z)\, E(z)\;}\tag{1}$$
(for λ = 1/2, so E^{2λ} = E). Equation (1) is the entire background model. It is ΛCDM with the
constant Λ replaced by the structure-gated, horizon-coupled term Ω_DE0 f_sat(z) E(z); at z = 0 with
f_sat = 1 it reproduces the observed Ω_DE0 = 1 − Ω_m.

**Replacement of Λ in the field equations.** SEDE is *holographic dark energy*, not modified gravity:
the Einstein tensor and Newton's constant are untouched, and the deformation enters only a smooth
(c_s² = 1; §3.4) perfect fluid on the right-hand side,
$$G_{\mu\nu} = 8\pi G\,(T^{\rm (matter)}_{\mu\nu} + T^{\rm (DE)}_{\mu\nu}),\qquad
\rho_{\rm DE} = \tfrac{3}{8\pi G}\,\Omega_{\rm DE0}\, H H_0\, f_{\rm sat},\quad p_{\rm DE} = w(z)\rho_{\rm DE}.$$
Equivalently, Λ becomes a *running* scalar Λ_eff(z) = 3 Ω_DE0 H H₀ f_sat that equals the observed Λ
today and is negligible in the early universe. The choice of right-hand side (fluid) over left-hand side
(geometry) is not cosmetic: a Barrow-*modified-gravity* Friedmann equation would rescale the
early-universe expansion by (M_P/H)^Δ ~ 10⁶¹, a Big-Bang-nucleosynthesis catastrophe (§8, Result 11).
Keeping the deformation in the dark-energy fluid leaves the early universe exactly standard.

**Energy conservation.** Because SEDE is general relativity with a fluid source, the Bianchi identity
$\nabla^\mu G_{\mu\nu}=0$ forces the **total** stress–energy to be covariantly conserved,
$\nabla^\mu(T^{(\rm matter)}_{\mu\nu} + T^{(\rm DE)}_{\mu\nu})=0$, identically — there is no violation. At
the background level matter dilutes standardly, ρ_m ∝ a^{−3}, and dark energy is a **self-conserved
fluid** whose density evolves entirely through its own pressure work,
$$\dot\rho_{\rm DE} + 3H\,(1+w_{\rm DE})\,\rho_{\rm DE} = 0, \qquad
w_{\rm DE}(a) = -1 - \tfrac13\,\frac{d\ln\rho_{\rm DE}}{d\ln a}$$
(the w(z) of §3.2 and Fig. 4). There is no literal energy transfer from matter to dark energy: such a
transfer of order ρ_crit over a Hubble time would spoil ρ_m ∝ a^{−3} and with it BBN and the CMB. This is
the crucial accounting point of the model. The **magnitude** of ρ_DE is anchored to the horizon
gravitational scale ∼ M_P² H² ∼ ρ_crit through the CKN bound and the flatness normalisation (§8.2),
not to energy borrowed from structure: the gravitational
binding energy released as halos virialise (E_bind ∝ M^{5/3}) is ∼ σ²/c² ∼ 10⁻⁶ of the collapsed
rest-mass energy, six orders of magnitude too small to *be* the dark energy. This makes the genuinely
new content precise: with the area law (Δ=0) the ratio ρ_DE/ρ_crit ∝ H^{−Δ} is a fixed fraction that
merely tracks — GR's own horizon energy repackaged (Jacobson) — whereas the Barrow deformation Δ>0
makes it grow as H falls, turning a tracking term into dark energy that comes to dominate; the new
physics is *entirely* the deformation ("is the energy new?" ⟺ "is Δ≠0?"). What structure formation
supplies is therefore not energy but **entropy bookkeeping**: by the apparent-horizon first law (A.1),
the deposited binding entropy ΔS_AH = E_bind/T_AH fixes the *weight* p = 5/3 that sets the shape γ of the
gate f_sat (§3.1) — it determines *when* and *how sharply* dark energy turns on, while the overall scale,
and the normalisation f_sat(0) = 1, come from the horizon (§3, §8). Read this way, SEDE's crossing of
w = −1 is **not a phantom field** — there is no negative kinetic term or gradient instability — but the
pressure work of a horizon-sourced fluid whose density rises as f_sat grows (Result 13: 1 + w is the horizon
entropy-production rate). (A single canonical or k-essence field cannot cross w = −1 without its kinetic
term ρ + p = 2X·P_X passing through zero — the Vikman no-go — so a fundamental-field description would
require a ghostly quintom; SEDE avoids this because its phantom is *effective*: ρ_DE = T_AH·s_grav is
built from covariant horizon scalars, not a fundamental scalar.) Finally, since f_sat is a functional of the *background* growth (c_s² = 1,
§3.4), the coupling lives at the background level only and introduces no local energy non-conservation in
the perturbations.

![Fig1](../output/paper/fig1_mechanism.png)
**Figure 1.** The mechanism. *Left:* the effective activation fraction (structure-growth gate) f_sat(z) rises toward 1
today as structure forms, activating the dark-energy gate — the effective running term Λ_eff = 3Ω_DE0 H H₀ f_sat
*fades in* (computed from the model's self-consistent linear growth factor D(z), §2). *Right:* the cosmic density inventory; SEDE's
ρ_DE is orders of magnitude below matter and radiation at high z and falls below ΛCDM's constant ρ_Λ,
so the early universe (BBN, recombination) is standard — the deformation acts only on the late-time
dark-energy fluid.

---

## 3. The fixed inputs — why the dark sector adds no fitted parameter

In flat ΛCDM the dark-energy amplitude Ω_Λ0 = 1 − Ω_m − Ω_r is itself fixed by the flatness closure,
not an independent fit; SEDE keeps exactly this — no independent dark-energy amplitude beyond flatness
— and differs only in that the *redshift dependence* is set by the horizon/growth prescription rather
than held constant. What that prescription fixes is the triplet (γ, w(z)-shape, λ) plus the sound
speed c_s² — and each is **fixed by a stated theoretical choice** (not fitted to cosmological data),
so SEDE adds no free dark-energy parameter. We are explicit about the status of each choice below
(some follow from flatness or the model construction; others, notably γ and c_s², are fixed
prescriptions). The choices separate along the two scales of §2. The **magnitude** of dark energy is
the horizon scale, ∼ M_P²H² ∼ ρ_crit (the CKN bound normalised by flatness, §8.2) — not a fitted
dark-sector parameter. What the choices control is the **shape**: *when* and *how sharply* that energy
turns on. We therefore lead with the structure
coupling γ — the shape, set by the binding energy — and then give w₀, λ, and c_s².

### 3.1 The structure coupling γ ≈ 1.5 from a halo-binding-energy derivation (Theorem 4C)
γ sets the sharpness of the structure gate, γ = d ln Σ_S / d ln σ8² with Σ_S the entropy-weighted
collapsed mass, Σ_S = ∫ M^p (dn/d ln M) d ln M — the derivative is with respect to the **variance**
σ8² because the gate's argument is x = D² = σ8²(z)/σ8²(0) (§3.4, A.2); equivalently γ = ½ d ln Σ_S/d ln σ8. The weight exponent p is fixed by *where* the
structure-formation entropy is deposited — and the reservoir is fixed *before any data enter*, by the
model's own §2 identity. Because SEDE identifies ρ_DE = T_AH·s_grav (§2), f_sat is **by construction**
the cosmic-*horizon* entropy activation fraction (A.1), so the binding energy released by virialisation
is accounted in the **horizon** ledger at its temperature T_AH ∝ H — *not* the halo's internal virial
temperature T_vir ∝ M^{2/3}. With T_AH mass-independent at a fixed epoch, the Clausius relation gives
ΔS_AH = E_bind/T_AH ∝ E_bind, and the virial binding energy E_bind ∝ M^{5/3} then fixes **p = 5/3**.
This reservoir choice is forced by the §2 identity; we flag explicitly that an earlier version used
T_vir (giving p = 1), which was a **reservoir error corrected for consistency with §2**, not a tuning —
though we acknowledge it is a discrete modelling revision, and that the corrected value happens to be
the one favoured in diagnostic fits.

γ has a transparent closed form. A self-similar reduction gives **γ = (p − 1)⟨1/α⟩_Σ** (α ≡ −d ln σ/d ln M
the effective spectral slope, ⟨·⟩_Σ the entropy-weighted average), which reproduces the full
mass-function integral to **0.1%** (A.3) — so the dark-sector input is *only* p = 5/3, with ⟨1/α⟩
ordinary halo physics. The *numerical* γ then follows from a standard, **pre-specified** halo model and
inherits its systematics: a Sheth–Tormen mass function over 10¹⁰–10¹⁶ M⊙ with an EH98 transfer gives **γ = 1.50**
(Result 4C), and varying the mass function (Press–Schechter, Tinker08, Despali16, Watson13), mass
range, and transfer function spans **γ ∈ [1.14, 1.60]** (median 1.50; `run_gamma_systematics.py`),
dominated by the mass-range low-cutoff and robust (1.48–1.57) to the mass-function and transfer choice.
We therefore claim precisely: the reservoir — hence the weight p = 5/3 — is fixed by the model's own
construction; the numerical γ ≈ 1.5 carries standard halo-model systematics of order ±0.2, **none
tuned to cosmological data**. That it lands near the diagnostic-fit value is then a *consistency check*,
not the motivation. (We do **not** claim a demonstrated microphysical transport channel depositing the
binding heat onto the horizon; the prescription is the accounting, and the smooth-DE closure of §3.4
ensures it induces no local dark-energy perturbation.) Binding energy thus sets the gate *shape*, not
the magnitude (§2, §10 ledger).

### 3.2 The equation of state w(z) from the fixed background, conservation, and flatness
The model has a single physical equation of state — that of the canonical fixed-point background
(eq. 1) — obtained from energy conservation,
$$w(z) = -1 - \tfrac13\,\frac{d\ln\rho_{\rm DE}}{d\ln a},$$
with ρ_DE = Ω_DE0 ρ_crit,0 E f_sat and the present-day normalisation fixed by spatial flatness
(Ω_DE0 = 1 − Ω_m − Ω_r, f_sat(0) = 1). Evaluated on the self-consistent background this gives
**w₀ ≈ −1.0**, **crossing −1** at low redshift (z ≈ 0.2) and dipping to ≈ −1.08 by z ≈ 2 — DESI DR2's
thawing/crossing quadrant. So w(z) is fixed by Ω_m (through the background) with no free EOS
parameter; this is the w(z) used throughout (§5, Fig. 4). Three w₀ values appear in this paper for *different* objects; to avoid confusion:

| object | scaling | role | w₀ |
|---|---|---|---|
| **canonical SEDE** | ρ_DE ∝ H f_sat | the actual model | ≈ −1.0, crosses −1 |
| λ=1 "temperature-factor" cousin | — | closed-form structural estimate (side note) | ≈ −0.85 |
| area-law branch | ρ_DE ∝ H² f_sat | disfavoured diagnostic alternative (§5.6, §8.1) | ≈ −0.49 |

Only the first row is the model's prediction; the −0.85 is a sanity check on the Ω_m-dependence
(closed form w₀ = (4Ω_m/3 − 1)/(1 − Ω_m)) and the −0.49 is the Δ=0 alternative the data disfavour.

### 3.3 The H-coupling λ = 0.5 from the horizon's volume entropy
λ fixes the *magnitude*'s scaling with H, with no free parameter: the gravitational entropy is taken
to be the coarse-grained entanglement entropy in the apparent-horizon **volume**, S_grav ∝ V_AH ∝ R_AH³,
i.e. a *constant* horizon entropy density. The conjugate horizon-fluid ansatz (ρ_DE = T_AH·s_grav,
T_AH ∝ H) then gives ρ_DE ∝ H ≡ H^{2λ} with **λ = 1/2** — no deformation parameter. Equivalently, in the
Barrow language S = (A/A₀)^{1+Δ/2} this is the maximal deformation Δ = 1 (s_AH ∝ H^{1−Δ}, λ = 1−Δ/2):
the volume law *is* Δ=1 (A^{3/2} = R³ = V). We keep Δ only as a **falsifiable test** of the volume
law (its data measurement, §5/§6), not as a model input; §8 motivates the volume/space-filling state
(capping Δ≤1 and postulating saturation). (Verified Δ-free ≡ Δ=1 bit-for-bit;
`friedmann.E_SEDE_volume`, `run_volume_equiv.py`.)

### 3.4 The sound speed c_s² = 1 from a background functional
SEDE's dark energy is *background-growth-gated* (its f_sat depends on the growth history), which might
suggest it clusters with structure. We adopt the **minimal smooth-dark-energy closure**: f_sat is a
functional of the **background** growth factor D(z) — a function of cosmic time only — so the structure
coupling is **background-level / nonlocal / coarse-grained, not local**. Hence ρ_DE is homogeneous by
construction, ∂f_sat/∂(local δ_m) = 0, and the dark energy has the standard metric-driven response with
**c_s² = 1**. We are explicit about the wavelength scope: this is a statement about the **sub-horizon**
response — there is no small-scale dark-energy clustering — and it is what the c_s² = 1 validation
below tests. A distinct *long-wavelength, separate-universe* modulation of the background gate by a
super-survey mode is a separate, super-horizon effect (the optional near-critical layer, deferred to a
companion paper, §6); it is not sub-horizon clustering and is not probed by the validation here.
Operationally SEDE is implemented as an **effective smooth dark-energy fluid using a PPF crossing
prescription** (which handles the w = −1 crossing; §2's Vikman-no-go argument is why a fundamental
single perfect fluid would *not* cross −1 without extra structure, and why SEDE's phantom is
*effective*). This is a modelling closure, not a theorem; a full CLASS Boltzmann treatment confirms the
resulting fσ8(z) matches this smooth-dark-energy growth to 0.1% (§4).

| input | value | status — fixed by |
|---|---|---|
| λ (H-coupling) | 0.5 | **postulate** — volume-law horizon entropy S_grav ∝ V_AH (≡ Barrow Δ=1), §3.3, §8 |
| w(z) (EOS) | w₀ ≈ −1.0, crosses −1 | **derived** — spatial flatness + eq. 1, §3.2 (no free EOS parameter) |
| γ (structure gate) | ≈ 1.50 (range [1.14,1.60]) | **fixed prescription** — halo-binding-entropy weighting (fixed, not fitted), §3.1 |
| c_s² (sound speed) | 1 | **closure** — minimal smooth-DE treatment (validated numerically), §3.4 |

![Fig2](../output/paper/fig2_inputs.png)
**Figure 2.** The two fixed dark-sector couplings. *Left:* the H-coupling λ = 1 − Δ/2 as a function
of the Barrow deformation; the volume-law endpoint Δ = 1 (the postulate of §8) fixes λ = 0.5. *Right:* the
structure coupling γ = d ln Σ_S/d ln σ8² (variance convention, §3.1) as a function of the entropy-weight exponent p; the
binding-energy-into-the-horizon argument (Result 4C / Prescription 4C, A.3) fixes p = 5/3, giving γ = 1.50 — close to the
value favoured in diagnostic fits. Neither is a fitted parameter.

---

## 4. Data & methodology

### 4.1 Data
We use a standard late-time-plus-compressed-CMB compilation, with both models fit to the identical
data vector:
- **BAO:** DESI DR2 (galaxy, QSO, and Lyα tracers, z = 0.3–2.33), as D_M/r_d and D_H/r_d.
- **Supernovae:** Pantheon+ as the baseline, with DES-5YR and Union3 used independently for the
  robustness cross-check of §5.4.
- **Cosmic chronometers:** the Moresco et al. H(z) compilation (model-independent expansion rates).
- **Redshift-space distortions:** the Gold-2018 fσ8(z) growth compilation (a standard, vetted set
  chosen for comparability with prior growth analyses; combining it with DESI DR2 BAO follows common
  late-time practice, and the audit below corrects a few mis-entered points), with the data audit of
  App. C (a small number of mis-entered points corrected to their published values).
- **CMB:** the compressed Planck distance priors (shift parameter R, acoustic scale l_A, ω_b), with
  the full plik_lite likelihood used as a cross-check.
- **Local distance ladder:** the SH0ES H₀ prior.
- **Out-of-sample (§5.3):** the pre-DESI eBOSS DR16 full-shape BAO, held entirely out of the main fit.

### 4.2 CAMB-in-the-loop: a physical sound horizon
The single most important methodological point is that the sound horizon is **not** a free parameter.
Earlier exploratory fits that floated r_d found large apparent preferences for SEDE; those did not
survive a physical sound horizon. Here, at every point in the chain, r_drag, r_star, and θ* are
computed self-consistently from the background by CAMB ("CAMB-in-the-loop"), so the BAO and CMB
distances are tied to the same early-universe physics for both models. This is what makes the ΔDIC of
§5 a like-for-like comparison rather than an artifact of a free standard ruler. The CMB enters as the
compressed (R, l_A) distances computed on the SEDE background through the same CAMB call. *Scope of
the compressed likelihood, and its full-likelihood validation:* the compressed priors capture the
background distance to last scattering (and hence the early-dark-energy fraction that drives the Δ
leverage of §5.6), but not the full acoustic-peak/polarization or lensing/ISW information. We therefore
**validate the headline against the full primary CMB likelihood**: a CAMB-in-the-loop fit on the full
Planck plik_lite TTTEEE + lowℓ TT/EE, with the five cosmological parameters and the Planck nuisance
*all free*, gives **Δχ²(SEDE−ΛCDM) = −0.43** — the preference sitting in the high-ℓ peak structure
(−0.37), *negative* (mildly favouring SEDE), **not** the large positive Δχ² that the holographic-DE
early–late tension would produce (§5.1, §9; `run_full_cmb_mcmc.py`). This directly tests the
acoustic-peak channel the compression omits, and we have confirmed it against the **full non-lite plik
with all 15 foregrounds sampled** (Δχ² = −0.33, consistent with the lite −0.43; §9). *Stated scope of
what remains:* only the ACT DR6 *primary* multifrequency likelihood is not run (ACT *lensing*, the
DE-sensitive channel, is, evaluated at the primary best-fit; §5.1). It is **formally deferred** under a
pre-registered trigger (PREREGISTRATION.md §E): the `act_dr6_mflike` likelihood and its multi-GB data
products are not installed, and with three concordant CMB datasets — plik_lite, full non-lite plik, and
ACT DR6 lensing — all at |Δχ²| < 1, a fourth primary likelihood is not expected to overturn them; the
committed (falsifiable) expectation is |Δχ²(SEDE−ΛCDM)| < 1, a large positive value being a genuine
failure. The %-level metric effects (§4.4) are reported but not used in the headline fit.

### 4.3 Inference and model comparison
Parameters are sampled by MCMC, marginalising over Ω_m, H₀ (with the SH0ES prior), ω_b, the
amplitude/σ8, and r_d-determining physics — the **same five** for SEDE and ΛCDM, since SEDE's dark
sector is fixed (§3). Model comparison uses the deviance information criterion (DIC), which penalises
the effective number of parameters p_D; because p_D is nearly equal (numerically close) for the two models, ΔDIC reduces
essentially to the χ² difference. Since DIC alone is distrusted for non-Gaussian posteriors, we report
the **AIC and BIC as cross-checks** (`run_info_criteria.py`): with k = 5 fitted parameters for *both*
models (Δk = 0), AIC = χ² + 2k and BIC = χ² + k ln N (N = 1628) give **ΔAIC = ΔBIC = ΔDIC = Δχ² = −2.95**
(−3.17 for the full-CMB joint) — the criteria coincide because the comparison is **penalty-free**, the
preference is not bought with complexity. We also compute a genuine **Bayesian evidence** by nested
sampling (dynesty, `run_bayesian_evidence.py`), integrating the joint likelihood over the shared prior:
ln Z = −737.6 (ΛCDM) and −735.7 (SEDE), giving **ln B(SEDE−ΛCDM) = +1.89 ± 0.33** — *positive but weak*
on the Jeffreys/Kass–Raftery scale, favouring SEDE, with the Occam factor explicit (and ≈ 0, since the
prior volume is identical). The evidence slightly *exceeds* ½Δχ²_min (1.48), i.e. SEDE is well-supported
across the prior, not only at the best fit. With the contested SH0ES prior dropped the evidence falls to
**ln B = +0.95 ± 0.41** — still favouring SEDE but formally *inconclusive* (|ln B| < 1) — so about half
the evidence is carried by SH0ES and the rest by the CMB distance (§5.4 ii′). Significance is then
assessed not by any single criterion but by the false-preference calibration of §5.2.

### 4.4 Perturbations: smooth dark energy, validated by CLASS
SEDE's dark energy has c_s² = 1 and does not cluster (§3.4), so its growth is captured analytically by
a smooth-dark-energy treatment. We validated this against a full Boltzmann computation: implemented in
CLASS as a PPF fluid (which handles the w = −1 crossing), the full-perturbation fσ8(z) matches the
analytic background-only growth to **0.1%** over z = 0.3–1.5. The smooth-DE closure is therefore
numerically validated for the observables used here, and the metric-level effects (CMB lensing
C_ℓ^φφ, low-ℓ ISW, S8) are all at the %-level relative to ΛCDM.

## 5. Results

### 5.1 The joint fit and the equation of state
In the marginalised joint analysis, SEDE is **mildly preferred** over ΛCDM, with a deviance
information criterion difference ΔDIC ≈ −2.9 in SEDE's favour. Because the dark sector adds no fitted
parameter (§3), the two models share the same five parameters and the same effective complexity,
p_D(SEDE) ≈ p_D(ΛCDM): the preference is in fit, not bought with parameters. **The headline survives the
full primary CMB likelihood**, not only the compressed priors: a marginalised CAMB-in-the-loop fit on
the full Planck plik_lite TTTEEE + lowℓ TT/EE (all five cosmological parameters and the Planck nuisance
free) gives Δχ²(SEDE−ΛCDM) = −0.43, the preference carried by the high-ℓ acoustic peaks (−0.37) with the
low-ℓ ISW-sensitive channel flat (≈ −0.03) — i.e. *mildly favouring* SEDE rather than the large positive
Δχ² the standard holographic-DE early–late tension would produce (§4.2, §9). More decisively, **folding
the full primary CMB *into* the joint** (replacing the compressed R, ℓ_A priors with plik_lite TTTEEE +
lowℓ, profiled over all parameters; `run_joint_fullcmb.py`) leaves the preference intact and marginally
stronger: joint Δχ²(SEDE−ΛCDM) = **−3.17** (ΛCDM 2443.99 → SEDE 2440.83), versus ≈ −2.9 with compressed
CMB — so the headline is **not an artifact of the CMB compression**. The physical reason the CMB poses
no obstacle is that SEDE carries *less* early dark energy than ΛCDM and deforms only the late-time
expansion (Figure 5): there is no early–late wall of the kind that disfavours standard holographic DE.
The recovered background parameters are physical and close to ΛCDM (Ω_m ≈ 0.30, H₀ ≈ 68.8, σ8 ≈ 0.76). The absolute
χ²/dof ≈ 3 (§5.5) reflects the heterogeneous, partly-compressed compilation (and the standard-but-not-
goodness-of-fit dof definition); it should *not* be read as a goodness-of-fit probability — we use
only matched-model differences (ΔDIC, Δχ²) on the identical data vector for inference. The canonical (λ = 1/2)
equation of state **crosses w = −1**, with a present-day value near −1 and a dip to ≈ −1.08 by z ≈ 2,
i.e. (w₀, w_a) ≈ (−0.98, −0.11) — squarely in DESI DR2's evolving-dark-energy (thawing/crossing)
quadrant, ~2.7σ from the ΛCDM point. The crossing is not imposed; it is where the expansion term and
the horizon-entropy-production term in 1 + w balance (Result 13, §6).

![Fig3](../output/paper/fig3_datafit.png)
**Figure 3.** Matched residual comparison (relative to ΛCDM). *Left:* DESI DR2 BAO residuals relative to the ΛCDM best fit, in
units of σ — the data (slate) with the SEDE prediction (red); SEDE is consistent with the measurements
and sits between ΛCDM and the points that pull away from it. *Right:* Pantheon+ supernova residuals
(binned) versus ΛCDM, with the SEDE shape overplotted. SEDE and ΛCDM are close at present precision; the
preference is a multi-probe aggregate (§5.4), not one decisive measurement.

![Fig4](../output/paper/fig4_eos.png)
**Figure 4.** The equation of state. *Left:* the canonical SEDE-H fluid w(z) (red) crosses −1 near z ≈ 0
and dips into the phantom regime (≈ −1.08 by z ≈ 2), the same sense as DESI DR2's CPL fit (green) and
distinct from ΛCDM (dashed). *Right:* the (w₀, w_a) plane — SEDE sits between ΛCDM and the DESI DR2
w0waCDM contour, in the evolving-dark-energy quadrant. The crossing is not imposed; it is where the
expansion and horizon-entropy-production terms in 1 + w balance (Result 13).

![Fig5](../output/fig_cmb_earlyde.png)
**Figure 5.** Why SEDE survives the CMB where standard holographic dark energy does not — the physical
content of the full-CMB result (full non-lite plik Δχ² = −0.33, plik_lite −0.43, ACT DR6 lensing −0.31;
all |Δχ²| < 1). *Left:* the early-dark-energy fraction Ω_DE(z). SEDE's structure gate (f_sat → 0 at high
z) holds it *below* ΛCDM (Ω_DE(z=3) = 0.029 < 0.035 at Ω_m = 0.30), whereas un-gated standard Barrow-HDE
carries substantially more early DE — the runaway that drives the early–late tension of refs. [Wu et al.
2509.02945]. *Right:* the fractional expansion deviation H/H_ΛCDM − 1. SEDE's deviation is sub-percent and
confined to z ≲ 3, decaying to zero through BBN and recombination, so the sound horizon, the acoustic-peak
positions, and the distance to last scattering are all essentially unchanged; standard Barrow-HDE deviates
over a far wider range. SEDE has *less* early dark energy than ΛCDM, not more — hence no early–late wall.
(`run_cmb_earlyde.py`.)

### 5.2 Is the preference real, or model flexibility?
A ΔDIC ≈ −3 is suggestive but not, on its own, decisive — and a structure-gated model could in
principle be winning by flexibility rather than by capturing a real signal. We test this directly with
a **false-preference calibration** (a parametric bootstrap). We generate mock datasets from a ΛCDM
truth — re-drawing *every* probe self-consistently, crucially including the compressed CMB distances
(R, l_A) and the SH0ES prior, not just the late-time data — and refit both models to each mock. Under
ΛCDM truth a flexible model would still tend to "win"; SEDE does not: the null distribution of
Δχ²(SEDE − ΛCDM) has mean **+0.42** (SEDE if anything fits ΛCDM-truth mocks slightly *worse*,
consistent with having no spare flexibility). The real data give Δχ² = −2.96, beyond essentially all
of the null draws: in a 2000-mock run (null mean +0.37, s.d. 1.25), **8 mocks in 2000 are as
SEDE-preferring as the real data, p ≈ 0.004 (~2.6σ), 95% upper limit 0.006**. (An earlier 500-mock run
gave p ≈ 0.002 but with only one mock at the threshold — resolution-limited; the 2000-mock value is the
resolved, slightly more conservative number.) Rerun with the contested SH0ES prior removed from both the
real fit and the mocks, the real value moves to Δχ² = −1.75 and the calibration gives **p = 0.030
(~1.9σ, 95% UL 0.050)** — the preference degrades but stays beyond the null, now carried by the
CMB-distance channel alone (§5.4 ii′).

Two caveats. First, an earlier version of this test held the CMB distances *fixed*
across mocks; that biased the null mean to −2.87 and gave p ≈ 0.5 — the artifact that made the
preference look like pure flexibility. Re-drawing the CMB self-consistently is the correct procedure
and is what yields p ≈ 0.004. Second, the significance is **CMB-distance-treatment-dependent**: an
independent sibling pipeline, with a different l_A/CAMB handling, obtains a milder real Δχ² ≈ −1.98
and p ≈ 0.025 (~2σ). Both analyses agree the preference is **not obviously attributable to model
flexibility**, and trace the 2σ–3σ spread to the CMB-distance engine rather than to SH0ES. We
therefore quote the preference as a **suggestive ~2–3σ**, pending independent reproduction.

![Fig6](../output/fig_F5_calibration.png)
**Figure 6.** The false-preference calibration. Histogram: the null distribution of Δχ²(SEDE − ΛCDM)
from ΛCDM-truth mocks (all probes, including the compressed CMB and SH0ES, re-drawn self-consistently),
centred near zero (mean +0.37) — SEDE has no flexibility advantage under ΛCDM truth. Red line: the real
data, Δχ² = −2.96, in the far tail. In the 2000-mock run only 8 mocks are as SEDE-preferring as
the real data, p ≈ 0.004 (~2.6σ, 95% UL 0.006); the preference is not obviously attributable to flexibility.

### 5.3 Out-of-sample check
A further check against a flexibility explanation is out-of-sample prediction. We trained SEDE on
**real pre-DESI data only** (eBOSS DR16 full-shape BAO, LRG/QSO) and used it to predict the later,
held-out DESI DR2 BAO: in this split SEDE matches the held-out DESI DR2 trend **better than ΛCDM by Δχ² = −8.2**; a
low-z→high-z split of the same kind gives −9.9. A model winning purely by flexibility would tend to do
*worse* out of sample, so this is encouraging — i.e. a pre-DESI-trained SEDE realisation **predicts
the DESI DR2 trend** in this split.
*Two clarifications so this number is read correctly.* (a) **Not comparable to the in-sample ΔDIC.** The
−8.2 is the Δχ² on the **held-out 13-point DESI BAO subset alone**, whereas the in-sample ΔDIC ≈ −2.95 is
the net over the **full 1628-point compilation** (most probes of which do not discriminate; §5.4); the
two are different objects on different data, so "OOS −8.2 exceeds in-sample −3" is *expected*, not
anomalous (a single-probe Δχ² is naturally larger than a diluted whole-compilation total) and is **not**
evidence of over-fitting (which would push OOS the other way). (b) **Direction-dependent.** The signal
lives in the low-z → high-z direction (the DESI evolving-DE direction, Δχ² = −9.9); the reverse split
(high-z → low-z) is flat and marginally favours ΛCDM (+0.3), so part of the OOS strength is that
eBOSS-trained parameters land where DESI's specific tracers sit — partly favourable, to be independently
reproduced before strong weight is placed on it.

![Fig7](../output/paper/fig5_oos.png)
**Figure 7.** Out-of-sample test. BAO residuals relative to ΛCDM (in σ) for the **training** set
(pre-DESI eBOSS DR16, purple) and the **held-out** DESI DR2 data (slate), with the SEDE prediction
(red) from parameters fit to eBOSS alone. SEDE, trained without ever seeing DESI, predicts the held-out
DESI BAO better than ΛCDM (Δχ² = −8.2) in this split — encouraging, and to be independently reproduced.

### 5.4 Robustness
The preference is not a single-dataset artifact. (i) **Supernova-set independence:** Δχ²(SEDE − ΛCDM)
is −2.95 / −2.87 / −2.97 with Pantheon+ / DES-5YR / Union3 respectively — not a Pantheon+ effect.
(ii) **Per-probe decomposition (where the preference comes from).** Decomposing the joint
Δχ²(SEDE − ΛCDM) = −2.95 into per-probe contributions *at the joint best-fit*
(`run_probe_decomposition.py`) is more informative than the leave-one-out, and we report it plainly:
the preference is carried by the **geometry channels** — the compressed CMB distance (R, ℓ_A)
contributes **−1.73** and the SH0ES H₀ prior **−1.25** (together ≈ 100% of the total) — while DESI BAO
(−0.14), fσ8 (−0.03) and ω_b (−0.07) are nearly flat and Pantheon+ SN (+0.16) and cosmic chronometers
(+0.11) marginally favour ΛCDM. So this is **not a broad multi-probe consensus**: SEDE's evolving w(z)
fits the distance to last scattering better and accommodates a marginally higher H₀ (68.9 vs 68.7) that
eases SH0ES — a *paired* (CMB-distance + H₀) effect. We show the decomposition precisely because the
**leave-one-probe-out is blind to this paired structure** (dropping either member leaves the other, so
every single-probe holdout stays negative: range [−3.06, −2.16], mean −2.64) — the holdout robustness is
real but does not by itself establish a consensus. Two honest qualifications follow. First, the
CMB-distance contribution is *not* a compressed-prior artifact: folding the full primary CMB into the
joint **strengthens** it (Δχ² = −3.17, §5.1). Second, SEDE does **not** resolve the Hubble tension — its
best-fit H₀ ≈ 69 still sits ~4σ from SH0ES (§7, §9) — but its geometry does ease the SH0ES prior by
Δχ² ≈ 1.2, so part of the preference is a mild H₀ accommodation, which we flag rather than hide.
(ii′) **No-SH0ES robustness (the contested prior removed).** Because the SH0ES H₀ prior is the field's
most disputed input — and DESI's own baseline excludes it — we recompute all three headline statistics
with SH0ES dropped (`run_no_shoes_robustness.py`). The preference **weakens but survives, and what
remains sits in the harder channel.** The joint Δχ²(SEDE − ΛCDM) = **−1.75** (from −2.95), now carried by
the compressed CMB distance (−1.37) and DESI BAO (−0.45), with the best-fit H₀ falling to **68.2** once it
is no longer pulled toward 73 — i.e. the CMB-distance improvement is genuine, not an H₀ accommodation. The
nested-sampling evidence falls to **ln B = +0.95 ± 0.41** — still favouring SEDE but now *inconclusive*
(|ln B| < 1) on the Jeffreys scale — and the false-preference calibration (§5.2) gives **p = 0.030**
(6/200 ΛCDM-truth mocks beyond the real value; ~1.9σ, 95% UL 0.050) versus p ≈ 0.004 (~2.6σ) with SH0ES.
So roughly 40% of the χ² preference and about half the Bayesian evidence are carried by SH0ES; the
remainder is **modest but robust and lives in the CMB-distance channel rather than the contested local
prior.** We read this as sharpening, not weakening, the paper's "suggestive, not established" posture:
with the one disputed prior removed, SEDE remains *mildly* preferred (~1.9σ) by hard geometric data, and
no part of that residual depends on resolving the Hubble tension.
(iii) **Code independence:**
the SEDE expansion history reproduced with CLASS matches CAMB to 5×10⁻⁵ in r_drag (≪ 0.2%) — not a
CAMB artifact. (iv) **Growth:** the full-Boltzmann γ_growth ≈ 0.55 confirms standard-GR dark energy,
not modified gravity. (v) **Tracer-level robustness:** recent analyses argue the DESI DR2 evolving-DE
signal is driven mainly by the LRG1 (z≈0.51) and LRG2 (z≈0.71) BAO bins [arXiv:2508.10514, 2504.15222].
Dropping these bins individually and together (`run_tracer_loo.py`) leaves SEDE's preference intact —
Δχ²(SEDE − ΛCDM) = −3.27 (−LRG1), −2.76 (−LRG2), and **−3.00 with both removed** (versus −2.95 for the
full set), negative in every tracer holdout (range [−3.27, −2.76]). SEDE's mild preference is therefore
*not* the LRG1–2 artefact the steep-CPL signal is attributed to: it is a geometry-channel preference
(CMB-distance + H₀, §5.4(ii)) that does not rest on the contested bins. The dissociation is quantitative: the ΛCDM-deviation measured by the
ratio(ω_m) early-vs-late consistency metric drops from 2.6σ to 1.2σ when LRG1–2 are removed
[arXiv:2407.04385], whereas SEDE's Δχ² is essentially unchanged (−2.95 → −3.00) — SEDE's preference and
the LRG-driven steep-CPL deviation are *different objects*. This is the natural answer to the post-DR2
robustness critiques (§5.7).

### 5.5 Evidential summary
SEDE is mildly preferred over ΛCDM (ΔDIC ≈ −2.9; calibration ~2–3σ); the preference is suggestive — not
obviously attributable to model flexibility, and supported by a pre-DESI → DESI out-of-sample check —
and the model adds no continuously-sampled dark-energy parameter. But it is modest: both models have
χ²/dof ≈ 3, so SEDE is "less wrong," not "right." This is stronger than indistinguishable, short of
established. DESI DR3 and Euclid provide the sharp future test (§6).

### 5.6 The deformation Δ from current data — a conditional diagnostic
The Barrow parameter Δ is a useful *test* axis: the volume postulate corresponds to Δ = 1, the area
law to Δ = 0. SEDE's *equation of state* is CPL-degenerate in shape (max|w − w_CPL| = 0.027), so a
diagnostic must use channels other than the EOS curvature. We stress at the outset that the numbers
below are **conditional diagnostics** (computed with the other cosmological and nuisance parameters
held at or near their best-fit values, not a full marginalised likelihood with Δ varied); they
indicate where current data sit, not a marginalised detection. Four channels:

| channel | what it uses | Δ̂ (conditional profile — **not a posterior**) |
|---|---|---|
| CMB shift parameter R (Planck) | early-DE fraction, Ω_DE(z_rec) ∝ H^{−Δ} | 0.98 ± 0.11 |
| lever arm: DESI BAO + Lyα(z=2.33) + R | ρ_DE/ρ_crit ∝ H^{−Δ} across the H-range | 1.03 ± 0.11 |
| DESI DR2 (w₀,w_a) quadrant (official chains) | EOS values + crossing | 0.93 / 0.83 / 0.90 (PantheonPlus/DESY5/Union3) |
| growth fσ8 (Gold-2018) — **orthogonal** | growth amplitude, not distance | 0.0 ± 0.61 |

> **These are not posterior constraints on Δ; they are profile/conditional diagnostics** — each Δ̂ is
> obtained with the other cosmological and nuisance parameters held at (or near) their best-fit values,
> *not* marginalised. The quoted ± is the conditional curvature width, **not** a marginalised error bar;
> read each row as "the data lean to Δ ≈ 1 in this channel," not as a measurement. The robust,
> marginalised Δ measurement awaits DESI DR3 + Euclid (§6).

These diagnostics **prefer the volume endpoint (Δ ≈ 1) and disfavour the area-law branch (Δ = 0)**,
consistent with the EOS argument of §8.1 (the area-law branch has w₀ = −0.49 and an outsized early-DE
fraction). The growth channel is weak but **orthogonal** (growth amplitude, not expansion) and is
consistent with the geometric value at 1.6σ — no growth–geometry tension.

The high apparent significance of the area-law disfavouring is a **conditional** statement and must
not be read as a marginalised detection: the CMB and lever-arm channels share the high-z expansion,
so their combination is an upper bound on the information, and it is conditional on Ω_m, H₀, r_d —
marginalising, the area-law signal is partly reabsorbed by shifts in those parameters. This is exactly
why the *joint, marginalised* model preference is only modest (ΔDIC ≈ −3, §5.1–§5.5), not
overwhelming; the two statements are not in tension once the conditional-vs-marginal distinction is
kept. The (w₀,w_a) channel also pulls slightly low because SEDE's gentle w_a under-matches DESI's
steeper preferred w_a — the same mild ~2.7σ EOS tension visible in the joint fit (Fig. 4). **Honest
summary: current late-time data prefer the volume endpoint and disfavour the area-law branch, but a
robust *marginalised* measurement of Δ awaits DESI DR3 and Euclid (§6).**

### 5.7 Relation to the post-DESI DR2 literature

The DESI DR2 BAO release [arXiv:2503.14738] sharpened interest in evolving dark energy, with a
preference for the (w₀ > −1, w_a < 0) quadrant; the subsequent literature is actively debating whether
that signal is physical or an artefact of dataset combination, tracer selection, or supernova
calibration [reviewed in arXiv:2602.05368]. Two threads bear directly on SEDE. First, **robustness
critiques** argue the steep signal reflects CMB/BAO/SN combination tension [arXiv:2504.15222] and is
driven mainly by the LRG1–2 BAO bins [arXiv:2508.10514, 2407.04385]. This *helps* SEDE rather than
threatening it: SEDE does not chase the steep w_a (it sits ~2.7σ from the DESI CPL point, §5.6), and
its preference survives a **tracer-level** leave-one-out that removes exactly those LRG1–2 bins (§5.4,
Δχ² = −3.00 with both dropped) — so SEDE realises the *conservative* reading the critiques call for,
rather than the fragile steep-CPL one. Second, **reconstruction and mechanism papers**: model-agnostic
Gaussian-process reconstructions place the w = −1 crossing at z_wt ≈ 0.46 (+0.24/−0.12)
[arXiv:2511.02220], versus SEDE's canonical z ≈ 0.19 — a ∼2.3σ offset (SEDE crosses more recently)
that is a genuine future discriminator, not yet a tension given the broad reconstruction errors. Among
mechanism competitors, data-driven analyses read the evolving signal as a *coupled* dark sector
(∼2.2σ DE–DM interaction at low z) [arXiv:2504.00985]; SEDE's CHR ledger is exactly such a coupling but
with the interaction *fixed* to the structure-growth rate, c_s² = 1, and instability-free (§7), rather
than a free coupling. And the **quintom** programme [arXiv:2505.24732, 2601.02284] rests on a *no-go
theorem* — no single canonical field or perfect fluid crosses w = −1, so viable models require two
fields, higher derivatives, modified gravity, or an effective-field-theory description. SEDE is the
last of these: it crosses w = −1 *effectively*, as a smooth horizon fluid in the GW-safe EFT corner
(α_T = α_M = 0, c_s² = 1; §7), avoiding the ghost/two-field machinery, and — unlike the two-field
quintom whose late attractor is phantom-dominated (a Big-Rip tendency) [arXiv:2601.02284] — SEDE's
future attractor is de Sitter (w → −1, with f_sat *saturating* to a finite ≈1.23 as growth freezes, not
diverging; §2, App. A.7, Result 12), so it does not predict a Big Rip.

Two further consistency points sharpen the post-DR2 placement. (a) **Turyshev's model filter.** The
DR2 status review [arXiv:2602.05368] screens evolving-DE models by *gravitational-wave propagation* and
*perturbation stability*, and stresses that the CPL preference is sensitive to redshift-dependent SN
calibration at the few×10⁻² mag level. SEDE passes the filter by construction (α_T = 0 ⟹ c_GW = c;
c_s² = 1 ⟹ gradient/ghost-stable; §7), and its SN-set independence (Δχ² = −2.95/−2.87/−2.97 across
Pantheon+/DES-5YR/Union3; §5.4) limits the calibration exposure. (b) **The r_d-independent shape test.**
Turyshev's calibration-free observable F_AP(z) ≡ D_M/D_H = (D_M/r_d)/(D_H/r_d) cancels the sound horizon
and is also H₀-independent, isolating the pure expansion shape E(z). On the six DESI DR2 tracer bins
that report both ratios (`run_fap_test.py`), SEDE's evolving-w shape gives χ²/n = 0.96 — consistent with
the data and indistinguishable from ΛCDM (Δχ² = −0.14) — with **no** reliance on r_d, H₀, the CMB sound
horizon, or SN calibration, the very systematics the robustness debate targets. (The LRG1 bin sits
∼1.8σ low for *both* SEDE and ΛCDM, i.e. it is a data feature, not a model failure.) SEDE is thus a
*physical* member of the post-DR2 evolving-DE family, distinguished by its mechanism (the
growth–expansion lock, §6), by passing the stability/GW filter, and by surviving the robustness and
calibration-free tests the steep signal may not.

## 6. Predictions & falsifiability

Because SEDE adds no fitted dark-energy parameter, it makes sharp predictions rather than fits. We
**pre-registered** these (a sha256-locked statement deposited ahead of the future data; App. D) to
forbid post-hoc tuning:

| quantity | SEDE prediction | discriminating from |
|---|---|---|
| Barrow deformation Δ | **1.0 ± 0.09** (forecast) | a smooth horizon Δ = 0 |
| equation of state (w₀, w_a) | ≈ (−0.98, −0.11), crossing −1 | the ΛCDM point (−1, 0) at ~2.7σ |
| w = −1 crossing redshift (consistency check, *not* a discriminator) | z ≈ 0.19 (genuine volume-law value; ∼2.2σ mild tension) | GP reconstructions z_wt ≈ 0.46 (+0.24/−0.12) [arXiv:2511.02220] |
| structure amplitude σ8 | ≈ 0.76 | (consistent with weak-lensing S8) |
| growth index γ_growth | ≈ 0.55 | modified gravity (γ_growth ≠ 0.55) |

**An honest near-term tension — and why z ≈ 0.19 is the genuine prediction, not a knob.** Canonical SEDE
crosses w = −1 at **z_cross = 0.195** (w₀ = −0.996), using the Gibbons–Hawking horizon temperature
T_AH = H/2π. This is **not** an arbitrary convention: T = H/2π is the de Sitter temperature appropriate
for a de Sitter-attractor dark-energy fluid, and we *tested* the alternative — applying the full
dynamical Cai–Kim factor (1−ε/2) to the volume-law fluid (`run_dynT_crossing.py`, the self-consistent ρ_DE ∝ H(1−ε/2)f
model). It **fails**: it has no stable self-consistent background, and in the stable perturbative limit
it gives w₀ = −1.15 and stays phantom (w < −1) at all z — **no crossing at all**, and in the *wrong*
DESI quadrant (w₀ < −1). The dynamical correction is largest in the matter era, exactly where f_sat → 0
and dark energy is negligible; forcing it where it does not belong destroys the late-time crossing. So
T = H/2π is the physically-correct choice for this sector, and z ≈ 0.19 is the genuine volume-law
prediction. (The unrelated z ≈ 0.59 is the *area-law* Δ = 0 model — ρ_DE ∝ H²(1−ε/2)f — a different
model, not a convention of this one.) Post-DR2 Gaussian-process reconstructions favour
z_wt ≈ 0.46 (+0.24/−0.12) [arXiv:2511.02220, 2504.00985]; canonical SEDE sits **∼2.2σ more recent** — a
*real, mild* tension (the reconstruction's lower error is tight), the one place current data pull against
canonical SEDE. We flag it as such, **not** as a discriminator: the crossing redshift is CPL-degenerate
(§6, below) and, read as a Δ-probe, would point the *wrong* way (toward Δ = 0), conflicting with the
early-DE diagnostics that favour Δ = 1 (§5.6). The model's actual test is the deformation amplitude Δ,
not the crossing; a DR3/Euclid measurement of the crossing is a useful consistency check, not the decider.

**The sharp future test.** The discriminating observable is the *amplitude* of the deformation Δ, not
the detailed shape of w(z) — indeed the SEDE w(z) is CPL-degenerate in shape to better than the DR3
precision (max |w − CPL| ≈ 0.027), so a CPL fit cannot distinguish the *mechanism*; Δ can. A Fisher
forecast in (Ω_m, Δ, ln A r_d) gives σ(Δ) ≈ 0.115 from DESI DR3 alone, 0.133 from Euclid, and
**σ(Δ) ≈ 0.087 combined** — an ~11σ separation of the volume-law horizon (Δ = 1) from a smooth one
(Δ = 0), under the Fisher assumptions. DESI DR3 and Euclid could therefore move SEDE substantially
toward confirmation or exclusion under those assumptions.

**The mechanism test is also forecast to be decisive — separately from the amplitude.** The σ(Δ) above
is the *geometry-side* deformation; the *distinctive* SEDE claim (dark energy sourced by structure) is
tested by measuring Δ on the **growth** side and checking it against the geometry side — the E1/P2
growth–expansion lock. A Fisher forecast for the growth-side deformation (`run_mechanism_forecast.py`)
gives σ(Δ)_growth = 0.25 from Euclid/LSST weak lensing alone, tightening to **0.13 (∼4σ)** when combined
with CMB-lensing, cluster abundance, RSD, and kSZ. Reframed as the parameter-free **null test**
r(z) ≡ ρ_DE^geom(z) − G[D^growth(z)] = 0 — dark energy a fixed function of the linear growth factor —
the test uses the full redshift shape with no dark-energy parameters; ΛCDM (ρ_DE = const, no
D-dependence) and Gough's information-DE (ρ_DE = G′[SMD], baryonic) violate it distinctly. So
Euclid/LSST-era growth data are forecast to **confirm or refute the structure-sourcing mechanism at
∼4σ**, not merely the EOS shape — closing the gap between "a w(z) with a thermodynamic story" and a
tested mechanism (§9). The exact data vector (the geometry leg DESI DR3 D_H/r_d→H(z) vs the growth leg
Euclid/LSST 3×2pt + CMB-lensing + clusters + RSD) and the **pre-registered decision rule** — CONFIRM if
\|Δ_grow − Δ_geom\| < 2σ with the null r(z) = ρ_DE^geom − G[D^grow] = 0 at χ²/dof ≲ 1; REFUTE if Δ_grow
is consistent with 0 at >3σ while Δ_geom = 1 — are locked in the pre-registration (PREREGISTRATION.md §F;
App. D), so the staged test carries no post-hoc freedom.

![Fig8](../output/paper/fig6_forecast.png)
**Figure 8.** The sharp future test. A Fisher forecast of the separation of the volume-law horizon
(Δ = 1) from a smooth one (Δ = 0): DESI DR3 reaches 8.7σ, Euclid 7.5σ, and the combination
σ(Δ) ≈ 0.087, an ~11σ separation. The mild present preference would become a strong test with the
next generation of surveys.

![Fig9](../output/fig_e1_mechanism.png)
**Figure 9.** The structure-sourcing mechanism on *current* data (the companion to the forecast of
Figure 8, and the paper's honest open frontier). SEDE's founding claim is that the activation fraction
f_sat read off the **expansion** (geometry leg: DESI D_H/r_d → H(z), blue) coincides with the one read
off **structure growth** (RSD fσ8 → σ8(z) → D(z), red); the black curve is the SEDE prediction both
should follow. The two legs **track at z < 1** but **diverge at z > 1** (the structure leg flattens while
the geometry/model declines). The significance is error-budget-dependent — χ²/dof = 0.40 in the
RSD-error budget (consistent) versus 1.40 in the precise-geometry budget (marginally disfavoured) — so on
present data the mechanism is **inconclusive and RSD-limited**: SEDE's EOS is favoured, but its
distinctive mechanism is not yet confirmed. Euclid/LSST weak-lensing tomography sharpens the structure leg
and decides this at ∼4σ (Figure 8; pre-registered decision rule, PREREGISTRATION.md §F).
(`run_e1_figure.py` / `run_e1_mechanism.py`.)

**The EOS as a meter of entropy production.** SEDE gives the equation of state a thermodynamic meaning
(Result 13): within the SEDE background ansatz, an exact identity 1 + w_DE = (1/3)[2λε − d ln f_sat/d ln a]
splits 1 + w into an expansion term (which pushes w > −1) and a horizon-entropy-production term (which
pushes w < −1). The phantom
crossing at z ≈ 0.2 is exactly where the two balance. Measuring w(z) thus probes the effective
activation rate of the horizon-fluid sector — the structure-growth-gated arrow of time made
observable in the dark-energy equation of state.

**Structure sourcing without an entropy-amplitude problem.** The identity above gives w(z) a
thermodynamic reading, and it invites the natural objection that a tiny structure-sourced entropy cannot
gate an order-unity dark-energy term — the deposited binding entropy is only
ε_dep ≈ Ω_m f_coll ⟨σ_v²⟩/c² ≈ 1.5 × 10⁻⁷ of the horizon entropy (the same six orders by which the
binding *energy* falls short of ρ_DE; §2). **The objection is a category error, and the amplitude
problem is *avoided*, not solved by amplification.** The activation coordinate is **not** the
deposited-entropy fraction ΔS_bind/S_AH; it is the structure variance x ≡ σ²(z)/σ²(0) = D²(z), an
**order-unity** quantity, which is the *argument* of the gate f_sat (the Press–Schechter collapsed
fraction; §3.1, App. A.2). The deposited binding entropy enters in a *different* place: it fixes the
entropy *weight* γ ≈ 1.5 — the *shape* of the kernel — through the heat ΔS = E_bind/T_AH delivered to
the horizon (§3.1, App. A.3). Hence ΔS_bind/S_AH ∼ 10⁻⁷ does **not** imply f_sat ∼ 10⁻⁷: they are
different objects. With the *argument* (variance, O(1)) and the *weight* (entropy, 10⁻⁷) kept distinct,
the gate activates simply because the structure variance grows by order unity — nothing is amplified and
nothing is tuned. This is a consequence of the fixed inputs of §3, **derived, not assumed.**

**The exact growth–expansion lock.** Write the gate as φ ≡ f_sat(x), x = D²(z), and define its derived
response to the activation coordinate,
χ_x ≡ ∂φ/∂x = d f_sat/dx = γ e^(−γx)/(1 − e^(−γ)),
which, once γ is fixed (§3.1), is an **order-unity** function with no free parameter — χ_x(x=0) ≈ 1.93,
χ_x(x=1) ≈ 0.43. The fluid equation of state then follows exactly from continuity (the
entropy-production identity above), using d ln φ/d ln a = (χ_x/φ)·x·(d ln x/d ln a) and
d ln x/d ln a = 2g with g ≡ d ln D/d ln a:
1 + w_DE = (1/3)[ 2λε − (χ_x/φ)·2x·g ],
verified against the direct fluid EOS to 3 × 10⁻³ (`run_chr_experiments.py`; suite tests V52–V53).
Because χ_x, x, and g are each derived or directly measured, the phantom term is **growth data**, not a
free function. This is the **w(z) ↔ σ²(z) lock**: given the measured expansion history, the dark-energy
equation of state and the structure-growth history are not independent; inverting the relation predicts
g(z) from w(z) to 0.6% within the model. A kinematic w(z) — including a CPL fit that reproduces SEDE's
w(z) to within the DR3 precision — obeys no such relation. **Expansion-only mimicry is therefore not
enough: SEDE requires growth–expansion phase-locking.** This is the test (P2 below) that separates the
structure-sourcing mechanism from a curve fit — w from BAO/SNe/CMB distances versus g, D, fσ8 from
RSD/weak lensing. It is parameter-free and uses nothing beyond §3; no criticality assumption enters.

**The observable form of the lock (P2), and a GR consistency null.** The lock derived above — the
phantom part of w(z) equals (1/3)(χ_x/φ)·2x·g reconstructed from growth data, obeyed by *no* kinematic
dark energy — is the model's one parameter-free, degeneracy-breaking prediction. Its established
observational realisation is the E_G estimator (the lensing-to-velocity amplitude ratio): SEDE is GR
with smooth dark energy, so it predicts no gravitational slip, E_G(z) = Ω_m,0/f(z), within <1% of ΛCDM
and consistent with the DESI-DR1 + weak-lensing measurement to z ≈ 1 [arXiv:2507.16098]
(`run_eg_test.py`). E_G is a *shared* GR null (with the GW-speed c_T = c and the
zero-slip nulls SEDE passes), not a SEDE-vs-ΛCDM discriminator — at 0.3% it sits far below the ∼5–20%
E_G errors; the model's only discriminator remains the deformation amplitude Δ (§5.6, §6). **Failure of
the lock (P2) would falsify the structure-sourcing mechanism itself** — SEDE would then be an ordinary
kinematic w(z) — so P2 is the test that separates the mechanism from a curve fit.

**A near-critical layer entailed by the volume law (predictions deferred to a companion paper).** The
*bare* gate is non-critical: it is the J = 0 limit of a Bragg–Williams area↔volume landscape, with a
finite, monotone susceptibility (χ_x ≈ 1.93 → 0.43), i.e. a smooth crossover. A genuine spinodal needs a
cooperative coupling J ≥ J_c = 4T between the horizon degrees of freedom — and that cooperativity is
**not** a new assumption. Volume-law entanglement (S ∝ V, the Δ = 1 postulate of §8.4) is realised by
maximally nonlocal, all-to-all connectivity (the SYK class / the nonlocal-connectivity requirement of
§8.5 and Appendix G), which is precisely the regime where mean
field is exact and the spinodal is a real, robust feature: the ordered phase exists for any coupling above
threshold and the critical susceptibility diverges (χ_peak ∝ √N → ∞), with no J-tuning, whereas the
area-law (local, short-range) class has no finite-temperature transition at all (`run_chr_soc.py`). Thus
the *existence* of the spinodal tracks Δ, and the same single postulate that fixes the background (horizon
volume law ⇒ Δ = 1) also supplies the cooperative criticality of the response sector. Self-organised
criticality fixes the operating point: the structure-growth drive is slower than the horizon scrambler
(SYK-saturated, ∼ Planckian) by ∼10⁵⁸, the slow-drive/fast-relaxation condition that holds the order
parameter on its spinodal, so f_sat = ½ at z* ≈ 1.2 (m = ½, the J_c point) is an attractor reached by the
GSL Δ → 1 flow (§8), not a tuned coincidence. At that operating point the response sector becomes
*near-critical* (the **Critical Horizon Response** regime) and acquires fluctuation-level signatures: a
localised enhancement of large-scale activity at z* (P3); a small, **nonzero**, redshift-localised
*separate-universe* (long-wavelength) modulation of growth/ISW peaked at z* (P4 — the kill test, distinct
from the *broad, sub-horizon* signal that clustering dark energy with c_s² < 1 would give, and from ΛCDM's
exact zero); and a critical-slowing transient in γ_growth (P5). These add **no fitted parameter** and now
follow from the volume-law postulate rather than a separate proximity-to-criticality hypothesis. Because
the horizon is driven *through* its spinodal at z* (not parked at a tuned critical point) with an
all-to-all coupling, the companion paper sharpens these into a specific signature: the z*-localised
fluctuations are scale-free in the **mean-field self-organised-criticality** class — a power-law amplitude
spectrum with avalanche exponent τ ≈ 3/2 and a 1/f temporal spectrum — which a single-scale clustering
dark energy cannot mimic. Because they are not required for — and do not affect — the core background
result, we develop them, with the separate-universe/super-sample-covariance formalism [arXiv:1711.07467]
and the supporting figures, in a **companion paper (in preparation)**. The corresponding kill test is
within reach of Euclid/Rubin weak-lensing tomography.

**Cross-horizon universality (speculative; companion).** Whether Δ = 1 extends to black-hole horizons —
giving S_BH ∝ A^{3/2} ∝ M³ and strongly altered primordial-black-hole evaporation — is a speculative
extension logically separate from the core model, which requires the deformation only in the dark-energy
sector. We do **not** adopt it: SEDE is consistent with the standard area law for astrophysical black
holes, and **on the evidence of this paper alone the black-hole-vs-cosmic-horizon asymmetry remains a
genuine conceptual cost (§8.5)**. The companion paper *offers a resolution* — a *state-dependent* count
in which a horizon is area-law when undriven (a black hole, in equilibrium) and volume-law when
sustainedly structure-driven (the cosmic horizon), which would make the asymmetry a predicted feature
rather than a defect — but that reframing **rests on a conjecture** (developed there); absent it, the
cost stands. The black-hole side is at least derivable: the structure entropy deposited at a black hole
is negligible against the entropy it *creates* (S_rad/S_BH ∼ (M/M_P)^{−1/2}), so no black hole — not even
a primordial one built entirely by collapse at formation — activates volume. Either way the asymmetry
sharpens a falsifier: an isolated *volume-law* black hole would imply a universal deformation. We pursue
this, and the Verlinde dark-matter connection (§7), in the companion paper.

## 7. Relation to other sectors

**Dark matter.** SEDE **does not derive dark matter**; it assumes the standard cold component. What it
*does* provide are falsifiable *relationships*: the equation of state is locked to Ω_m through the
fixed-point background — the canonical model giving w₀ ≈ −1 with a crossing set by the growth gate
(§3.2; the closed-form w₀ ≈ −0.85 enters only as a structural diagnostic, not the relation) — and the
coincidence problem is reframed: dark energy's activation gate f_sat(z) tracks the
(dark-matter-dominated) growth of structure, so "why now" corresponds to f_sat crossing ≈ 0.5 near
z ≈ 1, shortly before matter–DE equality. These are testable links, not a unification, and we state
them as such.

There is, however, a *principled* route to more, which we flag as a future direction and pursue in a
companion paper, not here. The same volume-law dark-sector entropy SEDE uses for dark *energy* is the
object from which Verlinde's emergent gravity (2016) derives apparent dark *matter* (baryonic
displacement of the de Sitter volume-law entropy → the MOND scale a₀ ∼ cH₀); if SEDE's Δ = 1 horizon
entropy and Verlinde's volume-law de Sitter entropy are the same object, **one volume-law postulate
could in principle source both dark energy and apparent dark matter**. We do **not** adopt this here:
Verlinde's proposal is contested (it struggles with clusters and lensing and lacks a covariant
CMB/structure completion — the regime SEDE *is* tested in), and our own cross-check
(`run_verlinde_crosscheck.py`) finds only an *order-of-magnitude* consistency — the CKN/flatness
normalisation sets a₀ within an O(1) factor of observed ((c√Ω_DE0 H₀)/2π ≈ 0.73 a₀, ∼25% low), not a
derivation. We record it as a motivated, partially-consistent target for the companion; SEDE's data
fits use standard cold dark matter throughout.

**Quantum and cross-field.** SEDE's horizon is fractal/volume-law while its *bulk* is standard GR
(holographic-DE scope), so it predicts a deformed horizon **without** a Lorentz-violating, foamy bulk.
Two real datasets are consistent with this restricted-scope picture (they constrain bulk effects SEDE
does not predict, so they do not contradict it): the Fermi GRB 090510 31 GeV photon excludes
Planck-scale linear Lorentz violation (E_QG > several M_P), and the Fermilab Holometer's null result
excludes bulk Planck-amplitude "holographic noise." The CKN reading gives a gravitating-vacuum cutoff
Λ_UV ≈ 4 meV at the dark-energy scale (§8.2). We tested one cross-field extrapolation and it **failed**:
the cosmological δ = 3/2 does *not* transfer to galaxy-cluster velocity statistics — a real SDSS/Coma
fit gives a Gaussian (Tsallis q = 1.01, not the q ≈ 1.5 a naïve transfer would predict). We report
this negative result plainly; the galaxy-dynamics link is the weakest cross-field claim and the data
do not support it.

**A contingent tension: cosmic birefringence.** Minimal SEDE predicts **β = 0** by symmetry (its
parity-even thermodynamic scalar has no φFF̃ coupling — the same structural fact as α_T = 0), whereas a
nonzero isotropic rotation is reported (ACT DR6 0.215°±0.074°, 2.9σ [arXiv:2509.13654]). We treat this
as a *real-but-contingent* threat, not a present falsification, and assess it in the §9 red-team
(threat 2): the signal's cosmic, dark-energy-sourced nature is not established (miscalibration-
degenerate α+β; the cleanest β×galaxy channel is null at SEDE's zero), but a confirmed cosmic β would
favour a single-axion rival that does SEDE's crossing and more.

**Where SEDE sits among dark-energy models — and what it borrows.** SEDE is a thermodynamic model, but
it is not outside the standard effective-field-theory description of dark energy, and locating it there
removes its main field-theoretic disadvantage. In the EFT of dark energy (Gleyzes–Gubitosi–Piazza–
Vernizzi), SEDE occupies the **safe corner**: the bulk is general relativity (holographic scope), so the
Planck-mass run rate **α_M = 0** (Ġ/G = 0, §6/V7) and the tensor speed excess **α_T = 0** — gravitational
waves travel at c, so SEDE is **GW170817-safe by construction** (|c_T/c − 1| = 0, against the ∼10⁻¹⁵
bound), having never switched on the higher Horndeski operators that GW170817 excludes. The dark energy
is smooth with **c_s² = 1** (§3.4), so it is gradient-stable (c_s² > 0) and does not cluster. SEDE is
thus the (α_T = α_M = 0, c_s² = 1) point — the *least*-constrained region of EFT-of-DE. **Field-theory
representation and stability through the crossing** (`run_eft_lagrangian.py`): a *single*
minimally-coupled scalar cannot reproduce the w = −1 crossing — for any P(X,φ), ρ + p = 2X P_X, so the
crossing forces P_X → 0 (verified: 1 + w → 0 at z ≈ 0.18), at which point the perturbation kinetic term
α_K ∝ P_X + 2X P_XX passes through zero (a ghost) and c_s² is ill-defined (the Vikman no-go). A
*fundamental* ghost-free crossing nonetheless exists with **braiding** (Kinetic Gravity Braiding /
Horndeski G₃(φ,X)□φ): the no-ghost determinant Q_S = α_K + (3/2)α_B² stays positive even when α_K → 0,
because a nonzero α_B carries it (we verify Q_S = 1.15 > 0 at the crossing for a representative
α_B = 1.5 Ω_DE realisation, with c_s² > 0 throughout) [Deffayet et al. 2010; Kimura & Yamamoto 2010].
SEDE, however, needs neither: ρ_DE = T_AH·s_grav is built from covariant horizon scalars, not a
fundamental field, and is treated as an effective smooth fluid (PPF, c_s² = 1) — gradient-stable
(c_s² = 1 > 0), with no propagating scalar to ghost, so the no-go simply does not apply (the
thermodynamic content, not a Lagrangian, carries the prediction).

The structure-coupling — SEDE's most distinctive and least-confirmed ingredient (the structure-sourced-
dark-energy idea itself traces to Gough's information dark energy, §7.6; what is specific here is the
horizon-thermodynamic realisation) — is also **free of the two instability classes that afflict coupled
dark energy** (`run_eft_stability.py`, V55). There is no
gradient or ghost instability because c_s² = 1 > 0; and there is no early-time large-scale instability
because the interaction fraction Ω_DE(z)·|d ln φ/d ln a| → 0 at high z (it peaks at ≈ 0.38 near z ≈ 0.5
and falls to ∼ 10⁻¹⁰ by recombination, since the dark-energy fraction itself vanishes as f_sat → 0). The
coupling switches *off* in precisely the radiation/early-matter era where coupled-DE instabilities are
seeded; the full perturbation evolution is validated end-to-end by CLASS-PPF (fσ8 to 0.1%, §4.4). So two
of the apparent disadvantages relative to other programmes in fact *invert*: against Horndeski SEDE is
less constrained (GW-safe), and against interacting dark energy it is more stable.

| framework | what SEDE shares or borrows | net position |
|---|---|---|
| ΛCDM (Planck) | empirical benchmark | **honest limit:** far less established — only DR3/Euclid closes this |
| quintessence (Ratra–Peebles) | dynamical DE; tracker/attractor crossing (§8.1, Result 12–13) | structure-locked w(z) with a non-tuned −1 crossing |
| CKN bound | the UV–IR energy bound fixes the DE scale (§8.2) | SEDE turns it into a cosmology; **limit:** needs the volume-law form on top |
| Li holographic DE | holographic DE with a horizon cutoff | uses the **apparent** horizon (no future-horizon causality/circularity); **limit:** volume-law is radical |
| standard Barrow-HDE (Saridakis 2020) | the Barrow entropy S ∝ A^{1+Δ/2} | SEDE's Δ enters through the H-coupling λ = 1−Δ/2 *and* the structure gate, not the HDE density directly — so neighbouring Barrow-HDE fits that prefer Δ ≪ 1 [arXiv:2506.03019] constrain a *different* model and do not exclude SEDE's Δ = 1 (§9); **limit:** the volume-law endpoint is the radical input |
| Hu–Sawicki f(R) | viable late-time acceleration | keeps GR (α_M = α_T = 0): no screening/GW tension; **limit:** less mature |
| CPL | the data language; w(z) is CPL-degenerate in shape (§6) | adds the growth–expansion lock CPL lacks; **limit:** model-dependent, not a neutral fit |
| early dark energy | — | **honest limit:** does not address H₀ (out of scope, r_d-pinned, §9) |
| Jacobson (1995) | Clausius → Friedmann; ρ_DE = T·s (App. A.1) | applies horizon thermodynamics to DE; entropy-law change **confined to the DE fluid** (§8.3), not gravity |
| Horndeski / EFT-of-DE | the α-function language | sits at α_T = α_M = 0, c_s² = 1 — GW170817-safe corner; **less** constrained, not less fundamental |
| Verlinde (2016) | volume-law de Sitter entropy; elastic/memory response (§6, §8.1) | independent motivation for the postulate; a₀ ∼ cH₀ consistent at O(1) (§7) |
| interacting DE | the dark-sector energy-exchange ledger (CHR, §6) | coupling fixed to growth, c_s² = 1, off at high z ⟹ **instability-free** |

The genuinely irreducible disadvantages, after these borrowings, reduce to two, and we state them
plainly: SEDE is **empirically younger** than ΛCDM (the decisive tests are DR3/Euclid, §6), and it rests
on the **one volume-law postulate** (§8.4) that motivation — Barrow, the GSL, CKN, Verlinde — supports
but does not derive. Everything else in the table is either shared, borrowed, or an honest matter of
scope.

**Symmetry origin of the null predictions.** SEDE's many "null" results are not independent assumptions:
they follow from a *single* structural condition — *the dark sector adds no symmetry-breaking operator
to the GR + Standard-Model effective action; it is a parity-even scalar (ρ_DE = T_AH·s_grav·f_sat)
minimally coupled through T_μν.* Each precision null is then forced by a specific symmetry:

| null prediction | symmetry that forces it |
|---|---|
| β = 0 (no cosmic birefringence) | parity in electromagnetism → no φ F F̃ coupling |
| α_T = 0, c_T = c (GW170817-safe) | 2-derivative diffeomorphism-invariant gravity → no Chern–Simons / Galileon term |
| α_M = 0, Ġ/G = 0, no fifth force, equivalence principle | minimal (metric-only) coupling → no non-minimal ξφ²R |
| no Lorentz violation (Fermi GRB, Holometer) | bulk Lorentz invariance |

So β = 0 and c_T = c are *the same fact* (§7 birefringence; absence of a parity-violating /
modified-propagation sector), and the whole cluster — the reason SEDE passes every GW, fifth-force,
birefringence, and Lorentz null *at once* — is one principle, not a list. The de Sitter endpoints add a
fifth: f_sat → const (with H → const) ⟹ a maximally-symmetric de Sitter phase ⟹ w = −1 by symmetry (the
dark-energy bracket of Result 12; App. A.7), with the crossing in between dynamical. We note the limit
sharply: this symmetry argument protects SEDE's *safety*, not its *content*. The distinctive,
load-bearing inputs — the volume law Δ = 1 (a geometric/thermodynamic postulate, §8.4), the coupling
γ ≈ 1.5 (halo binding physics, §3.1), and the crossing redshift (a dynamical balance, Result 13) — are
**not** symmetry-derivable; symmetry buys SEDE nothing there, and the one open postulate is untouched.
The symmetry-protected cluster is precisely the part SEDE shares with any minimally-coupled smooth dark
energy; what makes it SEDE remains thermodynamic and postulate-based.

### 7.6 Relation to information dark energy and other structure-sourced models
The idea that dark energy is *sourced by the entropy/information produced as structure forms* is **not
original to SEDE**, and we state the priority plainly. It is due to **M. P. Gough**, in a sustained but
under-cited programme that we date — from a NASA-ADS author search — to **2006–2008** ("On the Similarity
of Information Energy to Dark Energy," 2006; "Information Equation of State," 2008), developed through
*Holographic Dark Information Energy* [*Entropy* **13**, 924 (2011), arXiv:1105.4461] to a recent
observational analysis [*Entropy* 2025] that fits the stellar-mass-density history to the dark-energy
equation of state from Planck/Pantheon+/DES/DESI, finding ρ_DE ∝ SMD(a)^{p} with p ≈ 0.5 (R² = 0.93).
In Gough's model dark energy is the Landauer energy of the information/entropy of star-heated gas and
plasma. SEDE shares his central physics — structure's entropy *is* dark energy, with a self-regulating
feedback (acceleration suppresses growth) that addresses the coincidence problem and predicts
time-limited acceleration — and we credit it as the priority for the structure-sourced-dark-energy idea.

What distinguishes SEDE is the **mechanism and a testable consequence**, not the founding idea.
(i) *Mechanism:* Gough uses the thermal/information entropy of hot **baryons** (∝ T, residing in the
bulk); SEDE uses the **gravitational** binding entropy of bound structure deposited on the **apparent
horizon** (ρ_DE = T_AH·s_grav, Barrow volume-law Δ=1). These are physically different entropies in
different places. (ii) *Prediction:* the two are **not** degenerate. SEDE's driver is the gravitational
growth factor D(z) (dark matter), Gough's is the baryonic stellar-mass-density SMD(z) (star formation);
computing both (`run` comparison, `output/fig_gough_vs_sede.png`) gives dark-energy histories that agree
to ≲ 6% for z ≲ 0.6 — where supernova data are strongest, explaining why both fit current data — but
**diverge by up to ~50% by z = 3**, SEDE staying flatter (more high-z dark energy) because gravitational
growth persists where star formation declines. The coincidence that both land on a "0.5" exponent is
superficial: Gough's is a power of stellar mass density, SEDE's is the H-exponent λ = 1−Δ/2 from
Barrow Δ=1. (iii) SEDE adds a **parameter-free** dark sector and a w-crossing rather than a fitted CPL.
So a clean way to separate them observationally is to ask whether dark energy tracks **D(z)** (gravity)
or **SMD(z)** (stars), with the high-z dark-energy fraction (DESI Lyα, future surveys) as the lever.

Other structure-/entropy-linked programmes differ more sharply and are catalogued in `NOVELTY_CHECK.md`:
entropic-force and Barrow/Tsallis holographic dark energy source the term from the **horizon** (structure
downstream — and standard entropic-force DE in fact *struggles* with structure growth, which SEDE's
design inverts); backreaction (Buchert, Räsänen) and timescape (Wiltshire) make acceleration an
**apparent** inhomogeneity effect with no real dark energy; matter-creation/CCDM (Lima) replaces dark
energy with a particle-creation pressure; and cosmologically-coupled black holes (Farrah et al.) use one
special structure via vacuum interiors. A systematic priority pass over **both INSPIRE-HEP and NASA ADS**
(citation-ranked, multiple phrasings; documented in `NOVELTY_CHECK.md`) returns **no precedent** for
SEDE's specific construction — "dark energy = gravitational binding entropy of bound structure on the
apparent horizon" (`abs:("dark energy" "binding energy" horizon)` yields zero relevant hits in either
database) — nor for its terminology. So while the *general* structure-sourced-dark-energy idea belongs to
Gough (since ~2006–2008), SEDE's *gravitational-horizon mechanism, its gravity-driven (D(z), not SMD)
signature, and its parameter-free first-principles dark sector* appear to be new. We also ran Gough's
model through our identical CAMB-marginalised pipeline (its ρ_DE ∝ SMD^{0.5} expressed as a w(a) table):
both are viable, with SEDE modestly preferred (ΔDIC ≈ −3.3 vs Gough's ≈ −0.45 relative to ΛCDM), and
Gough's fit favouring a higher H₀ — consistent with his Hubble-tension claim. In summary, SEDE is
a distinct, more-derived realisation of Gough's idea, decided from ΛCDM and from Gough alike only by the
high-z dark-energy and growth data to come (§6).

The closest *entropic-acceleration* relative is the General Relativistic Entropic Acceleration (GREA)
programme of García-Bellido and Espinosa-Portalés [arXiv:2106.16014, 2106.16012], which — like SEDE —
removes Λ and sources late-time acceleration from the entropy of the cosmic horizon, with the Helmholtz
free energy F = U − TS (not matter alone) as the gravitational source and entropy *production* as the
engine. The differences are sharp and observable: GREA uses the **boundary (area)** horizon entropy and a
single parameter (the horizon-to-curvature ratio), with the growth driven by the horizon's own
*expansion*, giving a quintessence-like w(a); SEDE uses the **volume** entropy (Δ = 1) **gated by
structure** (f_sat), which locks w(z) to the growth history and yields a phantom crossing, with a
parameter-free dark sector. The connection is quantitative — SEDE's structure-driven entropy production
maps to a positive (second-law) effective bulk viscosity that peaks at the gate-activation epoch, casting
SEDE as a *structure-gated* member of the GREA family — and it is the deformation Δ that discriminates
them: were SEDE's residue to resolve to the area branch (Δ = 0), its phenomenology would collapse toward
GREA's. So SEDE is not "GREA plus a volume postulate": the structure-gating is itself a falsifiable
signature — the w(z)↔growth phase-lock (§6) that neither GREA's expansion-driven w(a) nor a kinematic CPL
fit shares. We develop this relationship, and SEDE's foundational status, in the companion paper.

## 8. The quantum-gravity underpinning

The energy-conservation account of §2–§3 took two horizon properties as given: the **magnitude**,
ρ_DE ~ ρ_crit rather than the QFT vacuum's M_P⁴, and the **H-scaling** λ = 1/2, which follows from
the horizon entropy being **volume-law** (S_grav ∝ V_AH; §3.3) with *no deformation parameter*. This
section supplies their quantum-gravitational basis — it is the support beam under the construction,
and a reader content to take the two scales as inputs may treat it as optional. (Its full development —
the reduction route map, the driven non-equilibrium-steady-state account, and the closure analysis — is
the subject of a forthcoming companion paper; **none of the cosmological results of §§1–7 depends on it**,
and the companion's mechanisms, where invoked above, are flagged as conjectural.) We argue why the
horizon entropy is volume-law rather than area-law, anchor the magnitude to the holographic bound and
the flatness normalisation, and isolate the single statement that remains a postulate. Throughout, the Barrow parameter Δ appears
only to *label* the area-to-volume interpolation (S ∝ A^{1+Δ/2}: Δ=0 area-law, Δ=1 ⟺ volume-law) and as
the **falsifiable test** of the volume postulate (measured in §5.6, forecast in §6); it is not an
input to the model.

### 8.1 Why volume-law, not area-law — consistency motivations (none a derivation)
The dark-sector entropy is the coarse-grained entanglement entropy contained in the apparent-horizon
**volume**, S_grav ∝ V_AH ∝ R_AH³ — equivalently the Barrow endpoint Δ=1 (A^{1+Δ/2} = A^{3/2} = R³).
The open physical question is sharp: why the volume law, rather than the Bekenstein–Hawking *area*
law? We give three considerations in the main line and collect three further, more speculative ones in
Appendix G. **None is a microscopic proof**; together they motivate the volume endpoint and make the
area-law branch disfavoured in the diagnostics.

1. **Geometric ceiling caps the deformation** (Result 9). The horizon's fractal dimension is
   d_H = 2 + Δ, capped at the space-filling value d_H = 3 for a two-surface in three-space. Since
   (for a volume law in embedding dimension d) Δ = 2/(d−1), the case d = 3 gives Δ = 1. We are precise
   about what this delivers: the ceiling bounds the deformation, **Δ ≤ 1** — it does *not* select Δ = 1.
   Landing *on* the endpoint requires, in addition, that the entropy is monotone in Δ (so the cap is
   saturated); since ρ_DE = T·s makes the free energy F = E − TS vanish identically, entropy is the
   governing potential and the generalised second law is consistent with — is a **fixed point of** — the
   saturated (volume-law) state. We stress that this is a **consistency check, not a selection
   principle**: "total entropy does not decrease" is not "the parameter takes its entropy-maximising
   value," and the saturation of the cap *is* the volume-law postulate (§8.4) in effective form, not an
   independent derivation of it. Volume-law scaling is the signature of a thermalised, maximally-scrambled
   state (Page curve, eigenstate thermalisation).
2. **The holographic bound fixes the scale, not the form.** The CKN energy bound (§8.2) pins the
   *magnitude* but not the exponent — its own state count gives a *sub*-area-law (δ = 3/4), the
   opposite of an enhancement. The volume law is therefore a statement about the horizon's quantum
   *state* (typical / thermalised), which the energy bound cannot supply. This isolates the open
   content to one sentence: *the horizon is in a typical, volume-law-entangled state.*
3. **The area-law branch is disfavoured directly by current data.** *Saturating* the CKN bound — the
   natural "no-deformation" alternative — is algebraically the area-law case Δ = 0, whose equation of
   state (w₀ = −0.49) and large early-dark-energy fraction are disfavoured in the conditional
   diagnostics of §5.6. On current data the volume branch is **favoured over the area-law
   alternative**, not merely an optional choice — though, as §5.6 stresses, the robust marginalised
   measurement is still to come. The push is *monotone* in Δ: the exact identity
   Ω_DE(z) = Ω_DE0 · f_sat(z) · E(z)^{−Δ} (from ρ_DE ∝ H^{2−Δ}f_sat) and E > 1 in the past make a larger
   Δ suppress the early-dark-energy fraction more strongly — Ω_DE(z = 3) falls from 0.127 (Δ = 0) through
   0.060 (Δ = 0.5) to 0.029 (Δ = 1, below ΛCDM's 0.035; `run_cmb_earlyde.py`). So the geometric ceiling
   of point 1 (Δ ≤ 1) and the early-dark-energy lever here *cooperate*: the ceiling bounds Δ from above
   and the CMB pushes it to that bound, the volume endpoint. This remains a **diagnostic, not a
   derivation** — but it is a third, independent reason the data sit at Δ = 1.

Three further consistency motivations — a *driven non-equilibrium steady state* (the structure gate
maintains the volume-law state despite a long scrambling time), a *roughening / Edwards–Wilkinson
universality-class* analogy for the horizon surface, and the *independent volume-law entanglement
entropy* of Verlinde's emergent-gravity programme — are collected in **Appendix G**. Each makes the
postulate less *ad hoc*; none derives it.

### 8.2 Why ρ_DE ~ ρ_crit — the CKN bound and the vacuum-energy scale
The deformation fixes the *scaling* of dark energy with H but not its *magnitude*. The magnitude
comes from the Cohen–Kaplan–Nelson (CKN) holographic energy bound (Principle 0): an effective field
theory in a region of size L cannot hold more energy than a black hole of that size, ρ L³ ≲ M_P² L,
i.e. ρ ≲ M_P²/L². With the IR cutoff at the cosmic horizon, L = c/H, this gives
$$\rho_{\rm DE} \lesssim M_P^2 H^2 \sim \rho_{\rm crit} \;\Rightarrow\; \Omega_{\rm DE} \sim \mathcal{O}(1).$$
The naïve QFT vacuum, ρ ~ M_P⁴, overshoots this by ρ_P/ρ_crit ~ 10¹²² — the cosmological-constant
problem. SEDE never invokes the QFT vacuum: dark energy is the horizon's conjugate thermal energy
density, and the CKN relation supplies the **horizon-scale upper bound** ρ_DE ≲ M_P²H². This is a bound, not an equality: the
present amplitude is fixed by spatial flatness (Ω_DE0 = 1 − Ω_m), in the same bookkeeping sense that
ΛCDM fixes Ω_Λ0 — once the volume-law branch is chosen, the bound is saturated to O(1) by that
normalisation **today**. We note explicitly that the saturation is a present-day statement: since
ρ_DE ∝ H (volume law) while the bound ∝ H², the ratio ρ_DE/(M_P²H²) = Ω_DE(z) falls at higher z
(0.70 → 0.22 → 0.03 at z = 0, 1, 3), so dark energy sits *below* the CKN bound in the past — consistent,
as CKN is an upper bound. So CKN *addresses* the vacuum-energy-scale problem (it explains why the natural
gravitational horizon scale is M_P²H² ~ ρ_crit, not M_P⁴) rather than deriving the magnitude from
first principles. The same bound, read as a UV–IR relation, gives a gravitating-vacuum cutoff
Λ_UV = (M_P H₀)^{1/2} ≈ 4 meV — numerically the observed dark-energy scale (§7).

### 8.3 The seam, BBN-safety, and the holographic scope
A naïve *microscopic* volume-law entropy density would be Planckian, so a naïve ρ_DE = T·(S_vol/V) would
overshoot ρ_crit by (M_P/H)^Δ ~ 10⁶¹ — the square root of the 10¹²² hierarchy (Result 11). The same
factor appears as a clean *no-go* (`run_deriv_E_nogo.py`): the Bekenstein bound on the horizon's energy
E = ρ_crit V_H gives S ≤ 2πRE/ħc = ¼(A/ℓ_P²) — exactly the Gibbons–Hawking *area* law (verified to within
the O(1) geometric factor) — so the volume law cannot be obtained from any maximum-entropy or energy-bound
argument; it *exceeds* the bound by precisely R/ℓ_P ~ 10⁶¹. This sharpens rather than weakens the
postulate: it proves the missing ingredient is a *state/counting* input (a Planckian-density, volume-law
horizon), not an energy bound, and fixes its exact size. The *same*
factor controls three apparently separate issues — this seam, the BBN bound on a Barrow-modified
Friedmann equation, and the modified-gravity-vs-holographic fork — which are one object. The resolution
is the holographic-DE scope adopted in §2: the deformation acts **only on the dark-energy fluid**
(magnitude from CKN, scaling from Barrow), **not** on the gravitational sector. The early-universe
expansion is then exactly standard (ΔN_eff = ΔY_p = 0, verified; BBN-safe), and the naïve T·(S/V) is
precisely the modified-gravity path the scope forbids.

This is the **load-bearing assumption of the model**. In the *modified-gravity*
reading — Barrow entropy deforming the Friedmann equations [arXiv:2110.00059] — BBN forces
**Δ ≲ 1.4×10⁻⁴** [Barrow, Basilakos & Saridakis, arXiv:2010.00986] and SN+BAO+H(z) give Δ ~ 10⁻⁴ (the
area law at 2σ [arXiv:2108.10998]), which read literally would exclude SEDE's Δ = 1. Those bounds do
*not* apply to SEDE: it keeps standard gravity and lets the volume-law entropy act only on the
dark-energy fluid via ρ_DE = T·s (§2, App. A.1), so there is no extra Friedmann term, no BBN deviation,
and Δ is unconstrained by them. The headline Δ = 1 thus rests on this single scope choice — legitimate
and standard (it *is* the holographic-dark-energy framework), but not free: a reader who rejects it
recovers the tightly-bounded modified-gravity Barrow cosmology instead.

### 8.4 The one irreducible postulate

> **Postulate (volume-law horizon).** The cosmic apparent horizon carries the coarse-grained
> *entanglement* entropy of its **volume**, S_grav ∝ V_AH ∝ R³ — a constant horizon entropy density,
> i.e. a thermalised, volume-law (non-extensive) state — rather than the Bekenstein–Hawking area law.

The results in §8.1–§8.3 are bounds, consistency arguments, or motivations *given* this postulate; the postulate
itself is not derived from a microscopic theory of quantum gravity. §8.1 *reduces* it (the CKN bound
fixes the scale, so the only open content is the horizon's entanglement *state*), *motivates* it
(non-extensivity, the GSL, the driven steady state, the Edwards–Wilkinson universality class, and the
independent volume-law entanglement entropy of Verlinde's emergent gravity), and
shows the *area-law branch is disfavoured by current diagnostics* — but a microscopic state count
that delivers the volume law remains the model's single open foundational item (§9). Two points of
emphasis. First, the model carries **no fitted deformation parameter**: Δ is not an input but the
falsifiable *test* of this postulate — current data already prefer the volume endpoint and disfavour
the area-law branch (§5.6, a conditional diagnostic), with the robust marginalised measurement to
come from DESI DR3 + Euclid (§6), alongside cross-horizon predictions (a deformed black-hole entropy
S ∝ A^{3/2}; §7, presented as speculative). Second, the open question is now stated plainly — *why is
the horizon entropy volume-law (coarse-grained entanglement) rather than area-law?* — rather than
hidden in a value of a tunable parameter. We do not claim to have derived quantum gravity; we claim
to have reduced its dark-energy content to one testable statement and removed it as a fitted
parameter.

**On the plausibility of a volume law (genre, not derivation).** A volume-law horizon entropy is a
strong hypothesis, but volume-law *entanglement* entropy is at least not exotic in cosmological field
theory. The entanglement entropy of a field subregion — area-law in the de Sitter adiabatic vacuum —
grows and saturates to a *volume law* after a de Sitter→flat transition, ending in a partially
thermalised (generalised-Gibbs) state [Khlebnikov & Sheoran, arXiv:1907.00487]; volume-law entanglement
also arises from field–curvature coupling [arXiv:2306.08357]. We cite these only to establish *genre and
plausibility* — that volume-law entropy is a real feature of cosmological QFT — and we explicitly resist
over-reading them. Two caveats: (i) these results concern the entanglement entropy of field
*subregions* under specific transitions, not the cosmic-*horizon* entropy SEDE invokes; and (ii) their
assignment (de Sitter = area-law, flat = volume-law) is *not* SEDE's and in fact runs opposite to it
(SEDE places the volume law in the de Sitter brackets). So the literature makes a volume-law hypothesis
less exotic, but it neither *anticipates* nor *derives* the specific postulate — which remains the
model's one open foundational item (§9). (We note for completeness that Padmanabhan's
holographic-equipartition route, dV/dt = N_sur − N_bulk, does *not* help here: its surface count is
area-law, N_sur = A/ℓ_P², and it reproduces the *standard* area-law Friedmann equation — it is the
baseline SEDE departs from, not independent support for the volume law.)

### 8.5 The postulate is a dof-*counting* claim, not a thermalisation claim

The phrase "volume-law" conflates two claims of very different status. **(i) The state:** are the
horizon degrees of freedom maximally entangled (thermal) or in a low-entanglement ground state?
**(ii) The counting:** does the *number* of horizon dof grow as the area, N ∝ R^{d-2}, or the spatial
volume, N ∝ R^{d-1}? Since S ∼ (state factor) × N, the deformation Δ — through S ∝ A^{1+Δ/2} ∝
R^{(d-2)(1+Δ/2)} — is fixed by the **counting** (area→Δ=0, volume→Δ=1 in d=4); the state only sets
whether the prefactor is maximal. SEDE's postulate is, precisely, the **counting** claim.

*The state half is settled, and is not peculiar to SEDE.* Maximal scrambling has a precise meaning in
the Sachdev–Ye–Kitaev model [Sachdev & Ye 1993; Kitaev 2015], the canonical chaotic, black-hole-dual
system. But maximal scrambling **does not fix Δ**: a black hole is a *maximal* scrambler yet
**area-law** (S = A/4, Δ=0, reproduced by string microstate counts). Scrambling is therefore shared
with the Δ=0 black hole; and SYK, being all-to-all (geometry-free), cannot even pose the counting
question. (An explicit Majorana-SYK diagonalisation confirming the state is maximally entangled —
Wigner–Dyson statistics, Page-value saturation, complete OTOC scrambling — is given in Appendix G;
it settles the *state*, not the counting.)

*The counting half is where SEDE lives, and area is the default.* In causal-set quantum gravity — the
setting where the counting question is natural — both scalings appear, but the *canonical*
horizon-entropy object, the causal-**link** count that recovers the Bekenstein–Hawking law [Sorkin and
collaborators], is the **area** one (R^{d-2}), in flat space and identically in de Sitter. So even in
the framework most hospitable to volume-counting, **area-law is the default**, and SEDE's Δ=1 is the
specific, non-default identification of horizon entropy with the **bulk** count. (The causal-set
sprinkle measurements, the entropy-bound check, and a Ryu–Takayanagi holographic-code realisation —
volume-law requires *nonlocal* horizon connectivity — are in Appendix G; their verdict is that the
volume count is **allowed and realisable but not selected**: no consistency argument available to us
derives it, and none excludes it.)

**The sharpened postulate, properly located.** The one open item is therefore *not* about
thermalisation or chaos (settled, and shared with area-law black holes) but is purely a
**dof-counting** statement:

> the cosmic horizon's entropy counts its **bulk** (volume-scaling, N ∝ R³ ⟺ d_H = 3) degrees of
> freedom rather than its **boundary** (area-scaling, links) ones — a count that exceeds the
> holographic area bound by R/ℓ_P (the seam of §8.3).

It reduces to two deciders. *Theoretically* it is precisely *"is the Hilbert-space dimension of the de
Sitter static patch e^{Area/4} or e^{Volume}?"* — a frontier problem a complete dS holography would
settle. *Empirically* it is the deformation Δ: Δ=0 (area) is already disfavoured by present data
(§5.6), and DESI DR3 + Euclid will measure Δ to ~0.09 (§6). SEDE's postulate is thus not an arbitrary
assumption but a *consistent, realisable bet on a well-posed open question of quantum gravity*, with a
decisive test imminent. Tellingly, the one horizon we can check today — the black hole — is area-law; on
the evidence here this asymmetry is a genuine conceptual cost. The companion paper *offers* a resolution
— a state-dependent count in which a black hole is area-law because *undriven* and only the
structure-driven cosmic horizon activates the volume law — which, if its (conjectural) mechanism holds,
would turn the cost into a predicted feature, with the black-hole side derivable and a sharp falsifier (an
isolated volume-law black hole).

### 8.6 Is the postulate derivable? Summary of the reduction (full treatment in companion paper)

The volume-law postulate is not monolithic: it splits into a *state* (the horizon dof are maximally
entangled/thermal), a *form* (the entropy scales as the volume), a *scale* (the magnitude is ρ_crit, not
ρ_Planck), and a *count* (the number of dof grows as the bulk, N ∝ R³, not the boundary, N ∝ R²). The
deformation Δ is fixed by the **count** alone. A systematic study of how derivable each piece is — five
routes (de Sitter holography/SYK, entanglement first law in a thermal state, Verlinde emergent gravity,
gravitational non-additivity, and a Bekenstein no-go), the reduction of the count to a *driven
non-equilibrium steady state*, and three closure routes — is developed in a companion paper; we state the
conclusions here, as they bear on how *ad hoc* the postulate is.

The *state* is first-principles (a maximal scrambler is volume-law-entangled; Appendix G, §8.5); the
*form* reduces to thermalisation (volume-law entanglement is the generic behaviour of a thermal reference
state for regions larger than the thermal length); the *scale* is CKN (§8.2). What remains irreducible is
purely the **count**, and a Bekenstein no-go makes precise *why*: maximum-entropy/energy bounds give the
*area* law (§8.3), so the volume count cannot be obtained from any energy argument — it is a genuine
state/counting statement. The count is, in turn, downstream of *connectivity* (locally-connected horizon
dof give area-law entropy, nonlocal give volume-law) and connectivity of *interaction range*: gravity is
strongly long-range (1/r with α = 1 ≤ d), hence non-additive — the volume class is the *natural* one for a
self-gravitating horizon, with area-law arising only when the holographic bound intervenes for an
*isolated, equilibrium* horizon (a black hole). The cosmic horizon differs in being continuously *driven*
by structure formation, which converts the black-hole/cosmic asymmetry (§7) into a discriminator
(equilibrium → area, driven → volume) and reduces the residue to a single, falsifiable **driven-NESS
statement**: a strongly-long-range horizon, driven by structure, self-organises to its area↔volume
spinodal and locks the volume branch. The companion paper supplies both microscopic maps this requires
(the cooperative coupling J ≥ J_c is the gravitational binding itself; the drive is the
structure-deposited entropy, clearing the barrier at z ≈ z*) and verifies the volume-lock is robust over a
broad, untuned basin — so the static *counting* postulate is replaced by the weaker, more physical
statement that the horizon free energy is bistable with a sub-unity barrier.

**Closure** — deriving even that bistable free energy from first principles — would reduce the dark sector
to Ω_m plus established physics. It takes three precise forms: a de Sitter-holography computation of the
bulk count (open, and provably not shortcut-able through energy bounds); exhibiting the volume branch as a
real saddle/state (the bistable landscape follows from gravitational cooperativity given that branch); or
*measuring* Δ. The last is decisive and already in hand: DESI DR3 + Euclid pin Δ to ∼0.09 (§6), so Δ = 1
establishes the volume count as fact and Δ = 0 refutes it. The postulate is thus a single, sharply-posed,
falsifiable statement with an imminent test — not an unconstrained assumption.

## 9. Discussion: open problems & honest assessment

**The one open item.** Every dark-sector input is fixed by a stated choice (§3) given the single
postulate that the horizon is volume-law-entangled (§8.4). That postulate is motivated by gravity's non-extensivity
and the second law but is not derived from a microscopic state count; it is the model's one
irreducible input. The canonical-vs-microcanonical reframing (which gives area-law λ = 1 or
Planckian-volume λ = 1/2) does not give both for free, but it is not symmetric: the de Sitter
horizon has **negative specific heat** (C = −2S) and gravity is strongly long-range (α = 1 ≤ d,
§8.6), so by the standard non-additive-systems result [Campa, Dauxois & Ruffo 2009] the two
ensembles are *inequivalent* and the canonical one is ill-defined for the isolated self-gravitating
horizon. This selects the *microcanonical accounting* (and removes the canonical/area-law horn as
unphysical for a long-range, negative-specific-heat horizon); it does not by itself derive the volume
count — within the microcanonical ensemble, whether the *state
count* is itself volume (the dS-holography question of §8.6): a genuine open problem, and a
*falsifiable* one (§6). The systematic route map of §8.6 narrows it: the
*state* (maximal entanglement) is first-principles, the *form* (S ∝ V) reduces to thermalisation, the
*scale* is CKN, and the irreducible residue is purely the bulk-vs-boundary **count**. That count, in
turn, reduces to one provable statement — a **driven-NESS entropy theorem** (a strongly-long-range
horizon, driven by structure formation, settles in the volume-law steady state rather than the area-law
equilibrium) — which, if established, would close the residue *and* resolve the black-hole/cosmic-horizon
asymmetry (§7) as the same equilibrium-vs-driven mechanism. Both its microscopic maps are now
supplied (§8.6): connectivity → J ≥ J_c (the cooperative coupling is the largest eigenvalue of the
gravitational coupling matrix = the per-site binding = the *same* super-extensive quantity behind the
non-additivity, so volume-counting and volume-locking are one number, generic at the horizon's N), and
deposition → drive (the accumulated structure entropy, ∝ f_sat, clears the locking barrier at z ≈ z* with
magnitude set by the near-critical CHR identity ε_dep·χ ≈ 1, and locks Δ = 1 permanently). What survives is
the shape of the horizon order-parameter free energy F(m). The monotone structure drive is a *ratchet*
that necessarily carries the horizon through the area-branch spinodal at z* and the second law locks it at
the volume well (S_vol ≫ S_area); we verify this is robust over a broad, untuned basin of the landscape
parameters, drive amplitude, and initial condition (`run_soc_attractor.py`) — self-organised, not a tuned
knife-edge. So the static volume-law *counting* postulate is reduced to one *dynamical* and falsifiable
statement: the horizon's area↔volume free energy is bistable with J ≥ J_c (gravity supplies this) and a
sub-unity barrier. We state the theorem as the well-posed successor to the postulate — both maps supplied,
the residue reduced to that one structural assumption and shown robust — not as a finished proof.

**Scope and tensions.** H₀ is out of SEDE's scope (it is r_d-pinned; the model does not address the
distance-ladder tension and we do not claim it does). The structure-sourcing *mechanism* itself — the
defining SEDE claim that f_sat tracks growth — is only weakly confirmed: a geometry-vs-structure
consistency test (E1; the null test ρ_DE = G[D(z)] with f_sat read independently off the *expansion*,
DESI D_H/r_d → H(z), and off *structure growth*, RSD fσ8 → σ8(z)) tracks at z < 1 but is RSD-limited
at z > 1 and remains inconclusive (current data: χ²/dof = 0.40 in the RSD-error budget, 1.40 in the
precise-geometry budget; Figure 9). **The pipeline is implemented and staged** (`run_e1_mechanism.py`) — it runs
now and is ready to apply the moment Euclid/LSST weak-lensing tomography sharpens the structure leg, at
which point the forecast of §6 makes it decisive (∼4σ). Its data vector and CONFIRM/REFUTE thresholds are
**pre-registered** (PREREGISTRATION.md §F) so the staged test cannot be tuned after the data arrive. We have separated this mechanism into a derived, parameter-free
part — the growth–expansion lock P2, which follows from f_sat = f(D²) with no extra assumption (§6),
testable with exactly that data — and an optional near-critical layer (deferred to a companion paper)
that supplies a further kill test. But a *derived consistency relation* is not yet a *confirmation*, and
the mechanism stays unconfirmed until the lock P2 is measured. We stress this is a *near-term* gap, not a
permanent one: a Fisher forecast (§6, `run_mechanism_forecast.py`) shows the growth-side mechanism test
reaches ∼4σ with combined Euclid/LSST weak lensing + CMB-lensing + clusters + RSD — so the distinctive
claim is decidable on the same timescale as the amplitude Δ, not deferred indefinitely. And the cosmic-
birefringence signal (§7) is a live tension. We
emphasise the honest framing: SEDE's *phenomenology* (w crossing −1, the volume-law endpoint) is
mildly favoured in the current compiled-data analysis with an encouraging out-of-sample check, but its
*distinctive mechanism* (the f_sat gate tracking structure growth) is the unconfirmed frontier.

**Threats from the recent literature (a red-team).** We list the strongest recent results that could
damage SEDE, with our honest assessment of each.

1. *Holographic DE disfavoured by full CMB / the early–late tension.* Wu et al. [arXiv:2509.02945] find
   that holographic and Ricci DE are disfavoured by full ACT DR6 + DESI through an early-vs-late tension
   that compressed (R, ℓ_A) priors hide — a direct challenge, since SEDE's headline ΔDIC used compressed
   CMB. We tested it. **Marginalised full primary likelihood:** a CAMB-in-the-loop fit on the full
   Planck plik_lite TTTEEE + lowℓ TT/EE, with the five cosmological parameters and the Planck nuisance
   *all free* (not θ\*-fixed), gives **Δχ²(SEDE−ΛCDM) = −0.43** (`run_full_cmb_mcmc.py`), the preference
   sitting in the high-ℓ peak structure (−0.37) — *negative*, i.e. mildly favouring SEDE, **not** the
   large *positive* Δχ² the HDE early–late tension would produce; the best-fits are physical
   (Ω_m = 0.295/0.315 for SEDE/ΛCDM). This is corroborated by two θ\*-matched single-point proxies — the
   same plik_lite (Δχ²_CMB = −0.50, `run_plik_check.py`) and the official ACT DR6 CMB-lensing likelihood
   (Δχ²_ACT-lens = −0.32, `run_act_lensing.py`, validated against the package's fiducial χ² = 14.06). ACT
   DR6 lensing re-evaluated **at the marginalised primary best-fit** (not a fixed fiducial) gives
   Δχ²_ACT-lens = −0.31 (ΛCDM 14.10 → SEDE 13.79), and the low-ℓ ISW-sensitive channel is flat (−0.03).
   All are |Δχ²| < 1: SEDE does **not** inherit the standard-HDE tension, because it has *less*
   early dark energy than ΛCDM (f_sat → 0, Ω_DE(z=3) ≈ 0.035 < 0.040), not the runaway high-z behaviour
   that breaks original HDE. **Full non-lite plik (foregrounds sampled):** we have since run the
   *complete* high-ℓ plik TTTEEE likelihood with all 15 foreground/calibration nuisance parameters
   sampled at their Planck priors (`run_full_cmb_mcmc.py --fullplik`), the gold-standard primary CMB —
   it gives **Δχ²(SEDE−ΛCDM) = −0.33** (TTTEEE −0.27, lowl.TT −0.23, lowl.EE +0.17), confirming the
   foreground-marginalised plik_lite result (−0.43) to within the minimiser noise. So **letting the
   foregrounds float does not reveal a hidden tension** — the standard-HDE early–late wall is simply
   absent for SEDE. *Residual scope:* only the ACT DR6 *primary* multifrequency likelihood remains
   not run (ACT *lensing*, the DE-sensitive channel, is; §above) — **formally deferred** under the
   pre-registered trigger (PREREGISTRATION.md §E; `act_dr6_mflike` + multi-GB products not installed),
   with the falsifiable committed expectation |Δχ²| < 1 given three concordant CMB datasets already at
   that level; a full MCMC for ΔDIC is deferred (expected ≈ Δχ², p_D shared). The full primary
   likelihood is also folded *into the joint* (§5.1).
2. *Cosmic birefringence — a contingent threat, with a rival waiting if it is confirmed.* SEDE predicts
   β = 0; a nonzero *rotation* is reported (ACT DR6 2.9σ [arXiv:2509.13654]; Planck ~0.3°), and *if* it
   is cosmic and dark-energy-sourced, a single axion explains β + acceleration [arXiv:2504.17638] and an
   interacting axion sector reproduces the phantom crossing + β + dark matter [arXiv:2506.12589] — a
   parsimony threat, since the rival does SEDE's crossing and more with one field. But the cosmic, DE-
   sourced nature is *not established*: the measurement is of α + β (miscalibration-degenerate), Planck
   and ACT disagree by >3σ, SPIDER is incompatible with the coupling [arXiv:2510.25489], and the
   β × LSS cross-correlation is null (p ≈ 37%, [arXiv:2509.22273]) — at SEDE's predicted zero, not the
   axion (§7). So this is a *contingent* threat, not a present falsification, and the cleanest current
   channels favour β = 0. The discriminator is the pair {w-crossing, β}; β = 0 is itself symmetry-derived
   (§7), tied to GW-safety, so a confirmed cosmic β would also predict GW parity violation.
3. *The DDE signal may be a low-z SN systematic.* Huang, Cai & Wang [arXiv:2502.04212] show the DESI
   evolving-DE preference is biased by an external low-z SN sample in DESY5 and drops to <2σ once
   corrected. *Defence (partial):* the bias is DESY5-specific ("in contrast to the uniform behaviour of
   PantheonPlus"), and SEDE's headline uses Pantheon+; more decisively, SEDE's leave-one-probe-out
   **drops supernovae entirely and stays preferred (Δχ² = −2.47)**, and F_AP / tracer-LOO are SN-free.
   SEDE is not SN-driven — but if the field's evolving-DE evidence softens, SEDE's *context* erodes, and
   its "suggestive, not established" framing must hold.
4. *The nearest Barrow-HDE fit prefers Δ ≪ 1.* Luciano, Paliathanasis & Saridakis [arXiv:2506.03019] fit
   standard Barrow-HDE to DESI DR2 + SN + CC and obtain **Δ < 0.47–0.54** (Δ = 0 within 1σ, mild positive
   tendency) — which, read literally, disfavours SEDE's Δ = 1. *Defence (shown, not asserted; Fig. 10,
   `run_delta_orthogonality.py`):* both models share the scaling ρ_DE ∝ H^{2−Δ}; the *only* difference is
   SEDE's structure gate f_sat(z), and the gate is exactly what decouples Δ from the observable the
   Barrow-HDE bound rests on. At the **same** Δ = 1 the two models differ sharply: SEDE's gate
   (f_sat → 0 at high z) holds the early-dark-energy fraction to Ω_DE(z=3) = 0.029 — CMB-safe — whereas
   un-gated Barrow-HDE has Ω_DE(z=3) = 0.147 (5× larger) and w₀ = −0.77 with no w = −1 crossing, versus
   SEDE's w₀ = −1.0 in DESI's crossing quadrant. So the Δ-constraint derived from *Barrow-HDE's*
   early-DE/EOS does not transfer: SEDE's Δ = 1 produces viable observables precisely where Barrow-HDE's
   Δ = 1 does not, and in SEDE's own construction Δ = 1 is data-preferred over Δ = 0 (§5.6). The bound
   constrains a sibling model, not SEDE — and Fig. 10 shows the orthogonality rather than asserting it.

![Fig10](../output/fig_delta_orthogonality.png)

**Figure 10.** Why the standard Barrow-HDE bound Δ ≪ 1 does not constrain SEDE's Δ = 1. Both models share
the holographic scaling ρ_DE ∝ H^{2−Δ}; the only difference is SEDE's structure gate f_sat(z). *Left:*
the early-dark-energy fraction Ω_DE(z=3) — the quantity the CMB and the standard-HDE EOS fit actually
constrain — versus Δ. Un-gated Barrow-HDE (red) carries large, strongly Δ-dependent early DE
(Ω_DE(z=3) = 0.15 even at Δ = 1, rising to 0.70 at Δ = 0), so its fits pin Δ < 0.5; SEDE's gate (blue)
drives f_sat → 0 at high z, holding Ω_DE(z=3) ≈ 0.03 — at or below the ΛCDM/CMB-safe level (grey dashed)
— for *all* Δ. At Δ = 1 the two models differ five-fold (0.029 vs 0.147), and SEDE additionally crosses
w = −1 (w₀ = −1.0) where un-gated Barrow-HDE does not (w₀ = −0.77). *Right:* consequently the data prefer
**orthogonal** Δ in the two parametrisations — Δ̂ < 0.5 for Barrow-HDE (EOS/early-DE; Luciano+ 2025) and
Δ̂ ≈ 1.0 for SEDE (amplitude/early-DE; §5.6) — with no contradiction. The gate is the difference.
(`run_delta_orthogonality.py`.)

Net: the two threats that targeted the model itself (full-CMB tension; the Barrow-Δ bound) are answered
by direct tests or a stated model distinction; the two that target the *context* (birefringence; the
low-z SN systematic) are real and partially outside SEDE's control — birefringence especially.

**Evidential status.** SEDE is a falsifiable, fixed-parameter, Λ-free alternative to ΛCDM, mildly
preferred (suggestive ~2–3σ) with an encouraging out-of-sample check, resting on one open postulate
and a few stated theoretical choices. It is stronger than "indistinguishable from ΛCDM" and weaker
than "established." Its value is twofold: a *physical story* (no fine-tuned Λ, w(z) from Ω_m, the
coincidence reframed) and a *decidable prediction* (the volume law, Δ = 1). The data, not the theory,
will settle it.

## 10. Conclusions

We have presented Structural Entropy Dark Energy: a fixed-parameter dark-energy model in which dark
energy is identified with the conjugate thermal energy density of a volume-law-entangled cosmic
horizon and activated through a structure-growth gate. The model is a single Friedmann fixed-point equation with **no fitted dark-sector
parameter** — λ, w(z), γ, and c_s² are each fixed by a stated theoretical choice rather than tuned —
that replaces the cosmological constant by a structure-gated, horizon-coupled term. It **addresses**
the vacuum-energy-scale problem (the scale is the CKN holographic bound, normalised by flatness, not
the QFT vacuum) and reframes the coincidence problem (the activation gate tracks the late-time
maturation of the structure-growth history). In
a marginalised, CAMB-in-the-loop joint analysis it is mildly preferred over ΛCDM (ΔDIC ≈ −2.9); the
preference is suggestive (calibration ~2–3σ, not obviously flexibility) with an encouraging pre-DESI →
DESI out-of-sample check, and its equation of state crosses −1 in DESI's evolving-DE quadrant. Its
quantum-gravitational content reduces to one falsifiable postulate — a volume-law horizon entropy.
The sharp future test is near: under the stated Fisher assumptions, DESI DR3 and Euclid could
constrain the effective deformation Δ to σ ≈ 0.09, separating SEDE's volume-law horizon from a smooth
one at ~11σ. SEDE is suggestive, not established — and, crucially, it is decidable.

## Code and data availability

All code, data vectors, and analysis pipelines are publicly available at
`github.com/spsingularity/SEDE-release`. A single entry point, `reproduce_all.py`, regenerates every number and
figure in this paper from staged drivers; the model and likelihood live in the `sede/` package, and the
verification suite (`sede/verification.py`) checks the quoted results numerically. The observational
inputs are standard public datasets (DESI DR2, Pantheon+/DES-5YR/Union3, Moresco cosmic chronometers,
Gold-2018 fσ8, compressed Planck, eBOSS DR16), retrieved via the documented loaders. A tagged release
is archived at Zenodo, DOI [10.5281/zenodo.21050314](https://doi.org/10.5281/zenodo.21050314).

## Use of AI tools

Artificial-intelligence tools (large language models) were used to assist with drafting and editing the
manuscript, with developing and cross-checking the analysis code, and with literature searches. The
author conceived the model, directed and independently verified all analyses, and takes full
responsibility for the content, including all derivations, results, and conclusions. No AI tool is an
author.

## References

*Cited inline by author and year. All 45 arXiv identifiers were verified to resolve on the live arXiv
API, with the returned titles checked against the cited content (`verify_refs.py`); the 19 post-2024
identifiers — the highest hallucination risk — all resolve to real papers with matching titles.
Journal volume/page details should be confirmed against the journals of record at proof.*

**Foundations — horizon thermodynamics, holography, entropy.**

- Jacobson, T. (1995). Thermodynamics of spacetime: the Einstein equation of state. *Phys. Rev. Lett.* **75**, 1260. arXiv:gr-qc/9504004.
- Cai, R.-G. & Kim, S. P. (2005). First law of thermodynamics and Friedmann equations of Friedmann–Robertson–Walker universe. *JHEP* **02**, 050. arXiv:hep-th/0501055.
- Gibbons, G. W. & Hawking, S. W. (1977). Cosmological event horizons, thermodynamics, and particle creation. *Phys. Rev. D* **15**, 2738.
- Bousso, R. (2002). The holographic principle. *Rev. Mod. Phys.* **74**, 825. arXiv:hep-th/0203101.
- Cohen, A. G., Kaplan, D. B. & Nelson, A. E. (1999). Effective field theory, black holes, and the cosmological constant. *Phys. Rev. Lett.* **82**, 4971. arXiv:hep-th/9803132.
- Page, D. N. (1993). Information in black hole radiation. *Phys. Rev. Lett.* **71**, 3743. arXiv:hep-th/9306083.
- Padmanabhan, T. (2012). Emergence and expansion of cosmic space as due to the quest for holographic equipartition. arXiv:1206.4916.
- Gough, M. P. (2008). Information equation of state. *Entropy* **10**, 150. *(Earliest of the priority lineage — dark energy from the information/entropy of structure; §7.6.)*
- Gough, M. P. (2011). Holographic dark information energy. *Entropy* **13**, 924. arXiv:1105.4461. *(Priority: dark energy as the entropy/information energy of star formation — the structure-sourced-dark-energy idea SEDE shares; §7.6.)*
- Gough, M. P. (2025). Evidence for dark energy driven by star formation: information dark energy? *Entropy* (2025). *(Fits stellar-mass-density history to the DE equation of state; ρ_DE ∝ SMD^{0.5}.)*
- Khlebnikov, S. & Sheoran, A. (2019). The first heat: production of entanglement entropy in the early universe. *Phys. Rev. D* **100**, 103515. arXiv:1907.00487. *(de Sitter → flat: entanglement entropy saturates to a volume law.)*
- Belfiglio, A. & Luongo, O. (2023). Entanglement area-law violation from field–curvature coupling. arXiv:2306.08357.

**Non-extensive / deformed entropy and holographic dark energy.**

- Barrow, J. D. (2020). The area of a rough black hole. *Phys. Lett. B* **808**, 135643. arXiv:2004.09444.
- Tsallis, C. & Cirto, L. J. L. (2013). Black hole thermodynamical entropy. *Eur. Phys. J. C* **73**, 2487. arXiv:1202.2154.
- Saridakis, E. N. (2020). Barrow holographic dark energy. *Phys. Rev. D* **102**, 123525. arXiv:2005.04115.
- Barrow, J. D., Basilakos, S. & Saridakis, E. N. (2021). Big Bang Nucleosynthesis constraints on Barrow entropy. *Phys. Lett. B* **815**, 136134. arXiv:2010.00986. *(Modified-Friedmann route: Δ ≲ 1.4×10⁻⁴.)*
- Asghari, M. & Sheykhi, A. (2022). Observational constraints of the modified cosmology through Barrow entropy. arXiv:2110.00059.
- Leon, G. *et al.* (2021). Barrow entropy cosmology: an observational approach. arXiv:2108.10998.

**Emergent gravity and the dark universe.**

- Verlinde, E. P. (2017). Emergent gravity and the dark universe. *SciPost Phys.* **2**, 016. arXiv:1611.02269 (2016).
- Verlinde, E. P. (2011). On the origin of gravity and the laws of Newton. *JHEP* **04**, 029. arXiv:1001.0785.
- Milgrom, M. (1983). A modification of the Newtonian dynamics as a possible alternative to the hidden mass hypothesis. *Astrophys. J.* **270**, 365.
- McGaugh, S. S., Lelli, F. & Schombert, J. M. (2016). Radial acceleration relation in rotationally supported galaxies. *Phys. Rev. Lett.* **117**, 201101. arXiv:1609.05917.

**Structure formation.**

- Press, W. H. & Schechter, P. (1974). Formation of galaxies and clusters of galaxies by self-similar gravitational condensation. *Astrophys. J.* **187**, 425.
- Sheth, R. K. & Tormen, G. (1999). Large-scale bias and the peak-background split. *Mon. Not. R. Astron. Soc.* **308**, 119. arXiv:astro-ph/9901122.

**Data.**

- DESI Collaboration (Abdul-Karim, M. *et al.*) (2025). DESI DR2 results II: measurements of baryon acoustic oscillations and cosmological constraints. *Phys. Rev. D* **112**, 083515. arXiv:2503.14738.
- Planck Collaboration (2020). Planck 2018 results. VI. Cosmological parameters. *Astron. Astrophys.* **641**, A6. arXiv:1807.06209.
- Brout, D. *et al.* (2022). The Pantheon+ analysis: cosmological constraints. *Astrophys. J.* **938**, 110. arXiv:2202.04077.
- DES Collaboration (2024). The Dark Energy Survey: cosmology results with ∼1500 type Ia supernovae (DES-SN5YR). arXiv:2401.02929.
- Rubin, D. *et al.* (2023). Union3 / UNITY1.5: cosmological constraints from a unified type Ia supernova compilation. arXiv:2311.12098.
- Alam, S. *et al.* (eBOSS Collaboration) (2021). Completed SDSS-IV eBOSS: cosmological implications. *Phys. Rev. D* **103**, 083533. arXiv:2007.08991.
- Eskilt, J. R. & Komatsu, E. (2022). Improved constraints on cosmic birefringence from the WMAP and Planck data. *Phys. Rev. D* **106**, 063503. arXiv:2205.13962.

**Post-DESI DR2 dark-energy literature** (arXiv identifiers verified 2026-06; confirm journal details before submission).

- Turyshev, S. G. (2026). Dark energy after DESI DR2: observational status, reconstructions, and physical models. *Phys. Rev. D* (2026). arXiv:2602.05368.
- Wang, D. & Mota, D. F. (2025). Did DESI DR2 truly reveal dynamical dark energy? *Eur. Phys. J. C* (2025). arXiv:2504.15222.
- Chaudhary, H., Capozziello, S. *et al.* (2025). Evidence for evolving dark energy from DESI DR2 BAO and Pantheon+, DES-Dovekie, and Union3. arXiv:2508.10514. *(LRG1–2-driven preference.)*
- Liu, G., Wang, Y. *et al.* (2024). Impact of LRG1 and LRG2 in DESI 2024 BAO data on dark-energy evolution. arXiv:2407.04385.
- You, C., Wang, D. *et al.* (2025). Dynamical dark energy implies a coupled dark sector: insights from DESI DR2 via a data-driven approach. arXiv:2504.00985.
- Zhang, X., Xu, Y.-H. *et al.* (2025). Reconstruction of dark energy using DESI DR2. arXiv:2511.02220. *(GP w(z); crossing z_wt ≈ 0.46.)*
- Cai, Y., Ren, X., Qiu, T., Li, M. & Zhang, X. (2025). The quintom theory of dark energy after DESI DR2. arXiv:2505.24732.
- Thanankullaphong, P., Sahoo, P. *et al.* (2026). Quintom dark energy: future attractor and phantom crossing in light of DESI DR2. arXiv:2601.02284.

**Critical analyses and tensions (red-team).**

- Wu, P.-J., Li, T.-N., Du, G.-H. & Zhang, X. (2025). Observational challenges to holographic and Ricci dark energy paradigms: insights from ACT DR6 and DESI DR2. arXiv:2509.02945. *(Full-CMB HDE tension; addressed in §9 via the ACT-lensing and plik-lite tests.)*
- Huang, L., Cai, R.-G. & Wang, S.-J. (2025). The DESI DR1/DR2 evidence for dynamical dark energy is biased by low-redshift supernovae. arXiv:2502.04212.
- Luciano, G. G., Paliathanasis, A. & Saridakis, E. N. (2025). Constraints on Barrow and Tsallis holographic dark energy from DESI DR2 BAO data. arXiv:2506.03019. *(Standard Barrow-HDE: Δ < 0.47–0.54.)*
- Diego-Palazuelos, P. & Komatsu, E. (2026). Cosmic birefringence from the Atacama Cosmology Telescope Data Release 6. arXiv:2509.13654. *(β = 0.215°±0.074°, 2.9σ.)*
- Barman, B. & Girmohanta, S. (2026). Implications of DESI for dark matter & cosmic birefringence. arXiv:2506.12589. *(One axion sector: phantom crossing + β + DM — the unified rival to SEDE.)*
- Lin, W., Visinelli, L. *et al.* (2025). Testing quintessence axion dark energy with recent cosmological results. arXiv:2504.17638. *(Axion: β + acceleration at 1σ.)*
- Nakagawa, S., Nakai, Y. *et al.* (2025). Interpreting cosmic birefringence and DESI data with an evolving axion in ΛCDM. arXiv:2503.18924.
- Nakagawa, S., Naokawa, F., Namikawa, T. & Komatsu, E. (2023). Redshift evolution of cosmic birefringence in CMB anisotropies. *Phys. Rev. D* **107**, 083529. arXiv:2301.07971. *(Universal late-time β(z) profile.)*
- Arcari, S., Bartolo, N. *et al.* (2025). Stairway to axions: cross-correlation of birefringence and galaxies from NPIPE and Quaia. arXiv:2509.22273. *(Birefringence tomography.)*
- BICEP/Keck Collaboration (2026). BICEP/Keck XXI: constraints on early-universe parity violation from multipole-dependent birefringence. arXiv:2603.06812.
- Yin, L., Xiong, S., Kochappan, J., Lee, B.-H. & Ghosh, T. (2025). Constraints on cosmic birefringence from SPIDER, Planck, and ACT observations. arXiv:2510.25489. *(α+β degeneracy; Planck–ACT differ >3σ; SPIDER incompatible with the coupling.)*

**Growth–geometry and response probes.**

- Rauhut, S. J., Blake, C. *et al.* (DESI Collaboration) (2025). Testing gravitational physics by combining DESI DR1 and weak lensing datasets using the E_G estimator. arXiv:2507.16098.
- Barreira, A. & Schmidt, F. (2017). Complete super-sample lensing covariance in the response approach. *JCAP* **06**, 053. arXiv:1711.07467.

**Codes and inference.**

- Lewis, A., Challinor, A. & Lasenby, A. (2000). Efficient computation of CMB anisotropies in closed FRW models (CAMB). *Astrophys. J.* **538**, 473. arXiv:astro-ph/9911177.
- Blas, D., Lesgourgues, J. & Tram, T. (2011). The Cosmic Linear Anisotropy Solving System (CLASS) II. *JCAP* **07**, 034. arXiv:1104.2933.
- Foreman-Mackey, D., Hogg, D. W., Lang, D. & Goodman, J. (2013). emcee: the MCMC hammer. *Publ. Astron. Soc. Pac.* **125**, 306. arXiv:1202.3665.
- Torrado, J. & Lewis, A. (2021). Cobaya: code for Bayesian analysis of hierarchical physical models. *JCAP* **05**, 057. arXiv:2005.05290.

## Appendix A — Key results, prescriptions, and assumptions

This appendix documents the six load-bearing items, labelled by their actual epistemic status — one
defining ansatz (A.1), supporting derivations (A.2, A.4, A.5), a fixed prescription (A.3), and a
motivated postulate (A.6). The headings give that status; the parenthetical "Result/Prescription/
Motivation N" numbering corresponds to the repository code labels `THEOREM_N` in `sede/theory.py`
(retained there for code stability, *not* a claim that each is a mathematical theorem). We work in
units c = ℏ = k_B = 1 and write 8πG = 1/M_P². Source statements `THEOREM_1`, `THEOREM_3`,
`THEOREM_4C`, `THEOREM_5`, `THEOREM_8`, `THEOREM_9`; the remaining items (2, 4, 4B, 6, 6B, 10–13, 5D)
are catalogued at the end with one-line statements.

### A.1 The conjugate horizon-fluid ansatz: ρ_DE = T_AH·s_grav

**Claim.** SEDE *identifies* the dark-energy density with the conjugate energy of the horizon
entropy, ρ_DE = T_AH·s_grav·f_sat, with s_grav the (volume-law) horizon entropy density and
f_sat ∈ [0,1]. The Clausius relation supplies the *dynamics*; the absolute density ρ = T·s is the
model's defining **ansatz**, motivated by horizon thermodynamics — not a theorem.

**Setup.** At the flat-FRW apparent horizon R_AH = 1/H, the dynamical Cai–Kim temperature is
T_AH = (H/2π)(1 − ε/2), ε ≡ −Ḣ/H². We adopt the Gibbons–Hawking limit T_AH = H/2π (ε → 0), the leading
temperature appropriate for the de Sitter-attractor dark-energy sector; the dynamical (1 − ε/2)
correction is sub-leading where dark energy matters and is disfavoured in full (§2, §6), so it does not
enter the calibrated E(z). The gravitational entropy is the coarse-grained **entanglement
entropy in the horizon volume**, S_grav ∝ V_AH ∝ R_AH³ (§3.3, §8.1), i.e. a *constant* entropy
density s_grav = s_0 (in contrast to the Bekenstein–Hawking area-law density s_AH = ¾M_P²H ∝ H, which
would instead give the area-law model ρ_DE ∝ H² — the branch disfavoured in the conditional
diagnostics of §5.6).

**Construction.**
*Step 1 (Clausius — the rigorous, dynamical content).* Hayward's first law for a dynamical horizon,
δE = T_AH dS_AH + W δV with work density W = (ρ − p)/2, applied to a shell of matter crossing the
horizon in time dt with flux −dE = (ρ + p)H V dt, gives the Clausius relation (ρ + p)H V = T_AH Ṡ_AH.
This relates *changes* and is standard horizon thermodynamics.

*Step 2 (Bousso bound as motivation for a bounded f_sat).* The covariant entropy bound gives S_matter ≤ A/(4l_P²) = S_AH on
any light-sheet bounded by the horizon. This *motivates* interpreting f_sat as a bounded effective
saturation fraction; in SEDE we impose 0 ≤ f_sat ≤ 1 and f_sat(0) = 1 with the redshift dependence
supplied by the halo prescription (A.2/A.3), and we do **not** identify f_sat literally with
S_matter/S_AH (consistent with the small binding-entropy budget, App. F).

*Step 3 (the SEDE identification — the ansatz).* SEDE posits the dark-energy density to be the
conjugate energy of the volume-law horizon entropy,
$$\rho_{\rm DE} \equiv T_{AH}\,s_{\rm grav}\,f_{\rm sat} = \frac{H}{2\pi}\,s_0\,f_{\rm sat} \;\propto\; H\,f_{\rm sat}.$$
Fixing the one constant s_0 by spatial flatness [ρ_DE(0) = Ω_DE0 ρ_crit,0, f_sat(0) = 1] gives
$$\rho_{\rm DE}(z) = \Omega_{\rm DE0}\,\rho_{\rm crit,0}\,\frac{H}{H_0}\,f_{\rm sat}(z)
= \Omega_{\rm DE0}\,\rho_{\rm crit,0}\,E(z)\,f_{\rm sat}(z),$$
i.e. the fixed-point term Ω_DE0·E·f_sat of eq. (1) — the **volume-law (λ = 1/2)** scaling, *not* the
area-law H². This is an identification motivated by horizon thermodynamics, not a derivation: the
Clausius relation (Step 1) fixes the dynamics, while the absolute density is the postulate (§8.4).

*Step 4 (scale / coincidence).* ρ_DE ~ M_P²H₀² ~ (10⁻³ eV)⁴ is small because H₀/M_P ~ 10⁻⁶¹. The CKN
bound (§8.2) *explains* why this horizon scale is ~ρ_crit — as an upper bound, taken here to be
saturated to O(1) by the flatness normalisation above — replacing the QFT-vacuum conflation behind the
10¹²⁰ problem rather than deriving the magnitude from scratch. ∎

### A.2 Result 3 — the structure gate f_sat(z)

**Claim.** f_sat(z) = [1 − e^{−γ x(z)}]/[1 − e^{−γ}], where x ≡ σ8²(z)/σ8²(0) = D²(z) (D the linear
growth factor, D(0) = 1).

**Setup.** We model the cumulative activation gate directly (dimensionless, avoiding any mix of area-
and volume-law densities) through a normalised growth-weighted deposition kernel
K(t) ∝ (dσ8²/dt) e^{−γσ8²/σ8²(0)} — the structure-growth rate dσ8²/dt > 0 times a Boltzmann factor
e^{−γx} for the exponential rarity of massive halos (γ from A.3) — taking f_sat to be its normalised
cumulative from early times to epoch z (cosmic time t),
$$f_{\rm sat}(z) = \frac{\int_0^{t(z)} K\,dt}{\int_0^{t_0} K\,dt},$$
so that f_sat(∞) = 0 (early, t→0) and f_sat(0) = 1 (today, t = t_0) by construction (the closed-form
gate below is this cumulative-kernel form, not a relaxation ODE).

**Proof.** The kernel integrated over cosmic time is K dt = (dσ8²/dt) e^{−γσ8²/σ8²(0)} dt =
e^{−γu} dσ8² with u ≡ σ8²/σ8²(0); the Hubble factors cancel and the integral reduces to one over u.
With u running from 0 (z′→∞, no early structure) to x(z) ≡ σ8²(z)/σ8²(0) = D²(z),
$$f_{\rm sat}(z) = \frac{\int_0^{x(z)} e^{-\gamma u}\,du}{\int_0^{1} e^{-\gamma u}\,du}
= \frac{1 - e^{-\gamma x(z)}}{1 - e^{-\gamma}}, \qquad x(z) = D^2(z),$$
where the denominator (x = 1, i.e. z = 0) enforces f_sat(0) = 1.
Checks: f_sat(0) = 1; f_sat(z ≫ 1) → 0 as x → 0; f_sat ∈ [0,1] for x ∈ [0,1]. The form is the CDF of a
truncated exponential of rate γ on x ∈ [0,1] — the normalised activation fraction reached
by epoch z. ∎

### A.3 Prescription 4C — the fixed structure coupling γ ≈ 1.496

**Claim.** The entropy weight exponent is p = 5/3, giving γ = d ln Σ_S^{(5/3)}/d ln σ8² = 1.496
(the derivative is w.r.t. the variance σ8², the gate argument x = D² = σ8²/σ8²(0) of A.2; equivalently
½ d ln Σ_S/d ln σ8 = 1.496, i.e. d ln Σ_S/d ln σ8 = 2.993). A self-similar reduction gives the closed
form γ = (p−1)⟨1/α⟩_Σ with α ≡ −d ln σ/d ln M the effective spectral slope over the entropy-weighted
mass range; it reproduces the exact mass-function integral to 0.1% (`run_gamma_systematics.py`,
`BINDING_ENERGY_MECHANISM.md` §7). Entropy is deposited by the ν_S ≈ 2.4 (~2σ) clusters at
M_S ~ 10¹⁴·⁵ M⊙/h, giving an independent cluster-abundance handle; p = 1 gives ≈ 0.27 (∼5× below
p = 5/3), a weak gate versus the real one.

**Proof.** In the SEDE horizon-fluid ansatz (A.1), f_sat is the effective fraction of the
**cosmic-horizon** entropy capacity that is gravitationally activated, and ρ_DE = T_AH·s_grav·f_sat —
so the entropy that defines the f_sat weighting is accounted in the *horizon* reservoir at the *horizon* temperature
T_AH = H/2π. By Clausius, heat Q deposited into a
reservoir at temperature T adds entropy Q/T. When a halo of mass M virialises it releases gravitational
binding energy E_bind ∝ M^{5/3} (NFW; R_vir ∝ M^{1/3}; verified slope 1.64). Under the SEDE
prescription this released binding energy is accounted in the horizon entropy ledger at the horizon
temperature T_AH (we do not claim a demonstrated transport mechanism), so the entropy assigned in this ledger is
$$\Delta S_{AH} = \frac{E_{\rm bind}}{T_{AH}} \propto \frac{M^{5/3}}{H}.$$
Because T_AH ∝ H is identical for every halo at a given epoch (mass-independent), the per-halo mass
weight is M^{5/3}: **p = 5/3**. With Σ_S = ∫ M^{5/3}(dn/d ln M)d ln M and a Sheth–Tormen mass function
(COLOSSUS, EH98 transfer, M ∈ [10¹⁰, 10¹⁶] M_⊙),
$$\gamma = \frac{d\ln\Sigma_S^{(5/3)}}{d\ln\sigma_8^2} = \tfrac{1}{2}\,\frac{d\ln\Sigma_S^{(5/3)}}{d\ln\sigma_8} = 1.4964 \approx 1.50$$
(the derivative is w.r.t. the variance σ8² because the gate's argument is x = D² = σ8²/σ8²(0); equivalently d ln Σ_S/d ln σ8 = 2.993).
The earlier p = 1 (Result 4B) divided E_bind by the *halo's* virial temperature T_vir ∝ M^{2/3}, which is
the entropy of the halo's own internal state — correct for the halo, wrong for the entropy added to the
cosmic horizon (where f_sat lives). The heat is the same; the reservoir, and hence the temperature, is
the horizon's. ∎

### A.4 Result 5 — the structural-EOS diagnostic (the w₀ side-note, not canonical w(z))

**Claim (structural EOS).** w_struct(0) = −1 + γ/[3(e^γ − 1)], well-approximated by the closed-form
algebraic expression w_alg = (4Ω_m/3 − 1)/(1 − Ω_m). (This is the side-note structural estimate of
§3.2, not the model's prediction, which is the canonical crossing w(z).)

**Proof.** From ρ_DE(z) = ρ_crit,0 Ω_DE0 f_sat(z) (A.1) with f_sat(x), x = D², the structural EOS in
the variance variable is w_struct(0) = −1 + (1/3) df_sat/dx|_{x=1}. Differentiating the Theorem-3 form,
$$\left.\frac{df_{\rm sat}}{dx}\right|_{x=1} = \frac{\gamma e^{-\gamma}}{1 - e^{-\gamma}} = \frac{\gamma}{e^\gamma - 1},$$
so w_struct(0) = −1 + γ/[3(e^γ − 1)]. At γ = 1.4964 this is 0.432, giving **w_struct(0) = −0.856**
(0.33σ from DESI DR2). Numerically γ/(e^γ − 1) = 0.432 ≈ Ω_m/(1 − Ω_m) = 0.451 at the observed
parameters, so the structural value is approximated to ~4% by the closed-form expression
w_alg = (4Ω_m/3 − 1)/(1 − Ω_m) = −0.849 (0.21σ from DESI). This is a numerical approximation in the
x = D² convention, **not** an exact identity. (The interacting-fluid effective EOS, from
Raychaudhuri + Friedmann at z = 0, is phantom, w_fluid(0) ≈ −1.15, because f_sat is still rising today;
the canonical λ = 1/2 background realises this as a crossing of −1, §5.1.) ∎

### A.5 Result 8 — the H-coupling λ = 1 − Δ/2

**Claim.** With the Barrow entropy S = (A/A₀)^{1+Δ/2}, the conjugate horizon-fluid ansatz gives
ρ_DE ∝ H^{2λ} f_sat with λ = 1 − Δ/2.

**Proof.** A quantum-gravitationally rough horizon carries the Barrow entropy S = (A/A₀)^{1+Δ/2},
Δ ∈ [0,1], in place of the area law. With A = 4π/H² and V_AH ∝ H^{−3}, the **Barrow-deformed** entropy
density (written s_Δ to distinguish it from the Bekenstein–Hawking area-law density s_AH ∝ H) is
$$s_{\Delta} = \frac{S}{V_{AH}} \propto \frac{(H^{-2})^{1+\Delta/2}}{H^{-3}} = H^{1-\Delta},$$
$$\rho_{\rm DE} = T_{AH}\,s_{\Delta}\,f_{\rm sat} \propto H\cdot H^{1-\Delta}\cdot f_{\rm sat}
= H^{2-\Delta}f_{\rm sat} \equiv H^{2\lambda}f_{\rm sat}, \qquad \boxed{\lambda = 1 - \tfrac{\Delta}{2}.}$$
This is Barrow holographic dark energy with the Hubble IR cutoff (ρ_DE ∝ L^{Δ−2}, L = 1/H). The
endpoints are Δ = 0 ⟹ λ = 1 (the area-law branch, disfavoured) and Δ = 1 ⟹ λ = 1/2 (the volume-law
postulate, A.6). λ is therefore not a free parameter; it is set by the Barrow deformation. **Scope:** the
Barrow exponent multiplies only the f_sat-gated ρ_DE, *not* the gravitational sector, so H(z) at BBN
equals the standard matter+radiation value to ~8 figures (Ω_DE/ρ_tot ~ 10⁻¹⁶ at z ~ 4×10⁹); the
modified-gravity BBN bound Δ ≲ 10⁻⁴ does not apply, and constant Δ = 1 is BBN-safe. ∎

### A.6 Motivation 9 — Δ = 1 from a geometric ceiling and the GSL

**Claim.** Δ = 1 (the volume-law endpoint) is *strongly motivated* — not fitted to cosmological data —
by two independent principles; it is the central postulate (§8.4), not a microscopic derivation.

**Pillar 1 — geometric ceiling.** The Barrow deformation is the fractal (Hausdorff) dimension of the
horizon surface, d_H = 2 + Δ. A 2-surface embedded in 3-space cannot exceed the space-filling dimension
d_H = 3, so **Δ ≤ 1**, with Δ = 1 the extremum where the area-measure becomes a volume-measure
(S ∝ A^{3/2} ∝ R³). This is Barrow's defining range.

**Pillar 2 — thermodynamic attractor.** Three premises: (P1) the DE sector is purely entropic — SEDE
posits ρ_DE = T_AH s_grav, so the free-energy density F = ρ_DE − T_AH s_grav = 0 identically; the
governing potential is therefore the *entropy*, not F = E − TS (this is the load-bearing step: with a
nonzero Δ-dependent free energy, equilibrium would extremise F and Δ = 1 would not follow). (P2) the
generalised second law: S_total = S_Δ + S_matter is non-decreasing, S_matter independent of Δ. (P3) Δ is
a thermalising horizon degree of freedom (the standard Jacobson–Padmanabhan assumption). By (P1)–(P2),
maximising S_total over Δ reduces to maximising S_Δ = (A/A₀)^{1+Δ/2}, whose derivative is
$$\frac{dS_\Delta}{d\Delta} = \frac{S_\Delta}{2}\ln\!\frac{A}{A_0} + [\text{back-reaction through }A(\Delta)].$$
The explicit term carries ln(A/A₀) ≈ ln(10¹²²) ≈ 280 > 0 (a vastly super-Planckian horizon). The
back-reaction (Δ shifts ρ_DE ∝ H^{2−Δ}, hence H, hence A) cannot flip the sign: it *vanishes* at z = 0
(H = H₀ ⟹ ln H = 0 ⟹ ∂_Δ H^{2−Δ} = 0) and is ~1/280-suppressed elsewhere. Hence dS_Δ/dΔ > 0 on all of
[0,1]; a strictly increasing function on a closed interval attains its maximum at the right endpoint.
With the ceiling Δ ≤ 1 (Pillar 1) the constrained maximum is uniquely **Δ = 1**, and by (P3) a
thermalising horizon relaxes to it.

**Conclusion.** Under the additional assumptions that Δ is a thermalising horizon degree of freedom
(P3) and that entropy maximisation over the Barrow family is the correct variational principle (P1),
the GSL argument is a **consistency check** — the endpoint Δ = 1 is a GSL *fixed point* — rather than a
selection principle: the geometric ceiling d_H = 3 supplies the bound Δ ≤ 1, and saturation of that cap
*is* the volume-law postulate in effective form (via λ = 1 − Δ/2 this gives λ = 1/2). These assumptions
**motivate but do not microscopically derive** the endpoint (§8.4). Because S_Δ is monotone
in Δ, the same argument favours Δ = 1 at all epochs (the adopted constant Δ). The current conditional
diagnostics prefer values near Δ(0) ≈ 1.0–1.1 (§5.6), consistent with the endpoint. This is a strong
motivation, not a derivation: it reduces Δ = 1 from a fitted coupling to a GSL-equilibrium extremum
under the stated premises, leaving the microscopic origin open. ∎

### A.7 The remaining repository items (catalogue)
**Result 2** (reheating master equation): the same logistic df_sat/dt = H f_sat(1 − f_sat) − Γ f_sat governs
the *early* f_sat (inflaton→radiation), giving the U-shape of Result 12. **Result 4 / 4B** [superseded by 4C]:
the M^{5/3} energy weight, and the p = 1 virial-temperature mis-step. **Result 5D**: the (1 − ε/2) "λ = 1
cousin" background (w₀ = −0.85, no crossing) vs the canonical λ = 1/2 crossing background. **Result 6 / 6B**:
the smooth-DE growth validation and the c_s² = 1 closure (f_sat a background functional). **Result 10**: the QG
chain and the volume-law postulate. **Result 11**: the seam = BBN bound = modgrav/holographic fork.
**Result 12**: the f_sat U-shape (inflation and dark energy as the two de Sitter brackets). **Result 13**:
1 + w = (1/3)[2λε − d ln f_sat/d ln a], the EOS as horizon-entropy production. Full statements in
`sede/theory.py`; each is exercised by a numbered check in App. B.

## Appendix B — The verification suite
The repository ships a `sede/verification.py` suite of 63 numerical checks (V1–V63) covering: the
fixed-point solution (closed-loop H(z) to 10⁻¹³ over z ∈ [0, 10⁶]), BBN-safety (ΔN_eff = ΔY_p = 0),
the CAMB≡CLASS r_drag agreement (5×10⁻⁵), the CLASS perturbation validation (0.1%), the GSL along
cosmic history, the fixed (derived, Thm 4C) γ and the smooth-DE c_s² = 1 closure, the horizon-dimension
scaling λ(d_H) = 2 − d_H/2 (V59), the §8.5 state/counting diagnostics (SYK scrambling V60; causal-set
dof counting in Minkowski + de Sitter V61; the entropy-bound audit V62; the tensor-network
realisability test V63), and the cross-horizon and quantum-sector consistency checks.
`reproduce_all.py` regenerates every number and figure in this paper from the staged drivers.

## Appendix C — Data vectors, covariances, and the fσ8 audit
The full data vectors and covariance matrices, the positive-definiteness checks, the benign
Pantheon+ duplicate-CID treatment, and the fσ8 data audit (the small set of mis-entered growth points
corrected to their published Gold-2018 values, and the SH0ES value correction) — all reproduced from
`data/` and documented in `DATA_AUDIT.md`.

## Appendix D — False-preference calibration methodology
The parametric-bootstrap procedure of §5.2: generation of ΛCDM-truth mocks with **all** probes
(including the compressed CMB R/l_A and SH0ES) re-drawn self-consistently; the refit of both models to
each mock; the null distribution of Δχ²(SEDE − ΛCDM) (mean +0.37, 2000 mocks); the p ≈ 0.004 (95% UL 0.006); the
documented bug (CMB held fixed → biased null −2.87, p ≈ 0.5) and its fix; and the sibling-pipeline
cross-check (~2σ, CMB-distance-engine-dependent). The pre-registration sha256 and contents are listed.
The standalone `PREREGISTRATION.md` additionally records the triggered CMB follow-ups and their status
(§E: full non-lite plik **run**, Δχ²=−0.33; ACT DR6 primary formally deferred) and the locked E1
structure-sourcing decision rule (§F: data vector + CONFIRM/REFUTE thresholds). The no-SH0ES robustness
(§5.4 ii′; Δχ²=−1.75, ln B=+0.95, p=0.030) is reproduced by `run_no_shoes_robustness.py`.

## Appendix E — The QG derivation ledger
The full chain from §8 with each link labelled by status:
- **Jacobson–Cai–Kim (MOTIVATED / ANSATZ):** Clausius horizon thermodynamics motivates the
  horizon-fluid construction and supplies the conservation/dynamical logic (δQ = T dS ⟹ Friedmann);
  the *absolute* identification ρ_DE = T·s is the SEDE ansatz (A.1), not a theorem.
- **Non-extensivity ⟹ Tsallis–Cirto/Barrow power law (MOTIVATED):** the weakest link.
- **Second-law / GSL ⟹ volume-law endpoint Δ = 1 (MOTIVATED, given the volume-law thermalisation
  premise):** the geometric ceiling caps Δ ≤ 1 and the GSL is *consistent with* (a fixed point at) the
  saturated endpoint — a consistency check, not a selection principle; a microscopic derivation
  remains **OPEN**.
- **CKN scale (RIGOROUS BOUND, normalised by flatness):** the holographic energy bound gives the horizon-scale upper
  bound; the amplitude is normalised by flatness (§8.2).
- **The one irreducible postulate (OPEN):** a microscopic state count delivering the volume law.

The perturbative-QG (log-correction) vs thermalised (volume-law) regimes are reconciled as different
*states*, not a contradiction. This ledger matches §8.4 and Appendix F.

## Appendix F — The energy and entropy budget
This appendix records the back-of-envelope that the energy-conservation framing of §1–§2 turns on, and
the resulting ledger of what is derived, normalised, and open.

**Energy budget.** Dark energy today is ρ_DE c² ≈ 0.7 ρ_crit c² ≈ 6×10⁻¹⁰ J/m³. The gravitational
binding energy released by structure is bounded by |E_bind|/Mc² ~ σ²/c² per virialised system, with
σ ~ 10³ km/s for clusters (σ²/c² ~ 10⁻⁵) and ~200 km/s for galaxies (~10⁻⁷); even taking all of Ω_m
collapsed, ρ_bind c² ≲ 10⁻⁵ · 0.3 · 8.5×10⁻¹⁰ ~ 10⁻¹⁵ J/m³. Hence ρ_bind/ρ_DE ~ 10⁻⁶: **binding
energy is six orders of magnitude too small to *be* the dark energy.** A literal matter→dark-energy
transfer is excluded independently — moving O(ρ_crit) out of matter over a Hubble time would spoil
ρ_m ∝ a⁻³ and with it BBN and the CMB.

**Entropy budget.** Deposited at the cold horizon temperature, the binding entropy is ΔS_AH =
E_bind/T_AH; with T_AH = ℏH/2π this comes to ΔS_AH ~ 10⁻⁷ S_AH, i.e. ~10⁻⁷ of the horizon entropy.
So binding-entropy deposition does not, by itself, fill f_sat to O(1); it sets the *relative*
z-dependence (the weight p = 5/3 → γ), while the normalisation f_sat(0) = 1 is set by the horizon
capacity. This is the same family as the (M_P/H)^Δ ~ 10⁶¹ "seam" (Result 11), resolved by the
holographic-DE scope (Barrow acts only on the f_sat-gated dark-energy fluid).

**The ledger.**

| quantity | status | set by |
|---|---|---|
| magnitude ρ_DE ~ ρ_crit | **bounded + normalised** | CKN holographic upper bound (§8.2), amplitude fixed by flatness (O(1) saturation of the bound assumed) |
| H-scaling λ = 1/2 | **postulate** | volume-law horizon entropy S_grav ∝ V_AH (§8.1) |
| shape γ ≈ 1.5 | **derived (Theorem 4C)** | halo-binding-energy weighting E_bind ∝ M^{5/3} → p = 5/3 (§3.1) |
| sound speed c_s² = 1 | **closure** | minimal smooth-DE treatment, f_sat a background functional (§3.4) |
| normalisation f_sat(0) = 1 | **imposed** | present-day horizon-entropy normalisation |
| volume-law horizon (Δ = 1 microscopically) | **open** | the one irreducible postulate (§8.4) |

The honest reading: none of these is a parameter *fitted to cosmological data*, but neither are they
all "derived from first principles" — the magnitude is a holographic bound normalised by flatness (as
ΛCDM normalises Ω_Λ0), the H-scaling is the central volume-law postulate, γ is a fixed halo
prescription, and c_s² = 1 is a modelling closure. Binding energy enters SEDE entropically, fixing the
*shape* of the gate; the dark energy's *magnitude* is the horizon's own thermal scale, not energy
borrowed from structure. The single irreducible foundational input is the microscopic volume-law
postulate.

## Appendix G — Supporting motivations and numerical probes for the volume-law postulate

This appendix collects the material demoted from §8 to keep the main line focused: three further
consistency *motivations* for the volume endpoint (§8.1) and the *numerical probes* of the state and
counting halves of the postulate (§8.5). **None is a derivation**; each is a consistency check that
makes the postulate less *ad hoc* without selecting it. (The CHR near-critical layer, the cross-horizon
black-hole reading, and the Verlinde dark-matter connection are developed in a separate companion paper.)

### G.1 Three further consistency motivations for the volume endpoint
*(supplementing the geometric ceiling, the CKN scale/form split, and the direct area-law disfavouring
of §8.1)*

- **A driven steady state.** The horizon need not self-thermalise — its scrambling time
  (∼ ln S / H ≈ 280/H) far exceeds its age. Instead the structure-growth gate (the rise of f_sat)
  *drives* it off the ground state; with a scrambling-limited relaxation rate the volume-law state is
  *maintained* as a non-equilibrium steady state throughout structure formation (the long scrambling
  time becomes an asset). The volume law is then a fixed point of SEDE's own dynamics, not an external
  input.
- **A roughening universality class** (developed with the sibling SEDE_V2 collaboration). Modelling the
  horizon as a stochastic growing surface, ∂h/∂t = ν∇²h + (λ_K/2)(∇h)² + η, the deformation is fixed by
  universality class: the thermal/linear Edwards–Wilkinson class (λ_K = 0) flows to the volume endpoint
  and corresponds to SEDE's purely-entropic F = 0 condition, whereas a generic nonlinear (KPZ) surface
  gives a sub-volume value. This is the most speculative of the motivations and is included only as such.
- **An independent emergent-gravity proposal** (Verlinde 2016). In Verlinde's programme the de Sitter
  vacuum carries a **volume-law** contribution to its entanglement entropy, S ∝ V_H, in addition to the
  horizon area law — the same super-area scaling SEDE adopts at Δ = 1. A second, independently-motivated
  route arriving at a volume-law entropy makes the postulate less *ad hoc*; it remains a motivation, not
  a proof, since Verlinde's volume-law derivation is itself heuristic.

### G.2 The state is maximally entangled (Majorana-SYK diagonalisation)
The *state* half of the postulate (§8.5) — that the horizon dof are maximally entangled (thermal),
not in a low-entanglement ground state — is settled and not peculiar to SEDE. Exact diagonalisation of
the Majorana Sachdev–Ye–Kitaev model [Sachdev & Ye 1993; Kitaev 2015], the canonical chaotic,
black-hole-dual system saturating the Maldacena–Shenker–Stanford bound λ_L ≤ 2πT [arXiv:1503.01409]
(`run_syk_scrambling.py`, V60), confirms: the gap ratio reproduces Wigner–Dyson statistics with the
correct N-mod-8 symmetry class [García-García & Verbaarschot 2016; Cotler et al 2017; Atas et al 2013],
eigenstates saturate the Page value [Page 1993; Vidmar & Rigol 2017], and the OTOC scrambles completely
— while an integrable control does none of these. But this **does not fix Δ**: a black hole is a
*maximal* scrambler yet **area-law** (S = A/4, reproduced by string microstate counts), and SYK, being
all-to-all (geometry-free), cannot even pose the counting question.

![FigG1](../output/fig_syk_scrambling.png)
**Figure G1.** SYK as the horizon scrambler (q=4 chaotic vs q=2 integrable control). *Left:* gap ratio
climbs through the N-mod-8 RMT classes (GOE→GUE→GSE) — the certificate of maximal chaos. *Centre:*
eigenstate entanglement saturates the Page value. *Right:* the OTOC scrambles completely (C→2). This
fixes the *state*, not the counting.

### G.3 The counting question — numerical probes (allowed, realisable, not selected)
The *counting* half (§8.5) — does the horizon dof number grow as area (N ∝ R^{d-2}, Δ=0) or bulk
volume (N ∝ R^{d-1}, Δ=1)? — is where SEDE's postulate lives. Three numerical probes
(`run_entropy_bounds.py`, `run_causal_set.py`, `run_tensor_network.py`; V62, V61, V63):

| question | probe | result |
|---|---|---|
| **Allowed?** by the established entropy bounds | volume entropy vs the smooth-area bounds | S ∝ R³ overshoots the holographic/Bousso/Bekenstein bound A/4 by ∝ R (~10⁶¹ at the cosmic horizon); consistent **only** as a genuine d_H=3 fractal horizon (true area ∝ R³ ⟹ S = A_fractal/4 *exactly* saturates the bound) |
| **Selected?** by quantum gravity | causal-set sprinkle, Minkowski **and** de Sitter | both scalings present; the canonical horizon object — causal *links* — is **area** (R^{d-2}, the Bekenstein–Hawking recovery), in flat space and, because dS is conformally flat, identically in de Sitter. **Volume is not the default.** (measured: bulk count 1.92≈2 ⟹ Δ=1; link count 0.98≈1 ⟹ Δ=0, in 2+1D) |
| **Realisable?** in a holographic code | Ryu–Takayanagi min-cut, local vs nonlocal | local (lattice) connectivity ⟹ area law (exponent 1.0); **nonlocal (expander) connectivity ⟹ volume law** (exponent 1.9). Δ=1 is realisable; its required ingredient is *nonlocal* horizon connectivity (the all-to-all structure of SYK). |

The pattern: volume-counting is **allowed** (not by evading the bounds but by the horizon being
genuinely space-filling — the postulate itself), **realisable** (a nonlocal holographic code produces
S ∝ R^{d-1} exactly), but **not selected** (the frameworks we can compute in default to area for the
canonical horizon object, and a black hole is area-law). No consistency argument available to us
*derives* the volume count; equally, none *excludes* it. The two deciders are stated in §8.5: the dS
static-patch Hilbert-space dimension (theoretical) and the deformation Δ (empirical, §5.6/§6).
