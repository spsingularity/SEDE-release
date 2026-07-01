# Is the volume-law horizon a postulate or a consequence? A driven non-equilibrium reduction of Structural Entropy Dark Energy

*(Companion to "Structural Entropy Dark Energy" (Paper I) — foundations/quantum-gravity
target, prepared for submission. All quantitative claims are reproduced by the scripts named
inline and collected in `reproduce_all.py`, each with validation assertions; epistemic tags —
DERIVED / REDUCED / NEGATIVE / OPEN — flag the status of each step. Figures 1–5 have sources in
`output/`, regenerable via `make_foundations_figures.py`; references are cited inline by author
and year. The cosmological results of Paper I stand independently of the reduction attempted here.)*

## Abstract

Structural Entropy Dark Energy (SEDE) reproduces the cosmological data with no fitted
dark-sector parameter, at the price of a single foundational input: the cosmic apparent
horizon carries the coarse-grained entanglement entropy of its **volume**, S ∝ V ∝ R³
(the Barrow endpoint Δ = 1), rather than the Bekenstein–Hawking area law. We ask how much
of that input is genuinely irreducible. We show the postulate decomposes into four
sub-claims — *state, form, scale, count* — of which the first three are first-principles or
reducible (a maximal scrambler is volume-law-entangled; volume-law entanglement is the
generic behaviour of a thermal reference state; the magnitude is the Cohen–Kaplan–Nelson
scale), leaving the **count** (does the horizon entropy count bulk or boundary degrees of
freedom?) as the sole residue. A Bekenstein no-go proves the count cannot follow from any
energy bound. We then reduce the count to a dynamical statement: the count is fixed by the
*connectivity* of the horizon degrees of freedom, connectivity by the *range* of the
interaction, and gravity (1/r, strongly long-range) places the horizon in the non-additive
(volume) class — with area-law arising only for an isolated, *equilibrium* horizon (a black
hole), while the cosmic horizon is continuously *driven* by structure formation. This
converts the black-hole/cosmic asymmetry into a discriminator and yields a **driven
non-equilibrium-steady-state (NESS) conjecture**: a strongly-long-range, structure-driven
horizon self-organises to its area↔volume spinodal and locks the volume branch. We supply
both microscopic maps the conjecture requires — the cooperative coupling is the gravitational
binding itself (J = λ_max(W_grav) ≥ J_c, generic at the horizon's N), and the drive is the
structure-deposited entropy (clearing the barrier at z ≈ z*) — and verify the volume-lock is
robust over a broad, untuned basin. Finally we make "closure" precise: three routes
(double-scaled SYK, Euclidean de Sitter saddles, a thermal free energy) converge on the same
irreducible statement, the de Sitter static-patch state count, with a decisive empirical
settler (Δ via DESI DR3 + Euclid) already in hand. The volume-law postulate is thereby
narrowed from an unconstrained assumption to a single, sharply-posed, falsifiable statement.

---

## 1. Introduction

Horizon thermodynamics is the most robust pointer we have toward the microscopic structure of
gravity. The first law on a causal horizon yields the Einstein equations [Jacobson 1995], the
apparent-horizon first law yields the Friedmann equations [Cai & Kim 2005], and the
Bekenstein–Hawking area law S = A/4 organises black-hole and cosmological thermodynamics
alike. Yet the area law is an assumption about a *state*: it is the entanglement entropy of a
particular (ground-like) horizon state, and generalisations — Tsallis–Cirto non-extensive
entropy, Barrow's fractal-horizon S ∝ A^{1+Δ/2} [Barrow 2020], and the holographic
dark-energy programme built on them [Saridakis 2020] — explore what happens when that
assumption is relaxed.

Structural Entropy Dark Energy (SEDE; Paper I) is one such model, distinguished by sourcing
dark energy from horizon entropy *gated by structure growth*, ρ_DE = T_AH s_grav f_sat(z).
Its appeal is that the entire dark sector is parameter-free once Ω_m and one input are fixed:
the H-coupling λ = 1−Δ/2 = ½, the equation of state, the structure coupling γ, and the sound
speed all follow. The one input is the **volume-law postulate**: the horizon entropy is the
volume-law (Δ = 1) entanglement entropy, not the area law. Paper I states plainly that this
is the model's single irreducible assumption, motivated but not derived, and that it is
falsifiable — Δ is measured by upcoming surveys.

This paper asks the prior question: *how irreducible is it, really?* We are not trying to
solve quantum gravity; we are trying to determine the minimal, sharply-posed statement that
the cosmology actually rests on, and to derive everything reducible to it. The result is a
chain that trades a static, untestable counting postulate for a dynamical, falsifiable
self-organisation statement, with a clean hand-off to an open — but well-posed — problem in
de Sitter holography. Several steps borrow established machinery: the Cohen–Kaplan–Nelson
holographic energy bound [Cohen, Kaplan & Nelson 1999] for the scale; the Page curve and
eigenstate thermalisation for the state; the Campa–Dauxois–Ruffo classification of long-range
systems [Campa, Dauxois & Ruffo 2009] for the connectivity; Ryu–Takayanagi min-cut for the
count↔connectivity link; and the double-scaled-SYK/de Sitter correspondence
[Susskind 2021; Narovlansky & Verlinde 2023] for the closure attempt. The contribution is the
chain, and the precise localisation of what remains open.

Throughout, every quantitative claim is a runnable experiment; the inline `run_*.py` names and
the ledger in §8 give the correspondence, and `reproduce_all.py` runs them all with validation
assertions.

![Fig1](../output/foundations_fig1_reduction.png)

**Figure 1.** The reduction at a glance. The volume-law postulate splits into *state*, *form*,
*scale* and *count*; the first three are derived or reduced (green/blue) and the count is the residue
(orange). The count reduces through connectivity and interaction range to a driven non-equilibrium
steady state, ratcheted through the spinodal at z* and locked at volume, leaving a single
dS-holography question that the Δ measurement decides empirically. Each box names the experiment that
establishes it.

### 1.1 Relation to prior work

SEDE sits at the confluence of several programmes, and its closest relative deserves explicit
positioning.

**Entropy-production cosmic acceleration (GREA).** The General Relativistic Entropic Acceleration
theory of García-Bellido and Espinosa-Portalés [arXiv:2106.16014, 2106.16012] is the nearest
existing model: it too abandons Λ, sources the late-time acceleration from cosmic-*horizon*
entropy, takes the Helmholtz free energy F = U − TS (not matter alone) as the gravitational
source, and makes the breaking of time-reversal by *entropy production* the engine of
acceleration. We regard GREA as prior art for "acceleration from horizon entropy production
without Λ", and our driven-NESS picture (§4–6) is conceptually continuous with it. The
differences are sharp and testable: (i) GREA uses the **boundary (area)** horizon entropy and a
single parameter α (the horizon-to-curvature ratio), whereas SEDE's defining claim is the
**volume** entropy (Δ = 1) — precisely the counting residue this paper isolates; (ii) GREA's
entropy growth is driven by the horizon's own *expansion*, whereas SEDE's is *gated by structure
formation* (the f_sat factor), which locks w(z) to the growth history and yields a phantom
crossing rather than GREA's quintessence-like w(a); (iii) SEDE's dark sector is parameter-free.
Notably, were SEDE's residue to resolve to the *area* branch (Δ = 0), its phenomenology would
collapse toward GREA's — so the Δ measurement (§7) that decides SEDE's residue also discriminates
the two models. A holographic reading of GREA has recently appeared [arXiv:2511.19546], and a
broader open-system view of gravitational non-equilibrium thermodynamics is developing
[e.g. arXiv:1111.6510, 2507.03103]; SEDE's contribution to that line is the reduction of the
volume-law *counting* input to a driven steady state. The connection is also quantitative: GREA's
entropy production maps to an effective bulk viscosity with negative pressure, and SEDE's
structure-driven entropy production d ln f_sat/d ln a > 0 maps to a *positive* (second-law-respecting)
effective bulk viscosity ζ(z) = ρ_DE (d ln f_sat/d ln a)/(9H²) that peaks at the structure-formation
epoch rather than tracking the horizon expansion (`run_lit_inferences.py`); the w = −1 crossing is the
GREA negative-pressure balance, structure-timed. SEDE is, in this sense, a *structure-gated* member of
the GREA/entropic-dark-energy family — the gate being exactly what makes its w(z) lock to the growth
history (the §5 phase-lock) rather than to the horizon's own size.

![Fig5](../output/foundations_fig5_inferences.png)

**Figure 2.** Two consequences drawn from the new literature. *(a)* The de Sitter horizon has negative
specific heat (C = −2S); with gravity strongly long-range this makes the ensembles inequivalent and the
canonical (area-law) accounting ill-defined, selecting the microcanonical (volume) branch (§4). *(b)*
SEDE's structure-driven entropy production maps to a positive (second-law) effective bulk viscosity
ζ(z) that peaks at the gate-activation epoch (z ≲ 0.5; it tracks d ln f_sat/d ln a, not the z≈2 star-formation peak), casting SEDE as a structure-gated member of the
GREA/entropic-dark-energy family; w(z) crosses −1 at the viscous–dilution balance.

**Holographic dark energy and generalised entropy.** SEDE is, at the background level, a
structure-gated Barrow holographic dark energy: the CKN bound [Cohen, Kaplan & Nelson 1999] fixes
the scale, Barrow's fractal-horizon entropy [Barrow 2020] and its holographic-DE use
[Saridakis 2020] fix the form, and the Tsallis–Cirto non-extensive entropy [Tsallis & Cirto 2013]
supplies the δ ↔ Δ map. Our addition is the *derivation* question — which of these inputs is
irreducible — and the answer that only the count is.

**Horizon thermodynamics and entanglement.** The first-law derivations of Einstein
[Jacobson 1995] and Friedmann [Cai & Kim 2005] equations underlie the construction; the
entanglement-equilibrium refinement [Jacobson 2016] is the basis of route B; the Ryu–Takayanagi
min-cut [Ryu & Takayanagi 2006] supplies the count ⟺ connectivity link; and the volume-law =
thermal identification rests on eigenstate thermalisation [Deutsch 1991; Srednicki 1994;
Rigol et al. 2008].

**de Sitter holography.** The residue — "dim 𝓗(static patch) = e^{Area} or e^{Vol}?" — is posed in
the language of finite-dimensional de Sitter holography [Banks 2001] and is the target of the
double-scaled-SYK/de Sitter correspondence [Susskind 2021; Narovlansky & Verlinde 2023;
Lin 2022], which we use in the closure attempt (§7). The maximal-scrambler state invokes the
chaos bound [Maldacena, Shenker & Stanford 2015] and SYK random-matrix universality
[Cotler et al. 2017].

**Many-body and complexity inputs.** The long-range/non-additivity classification is
[Campa, Dauxois & Ruffo 2009]; the self-organisation framing is Bak–Tang–Wiesenfeld
self-organised criticality [Bak, Tang & Wiesenfeld 1987], whose cosmological applications to date
concern the *primordial* universe and the dark-matter sector rather than the late-time horizon.
The prototype's free-fermion *negative* (volume-law needs interactions, not free transport) is
consistent with the monitored free-fermion entanglement literature [Fisher et al. 2023], and that
genuine **volume-law non-equilibrium steady states** exist is established in driven many-body
systems [Ippoliti, Rakovszky & Khemani 2022] — a condensed-matter precedent for, not a derivation
of, the driven-NESS we invoke.

**Structure-sourced dark energy.** The general idea that dark energy is sourced by cosmic structure
predates SEDE in several distinct mechanisms: Gough's information dark energy (≈2008); the
spacetime-fragmentation/"FRW-islands" proposal of Hossain (2007); and the cosmological-backreaction
and averaging programme (Buchert; Wiltshire's timescape). SEDE Paper I credits and distinguishes
these. What is specific to SEDE — and, to our knowledge, without precedent — is the *combination*
the present paper analyses: the source is the **volume-law horizon entropy** (Δ = 1, a counting
claim), it is **driven/self-organised by structure** through the f_sat gate, and the dark sector is
**parameter-free**. The closest entropic relative remains GREA (above), which differs in using the
boundary (area) entropy. A systematic priority search returned no work combining these three
elements; this companion concerns that horizon-entropy foundation.

## 2. The postulate, decomposed

"Volume-law" conflates four claims of very different status:

- **state** — the horizon degrees of freedom are maximally entangled (thermal), not in a
  low-entanglement ground state;
- **form** — the entropy scales as the volume, S ∝ V;
- **scale** — the magnitude is ρ_DE ∼ ρ_crit, not ρ_Planck;
- **count** — the *number* of horizon dof grows as the bulk (N ∝ R^{d−1}) rather than the
  boundary (N ∝ R^{d−2}).

Since S ∼ (state factor) × N and the Barrow deformation enters as S ∝ A^{1+Δ/2} ∝
R^{(d−2)(1+Δ/2)}, the exponent Δ is fixed by the **count** (area → Δ=0, volume → Δ=1 in d=4);
the state only sets whether the prefactor is maximal. So the postulate's empirical content —
Δ — is the counting claim, and the other three are separable. The rest of the paper derives or
reduces state, form and scale, and isolates count.

## 3. How derivable is the postulate? The route map

Five routes (`run_deriv_{A..E}_*.py`):

**A — de Sitter holography / maximal scrambler [DERIVED state; OPEN count].** A Majorana-SYK
diagonalisation confirms the horizon *state* is maximally entangled / volume-law: chaotic
level statistics (⟨r⟩ ≈ 0.71) and Page-value saturation (S/S_Page ≈ 0.99). But SYK is
all-to-all / geometry-free and cannot pose the bulk-vs-boundary count; the count is the genuine
open dS-holography sub-target.

**B — entanglement first law in a thermal state [REDUCED form].** Jacobson's derivation gets
the area-law Einstein equations from maximal entanglement of the *vacuum*. Redone around a
*thermal* (de Sitter Gibbs) reference, the entanglement entropy acquires a volume term that
dominates the area term once the region exceeds the thermal length λ_th = 1/T. Volume-law is
thus the generic thermalised behaviour, not exotic. At the de Sitter temperature alone the
horizon sits just on the area side (R_H/λ_th = 1/2π); volume-law dominance requires
thermalisation above the de Sitter bath (toward Planck = maximal scrambling), whereupon the
density is Planckian — the CKN scale. So the *form* reduces to thermalisation, leaving the
*scale* to CKN.

**C — Verlinde emergent gravity [REDUCED, unification].** A volume-law de Sitter entanglement
entropy is exactly what Verlinde's emergent-gravity programme invokes for apparent dark
matter; one volume-law scale gives both ρ_DE ∼ ρ_crit and a₀ = cH₀/2π (≈ 0.87× the
radial-acceleration-relation value). A unification motivation that ties the postulate to an
independent anomaly; it inherits that programme's open issues.

**D — gravitational non-additivity [NEGATIVE].** The Tsallis–Cirto map δ = 1 + Δ/2 makes
Δ=1 ⟺ δ=3/2 exact, but pinning δ from a measurable non-additivity fails: real cluster
kinematics are Gaussian (q ≈ 1.01, `run_cluster_tsallis.py`), not q ≈ 1.5 — the weakest route.

**E — no-go from energy bounds [DERIVED].** The Bekenstein bound on the horizon energy gives
S ≤ 2πRE/ħc = ¼(A/ℓ_P²) — *exactly* the area law (verified S_Bek/S_area = 1.000, the known
10¹²²). The volume law therefore cannot come from any maximum-entropy/energy argument; it
exceeds the bound by precisely R/ℓ_P ∼ 10⁶¹ (the "seam"). This proves the missing ingredient is
a **state/counting** input and fixes its size.

*Net:* state (A) and form (B) and scale (CKN) are first-principles or reduced; the irreducible
residue is the bulk-vs-boundary **count**.

## 4. Reducing the residue to a driven NESS

The count is not a free coin-flip — but here we must be careful to avoid a circularity, because the
entanglement class cannot simply be *posited*. It is fixed by the Hamiltonian and the state, and the
decisive fact is that **the entanglement area law is a theorem with a hypothesis**: it is proved for
*local* (finite-range) Hamiltonians — Hastings (2007) rigorously in 1D via Lieb–Robinson bounds,
Brandão–Horodecki (2015) from exponential decay of correlations, reviewed in Eisert–Cramer–Plenio
(2010) — and the proof *requires* exponential clustering. Gravity is 1/r (α = 1), strongly long-range
(α ≤ d) and non-additive [Campa, Dauxois & Ruffo 2009], and there the hypothesis fails: in an explicit
gapped free-fermion model we find the correlation function decays as a *power law* for α = 1 (slope
−1.5) versus *exponentially* for the short-range reference (slope −6.9), exactly as Koffel–Lewenstein–
Tagliacozzo (2012) found (`run_v1_locality.py`). So the area-law theorem **does not apply** to a
gravitationally-bound horizon: area-law is *not guaranteed*, and the count is genuinely open between
area and volume with neither the default. We are explicit about what this does and does not give. It
does *not* by itself derive volume — a noncritical long-range ground state can still be area-law
[Kuwahara–Saito 2020], and indeed the gapped ground-state block entropy stays bounded for every α we
test. The volume law additionally requires the horizon *state* to be thermal / maximally scrambled
(route A), where entanglement is volume-law; and even that fixes the *state*, not the *count*. The
super-extensive non-additivity (`run_residue_longrange.py`; u ∝ N^{1−α/d}, slopes 0.53 in d=2, 0.68 in
d=3 vs 0.50, 0.67) is the thermodynamic face of the same long-range property. The honest conclusion is
that volume is the *motivated, non-default* class for a self-gravitating thermal horizon — area-law
arises only when the holographic bound intervenes for an *isolated, equilibrium* horizon (a black hole
saturates the Bekenstein bound, route E, and relaxes to area-law) — while the bulk-vs-boundary **count**
itself remains the open dS-holography residue. One caveat keeps this honest: the area-law theorems are
statements about many-body systems with a specified Hilbert space and Hamiltonian, and what the horizon's
degrees of freedom *are* — and what Hamiltonian governs them — is precisely the open dS-holography question.
So this is a motivated *conditional* — *if* the horizon is a long-range many-body system, *then* area-law is
not forced and volume is favoured — not a theorem about the actual horizon. What it establishes is that the
entanglement class is read off the (modelled) Hamiltonian rather than assumed — which is what defeats the
circularity charge — not that the count is derived.

This same long-range property has a second, independent consequence that bears on the
canonical-vs-microcanonical λ ambiguity (whether ρ_DE ∝ H² f_sat, area, or ∝ H f_sat, volume). A
non-additive system has *inequivalent* statistical ensembles [Campa, Dauxois & Ruffo 2009], the
signature being negative specific heat — and the de Sitter horizon indeed has C = −2S < 0 throughout
its history (`run_lit_inferences.py`). For such a system the canonical ensemble is ill-defined; the
isolated self-gravitating horizon must be treated *microcanonically* (by its state count). The
canonical accounting is precisely the area-law (λ = 1) branch, so ensemble inequivalence removes it
as unphysical and selects the microcanonical (volume, λ = 1/2) branch the model uses. This does not
close the residue — whether the microcanonical state count is itself volume is the dS-holography
question — but it eliminates one of its two horns from first principles. The
cosmic horizon differs in one respect: it is continuously *driven* by structure formation (the
f_sat gate, dS_struct/dt > 0; present drive ≈ 0.45). This turns the black-hole/cosmic asymmetry
into a **discriminator** — equilibrium → area, driven NESS → volume — and reduces the residue
to a single, falsifiable statement:

> **Driven-NESS conjecture.** A strongly-long-range horizon, driven by structure formation,
> settles in the volume-law steady state rather than relaxing to the area-law equilibrium.

## 5. The driven-NESS conjecture and its two maps

**Prototype (`run_driven_ness.py`).** Generic driven-dissipative free fermions do *not* realise
the effect (uniform drive is trivially extensive and connectivity-independent; boundary drive
decays to area-law) — a NEGATIVE that sharpens the content to *cooperative bistable locking*.
With the area↔volume order parameter under cooperative (Wilson–Cowan) dynamics, a transient
structure-drive pulse locks the volume branch (m → 0.93, persistent) for long-range
cooperativity (J > J_c), while short-range tracks-then-relaxes to area and the undriven case
stays at the area ground. Only drive × long-range selects volume.

**Map 1 — connectivity → J ≥ J_c [DERIVED] (`run_gravitational_coupling.py`).** The cooperative
coupling is the largest eigenvalue of the gravitational coupling matrix, J_eff = λ_max(W_grav),
W_ij = g/r_ij^α. For a homogeneous horizon this equals the per-site gravitational binding
(λ_max/row-sum = 1.02) — the *same* super-extensive quantity behind the non-additivity of §4.
So volume-counting and volume-locking are one number. It scales as N^{1−α/d} (slopes 0.52, 0.67;
a short-range α>d control is bounded, N^{0.01}), so since α=1≤d it exceeds any fixed J_c=4T at
the horizon's N ∼ 10¹²² — bistability is generic, not tuned.

![Fig2](../output/criticality_soc.png)

**Figure 3.** Why the criticality is the volume law, not a separate assumption. *(a)* The area↔volume
free energy is a smooth crossover for the bare gate (J = 0) and develops a spinodal only for J ≥ J_c.
*(b)* The bare-gate susceptibility is finite; the spinodal diverges at the operating point m = ½. *(c)*
The cooperative coupling J = λ_max(W_grav) is the gravitational binding itself: all-to-all (volume-law)
connectivity gives χ_peak ∝ √N → ∞ (a real spinodal), while local (area-law) connectivity is bounded.

**Map 2 — deposition → drive [DERIVED] (`run_deposition_drive.py`).** The drive is the
accumulated structure-deposited horizon entropy (the collapsed fraction ∝ f_sat). Its timing is
derived — the drive clears the area-well barrier h_c when f_sat ≈ h_c, i.e. at z ≈ 1.45 ≈ z* — and
integrating the order-parameter dynamics with this physical drive flips the horizon area→volume
near z* and locks it (Δ = 1 permanent), whereas an un-amplified drive never leaves area.

![Fig3](../output/deposition_drive.png)

**Figure 4.** The deposition → drive map. The drive is the accumulated structure entropy (∝ f_sat);
it crosses the area-well barrier h_c at z ≈ z*, flipping the order parameter to volume and locking it
there permanently (red), whereas the un-amplified bare deposition stays on the area branch (blue).

**The drive is a violent-relaxation deposit, and the γ-mechanism is the same event**
(`run_violent_relaxation.py`). Structure forms by gravitational collapse, i.e. by *violent
relaxation* [Lynden-Bell 1967], which drives a long-range system on the dynamical time to a
quasi-stationary state rather than to thermal equilibrium. This single event does two things SEDE
otherwise treats separately. First, it dissipates the halo binding energy E_bind = Mσ_v² ∝ M^{5/3}
(virial: R_vir ∝ M^{1/3}, σ_v² ∝ M^{2/3}) into the horizon ledger — the mass weight that *is* the
structure coupling γ ≈ 1.5 (Paper I, Theorem 4C). Second, that same deposited binding energy
(∝ ⟨σ_v²⟩) is the amplitude of the NESS drive above. So the γ-weight and the driven-NESS drive are
one quantity from one process — the volume lock is the horizon's Lynden-Bell quasi-stationary state,
fed by violent relaxation. Two corollaries follow at no cost: collapse is prompt (t_dyn/t_Hubble =
Δ_vir^{−1/2} ≈ 0.07), so the horizon is driven by an ongoing sequence of such events whose envelope
is f_sat; and the quasi-stationary state's lifetime in a long-range system scales as N^δ
[Campa, Dauxois & Ruffo 2009], so at the horizon's N ∼ 10¹²² the lock is permanent — a second
permanence argument beside the second-law one (route E). Both the halo QSS and the de Sitter horizon
have negative specific heat, consistent with the microcanonical (volume) selection of §4.

## 6. Driven *through* the spinodal: ratchet and robustness

Developing the conjecture corrected its statement (`run_soc_attractor.py`): the system locks at
the *volume* well, not *at* the spinodal. Because the structure drive ∝ f_sat is monotone, it
*necessarily* climbs to the field at which the area well vanishes — the area-branch spinodal — at
z ≈ z*; the crossing is a forced **ratchet**, and the generalised second law then locks the
volume well (volume is the global entropy maximum, S_vol ≫ S_area by R/ℓ_P, route E). The
near-critical fluctuation signatures are therefore a *transient* at z*, not a permanent state.
The self-organisation content is a **robustness** claim, verified: volume-locking holds over a
broad contiguous basin of the landscape parameters (J ≥ J_c, θ), over drive amplitudes 0.6–3×,
and from all sub-barrier initial conditions — an attractor, not a knife-edge (it fails only at
strong holographic pull). What is *not* derived is the absolute horizon free energy F(m); only
that a broad, untuned range works. So the residue's sharpest form is "the horizon area↔volume
free energy is bistable with J ≥ J_c (gravity supplies this) and a sub-unity barrier" — strictly
weaker, more physical, and falsifiable than "assume volume counting".

**A new, sharp prediction: the crossing carries mean-field SOC statistics**
(`run_soc_signature.py`). Because the system is driven *through* its spinodal at z* (not parked at
a tuned critical point) and the coupling is all-to-all/mean-field (§5), the fluctuations at the
crossing are scale-free in the **mean-field self-organised-criticality** universality class
[Bak, Tang & Wiesenfeld 1987]. That class is parameter-free: the avalanche-size distribution is
P(s) ∝ s^{−3/2} (we recover τ = 1.51 from the critical branching process), and the susceptibility
and autocorrelation time diverge as the order parameter crosses the spinodal — a *transient*
critical slowing-down localised at z* (we find the peak at z ≈ 1.08 ≈ z*), reddening the temporal
spectrum toward 1/f. This turns the companion-paper CHR layer (Paper I §6, predictions P3–P5) from
"an enhancement at z*" into a specific, falsifiable statement: any z*-localised excess in growth/ISW
should be **scale-free with exponent τ ≈ 3/2 and a 1/f spectrum**, generic (self-organised, not
tuned), and *transient*. This discriminates SEDE from clustering dark energy / k-essence (which
carries a single characteristic scale, the dark-energy sound horizon, not a scale-free spectrum) and
from ΛCDM (no z*-localised feature), and is within reach of redshift-resolved weak-lensing
tomography and ISW cross-correlation (Euclid/Rubin). We flag the status honestly: this is a
*conditional* prediction. The exponent τ ≈ 3/2 is parameter-free, but the *amplitude* of the
z*-localised excess is not yet forecast — it is set by the deposited-entropy fraction ε_dep ∼ 10⁻⁷ and
the near-critical susceptibility, and a quantitative amplitude and survey signal-to-noise are deferred
to the companion CHR analysis. So τ ≈ 3/2 is a falsifiable shape *conditional on the CHR layer being
detectable at all*, not a claim that it is.

![Fig4](../output/foundations_fig4_socsignature.png)

**Figure 5.** The new prediction. *(a)* The near-spinodal fluctuations fall in the mean-field
self-organised-criticality class, whose avalanche-size law is the parameter-free P(s) ∝ s^{−3/2}
(recovered as τ = 1.51 from critical branching). *(b)* As the structure drive carries the order
parameter to the spinodal, the area-phase susceptibility and autocorrelation time diverge — a
*transient* critical slowing-down localised at z*, reddening the spectrum toward 1/f. Any z*-localised
growth/ISW excess should therefore be scale-free with τ ≈ 3/2, which single-scale clustering dark
energy cannot mimic.

## 7. What would make it a closure

A closure derives F(m) itself, reducing the dark sector to Ω_m + established physics. The three
routes (`run_closure_attempts.py`) converge:

- **A (DSSYK).** Entropy extensive in the dof count (S_max = (N/2)ln 2 ∝ N, slope 1.00); cannot
  be super-extensive (volume ∝ N^{3/2}). We state this plainly: **read in the standard
  Narovlansky–Verlinde dictionary, the leading concrete model of de Sitter holography maps to the
  de Sitter *area* (Gibbons–Hawking) — i.e. it currently favours the area branch, a real tension
  SEDE's volume postulate must eventually overcome.** The mitigating point — that DSSYK is
  all-to-all / geometry-free and so cannot, by itself, *pose* the bulk-vs-boundary counting question
  — is secondary, not a neutraliser. As things stand, route A is evidence to be answered, not used.
- **B (Euclidean dS saddles).** Einstein and Barrow horizon free energies each have a single
  saddle — no volume-law saddle of the bare path integral (route-E no-go in saddle form).
- **C (thermal F(m)).** The bistable landscape *is* derivable —
  F(m)/T = −(J/2)m² + θm + [m ln m + (1−m)ln(1−m)] from gravitational cooperativity
  (J = λ_max ≥ J_c), the holographic pull (route E), and configurational entropy: area well
  m ≈ 0.07, volume well m ≈ 0.93 whenever J ≥ J_c. Level-2 existence, conditional on the volume
  branch being accessible.

All three reduce to one statement — *the horizon's dof count is the bulk (volume) one* — which A
and B cannot supply and C uses. Closure has three precise forms: **Level 1**, a dS-holography
computation of the bulk count (open, not shortcut-able past route E); **Level 2**, exhibiting the
volume branch as a real saddle/state (C gives the landscape conditional on it); **Level 3**,
*measuring* Δ (DESI DR3 + Euclid to ∼0.09; Δ=1 establishes, Δ=0 refutes — decisive, in hand).

### 7.1 The count is state-dependent — and that resolves the "everything points to area" tension

Routes A and B, and the standard handles generally (the de Sitter modular Hamiltonian is a *local*
boost; DSSYK in the Narovlansky–Verlinde dictionary; the CKN holographic bound; black-hole
thermodynamics), all favour the **area** count. Read as universal statements about "the" horizon
count they look like a standing objection. The resolution is that they are not universal: **the count
is state-dependent**, and every one of those handles is an *equilibrium* horizon. The
equilibrium-vs-driven discriminator of §4 is exactly this — the f_sat gate *is* the fraction of the
volume capacity that structure-driving activates, N_eff = N_area + (N_vol − N_area)·f_sat, so the
equilibrium limit (f_sat → 0) is area and the driven, locked limit (f_sat → 1) is volume
(`run_residual_resolution.py`). (This is *not* the rejected running-Δ background of §3: the cosmic
horizon is driven past the spinodal early and *locks* at constant Δ = 1; the state-dependence is
*across horizons*, driven cosmic vs equilibrium black hole, not a running Δ over cosmic time.) So the
area results are SEDE's equilibrium baseline, not counter-evidence. One case needs care because it is
*not* obviously equilibrium: the DSSYK correspondence is a duality with the de Sitter static patch —
the cosmological horizon itself. The point is that it is a duality with *eternal* (matter-free)
de Sitter, an equilibrium spacetime with no structure to drive it (f_sat ≡ 0); the real universe's
apparent horizon is embedded in a matter-plus-structure cosmology and is driven. So "DSSYK = area" is
the area of *equilibrium* de Sitter — again the f_sat → 0 limit — and is consistent with a volume-law
*driven* horizon, not a contradiction of it. We are explicit that the state-dependent count is itself
an effective-model statement: writing N_eff as a function of f_sat, and reading a black hole as the
f_sat = 0 case, extends f_sat from its cosmological definition (a function of the growth D(z)) to a
putative property of any horizon's driving — physically motivated, but the same conditional as V1, one
layer out.

This makes the residue **doubly falsifiable, on two independent horizons**. (i) The cosmic horizon
(driven) should have Δ = 1 — DESI DR3 + Euclid. (ii) The black hole (equilibrium) should have Δ = 0 —
and the observable is primordial-black-hole evaporation: SEDE *predicts that PBHs evaporate normally*
(standard Hawking), because black holes are undriven, whereas a *universal*-volume theory forbids it
(T_B/T_H ∼ 7×10⁻²¹, evaporation time × 10⁸¹). The PBH channel therefore tests the *mechanism* — it
distinguishes driven-NESS SEDE (Δ_BH = 0) from universal-volume quantum gravity (Δ_BH = 1) — and the
black hole being area-law, far from a problem, is a SEDE prediction. This is a *future* discriminator,
not present-day evidence: SEDE's Δ_BH = 0 means standard Hawking evaporation, consistent with all
current data (including ~10¹⁵ g PBH abundance bounds) precisely because it predicts no new black-hole
effect; the channel discriminates only on a positive PBH-evaporation detection or a robust population
statement.

This also yields a new statement about black holes in general: **area-law is the *equilibrium* default,
not a fundamental geometric law** — a black hole is area-law *because it is undriven* (f_sat → 0), and the
cosmic horizon's uniqueness is a calculation, not an assumption. The most-driven black hole conceivable is
a primordial one *at formation*, where the entire horizon is built by a collapsing overdensity (violent
relaxation); yet even there the structure entropy deposited is negligible against the black-hole entropy
*created*, S_rad/S_BH ∼ (M/M_P)^{−1/2}, so the formation drive f_form ∼ (M_P/M)^{1/2} ≪ h_c even after the
near-critical amplification (`run_pbh_driving.py`, `run_bh_new_physics.py`). Accreting and merging black
holes are weaker still. So SEDE robustly predicts *all* black holes are area-law (Δ_BH = 0), and the
cosmic apparent horizon — uniquely driven *sustainedly* to f_sat → 1 over the whole structure-formation
history — is the only volume-law horizon. The falsifier is correspondingly sharp: an *isolated* volume-law
black hole (anomalously slow PBH evaporation, or a modified-area GW ringdown) would mean Δ is *universal*
(Barrow-type) and would refute the state-dependent, driven-vs-equilibrium mechanism — so SEDE is *more*
falsifiable on black holes than a universal-Barrow theory, not less.

The honest theoretical residue is
correspondingly narrowed: not "derive the count," but "prove that a *driven* horizon activates the
volume capacity" — a statement now *consistent with* the equilibrium area results (their f_sat → 0
limit), with the two-channel measurement deciding it either way regardless of dS holography.

Two final caveats, stated openly. First, the *magnitude*: the structure drive is ε_dep ∼ 10⁻⁷, yet the
count it flips spans area to volume — a factor ∼ R/ℓ_P ∼ 10⁶¹ in activated degrees of freedom. The
near-critical bistability is therefore doing heavy lifting: the drive does not *add* the dof, it
*nucleates* the area→volume transition past the spinodal, and the second law carries the system to the
volume well, where it locks (§5–6). That a 10⁻⁷ trigger produces a 10⁶¹ change in the count is the
V1/§7.1 conditional stated quantitatively, and it rests entirely on the (modelled) bistable landscape.
Second, the *pre-registered* case: a measured cosmic Δ significantly below 1 — say Δ ≈ 0.5 from DESI
DR3 + Euclid — would be evidence *against* the locked-Δ = 1 picture, pointing either to incomplete
volume activation or to a genuinely intermediate count, and would *not* be a free parameter to absorb;
the model predicts Δ = 1 for the driven cosmic horizon and is falsified by a robust intermediate value.

### 7.2 The count as a horizon Hausdorff dimension, and what remains of it (`run_ds_count.py`)

The Level-1 question — "is dim 𝓗(de Sitter static patch) = e^{Area/4} or e^{Volume}?" — has a clean
geometric form: dim 𝓗 = e^{Area} ⟺ the horizon is a smooth 2-surface (Hausdorff dimension d_H = 2), and
dim 𝓗 = e^{Volume} ⟺ it is a space-filling fractal (d_H = 3), so the count *is* the horizon's Hausdorff
dimension, Δ = d_H − 2. The standard frameworks all give a smooth horizon: the Gibbons–Hawking maximum,
the finite-dimensional proposal [Banks 2001], the Type II₁ observer algebra in which S_dS is the maximal
entropy [Chandrasekaran, Longo, Penington & Witten 2023], and DSSYK. There is, however, a concrete
mechanism by which the *driven* horizon could differ: a horizon driven by a stochastic source is a
roughening 2-surface, and a self-affine surface of roughness exponent α has d_H = 3 − α, hence Δ = 1 − α.
The identification we make is that the roughening *noise is the structure drive* — so an undriven horizon
is the smooth minimum (d_H = 2, area) and a structure-driven horizon roughens (d_H → 3, volume). This is
the microscopic, horizon-geometry form of the state-dependent count of §7.1, and it *reconciles* the area
frameworks (which describe the eternal, equilibrium horizon) with SEDE's volume (the driven horizon) as
two limits of one statement, rather than as a tension. It reduces the residue to two sharp, well-posed
questions. (i) Is the driven roughening *linear* (Edwards–Wilkinson, the 2+1D lower-critical class,
α → 0 ⇒ Δ = 1) or *nonlinear* (KPZ, α ≈ 0.39 ⇒ Δ ≈ 0.6)? This is settled by two textbook inputs rather
than a conjecture (`run_kpz_membrane.py`): the KPZ coefficient *is* the normal growth velocity
(λ = v, from the height-representation geometry ∂h/∂t = v√(1+(∇h)²) ≈ v + (v/2)(∇h)²), and a driven
interface's velocity is proportional to its driving free energy (v = M·F/V, interface kinetics — for a
horizon, the membrane-paradigm statement that it grows in response to the free-energy flux across it
[Damour 1978; Thorne, Price & Macdonald 1986]). So λ = M·F/V ∝ F; and SEDE's free energy vanishes
identically (F = E − T_AH S = 0, since ρ_DE = T_AH s_grav; §8), giving λ = 0 ⇒ EW ⇒ Δ = 1. The one
physical assumption is that the horizon grows by membrane-paradigm interface kinetics — the standard
effective description of horizon dynamics. (ii) Does the marginal 2+1D roughening count as a genuine
d_H = 3 fractal? Here we must be candid (`run_ds_marginal.py`): the horizon is a 2-surface, which is
*exactly* the lower critical dimension of EW roughening (α = (2−d)/2 = 0 at d = 2) — below it the same
dynamics is a genuine fractal, above it is smooth, and at d = 2 it is *marginal*, space-filling only
logarithmically. So the roughening picture does **not robustly force** Δ = 1; it explains why the count is
delicate — the horizon sits precisely on the area↔volume boundary — and leaves the tie to sub-leading
physics. This tempers "EW ⇒ Δ = 1" to "EW ⇒ Δ = 1 *marginally*," and it means a measured Δ ∈ (0, 1) (the
log-suppressed case) would be *within* the picture, not a refutation.

The honest bottom line, then, is that we have **not resolved de Sitter holography** — no one has computed
dim 𝓗(static patch) non-perturbatively, and we do not claim to (`run_ds_resolve.py`). What we have done
is turn the count from an unstructured "area or volume?" into a *geometric, empirically-decidable*
statement (Δ = d_H − 2), reduce its theoretical residue to two named obstructions — the marginal
dimension (O1) and the Planck-scale validity of the effective descriptions (O2) — and identify the single
open computation an actual resolution requires. The equilibrium value is Δ = 0 and the driven value is
Δ → 1, marginally; the empirical Δ (DESI DR3 + Euclid, §6) is the arbiter.

## 8. Discussion

The chain reduces SEDE's one foundational input from a static, untestable counting postulate to
a single dynamical, falsifiable statement: *the de Sitter static patch's degrees of freedom count
its bulk, not its boundary.* Everything else — state, form, scale, both microscopic maps, and the
bistable landscape conditional on the volume branch — is derived or reduced. Two honest limits
bound the claim. First, this is a *reduction*, not a closure: the volume branch is used, not
derived, and deriving it is the open dS-holography problem (route E forbids energy-bound
shortcuts). Second, the dynamical model is effective/mean-field; the connectivity→J and
deposition→drive maps are physical, but F(m)'s absolute parameters are not computed from the
microscopic horizon, only shown to lie in a broad untuned basin. What makes the reduction
worthwhile is that it is *falsifiable where the postulate was not*: it predicts a transient
near-critical epoch at z* and ties Δ to the structure-driving history, and Δ itself is pinned to
∼0.09 by DESI DR3 + Euclid. The residue is then not a weakness in the cosmology but a clean
hand-off to a well-posed quantum-gravity question, with a measurement that settles it either way.

### Ledger

| step | result | tag | script |
|---|---|---|---|
| state | maximal entanglement (Page saturation) | DERIVED | `run_deriv_A_dsholo.py` |
| form | volume = thermal-state entanglement | REDUCED | `run_deriv_B_thermal.py` |
| scale | ρ_DE ∼ ρ_crit (CKN) | DERIVED | — |
| no-go | energy bounds ⇒ area, not volume | DERIVED | `run_deriv_E_nogo.py` |
| residue→range | gravity long-range (α≤d) ⇒ volume class | REDUCED | `run_residue_longrange.py` |
| conjecture | drive × long-range locks volume | MODELLED | `run_driven_ness.py` |
| map 1 | J_eff = λ_max(W_grav) = binding ≥ J_c | DERIVED | `run_gravitational_coupling.py` |
| map 2 | drive ∝ f_sat clears barrier at z* | DERIVED | `run_deposition_drive.py` |
| spinodal | ratchet through spinodal; robust basin | DERIVED (robust) | `run_soc_attractor.py` |
| closure | three routes ⇒ the dS state count | analysis | `run_closure_attempts.py` |
| ensemble | C = −2S + long-range ⇒ microcanonical; area horn removed | REDUCED | `run_lit_inferences.py` |
| unify | γ-mechanism = driven-NESS (one Lynden-Bell QSS) | DERIVED | `run_violent_relaxation.py` |
| prediction | z*-fluctuations scale-free, τ ≈ 3/2 (mean-field SOC) | PREDICTION | `run_soc_signature.py` |
| V1 defence | area-law is a *local*-Hamiltonian theorem; gravity (power-law correlations) voids it | ARGUED | `run_v1_locality.py` |
| residual | count is state-dependent (eq.→area, driven→volume); PBH evaporation tests it | REFRAME / PREDICTION | `run_residual_resolution.py` |
| dS count | count = horizon Hausdorff dim; driving roughens d_H 2→3; residue → 2 sub-questions | ANALYSIS | `run_ds_count.py` |
| **residue** | **the m≈1 volume branch is the horizon's state** | **OPEN (L1/2) / measurable (L3)** | dS holography / DR3+Euclid |

### Reproducibility

```
python reproduce_all.py --only derivA,derivB,derivC,derivD,derivE,residue,ness,jcoup,\
deposdrv,socatt,closure,litinf,violrel,socsig,v1loc,resid,bhnew,pbhdrv,dscount,kpzmem,dsmarg,dsresolve
```

Figures 1–5 are generated by `paper/make_foundations_figures.py`: Figs 1 (reduction chain), 2
(inferences) and 5 (SOC signature) are built directly there, while Fig 3 = `criticality_soc.png` and
Fig 4 = `deposition_drive.png` are produced by the `soc`/`deposdrv` run-scripts. Each `run_*.py` prints
a verdict and runs validation assertions; all stages PASS.

## Code and data availability

All code and analysis pipelines are publicly available at `github.com/spsingularity/SEDE-release`.
Every quantitative claim in this paper is produced by a named `run_*.py` script wired into the single
entry point `reproduce_all.py` (stages `soc`, `deriv{A–E}`, `residue`, `ness`, `jcoup`, `deposdrv`,
`socatt`, `closure`, `dscount`, `kpzmem`, `dsmarg`, `dsresolve`, …), each carrying validation
assertions; the figures are regenerated by `make_foundations_figures.py`. The model and its cosmology
live in the `sede/` package with the companion (Paper I); no additional observational data beyond
Paper I's standard public inputs is used here. A tagged release is archived at Zenodo, DOI
[10.5281/zenodo.21050314](https://doi.org/10.5281/zenodo.21050314).

## Use of AI tools

Artificial-intelligence tools (large language models) were used to assist with drafting and editing the
manuscript, with developing and cross-checking the analysis code, and with literature searches. The
author conceived the model, directed and independently verified all analyses, and takes full
responsibility for the content, including all derivations, results, and conclusions. No AI tool is an
author.

## References

BibTeX in `paper/refs.bib` (cosmology, shared with Paper I) and `paper/refs_foundations.bib`
(foundational additions for this paper). Key entries, by theme:

*Horizon thermodynamics & entanglement* — Jacobson, *Thermodynamics of Spacetime*, PRL 75 (1995)
1260 [`Jacobson:1995ab`]; Cai & Kim, JHEP 02 (2005) 050 [`Cai:2005ra`]; Jacobson, *Entanglement
Equilibrium and the Einstein Equation*, PRL 116 (2016) 201101, arXiv:1505.04753
[`Jacobson:2015hqa`]; Ryu & Takayanagi, PRL 96 (2006) 181602, hep-th/0603001 [`Ryu:2006bv`];
Gibbons & Hawking, PRD 15 (1977) 2738 [`GibbonsHawking1977`].

*Generalised entropy & holographic dark energy* — Cohen, Kaplan & Nelson, PRL 82 (1999) 4971
[`Cohen:1998zx`]; Barrow, PLB 808 (2020) 135643 [`Barrow:2020tzx`]; Saridakis, PRD 102 (2020)
123525 [`Saridakis:2020zol`]; Tsallis & Cirto, EPJC 73 (2013) 2487 [`Tsallis:2012js`];
Padmanabhan, arXiv:1206.4916 [`Padmanabhan:2012ik`]; Bousso, RMP 74 (2002) 825 [`Bousso:2002ju`].

*Entropy-production acceleration & non-equilibrium gravity* — García-Bellido &
Espinosa-Portalés, *Cosmic acceleration from first principles*, PDU 34 (2021) 100892,
arXiv:2106.16014 [`Garcia-Bellido:2021idr`]; Espinosa-Portalés & García-Bellido, PDU 34 (2021)
100893, arXiv:2106.16012 [`Espinosa-Portales:2021cac`].

*Quantum chaos, ETH, de Sitter holography* — Maldacena, Shenker & Stanford, JHEP 08 (2016) 106,
arXiv:1503.01409 [`Maldacena:2015waa`]; Cotler et al., JHEP 05 (2017) 118 [`Cotler:2016fpe`];
Page, PRL 71 (1993) 1291 [`Page:1993df`]; Deutsch, PRA 43 (1991) 2046 [`Deutsch:1991msp`];
Srednicki, PRE 50 (1994) 888, cond-mat/9403051 [`Srednicki:1994mfb`]; Rigol, Dunjko & Olshanii,
Nature 452 (2008) 854 [`Rigol:2007juv`]; Banks, IJMPA 16 (2001) 910, hep-th/0007146
[`Banks:2000fe`]; Susskind, JHAP 1 (2021) 1, arXiv:2109.14104 [`Susskind:2021esx`]; Lin, JHEP 11
(2022) 060, arXiv:2208.07032 [`Lin:2022rbf`]; Narovlansky & Verlinde, JHEP 05 (2025) 032,
arXiv:2310.16994 [`Narovlansky:2023lfz`].

*Emergent gravity* — Verlinde, JHEP 04 (2011) 029 [`Verlinde:2010hp`]; Verlinde, SciPost Phys. 2
(2017) 016 [`Verlinde:2016toy`].

*Long-range systems, criticality, monitored dynamics* — Campa, Dauxois & Ruffo, Phys. Rept. 480
(2009) 57, arXiv:0907.0323 [`Campa:2009jxa`]; Bak, Tang & Wiesenfeld, PRL 59 (1987) 381
[`Bak:1987xua`]; Hawking & Page, CMP 87 (1983) 577 [`Hawking:1982dh`]; Fisher et al., Ann.
Rev. Cond. Mat. Phys. 14 (2023) 335, arXiv:2207.14280 [`Fisher:2022qey`].

*Volume-law entanglement in cosmology* — Khlebnikov & Sheoran, arXiv:1907.00487
[`Khlebnikov:2019wbk`].

*Structure-sourced dark energy* — Gough, *Information dark energy* (≈2008) [`Gough2008`];
Pandev, *Structural Entropy Dark Energy: a fixed-parameter, growth-gated holographic
dark-energy model without a cosmological constant* (Paper I, companion to this work; in
preparation, 2025) — the cosmology this paper's foundations underwrite.
