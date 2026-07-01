# The Driven-NESS Programme — reducing SEDE's one postulate, end to end

**Scope.** SEDE's dark sector has no fitted parameter: w₀, λ = ½, γ ≈ 1.5 and
c_s² = 1 are all fixed once Ω_m and a single foundational input are given. That
input is the **volume-law horizon postulate** (§8.4): the cosmic apparent horizon
carries the coarse-grained entanglement entropy of its *volume*, S_grav ∝ V ∝ R³
(the Barrow endpoint Δ = 1), rather than the Bekenstein–Hawking *area* law. This
note records a single connected line of work that takes that postulate apart:
what is derivable, what reduces, what is irreducible, and exactly what would close
it. Every step is a runnable experiment (`run_*.py`, staged in `reproduce_all.py`)
with validation assertions; numbers quoted here are reproduced by those scripts.

Epistemic tags: **[DERIVED]** follows from stated inputs by standard physics;
**[REDUCED]** an assumption replaced by a weaker/more physical one; **[NEGATIVE]**
a tested idea that fails (reported because the failure is informative);
**[OPEN]** an isolated assumption or open problem.

---

## 1. The postulate, decomposed

"Volume-law" conflates four sub-claims of very different status. Separating them
is what makes the rest possible.

| sub-claim | meaning | status | where |
|---|---|---|---|
| **state** | horizon dof maximally entangled (thermal), not a low-entanglement ground state | **[DERIVED]** | route A |
| **form** | the entropy scales as the volume, S ∝ V | **[REDUCED]** to thermalisation | route B |
| **scale** | the magnitude is ρ_DE ∼ ρ_crit, not ρ_Planck | **[DERIVED]** (CKN) | §8.2 |
| **count** | the *number* of horizon dof grows as the bulk (N ∝ R³), not the boundary (N ∝ R²) | **[OPEN]** → the residue | §8.5, this note |

The deformation Δ is fixed by the **count** (area → Δ=0, volume → Δ=1); the state
only sets whether the prefactor is maximal. So the residue, after everything
below, is purely the counting claim — a de Sitter-holography question,
"dim 𝓗(static patch) = e^{Area/4} or e^{Volume}?"

---

## 2. How derivable is the postulate? The route map (A–E)

Five independent routes were tried (`run_deriv_{A..E}_*.py`). None closes the
postulate; together they locate it.

- **A — dS holography / SYK** [DERIVED state, OPEN count]. A Majorana-SYK
  diagonalisation confirms the horizon *state* from first principles: chaotic
  level statistics (⟨r⟩ ≈ 0.71) and Page-value saturation (S/S_Page ≈ 0.99). But
  SYK is all-to-all / geometry-free and **cannot pose** the bulk-vs-boundary
  count. State settled; counting is the open dS-holography sub-target.
- **B — entanglement first law around a thermal state** [REDUCED form]. Redoing
  Jacobson's construction around a *thermal* (de Sitter Gibbs) reference instead
  of the vacuum gives an entanglement entropy whose volume term dominates for
  R > λ_th = 1/T. So the volume *form* is the generic thermalised behaviour, not
  exotic. At T_dS alone the horizon sits just on the area side (R_H/λ_th = 1/2π);
  volume requires thermalisation above the de Sitter bath (→ Planck = maximal
  scrambling), whereupon the density is Planckian = the CKN scale. The form is
  derived from thermalisation; the scale is CKN. **Strongest reduction.**
- **C — Verlinde emergent gravity** [REDUCED, unification]. One volume-law de
  Sitter entropy underlies both ρ_DE ∼ ρ_crit and the apparent-dark-matter scale
  a₀ = cH₀/2π (≈ 0.87× the RAR value). A motivation that ties the postulate to an
  independent anomaly; inherits Verlinde's open issues.
- **D — gravitational non-additivity** [NEGATIVE]. The Tsallis–Cirto map
  δ = 1 + Δ/2 makes Δ=1 ⟺ δ=3/2 exact, but pinning δ from a measurable
  non-additivity fails: real cluster kinematics are Gaussian (q ≈ 1.01,
  `run_cluster_tsallis.py`), not the q ≈ 1.5 a kinetic mirror would need.
- **E — no-go from energy bounds** [DERIVED negative]. The Bekenstein bound on the
  horizon energy gives S ≤ ¼(A/ℓ_P²) — *exactly* the area law (verified
  S_Bek/S_area = 1.000, the known 10¹²²). The volume law cannot come from any
  energy/maximum-entropy argument; it exceeds the bound by precisely R/ℓ_P ∼ 10⁶¹
  (the §8.3 "seam"). This proves the missing ingredient is a **state/counting**
  input, and fixes its size.

**Net of the map:** state (A) and form (B) and scale (CKN) are first-principles
or reduced; the irreducible residue is purely the **count**.

---

## 3. Reducing the residue: count ⟺ connectivity ⟸ interaction range

The counting is not a free coin-flip (`run_residue_longrange.py`).

1. **count ⟺ connectivity** (Ryu–Takayanagi min-cut): locally-connected horizon
   dof give area-law entropy; nonlocal/expander connectivity gives volume-law.
2. **connectivity ⟸ range**: Newtonian gravity is 1/r^α with α = 1, and the
   Campa–Dauxois–Ruffo classification is that α ≤ d (embedding dimension) is
   *strongly long-range* and **non-additive** — entropy not area-additive. On the
   horizon surface (d=2) or bulk (d=3), α = 1 ≤ d. Verified: the per-particle
   energy is super-extensive, u ∝ N^{1−α/d} (measured slopes 0.53 in d=2, 0.68 in
   d=3 vs predicted 0.50, 0.67).
3. **black hole vs cosmic horizon**: a black hole is gravitational too, yet
   area-law. Resolution — *equilibrium vs driven*. The holographic bound (route E)
   relaxes an isolated, equilibrium horizon to area-law; the cosmic horizon is
   continuously **driven** by structure formation (the f_sat gate, dS_struct/dt > 0,
   present drive ≈ 0.45), a non-equilibrium steady state held off the area-law
   equilibrium. This turns the paper's named conceptual cost into a *discriminator*.

So the residue **[REDUCED]** to one provable statement — a **driven-NESS entropy
theorem**: *a strongly-long-range horizon, driven by structure formation, settles
in the volume-law steady state rather than relaxing to the area-law equilibrium.*

---

## 4. The driven-NESS theorem, prototyped

`run_driven_ness.py` tests the theorem in a tractable model, and the negative
result is as informative as the positive one.

- **[NEGATIVE]** Generic driven-dissipative free fermions do *not* realise it:
  uniform drive gives a trivially extensive thermal state (connectivity-
  independent); boundary drive decays to area-law. So it is **not** generic driven
  transport.
- **[DERIVED, model]** The content is *cooperative bistable locking*. With the
  area↔volume order parameter under Wilson–Cowan dynamics (a fixed threshold θ =
  the holographic pull; cooperativity J = connectivity), a transient
  structure-drive pulse: **driven + long-range (J=6>J_c)** locks the volume branch
  (m → 0.93, persists after the drive); **short-range (J=2<J_c)** tracks then
  relaxes to area; **undriven** stays at the area ground. Only drive × long-range-
  cooperativity selects volume.

The theorem needs two microscopic maps. Both are now supplied.

### 4.1 Map 1 — connectivity → J ≥ J_c  [DERIVED]

`run_gravitational_coupling.py`. The cooperative coupling is the largest (Perron)
eigenvalue of the gravitational coupling matrix, J_eff = λ_max(W_grav),
W_ij = g/r_ij^α. For a homogeneous horizon this **equals the per-site
gravitational binding** (verified λ_max/row-sum = 1.02) — the *same*
super-extensive quantity behind the non-additivity (§3). So volume-*counting* and
volume-*locking* are one number, not two assumptions. It scales as N^{1−α/d}
(slopes 0.52 in d=2, 0.67 in d=3); a short-range α > d control is bounded
(N^{0.01}). Since α = 1 ≤ d, J_eff exceeds any fixed J_c = 4T at the horizon's
N ∼ 10¹²² — bistability is generic, not tuned. (The black hole shares this floor:
also bistable, but undriven ⇒ area. Bistability is necessary, not sufficient.)

### 4.2 Map 2 — deposition → drive  [DERIVED]

`run_deposition_drive.py`. The drive is the accumulated structure-deposited
horizon entropy = the collapsed fraction ∝ f_sat (the cumulative field, to which
the hysteretic lock responds — not the rate). **Timing:** the drive clears the
area-well barrier h_c ≈ 0.42 when f_sat ≈ h_c, i.e. at z ≈ 1.45 ≈ z* = 1.15, as
structure matures. **Magnitude:** the bare deposition ε_dep ≈ 1.5×10⁻⁷ is far
below h_c, but the cumulative activated fraction is O(1); integrating with this
physical drive flips the horizon area → volume near z* and **locks** it (Δ = 1
permanent), whereas an un-amplified drive never leaves the area branch.

---

## 5. Driven *through* the spinodal — and a correction

Developing "the horizon is driven to the spinodal" (`run_soc_attractor.py`)
corrected the picture: the system locks at the *volume* well, not *at* the
spinodal. The right statement:

- **[DERIVED] ratchet.** Because the structure drive ∝ f_sat is monotone, it
  *necessarily* climbs to the field at which the area well vanishes — the
  area-branch spinodal — at z ≈ z*. The crossing is forced, not tuned. The system
  is carried *through* the spinodal and then locked at the volume well.
- **[DERIVED] irreversibility.** Volume is the global entropy maximum
  (S_vol ≫ S_area by R/ℓ_P, route E), so the generalised second law forbids
  relaxation back; Δ = 1 is permanent. The near-critical P3–P5 signatures are a
  **transient at z*** — exactly as §6 frames them, not a permanent state. (This
  also removes any need for a 10⁷ "near-critical amplification": the cumulative
  drive is simply O(1).)
- **[DERIVED, robust] self-organisation = robustness.** Volume-locking holds over
  a broad contiguous basin of the landscape parameters (J ≥ J_c, θ), over drive
  amplitudes 0.6–3×, and from *all* sub-barrier initial conditions — an attractor,
  not a tuned knife-edge. It fails only at strong holographic pull (large θ):
  broad basin, not universal.

**[OPEN]** What is not derived: the absolute horizon free energy F(m) (the values
of J and θ) — only that a broad, untuned range works. So the residue's sharpest
form is *not* "assume volume counting" but "assume the horizon's area↔volume free
energy is bistable, with J ≥ J_c (gravity supplies this) and a sub-unity barrier"
— a strictly weaker, more physical, falsifiable assumption.

---

## 6. What would make it a closure

A closure derives F(m) itself, reducing the dark sector to Ω_m + established
physics with no SEDE-specific input. The three named routes were tried
(`run_closure_attempts.py`) and **converge**.

- **A (DSSYK).** Entropy extensive in the dof count, S_max = (N/2)ln 2 ∝ N
  (slope 1.00); cannot be super-extensive (volume ∝ N^{3/2}). Geometry-free ⇒
  silent on area-vs-volume; the standard dictionary picks area. Settles the state,
  not the count.
- **B (Euclidean dS saddles).** Einstein and Barrow horizon free energies each
  have a *single* saddle — no volume-law saddle of the bare path integral (route-E
  no-go in saddle form).
- **C (thermal F(m)).** The bistable landscape **is derivable** —
  F(m)/T = −(J/2)m² + θm + [m ln m + (1−m)ln(1−m)] from gravitational
  cooperativity (J = λ_max ≥ J_c), the holographic pull (route E), and the
  configurational entropy of the area/volume mixture: area well m ≈ 0.07, volume
  well m ≈ 0.93 whenever J ≥ J_c. **Level-2 existence**, conditional on the volume
  branch being an accessible state.

All three reduce to one statement — *the horizon's dof count is the bulk (volume)
one; the m ≈ 1 branch exists.* A and B cannot supply it; C uses it. So the residue
is **provably the de Sitter state count**, not removable by these shortcuts.
Closure has three precise forms:

- **Level 1 (full):** a dS-holography computation deriving the bulk/volume count.
  Open; route E shows it cannot be shortcut through energy bounds.
- **Level 2 (existence):** exhibit the volume branch as a real saddle/state.
  C gives the landscape conditional on it; genuine Level-2 closure constructs the
  branch (a volume-law saddle of a controlled dS partition function, or a DSSYK
  volume fixed point in a thermalised regime).
- **Level 3 (empirical):** measure Δ. DESI DR3 + Euclid pin Δ to ∼0.09 (§6);
  Δ = 1 establishes the volume branch as fact, Δ = 0 refutes it — settling the
  residue without a derivation. **Decisive, and already in hand.**

### 6.1 The dS-count sharpened: horizon roughening (and where it stops)

Working the Level-1 count further (`dscount`, `kpzmem`, `dsmarg`, `dsresolve`)
sharpens it into a *geometric* statement without resolving it:

- **The count is the horizon's Hausdorff dimension**, Δ = d_H − 2: d_H = 2 (smooth
  2-surface) ⟺ area, d_H = 3 (space-filling) ⟺ volume. Every *equilibrium*
  framework — Gibbons–Hawking, Banks–Fischler, CLPW Type II₁, DSSYK — gives d_H = 2
  (area); these describe the *eternal* horizon and are correct there (`dscount`).
- **A driven horizon roughens.** Its roughening class fixes d_H. The class is
  Edwards–Wilkinson (not KPZ) because the KPZ lateral coefficient is the normal
  growth velocity λ = v = M·(F/V) (geometry + membrane-paradigm interface kinetics,
  Damour/Thorne–Price–Macdonald), and SEDE's free energy F = E − T_AH·S = 0
  identically (ρ_DE = T_AH·s_grav). F = 0 ⇒ λ = 0 ⇒ EW (`kpzmem`).
- **But the horizon is a 2-surface = the EW lower critical dimension** (α = (2−d)/2 = 0
  at d = 2): the driven roughening is *marginal* (logarithmic), so it does **not
  robustly force** d_H = 3. It places the count on a knife-edge (`dsmarg`).
- **Honest status: NOT resolved** (`dsresolve`). Two irreducible obstructions remain —
  (O1) the marginal dimension (sub-leading log physics not fixed by the effective
  description) and (O2) the Planck-scale validity of the membrane/EW descriptions —
  both reducing to one open computation: dim 𝓗(dS static patch), non-perturbatively.
  A measured Δ ∈ (0,1] is *within* the picture (marginal/log-suppressed vs genuine
  space-filling); Level 3 breaks the tie. Claiming a resolution would be over-reach.

---

## 7. Ledger

| step | result | tag | script / stage |
|---|---|---|---|
| state | maximal entanglement (Page saturation) | DERIVED | `derivA` |
| form | volume = thermalised-state entanglement | REDUCED | `derivB` |
| scale | ρ_DE ∼ ρ_crit from CKN | DERIVED | §8.2 |
| no-go | energy bounds ⇒ area, not volume | DERIVED | `derivE` |
| residue ⟶ connectivity ⟶ range | gravity long-range (α≤d) ⇒ non-additive (volume class) | REDUCED | `residue` |
| BH vs cosmic | equilibrium (area) vs driven NESS (volume) | REDUCED | `residue` |
| theorem | drive × long-range locks volume | DERIVED (model) | `ness` |
| map 1 | J_eff = λ_max(W_grav) = binding ≥ J_c, generic | DERIVED | `jcoup` |
| map 2 | drive ∝ f_sat clears barrier at z*, locks | DERIVED | `deposdrv` |
| spinodal | ratchet through spinodal; robust basin; GSL lock | DERIVED (robust) | `socatt` |
| closure | three routes converge on the dS state count; C derives bistable F(m) | analysis | `closure` |
| dS count | count = horizon Hausdorff dim (Δ=d_H−2); equilibrium ⇒ area | REDUCED | `dscount` |
| λ ∝ F | KPZ λ = M·F/V (membrane); F = 0 ⇒ EW ⇒ Δ→1 | REDUCED | `kpzmem` |
| marginal | horizon 2-surface = EW lower critical dim ⇒ marginal, not robust | HONEST LIMIT | `dsmarg` |
| dS status | not resolved; two obstructions (O1 marginal, O2 Planck) + 1 open QG comp | OPEN | `dsresolve` |
| **residue** | **the m≈1 volume branch is the horizon's state** | **OPEN (L1/2) / measurable (L3)** | dS holography / DR3+Euclid |

**Bottom line.** The volume-law postulate — once a bare, untestable counting
assumption — is reduced to a single, sharply-posed, falsifiable statement: *the
de Sitter static patch's degrees of freedom count its bulk, not its boundary.*
Everything else is derived or reduced; the rest follows once that branch holds
(C); deriving it is dS holography (open, not shortcut-able); measuring it is Δ via
DESI DR3 + Euclid (imminent, decisive). The residue is not a weakness in the SEDE
argument but a clean hand-off to an open problem in quantum gravity, with a
measurement that settles it either way.

---

## 8. Reproducibility

```
python reproduce_all.py --only derivA,derivB,derivC,derivD,derivE   # route map §2
python reproduce_all.py --only residue                              # §3
python reproduce_all.py --only ness,jcoup,deposdrv,socatt           # §4–5 driven-NESS chain
python reproduce_all.py --only closure                              # §6
```

Figures: `output/criticality_soc.png` (criticality = volume-law),
`output/driven_ness.png` (the locking), `output/deposition_drive.png`
(deposition → drive at z*). Each script prints a verdict and runs validation
assertions; all stages PASS.
