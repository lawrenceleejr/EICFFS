# Fragmentation does not know about the laboratory: a beam-energy test at the EIC

*Study note for the EICFFS framework.  Motivated by L. Lee, C. Bell, J. Lawless,
C. Nash, E. Nibigira, "Experimental impact of jet fragmentation reference frames
at particle colliders", Phys. Lett. B 866 (2025) 139561,
[arXiv:2308.10951](https://arxiv.org/abs/2308.10951).*

---

## Summary

Jet fragmentation happens in the rest frame of the colour-connected system, not
in the laboratory.  In neutral-current DIS that system is the whole hadronic
final state, of invariant mass *W*, and its rest frame is the γ*p frame.  The
claim is testable at the EIC in a way that no e⁺e⁻ measurement can match,
because the accelerator itself provides the boost.

At fixed (*W*, *Q*) the colour-frame physics is fixed.  Changing the beam
energies then changes *only* the laboratory frame.  Comparing the 5 × 41,
10 × 100 and 18 × 275 GeV configurations moves the same colour-singlet system
through laboratory momenta differing by up to a factor of seven.  If
fragmentation is a property of the colour rest frame, ⟨n₉₀⟩ must not move.

It does not.  Measuring the lab-frame dependence as the exponent
d ln⟨n₉₀⟩ / d ln|p|_lab across the three configurations:

| what is measured | n₉₀ | n_SD |
|---|---|---|
| all hemispheres, one beam energy, no control | **+0.277** | non-monotonic |
| leading anti-kT *R* = 0.4 lab jet, fixed (*W*, *Q*) | +0.038 | −1.48 |
| whole current hemisphere, from lab momenta | −0.050 | −0.513 |
| γ*p-frame jet, from lab momenta | −0.001 | −0.430 |
| whole current hemisphere, from colour-frame momenta | −0.008 | −0.007 |
| γ*p-frame jet, from colour-frame momenta | +0.015 | +0.016 |

n_SD is the iterated soft-drop multiplicity, an IRC-safe counting observable
used here as a cross-check on n₉₀ (Sec. 3).  It gives the same answer where it
matters — no lab dependence once the observable is built in the colour frame —
but it is far more fragile when computed in the laboratory, because its
angular cut is a laboratory angle and a boost rescales angles.

The steep inclusive slope is not fragmentation responding to the laboratory.
It is the (*W*, *Q*) content of the sample changing along the axis.  Once the
colour-frame kinematics are held fixed and the observable is defined in that
frame, the residual dependence on a factor of seven in laboratory momentum is
consistent with zero (`figures/frame_ladder.pdf`).

Everything between those two extremes is a choice made in the laboratory: a
fixed cone keeps a boost-dependent share of the shower, ordering constituents
by laboratory momentum is not boost-invariant, and an angular cut applied in
the laboratory is rescaled by the boost.  All three are quantified below.

---

## 1. Why the EIC, and why beam energy is the only clean lever

DIS has two kinematic degrees of freedom.  Within one beam configuration the
laboratory frame is a deterministic function of them, so once (*W*, *Q*) is
fixed nothing is left to vary: the current system's laboratory momentum is then
pinned.  At *W* = 15–22 GeV the hemisphere's median |p|_lab runs 2.1, 4.2, 9.0
and 17.4 GeV across the *Q* = 2.2–3.3, 3.3–5, 5–7.5 and 7.5–11 GeV bins, with a
central-68 % spread of only about a factor of two inside each cell
(`hemisphere_p_vs_q.pdf`).  A laboratory-momentum scan at fixed colour-frame
energy is therefore a *Q* scan in disguise, and *Q* is a physical scale that the
shower is entitled to know about.

This is the essential difference from the e⁺e⁻ → ZZ example of
arXiv:2308.10951, where the Z boost varies independently of m_Z.  The EIC
recovers an independent boost knob by changing the beams:

| configuration | √s | γ*p-frame rapidity in the lab at *W* = 15–22 GeV |
|---|---|---|
| 5 × 41 GeV | 28.6 GeV | small |
| 10 × 100 GeV | 63.2 GeV | intermediate |
| 18 × 275 GeV | 140.7 GeV | large |

Every EIC run plan includes all three.  The comparison needs no new apparatus
and no unfolding to a theoretical frame: it is the same measurement repeated at
three beam settings.

---

## 2. The measurement

**Samples.**  Pythia 8.317, neutral-current DIS, *Q*² > 1 GeV², *W* > 10 GeV,
lepton-beam ISR off, DIS dipole-recoil shower on.  2.4 M events at 10 × 100 and
1.2 M each at 5 × 41 and 18 × 275 GeV.

**Cells.**  *W* ∈ {10–15, 15–22, 22–28} GeV × *Q* ∈ {2.2–3.3, 3.3–5, 5–7.5} GeV,
requiring at least 400 entries per configuration.

**Objects.**  The current hemisphere of the Breit frame taken whole; jets
clustered in the γ*p frame with an angular (e⁺e⁻-style) anti-kT algorithm at
*R* = 0.4 rad; and, for contrast, the leading anti-kT *R* = 0.4 jet clustered in
the laboratory.

**Observable.**  n₉₀, the interpolated number of constituents carrying 90 % of
the object's scalar momentum, as defined in arXiv:2308.10951 Sec. 2, computed
either from laboratory momenta or from colour-frame momenta.

### 2.1 The whole current hemisphere

⟨n₉₀⟩ per cell, with the hemisphere's median laboratory momentum in GeV
(`beam_energy_hemisphere.pdf`):

| cell | 5 × 41 | 10 × 100 | 18 × 275 | lever | exponent |
|---|---|---|---|---|---|
| *W* 10–15, *Q* 2.2–3.3 | 2.03 @ 1.9 | 1.91 @ 3.7 | 1.87 @ 9.4 | ×4.9 | −0.050 |
| *W* 10–15, *Q* 3.3–5 | 2.60 @ 3.5 | 2.53 @ 8.2 | 2.53 @ 22.5 | ×6.5 | −0.014 |
| *W* 10–15, *Q* 5–7.5 | 3.26 @ 6.5 | 3.24 @ 17.2 | 3.26 @ 48.5 | ×7.4 | −0.001 |
| *W* 15–22, *Q* 2.2–3.3 | 2.15 @ 1.8 | 1.99 @ 2.1 | 1.87 @ 4.5 | ×2.6 | −0.126 |
| *W* 15–22, *Q* 3.3–5 | 2.68 @ 2.5 | 2.58 @ 4.2 | 2.51 @ 10.7 | ×4.4 | −0.043 |
| *W* 15–22, *Q* 5–7.5 | 3.40 @ 3.8 | 3.29 @ 9.0 | 3.25 @ 24.6 | ×6.5 | −0.023 |
| *W* 22–28, *Q* 2.2–3.3 | 2.39 @ 2.3 | 2.05 @ 1.7 | 1.93 @ 2.7 | ×1.6 | −0.086 |
| *W* 22–28, *Q* 3.3–5 | 2.78 @ 3.0 | 2.62 @ 2.8 | 2.52 @ 6.1 | ×2.2 | −0.088 |
| *W* 22–28, *Q* 5–7.5 | 3.46 @ 3.6 | 3.36 @ 5.5 | 3.22 @ 14.0 | ×3.9 | −0.051 |

Statistical errors on each ⟨n₉₀⟩ are 0.005–0.02.  The cell with the longest
lever, *W* = 10–15 and *Q* = 5–7.5, moves the hemisphere from 6.5 to 48.5 GeV of
laboratory momentum and changes ⟨n₉₀⟩ by 0.02 out of 3.25.

The residual is not an artefact of wide cells.  Reweighting each configuration
to a common (*W*, *Q*) distribution on an 8 × 8 grid inside every cell moves the
median exponent from −0.050 to −0.045.

### 2.2 The residual is the observable's ordering frame, not the physics

At fixed (*W*, *Q*) the hemisphere contains the same particles with the same
colour-frame momenta whichever beams produced it.  What differs is that n₉₀
orders constituents by *laboratory* momentum, and ordering is not
boost-invariant: a large boost pushes the ordering towards light-cone momentum,
which concentrates the object in fewer particles.

Recomputing n₉₀ from colour-frame momenta for exactly the same hemispheres
removes most of the residual (`beam_energy_ordering.pdf`):

| cell | lab-ordered | frame-ordered |
|---|---|---|
| *W* 10–15, *Q* 2.2–3.3 | −0.050 | −0.031 |
| *W* 15–22, *Q* 2.2–3.3 | −0.126 | −0.073 |
| *W* 15–22, *Q* 5–7.5 | −0.023 | −0.003 |
| *W* 22–28, *Q* 2.2–3.3 | −0.086 | −0.006 |
| *W* 22–28, *Q* 5–7.5 | −0.051 | −0.008 |
| **median over nine cells** | **−0.050** | **−0.008** |

For narrow γ*p-frame jets the two orderings agree, −0.001 against +0.015, both
consistent with zero within the cell-to-cell scatter: a collimated object is
nearly unaffected by a longitudinal boost, so its ordering survives.  The
ordering effect is specific to wide objects.

### 2.3 What a laboratory cone does

Repeating the identical cell comparison with the leading anti-kT *R* = 0.4 jet
clustered in the laboratory gives a median exponent of +0.038, ranging to +0.14,
with the lowest-energy configuration breaking away in the higher-*Q* cells
(`beam_energy_labjet.pdf`).  The cause is direct: a fixed cone does not hold a
fixed share of the current system.  At *Q* = 5–7.5 GeV the leading current jet
carries 0.60 of the hemisphere's laboratory momentum at p_T = 2 GeV and 1.05 at
p_T = 8 GeV, and above unity it is sweeping in the target side as well
(`capture_fraction.pdf`).

Widening the cone does not fix this.  With *Q* and colour-frame energy both held
fixed within one beam configuration, the residual dependence on lab p_T is
+0.41, +0.40, +0.38, +0.35 and +0.27 for *R* = 0.4, 0.8, 1.2, 1.6 and 2.4, and
only +0.06 for the whole hemisphere (`slope_vs_radius.pdf`).  The geometry is
not the obstacle — the hemisphere is compact in the laboratory, the radius
containing 90 % of its momentum having median 0.39 and reaching 0.76 in the
highest *W* slice — the obstacle is that anti-kT selects and splits by
laboratory p_T while the shower is organised in another frame.

---

## 3. What the effect looks like with no control at all

The inclusive measurement is the one an experiment would make by default, and it
shows the frame-dependent shift at full strength.  For the current hemisphere at
10 × 100 GeV, ⟨n₉₀⟩ rises from 1.59 to 4.06 between |p|_lab ≈ 1.3 and 32 GeV,
an exponent of +0.277 (`hemisphere_vs_p.pdf`).  Slicing in colour-frame energy
does not flatten it (median spread 90 %), because along each slice the median
*Q* climbs — from 2.0 to 12 GeV in the *E*_cm = 4–6 GeV slice.

For laboratory jets the same story appears as a fan in *W*.  At fixed lab-frame
jet momentum ⟨n₉₀⟩ rises with *W* by +35 % to +83 % between *W* = 10–15 and
32–45 GeV (`ffs_fan.pdf`, `ffs_slopegraph.pdf`), and by construction those jets
also differ in *Q* and in the share of the shower their cone holds.

The point of Sec. 2 is that all of this structure is bookkeeping about what was
put in the bin, and none of it is fragmentation responding to the laboratory.

---

## 3. An IRC-safe alternative: iterated soft-drop multiplicity

n₉₀ has a defect worth stating plainly, because the original framework claimed
otherwise: **it is infrared safe but collinear unsafe.**  Splitting a
constituent into two collinear halves changes how many particles are needed to
reach 90 % of the momentum, and the change does not go away as the splitting
angle goes to zero.

The natural IRC-safe counting observable is the iterated soft-drop
multiplicity n_SD (Frye, Larkoski, Thaler, Zhou,
[arXiv:1704.06266](https://arxiv.org/abs/1704.06266)): recluster the
constituents with Cambridge/Aachen in opening angle, walk the hardest branch,
and count branchings passing z > z_cut (θ/R₀)^β with θ > θ_cut.  Here
z_cut = 0.1, β = 0, θ_cut = 0.1 rad, R₀ = 1 rad, in the e⁺e⁻ (energy, angle)
form so that the same definition applies in any frame.

### 3.1 The safety test

Taking real current hemispheres and deforming them
(`irc_safety_test.py`, `figures/irc_safety.pdf`):

| deformation | ⟨n₉₀⟩ | ⟨n_SD⟩ |
|---|---|---|
| unmodified | 2.843 | 1.738 |
| every constituent split in two at δ = 0.1 rad | 5.686 | 2.271 |
| … at δ = 0.03 | 5.686 | 1.738 |
| … at δ = 0.001 | 5.686 | 1.738 |
| three particles added at ε = 10⁻² | 3.074 | 1.752 |
| … at ε = 10⁻⁴ | 2.845 | 1.738 |

A democratic collinear split **doubles** n₉₀ and it stays doubled however small
the angle: n₉₀ is collinear unsafe, and the violation is order unity, not a
correction.  n_SD returns to its unmodified value once the splitting falls
below θ_cut.  Both are infrared safe.

### 3.2 What n_SD does to the frame picture

Repeating the beam-energy ladder with n_SD gives the striking result in the
summary table, and it cuts both ways.

*Where the observable is built in the colour frame, the two agree exactly:*
−0.007 against −0.008 for the hemisphere, +0.016 against +0.015 for the γ*p
jet.  That is as it must be — at fixed (*W*, *Q*) the colour-frame final state
is identical, so any frame-defined observable is invariant.

*Where the observable is built in the laboratory, n_SD is an order of magnitude
worse:* −0.513 against −0.050 for the hemisphere, −0.430 against −0.001 for the
γ*p jet, −1.48 against +0.038 for the lab cone.  In the longest-lever cell
(*W* = 10–15, *Q* = 5–7.5 GeV) the lab-computed n_SD falls 1.85 → 1.02 → 0.10
across the three configurations, while the frame-computed value sits at
2.12 → 2.12 → 2.13, constant to half a percent (`beam_energy_sd.pdf`).

The reason is structural.  n₉₀ is built only from momentum ordering and has no
angular scale, so a longitudinal boost disturbs it only mildly.  n_SD has an
explicit angular scale, and a boost rescales angles: at high beam energy the
hemisphere is collimated below θ_cut and the branchings simply stop being
counted.  **IRC safety and frame robustness are different axes**, and an
observable can be excellent on one and poor on the other.

### 3.3 A trap in the inclusive plot

Inclusively at 10 × 100 GeV, ⟨n_SD⟩ of the hemisphere against |p|_lab reads
0.82, 1.06, 1.32, 1.50, 1.53, 1.34, 1.13, 0.82 across bins from 1.3 to 32 GeV.
It rises and then falls: the *Q* growth that adds branchings is cancelled at
high momentum by the collimation that pushes them below θ_cut.  A power-law fit
returns +0.010, which would read as "perfectly frame independent" and is
meaningless — the curve is not a power law, and the controlled test shows this
is the *most* lab-sensitive of all the definitions tried.  A flat inclusive
curve is not evidence of frame independence.

### 3.4 Which to use

For a measurement, n_SD computed in the colour rest frame is the best of both:
IRC safe and frame independent to under two percent.  n₉₀ in the colour frame
is equally frame independent and simpler to construct, at the cost of collinear
unsafety, which matters for comparison to fixed-order or resummed calculations
but not for a Monte Carlo comparison at hadron level.  n_SD computed in the
laboratory should be avoided in this context.  Its mean is also small at EIC
energies — 1.14 for the hemisphere at 10 × 100 — so its statistical reach is
weaker than n₉₀'s.

---

## 4. Method and framework

Pipeline: `generate_events.py` (Pythia 8, particle-level, stores q, k′, the beam
proton and the struck parton), `analyze_events.py` (jet finding in the lab and
in the γ*p frame, frame boosts, n₉₀, per-jet and per-event trees),
`make_figures.py` (one panel per PDF).  `./run.sh` reproduces the 10 × 100
sample in about five minutes on four cores; the two other beam configurations
are the same command with different `--electron-energy` and `--proton-energy`.

Verified: the n₉₀ definition and interpolation reproduce arXiv:2308.10951
Sec. 2, and the vectorised implementation agrees with a scalar reference to
1 × 10⁻¹¹; the Breit construction gives exactly zero boson energy and |q| = *Q*;
*W* reconstructed from the stored hadronic final state agrees with the DIS
invariants to 0.03 GeV on average, the residual being neutrinos from
heavy-flavour decay.

### Audit of the original framework

The premise of the original code was right, but it could not have measured the
effect.  All items are fixed on this branch.

| severity | finding | fix |
|---|---|---|
| **blocking** | The hard-process lepton was found with `status == 23` in the event record. The shower copies it, so the record carries −23 in ~96 % of events; those were silently dropped and the surviving 4 % are a biased subset with no lepton radiation or recoil. | Read it from `pythia.process`, where the status is always +23. Efficiency 100 %. |
| major | Lepton-beam ISR left on, so the stored q = k_beam − k′ included the radiated photon and was not the exchanged boson; *W*, *Q*², *x*, *y* were mis-reconstructed for radiative events. | `PDF:lepton = off`, the standard EIC choice; a flag restores it. |
| major | Breit frame and γ*p frame conflated. They differ by a boost along the boson axis; the struck quark carries *Q*/2 in one and *W*/2 in the other. Only the γ*p frame is the colour rest frame. | Both boosts computed and stored per jet; text corrected. |
| major | The sign of the predicted effect was stated backwards. | Corrected, and confirmed by simulation. |
| major | All jets used, including proton-remnant fragments at forward rapidity. | Jets flagged by Breit hemisphere; figures use current jets. |
| moderate | FastJet installed but never used: the run script never passed the flag, so a greedy cone was the production algorithm. | Anti-kT via FastJet by default. |
| moderate | No dipole-recoil shower for DIS; an arbitrary rapidity-ordering change instead. | `SpaceShower:dipoleRecoil = on`. |
| moderate | Uncertainty on ⟨n₉₀⟩ taken as √(mean/N), a Poisson assumption invalid for an interpolated observable. | Standard error of the mean. |
| moderate | The default 200 k sample has median *W* ≈ 8 GeV; events with *W* > 30 GeV and a current jet are ~1 % of it. | Generation-time *W* cut, parallel seeds. |
| minor | All electrons and muons removed from the final state rather than the scattered lepton and its radiation. | Removal by ancestry of the hard-process lepton. |

---

## 5. Caveats

* Pythia only.  arXiv:2308.10951 compared Pythia, Vincia and Herwig; a Herwig 7
  or Sherpa cross-check of the ladder is the obvious next step and is a
  generator swap in this framework.
* Detector effects are absent.  The three configurations put the current system
  at different rapidities, so the acceptance and momentum resolution of the same
  detector differ between them.  A frame-independence test is only as good as
  the control of that difference, and it is the main experimental risk.
* The hemisphere requires reconstructing the Breit frame event by event, so it
  inherits the resolution of the scattered-lepton or hadronic-method kinematics.
* The *W* = 22–28 GeV cells sit at the edge of the 5 × 41 phase space
  (W_max ≈ 28.6 GeV, y → 1) and carry the largest residuals; a measurement would
  either drop them or narrow the binning.
* *Q*² > 1 GeV² admits jets from boson–gluon fusion and QCD Compton. The colour
  rest frame is still the γ*p frame, but the jet is not the Born struck quark;
  the stored parton match isolates the Born-like sample.
* Two counting observables are studied here. A continuous jet shape (thrust,
  angularities, energy–energy correlators) should be added before claiming
  generality, and the n_SD working point should be scanned: θ_cut and z_cut set
  how strongly the lab-frame version reacts to the boost.
* n_SD at EIC energies is small, mean 1.14 for the hemisphere at 10 × 100,
  because the current system typically holds only three to five particles. The
  ladder is limited by that, not by sample size.

---

## 6. Prior work

Checked on 2026-09-05.

* **Citations of arXiv:2308.10951** on INSPIRE: three, on heavy-ion jet
  background subtraction, top-jet classification and quark/gluon tagging in CMS
  open data.  None mentions DIS, ep, the EIC, HERA or the Breit frame.
* **arXiv searches** for "colour/color rest frame", "frame-dependent
  fragmentation" with jets, and EIC + jet + fragmentation + frame return the
  reference paper and unrelated hits.
* **EIC literature.**  The Yellow Report
  ([2103.05419](https://arxiv.org/abs/2103.05419)), the EIC jet overview
  ([1911.00657](https://arxiv.org/abs/1911.00657)) and the NC-DIS substructure
  study ([2302.06941](https://arxiv.org/abs/2302.06941)) compare laboratory and
  Breit-frame jet finding and yields.  None studies a fragmentation observable
  at fixed laboratory momentum as a function of *W*, and none proposes a
  beam-energy comparison at fixed colour-frame kinematics.
* **HERA.**  ZEUS ([0803.3878](https://arxiv.org/abs/0803.3878)) and H1
  ([hep-ex/9707005](https://arxiv.org/abs/hep-ex/9707005)) measured
  multiplicities in the Breit current hemisphere against *Q* and in the γ*p
  current region against *W*, testing universality against e⁺e⁻.  Those are
  frame-corrected measurements designed to remove the boost; the test proposed
  here uses the boost as the independent variable.  Breit-frame jet algorithms
  for the EIC (Centauro, [2006.10751](https://arxiv.org/abs/2006.10751)) are
  likewise a way to undo the boost.
* HERA ran at essentially one beam-energy configuration for most of its life,
  with a short low-energy run at the end; the EIC's three configurations with
  comparable luminosity are what makes this test practical.

No prior proposal of a beam-energy frame-independence test of fragmentation was
found.
