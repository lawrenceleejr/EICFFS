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

| what is measured | n₉₀ | n_SD, standard form |
|---|---|---|
| all hemispheres, one beam energy, no control | **+0.277** | non-monotonic |
| leading anti-kT *R* = 1.2 lab jet, fixed (*W*, *Q*) | +0.011 | +0.007 |
| whole current hemisphere, from lab momenta | −0.050 | −0.018 |
| γ*p-frame jet, from lab momenta | −0.001 | +0.016 |
| whole current hemisphere, from colour-frame momenta | −0.008 | −0.007 |
| γ*p-frame jet, from colour-frame momenta | +0.015 | +0.016 |

n_SD is the iterated soft-drop multiplicity, an IRC-safe counting observable
used here as a cross-check on n₉₀ (Sec. 3).  In its standard hadron-collider
form — transverse-momentum fractions and a rapidity–azimuth distance, which are
invariant under a boost along the axis — it can be evaluated on laboratory
momenta as they are, and it agrees with n₉₀ on every rung.  (The colour-frame
rows use the e⁺e⁻ form evaluated in that frame, where the two forms coincide.)
The e⁺e⁻ form with an absolute opening-angle cut, applied in the laboratory, is
a different matter: a boost rescales angles and the observable collapses,
−0.47 for the hemisphere.  Sec. 3.5 separates the two.

The steep inclusive slope is not fragmentation responding to the laboratory.
It is the (*W*, *Q*) content of the sample changing along the axis.  Once the
colour-frame kinematics are held fixed and the observable is defined in that
frame, the residual dependence on a factor of seven in laboratory momentum is
consistent with zero (`figures/frame_ladder.pdf`).

Everything between those two extremes is a choice made in the laboratory, and
none of it is a property of fragmentation.  `figures/frame_breakers.pdf` ranks
the choices by the laboratory dependence they introduce on the same test: an
opening-angle cut applied in the laboratory (−0.47), no (*W*, *Q*) control at
all (+0.28), a laboratory momentum threshold on every particle (+0.25),
ordering a wide object's constituents by laboratory momentum (−0.05), and a lab
cone at *R* = 0.4 rather than the *R* ≈ 1 the EIC community uses (+0.04).  A
transverse-momentum threshold, the boost-invariant version of the third,
introduces +0.01.  All are quantified below.

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
*R* = 0.4 rad; and, for contrast, the leading anti-kT jet clustered in the
laboratory at *R* = 1.2, the radius EIC studies use, with *R* = 0.4 to 2.4 in
the radius scan.

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

Repeating the identical cell comparison with the leading anti-kT jet clustered
in the laboratory at *R* = 1.2, the radius EIC studies use, gives a median
exponent of +0.011, against +0.22 for the inclusive lab-jet curve drawn behind
the cells (`beam_energy_labjet.pdf`): a plain lab cone of sensible size is
already frame independent.  The same panel for soft-drop multiplicity in its
standard form (`beam_energy_labjet_sd.pdf`) gives +0.007 against +0.13 for the
inclusive curve: an IRC-safe observable, on lab momenta, in a lab cone, with
nothing boosted, and it is flat.  At the *R* = 0.4 of the reference paper the
exponent is +0.038, ranging to +0.14, with the lowest-energy configuration
breaking away in the higher-*Q* cells (`ladder_vs_radius.pdf`).  The cause is
direct: a fixed small cone does not hold a fixed share of the current system.  At *Q* = 5–7.5 GeV the leading current jet
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

Repeating the beam-energy ladder with n_SD in the e⁺e⁻ form gives a striking
result, and it cuts both ways.

*Where the observable is built in the colour frame, the two agree exactly:*
−0.007 against −0.008 for the hemisphere, +0.016 against +0.015 for the γ*p
jet.  That is as it must be — at fixed (*W*, *Q*) the colour-frame final state
is identical, so any frame-defined observable is invariant.

*Where the observable is built in the laboratory with an absolute opening-angle
cut, n_SD is an order of magnitude worse:* −0.513 against −0.050 for the
hemisphere, −0.430 against −0.001 for the γ*p jet, −1.48 against +0.038 for the
lab cone.  This is specific to the e⁺e⁻ form of the condition; Sec. 3.5 shows
the standard pp form is flat in the lab.  In the longest-lever cell
(*W* = 10–15, *Q* = 5–7.5 GeV) the lab-computed n_SD falls 1.85 → 1.02 → 0.10
across the three configurations, while the frame-computed value sits at
2.12 → 2.12 → 2.13, constant to half a percent (`beam_energy_sd.pdf`).

The reason is structural.  n₉₀ is built only from momentum ordering and has no
angular scale, so a longitudinal boost disturbs it only mildly.  An absolute
opening-angle cut does have a scale, and a boost rescales angles: at high beam
energy the hemisphere is collimated below θ_cut and the branchings stop being
counted.  The fix is the one hadron colliders already adopted — measure the
angle as a rapidity–azimuth distance, which is invariant under boosts along the
axis (Sec. 3.5).  **Frame robustness is a property of the variables, not of the
grooming**, and the lesson is about how a cut is written rather than about
soft-drop multiplicity itself.

### 3.3 A trap in the inclusive plot

Inclusively at 10 × 100 GeV, ⟨n_SD⟩ of the hemisphere against |p|_lab reads
0.82, 1.06, 1.32, 1.50, 1.53, 1.34, 1.13, 0.82 across bins from 1.3 to 32 GeV.
It rises and then falls: the *Q* growth that adds branchings is cancelled at
high momentum by the collimation that pushes them below θ_cut.  A power-law fit
returns +0.010, which would read as "perfectly frame independent" and is
meaningless — the curve is not a power law, and the controlled test shows this
is the *most* lab-sensitive of all the definitions tried.  A flat inclusive
curve is not evidence of frame independence.

### 3.4 Why n_SD responds so weakly to the physics

Two reasons, and the second is the important one.

**It is a coarse counter at EIC multiplicities.**  The current hemisphere holds
3.09 constituents on average — 20 % have one, 26 % have two — so there is
almost no logarithmic phase space for a grooming counter.  ⟨n_SD⟩ = 0.92 with
39 % of hemispheres giving exactly zero and 77 % giving zero or one, against
⟨n₉₀⟩ = 2.21 with a standard deviation of 1.28 and a tail to 18.  n_SD is close
to a binary variable here, so its statistical reach per event is poor.

**What variation it does show is mostly the angular cut, not fragmentation.**
Hold *Q* fixed at 3.3–5 GeV and slice in *W*:

| *W* | ⟨N_const⟩ | ⟨n₉₀⟩ | ⟨n_SD⟩ γ*p frame | ⟨n_SD⟩ object rest frame |
|---|---|---|---|---|
| 10–15 | 3.85 | 2.685 | 1.651 | 1.878 |
| 15–22 | 3.87 | 2.697 | 1.489 | 1.885 |
| 22–32 | 3.84 | 2.678 | 1.235 | 1.873 |
| 32–45 | 3.83 | 2.671 | 0.933 | 1.864 |
| **change** | **−0.4 %** | **−0.5 %** | **−43.5 %** | **−0.7 %** |

The hemisphere's particle content is *unchanged* across this range — that is the
classic HERA result that current-region multiplicity is set by *Q*, not by *W* —
and n₉₀ correctly reports nothing happening.  n_SD in the γ*p frame falls by
43 % anyway.  The reason is that the γ*p frame is the rest frame of the *whole*
hadronic system, but the current hemisphere still moves within it, with energy
≈ *W*/2 along the boson axis.  As *W* grows the same particles are collimated
into smaller opening angles, they fall below θ_cut = 0.1 rad, and the branchings
stop being counted.

Boosting into the object's own rest frame before applying the cut removes this
entirely: −0.7 % across the same range, and −0.016 on the beam-energy ladder.
The lesson generalises: **an absolute angular cut inherits whatever boost the
object has in the frame where the cut is applied**, and identifying the colour
rest frame of the event is not sufficient — the object of interest must also be
at rest.

### 3.5 Which frame, and which variables

Soft drop is normally written for pp with transverse-momentum fractions and an
angular distance in rapidity–azimuth,
z = min(p_T1, p_T2)/(p_T1 + p_T2) > z_cut (ΔR₁₂/R₀)^β
([arXiv:1402.2657](https://arxiv.org/abs/1402.2657)), and iterated soft drop
inherits that convention.  Those variables are invariant under boosts along the
*beam*, so in pp the choice of frame is neutralised by the variables rather than
by boosting anything.  The e⁺e⁻ form uses energy fractions and opening angles,
equivalent in the small-angle limit, and is applied in the CM frame where there
is no ambiguity.

Neither is automatically safe in DIS, because the boost relating the laboratory
to the colour rest frame runs along *P + q*, not along the beam.  The variant
used above — e⁺e⁻ variables in the laboratory — is the worst case, and that is
why it gave −0.474.  Measuring the *standard pp variables about the P + q axis*
fixes it with no boosting at all (`sd_frame_test.py`,
`figures/sd_frame_choice.pdf`):

| prescription | exponent | ⟨n_SD⟩ |
|---|---|---|
| e⁺e⁻ variables (E, θ) in the laboratory | −0.474 | 1.37 |
| **pp variables (p_T, ΔR) about the beam axis, in the lab** | **−0.016** | 1.48 |
| **pp variables (p_T, ΔR) about the P + q axis, in the lab** | **−0.009** | 1.49 |
| e⁺e⁻ variables in the object's own rest frame | −0.011 | 1.67 |

The important row is the second.  **Written the standard way, soft-drop
multiplicity is frame independent in the laboratory with no modification at
all.**  The literal pp prescription — p_T fractions and a rapidity–azimuth
distance about the beam — already gives −0.016, because at EIC kinematics the
γ*p system's momentum is dominated by the proton, so the boost axis is close to
the beam.  Using the exact P + q axis improves it slightly, to −0.009, and is
the principled choice, but the difference is small.

Every large soft-drop number reported in this note therefore belongs to the
e⁺e⁻ form with an absolute opening-angle cut, which is what the implementation
here used.  That form is the outlier, not the observable.

**Is the object rest frame well defined?**  Mathematically yes whenever the
object has m² > 0 and at least two constituents, and it preserves IRC safety: a
collinear splitting leaves the total four-momentum, and hence the boost,
exactly unchanged, while a soft addition moves it continuously.  But it has
three drawbacks.  A two-body object is *exactly* back to back in its own rest
frame — θ = 180.0°, and z = 0.500 for equal masses — so all angular information
is destroyed, and 26 % of Breit current hemispheres hold exactly two particles.
The boost factor is E/m, set by the object's mass, which for a narrow jet is
dominated by its softest wide-angle structure, so the frame is formally safe but
numerically jumpy.  And it is not the standard convention, so results would not
be directly comparable with pp or e⁺e⁻ measurements.

The recommendation is therefore the middle row: keep the standard soft-drop
variables and measure them about the P + q axis.  It is the same prescription
the Breit-frame jet algorithms use for the same reason.

*(That two-body degeneracy also exposed an implementation bug: reclustering at
a radius of exactly π leaves a back-to-back pair tied against the beam distance,
so it is never merged and the declustering finds no branchings.  Every n_SD
number here was regenerated with the radius above π.)*

### 3.6 Which to use

For a measurement, use n_SD in its standard hadron-collider form — p_T
fractions and a rapidity–azimuth distance, ideally about the *P + q* axis.  It
is IRC safe and frame independent in the laboratory with no boosting, and it
stays comparable with pp and e⁺e⁻ results.  The object's own rest frame works
equally well but is non-standard and degenerate for two-body objects.  Only the
e⁺e⁻ opening-angle form applied in the laboratory should be avoided.  n₉₀ in the colour frame is equally
frame independent and simpler to construct, at the cost of collinear unsafety,
which matters for comparison to fixed-order or resummed calculations but not for
a Monte Carlo comparison at hadron level.  Whichever variant is chosen, n_SD is
small at EIC energies — 1.45 in the object rest frame, with 21 % of hemispheres
giving zero — so n₉₀ carries more information per event despite being collinear
unsafe.

---

## 3.7 Dropping cones altogether

If a fixed cone is the problem, the obvious move is to abandon it and always
measure the whole current region.  That works, and it is the recommendation,
but the current region has choices of its own.  `object_choice_test.py` runs the
beam-energy test over them (`figures/object_choice.pdf`):

| definition | exponent | ⟨N_const⟩ | events with ≤ 1 particle |
|---|---|---|---|
| Breit current hemisphere | −0.012 | 3.2 | 39 % |
| … with lab \|p\| > 0.5 GeV | **+0.166** | 2.1 | 39 % |
| … with lab p_T > 0.15 GeV | +0.009 | 2.6 | 26 % |
| γ*p current region, y* > 0 | +0.005 | 8.4 | 0 % |
| … with lab \|p\| > 0.5 GeV | **+0.245** | 4.7 | 2 % |
| … with lab p_T > 0.15 GeV | +0.014 | 5.9 | 0 % |

Three things follow.

**The cone pathology is genuinely gone.**  The radius scan ran from +0.41 at
*R* = 0.4 to +0.06 for the whole hemisphere; the current region carries no
residual of it.

**Detector angular acceptance does not spoil it.**  This was the obvious worry,
and it is unfounded at EIC energies: 99.7 % of the current hemisphere's momentum
lies inside |η| < 3.5, and the exponent is unchanged at |η| < 4 or 3.5
(−0.012, −0.018).  It is the *target* region that disappears down the
beampipe.

**A momentum threshold does spoil it, and this is the replacement for the cone
effect.**  Requiring |p|_lab > 0.5 GeV per particle gives +0.166 for the Breit
hemisphere and +0.245 for the γ*p region — worse than the *R* = 0.4 cone.  The
reason is elementary: |p| is not invariant under a boost along the axis, so the
same particle passes at high beam energy and fails at low.  A **transverse**
momentum threshold is invariant under that boost, and costs almost nothing:
p_T > 0.15 GeV gives +0.009 and +0.014.  Any per-particle cut used in such a
measurement should be transverse.

**Which current region.**  The γ*p region (y* > 0) is the better object of the
two.  It holds 8.4 particles against the Breit hemisphere's 3.2 and has no empty
events, against 39 % of Breit hemispheres with one particle or fewer — which is
why the grooming observable of Sec. 3 failed there.  The two are not the same
measurement, though: at fixed *Q* the Breit hemisphere is *W*-independent
(N_const −2 %, ⟨n₉₀⟩ −2 % from *W* = 10–15 to 32–45, the classic HERA
current-region result), while the γ*p region grows with the string
(N_const +54 %, ⟨n₉₀⟩ +32 %).  The Breit hemisphere is the *Q*-only object; the
γ*p region is the one that responds to the colour-frame energy, and for a
frame-dependence study that is the point.

What remains unaddressed is reconstruction: both definitions need the frame
built event by event, so both inherit the resolution of the scattered-lepton or
hadronic-method kinematics, and particles near the boundary migrate.  That is
not simulated here and is the leading experimental question.

---

## 3.8 Relation to EIC jet conventions, and a correction

**What the field does.**  Most EIC simulation studies cluster *inclusively in the
laboratory* with anti-kT at a large radius, **R ≈ 1.0**, chosen for the low
hadron multiplicity, with inputs typically p_T > 250 MeV and |η| < 4.5
(lepton–jet correlations [1812.08077](https://arxiv.org/abs/1812.08077),
jet-based Sivers and Collins asymmetries
[2007.07281](https://arxiv.org/abs/2007.07281), heavy-flavour and jet studies
[2007.14417](https://arxiv.org/abs/2007.14417)).  On the theory side the Breit
frame is standard and growing: Centauro
([2006.10751](https://arxiv.org/abs/2006.10751)) was built to cluster the Born
configuration there, and a recent inclusive generalised-kT algorithm for DIS
([2606.13077](https://arxiv.org/abs/2606.13077)) is defined in the Breit frame
and shipped in fjcontrib, noting that in DIS, unlike pp, "no single algorithm
dominates the landscape".  Exclusive algorithms — Durham and Cambridge adapted
from e⁺e⁻ by treating the remnant as a particle of infinite momentum — are the
HERA legacy and have largely given way to inclusive longitudinally-invariant kT
in the Breit frame.  Separately, and without any clustering, there is the HERA
current-region tradition: multiplicities and fragmentation functions measured in
the current hemisphere of the Breit frame.  The object recommended in Sec. 3.7
sits in that tradition, so it is conventional rather than exotic.

**The correction.**  The *R* = 0.4 used throughout this note was inherited from
arXiv:2308.10951, not from EIC practice.  Repeating the beam-energy ladder at
the radii the field actually uses weakens the claim that a lab cone spoils
frame independence:

| leading lab jet | n₉₀ exponent | n_SD exponent, standard form | n_SD exponent, e⁺e⁻ form | ⟨N_const⟩ |
|---|---|---|---|---|
| *R* = 0.4 | +0.038 | +0.051 | −1.480 | 2.6 |
| *R* = 0.8 | +0.035 | +0.029 | −0.545 | 4.1 |
| *R* = 1.2 | +0.011 | +0.007 | −0.329 | 5.2 |
| *R* = 1.6 | +0.003 | −0.004 | −0.302 | 6.1 |
| *R* = 2.4 | −0.002 | −0.024 | −0.346 | 7.4 |

At the EIC-standard *R* ≈ 1 a laboratory cone is already frame independent for
n₉₀ to about a percent, comparable to the γ*p-frame jet (+0.015) and the
hemisphere (−0.008), and soft-drop multiplicity in its standard form tracks n₉₀
radius for radius.  On the beam-energy test the cone is therefore *not* the
problem; the problems are running with no control at all (+0.277) and
computing an opening-angle-cut observable in the laboratory (−0.3 to −1.5, at
every radius).  The summary table and `frame_ladder.pdf` therefore quote the
lab cone at *R* = 1.2; `ladder_vs_radius.pdf` shows the full scan.

**What this means for the radius scan of Sec. 3.7.**  That scan held *Q* and
colour-frame energy fixed *within one beam configuration* and varied lab p_T,
finding +0.41 at *R* = 0.4 falling only to +0.06 for the hemisphere.  That test
is not purely a frame test: at fixed (E_cm, *Q*) the lab p_T of a jet is set by
its angle θ* to the boson axis in the colour frame, so part of the residual is
genuine θ* dependence of the fragmentation rather than cone bookkeeping.  The
capture-fraction measurement (0.60 → 1.05 across p_T at *R* = 0.4) is direct
evidence that cone bookkeeping contributes, but the split between the two was
not separated, and the beam-energy ladder above is the cleaner statement.

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

* **Are *W* and *Q* measurable cell by cell?**  Yes, but not equally well in
  every cell.  *Q*² and *y* come from the scattered electron alone
  (*Q*² = 2*E*_e*E*'(1 + cos θ), *y* = 1 − (*E*'/2*E*_e)(1 − cos θ)), or from
  the hadronic final state (Jacquet–Blondel), and *W*² = *M*² + *ys* − *Q*², so
  *W* is a *y* measurement in disguise.  The electron method resolves *Q*² to a
  few per cent everywhere but its *y* resolution degrades as 1/*y*; the
  hadronic and Σ methods take over below *y* ≈ 0.05.  At fixed (*W*, *Q*) the
  three beam configurations sit at very different *y* (medians over the cells
  used here):

  | *W* cell | 5 × 41 | 10 × 100 | 18 × 275 |
  |---|---|---|---|
  | 10–15 GeV | *y* = 0.19–0.23 | 0.039–0.048 | 0.008–0.010 |
  | 15–22 GeV | 0.35–0.44 | 0.084–0.093 | 0.017–0.019 |
  | 22–28 GeV | 0.63–0.77 | 0.16 | 0.031–0.033 |

  The two ends of the lever arm are the hard ones.  At 18 × 275 the *W* = 10–15
  cells lie at *y* < 0.01, below the *y* > 0.01 cut EIC studies normally apply,
  and *W* there rests entirely on the hadronic method at the edge of its
  resolution.  At 5 × 41 the *W* = 22–28 cells lie at *y* = 0.6–0.8, where the
  scattered electron carries only 1–2 GeV and must be separated from
  photoproduction background.  The middle of the table, *y* between about 0.02
  and 0.5, is comfortable in every configuration; the *W* = 15–22 row is
  measurable in all three beams with standard methods, and even the extreme
  cells trade only lever arm, not the test itself.  *x* runs from 0.01 to 0.18
  across the cells and is identical between beams by construction.
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
