# Analysis Design: Observing Frame-Dependent Fragmentation at the EIC

**Target:** A phenomenology-level study worthy of a PRL, demonstrating that the
Frame-dependent Fragmentation Shift (FFS) effect of arXiv:2308.10951
(*Phys. Lett. B* 866 (2025) 139561) produces a large, unambiguous, and
measurable breakdown of lab-frame jet universality in neutral-current DIS at
the Electron–Ion Collider — and that binning jet structure in the
color-rest-frame momentum *restores* universality.

---

## 1. Executive summary — the PRL claim

> **Claim.** Jet substructure at the EIC, at *fixed lab-frame jet momentum*,
> depends strongly on the photon–proton invariant mass *W*: jets with
> identical lab kinematics differ by up to O(50%) in ⟨n₉₀⟩ across the
> accessible *W* range. This breaks the lab-frame universality implicitly
> assumed by jet calibration and tagging. When the same jets are instead
> binned by their momentum in the γ\*p center-of-momentum frame — the rest
> frame of the color-connected system — all *W* slices collapse onto a single
> universal curve. The EIC, with its multiple beam-energy configurations and
> event-by-event reconstruction of the boost (via the scattered electron),
> is the *only* collider where this frame can be reconstructed exactly,
> event by event. With ≥10 fb⁻¹ the effect is observable at ≫5σ; the
> measurement is systematics-limited, and the design below is built around
> ratio observables in which dominant systematics cancel.

Three ingredients make this PRL-worthy rather than a routine MC study:

1. **A falsifiable two-sided prediction** — *splay in the lab frame, collapse
   in the color frame* — with an explicit frame-independent null hypothesis,
   not merely "a trend in MC."
2. **Decoupling from known physics** — the analysis explicitly disentangles
   the FFS boost effect from Q² evolution of fragmentation (DGLAP), phase-space
   selection sculpting, and target-fragmentation contamination, which are the
   objections any referee will raise first.
3. **A uniquely-EIC experimental handle** — the same lab-frame jet momentum
   is sampled at very different *W* using the three EIC beam configurations
   with the *same detector*, so detector systematics largely cancel in
   cross-configuration ratios.

---

## 2. Physics background

### 2.1 The FFS effect

Hadronization models (Lund string, cluster) perform fragmentation in the rest
frame of the color-connected system, not the lab frame. arXiv:2308.10951
demonstrated this with e⁺e⁻ → 3 jets (color frame = lab frame) vs.
e⁺e⁻ → ZZ → 4 jets (color frame = boosted Z rest frame): for the ZZ case the
jet constituent structure is set by m_Z/2 ≈ 45 GeV *regardless of lab
momentum*, giving up to ~50% shifts in ⟨n₉₀⟩ at fixed lab |p|. A jet is
therefore not factorizable from its color-connected siblings.

### 2.2 Translation to NC DIS

In ep → e′ + X, the struck quark is color-connected to the proton remnant.
The rest frame of the full color-connected system (struck quark ⊕ remnant —
the Lund string) is the **γ\*p CM frame**, whose invariant mass is *W*:

| Quantity | Definition | Role |
|---|---|---|
| Q² = −q² | photon virtuality | hard scale (controls DGLAP evolution) |
| x = Q²/(2P·q) | Bjorken x | struck-parton momentum fraction |
| y = P·q/P·k | inelasticity | lab↔CM boost lever |
| W = √((P+q)²) | γ\*p CM energy | **color-frame energy: the FFS variable** |

In the γ\*p frame the current jet carries momentum ≈ W/2 — the exact analogue
of the m_Z/2 anchor in the ZZ benchmark. A jet selected at fixed
**lab-frame** momentum |p|_lab therefore corresponds to a CM-frame momentum
that grows with *W*; its string fragments at the W-scale, so its constituent
structure (⟨n₉₀⟩, multiplicity) **rises with W at fixed |p|_lab**.
Frame-independent fragmentation predicts a *flat* dependence (up to slow
ln Q² evolution, which we control by binning).

Crucially, unlike at the LHC, the boost between the lab and the color frame
is **reconstructed exactly event-by-event** from the scattered electron.
This turns the FFS effect from a population-level statement into a per-event
kinematic mapping — the strongest possible experimental configuration.

---

## 3. Hypotheses and money plots

- **H₀ (frame-independent fragmentation):** at fixed (|p|_lab, Q²), ⟨n₉₀⟩ is
  independent of W. Jet structure is a universal function of lab momentum.
- **H₁ (FFS):** at fixed (|p|_lab, Q²), ⟨n₉₀⟩ is a function of the jet
  momentum in the γ\*p frame, p_CM(|p|_lab, W, η_jet); ⟨n₉₀⟩ rises with W.
  Universality is restored when binning in p_CM.

**Figure plan for the Letter (4 figures):**

| Fig | Content |
|---|---|
| **1** | Concept + reach: boost factor between lab and γ\*p frame across the EIC (W, η_jet) plane for the three beam configurations; arrows showing how a fixed-|p|_lab jet maps to different p_CM. |
| **2** | **The effect:** ⟨n₉₀⟩ vs W in fixed (|p|_lab, Q²) bins. Pythia/Herwig/Sherpa envelope as a band; H₀ shown as the flat line anchored at the lowest-W bin. Lower panel: ratio to H₀. |
| **3** | **Universality restoration (the signature plot):** left panel, ⟨n₉₀⟩ vs |p|_lab in W slices (curves splay); right panel, ⟨n₉₀⟩ vs p_CM in the same W slices (curves collapse). Collapse quantified by χ²/ndf in the panel. |
| **4** | **Measurability:** projected significance of H₁ vs H₀ as a function of integrated luminosity, after detector smearing, for 10×100 and 18×275; with and without cross-beam-energy combination. |

---

## 4. Simulation campaign

### 4.1 Beam configurations (ePIC reference)

| Config | E_e × E_p | √s | Accessible W (jet analysis) | Role |
|---|---|---|---|---|
| **A** | 5 × 41 GeV | 28.6 GeV | 8–25 GeV | low-W anchor |
| **B** | 10 × 100 GeV | 63.2 GeV | 10–55 GeV | core dataset (current default) |
| **C** | 18 × 275 GeV | 140.7 GeV | 20–130 GeV | high-W lever arm |

The overlap regions (same |p|_lab, same Q², different beam config) give
same-detector cross-checks in which luminosity, tracking, and JES
systematics cancel in ratios.

### 4.2 Generators and physics variations

| Sample | Purpose |
|---|---|
| **Pythia 8.3** (Lund string), Monash tune | baseline; string fragments in color frame → embodies FFS |
| Pythia 8.3, 2 alternative tunes (e.g. Monash variations of `StringZ`/`StringPT`) | hadronization-parameter envelope |
| **Herwig 7** (cluster hadronization), DIS matrix elements | independent hadronization paradigm — the FFS prediction must survive |
| **Sherpa 3** (cluster + CSS shower) | third paradigm; shower-model spread |
| Pythia, MPI **on** (default is off) | low-Q² photon-structure contamination check |
| Pythia, QED radiation on/off | radiative-correction robustness of kinematic reconstruction |

The **theory band in Fig. 2/3 is the envelope of these samples**. The PRL
statement is that the FFS trend (slope in W, collapse in p_CM) is common to
all of them while its magnitude varies — i.e., the EIC measurement is also a
*discriminator between hadronization models*, which strengthens the Letter.

**H₀ pseudo-data:** constructed from the baseline sample by freezing the
⟨n₉₀⟩(|p|_lab, Q²) map measured in the lowest-W slice and broadcasting it to
all W — a concrete, simulation-derived frame-independent null.

### 4.3 Statistics

NC DIS at Q² > 25 GeV² has σ ~ O(10 nb); jet selection efficiency makes
statistics a non-issue at fb⁻¹ luminosities. Generate **2M events per
(config, generator, variation)** ≈ 30M events total. Effect size
Δ⟨n₉₀⟩ ≈ 1–2 units across the W range vs per-jet spread σ(n₉₀) ≈ 1.5 →
per-bin precision of 0.05 needs only ~10³ jets/bin: MC statistics must
exceed projected data statistics, hence 2M/sample.

---

## 5. Event and object selection

### 5.1 DIS event selection

- Scattered electron found; reconstruct kinematics with the **electron
  method** (baseline) and the **double-angle / Σ methods** (cross-checks —
  these are insensitive to QED FSR and electron energy scale in different ways;
  agreement between methods is a systematic handle).
- **Q² > 25 GeV²** baseline (suppresses resolved-photon/MPI contamination and
  keeps DGLAP evolution slow across the sample); analysis repeated in Q² bins.
- 0.05 < y < 0.95 (reconstruction quality; radiative tails).
- W within the per-configuration ranges of §4.1.

### 5.2 Jet selection

- anti-k_T, **R = 1.0 in the γ\*p frame** for the collapse analysis and
  **R = 0.4 in the lab** for the lab-frame splay analysis. (Clustering in the
  frame where the observable is defined avoids R-sculpting of the comparison;
  the lab-frame R = 0.4 analysis is what an experiment would actually
  calibrate, and is where the universality breakdown is the practical
  message.)
- **Current-hemisphere requirement:** jet axis in the current hemisphere of
  the **Breit frame** (p_z^Breit < 0). This removes target/remnant
  fragmentation, which is *not* the FFS signal and would contaminate the
  multiplicity observables. Quoted as a fiducial cut of the proposed
  measurement.
- Lab acceptance: |η_jet| < 3.5, jet p_T^lab > 4 GeV, |p|_lab ∈ [5, 30] GeV.
- Leading jet only (per-event statistical independence; simplifies the
  statistical model).

### 5.3 Observables

| Observable | Definition | Role |
|---|---|---|
| **n₉₀** | fractional # of constituents (momentum-ordered) carrying 90% of jet |p| — as implemented in `analyze_events.py::compute_n_x` | **primary** |
| n₇₅, n₉₅ | threshold variations | shape robustness |
| n₉₀^charged | charged constituents only, p_T > 0.2 GeV | what tracking actually measures; used for all smeared projections |
| N_charged | charged multiplicity | secondary, historically familiar |
| p_T D, LHA λ_{0.5}¹ | momentum dispersion, Les Houches angularity | cross-check that FFS appears coherently across the substructure basis |

All observables are computed from the same constituent list in both frames
(constituent four-vectors boosted with `DISKinematics.boost_to_gamma_p_cm`);
n₉₀ is invariant under the ordering-preserving part of the boost only if
computed per frame — we compute it **in the frame of the binning variable**.

---

## 6. Analysis strategy — controlling the confounders

This section is the referee-proofing. Each known effect that could fake or
dilute an FFS signal gets an explicit control:

| Confounder | Why it's dangerous | Control |
|---|---|---|
| **Q² evolution (DGLAP)** | multiplicity rises with ln Q²; Q² correlates with W | All results in **narrow Q² bins** ([25,50], [50,100], [100,250], [250,1000] GeV²). FFS predicts a strong W slope *within* a fixed Q² bin; DGLAP predicts none. Additionally reweight W slices to a common Q² spectrum within each bin. |
| **η/phase-space sculpting** | at fixed |p|_lab, different W populate different jet η; fragmentation observables have η-dependent acceptance | (i) Reweight W slices to a common jet-η distribution within each (|p|_lab, Q²) bin; (ii) repeat in restricted |η_jet| < 1 barrel fiducial volume; (iii) the collapse test (Fig. 3) is itself insensitive — sculpting cannot manufacture a collapse in p_CM. |
| **Spectrum sculpting within bins** | steeply falling |p| spectrum shifts the effective ⟨|p|⟩ per bin differently at different W | fine |p|_lab binning (2 GeV) + reweight to common in-bin spectrum. |
| **Target fragmentation** | remnant-side particles fake multiplicity trends with W | Breit-frame current-hemisphere cut (§5.2); vary hemisphere boundary as a systematic. |
| **Resolved photon / MPI at low Q²** | extra activity correlated with W | Q² > 25 GeV² baseline; MPI-on sample quantifies residual; repeat at Q² > 100 GeV². |
| **QED radiation** | distorts reconstructed W event-by-event → smears the binning variable | compare electron/Σ/double-angle reconstruction; radiative on/off generator pair; migration treated in smearing study. |
| **Hadronization model dependence** | "the effect is just a Pythia artifact" | Herwig + Sherpa envelope (§4.2); the Letter claims the *trend*, common to all models, and presents the spread as the model-discrimination opportunity. |

**The decisive internal cross-check** is the universality-restoration test:
none of the above confounders can cause curves that splay in |p|_lab binning
to *collapse* when re-binned in p_CM. The collapse is the fingerprint of the
boost mechanism specifically.

---

## 7. Detector realism (fast smearing)

A pheno-level PRL needs credible — not full-sim — detector effects:

- **Tracking** (drives n₉₀^charged): ePIC-like σ(p)/p ≈ 0.05% · p ⊕ 0.5%
  (barrel), degraded ×2 in 1 < |η| < 3.5; track efficiency 95% (p_T > 0.2 GeV,
  |η| < 3.5).
- **Scattered electron / kinematics:** EM calorimeter σ(E)/E ≈ 2%/√E ⊕ 1%,
  electron-method W resolution propagated event-by-event → W-bin migration
  matrix; the analysis is performed in reconstructed W with bin widths ≥ 2×
  the W resolution.
- **Implementation:** a `smear.py` module applying parametric smearing to the
  Parquet event record, producing a parallel smeared dataset run through the
  identical analysis. All Fig. 4 projections use smeared quantities;
  Figs. 2–3 show particle level with smeared results overlaid.

No unfolding is performed (appropriate for a projection paper); instead we
demonstrate that the smeared effect remains ≫ the smearing-induced
distortion.

---

## 8. Statistical methodology

1. **Per-bin estimator:** ⟨n₉₀⟩ with bootstrap uncertainties (event-level
   resampling, 1000 replicas) — handles the non-Gaussian n₉₀ distribution.
2. **Effect significance (Fig. 2):** in each (|p|_lab, Q²) bin, fit
   ⟨n₉₀⟩(W) with (a) a constant (H₀) and (b) the FFS-predicted shape
   f(p_CM(W)) with one normalization parameter (H₁); Δχ² → significance.
3. **Collapse metric (Fig. 3):** χ²/ndf of all W slices against the global
   universal curve, computed in both binnings; report the ratio
   χ²_lab/χ²_CM as the universality-restoration factor.
4. **Luminosity projection (Fig. 4):** scale bootstrap uncertainties to
   N_jets(L) per bin from the generator cross sections (L = 1, 10, 100 fb⁻¹
   per configuration); recompute Δχ²; include the systematics floor from §9
   as nuisance parameters with the stated priors. Report L needed for 5σ
   rejection of H₀ (expected: ≪ 1 fb⁻¹ statistically; the result will be the
   systematics-floor-limited significance, which is the honest and more
   interesting number).

---

## 9. Systematic uncertainties (projection budget)

| Source | Evaluation | Expected size on ⟨n₉₀⟩ ratio |
|---|---|---|
| Hadronization model | Pythia vs Herwig vs Sherpa envelope | dominant on absolute ⟨n₉₀⟩; ~few % on W-trend ratios |
| Tune / string parameters | Pythia tune variations | few % |
| Track efficiency | ±2% absolute | ~1–2% on n₉₀^charged, cancels in W-ratios |
| W reconstruction / QED radiation | method comparison + radiative variations | bin migration, ~few % |
| Jet energy scale | ±2% shift on constituents | small, cancels in ratios |
| Hemisphere definition | vary Breit-frame cut | small |
| Q²-evolution residual | closure on H₀ sample | quantifies method floor |

The headline observables are **ratios across W at fixed lab kinematics**
(and across beam configurations), chosen precisely because the experimental
systematics above cancel to first order, leaving the hadronization envelope
— which is the physics — as the dominant band.

---

## 10. Implementation roadmap (mapped to this repository)

**Phase 1 — Generation upgrades** (`generate_events.py`)
- [ ] Beam-configuration presets (`--config 5x41|10x100|18x275`).
- [ ] Raise default Q² floor to 25 GeV² (keep override); store generator
      cross section in Parquet metadata for luminosity scaling.
- [ ] Generator switch: add Herwig 7 and Sherpa runcards (Docker services),
      common HepMC→Parquet converter so `analyze_events.py` is
      generator-agnostic.
- [ ] Pythia tune/MPI/QED-radiation variation flags.
- [ ] Store scattered-electron four-vector per event (needed for smearing
      and reconstruction-method studies).

**Phase 2 — Analysis upgrades** (`analyze_events.py`, `utils/`)
- [ ] Boost constituents to γ\*p frame (`boost_to_gamma_p_cm`, already in
      `utils/dis_kinematics.py`); cluster R = 1.0 jets in that frame.
- [ ] Breit-frame current-hemisphere selection (`boost_to_breit` exists).
- [ ] Add observables: n₇₅/n₉₅, n₉₀^charged, p_T D, λ_{0.5}; histogram
      axes extended with Q² and p_CM.
- [ ] `utils/smear.py`: parametric smearing per §7, applied as an optional
      pass producing `events_smeared.parquet`.
- [ ] Leading-jet-only option; η/spectrum reweighting machinery.

**Phase 3 — Statistics & figures** (`make_plots.py` → split into
`make_results.py` + `make_plots.py`)
- [ ] Bootstrap profile machinery; H₀ construction; Δχ² and collapse metrics.
- [ ] The four PRL figures of §3, PRL column-width styling (mplhep).

**Compute:** 30M Pythia-class events ≈ O(100) CPU-hours — overnight on a
small batch allocation; Docker services already orchestrate the pipeline.

---

## 11. Risks and mitigations

| Risk | Likelihood | Mitigation |
|---|---|---|
| Effect smaller at EIC than in the ZZ benchmark (boosts are milder at low √s) | medium | the 18×275 configuration extends W to 130 GeV; cross-configuration combination maximizes the lever arm; if the smeared significance lands at 3σ instead of 5σ for early luminosity, the Letter's claim shifts to the 100 fb⁻¹ dataset — still a clean PRL message |
| Herwig/Sherpa show qualitatively different trend | low (cluster models also hadronize in the color frame) | this would itself be a publishable result: the measurement discriminates hadronization paradigms |
| W-resolution washout at high W (electron method degrades at low y) | medium | Σ/double-angle methods; restrict to y > 0.05; widen W bins at high W |
| Referee objects "this is just multiplicity vs W" | certain | Fig. 3 collapse + fixed-(|p|_lab, Q²) binning is the rebuttal built into the design; emphasized in the text |

---

## 12. One-paragraph abstract draft

> Hadronization occurs in the rest frame of color-connected systems, not the
> laboratory frame, implying that jets with identical laboratory kinematics
> carry different internal structure depending on event-level kinematics —
> the frame-dependent fragmentation shift. We show that neutral-current DIS
> at the Electron–Ion Collider provides the first environment where the
> color-connected rest frame is reconstructed exactly, event by event. Using
> Pythia 8, Herwig 7, and Sherpa simulations across three EIC beam
> configurations with parametric detector smearing, we find that the jet
> core-multiplicity observable ⟨n₉₀⟩ at fixed laboratory momentum varies by
> up to ~50% across the accessible W range at fixed Q², while collapsing
> onto a single universal curve when expressed in the photon–proton
> center-of-momentum frame. The effect is observable with high significance
> in early EIC data and provides a new, frame-resolved discriminator of
> hadronization models.

---

*Reference implementation baseline: this repository (Pythia 8 NC-DIS →
anti-k_T jets → n₉₀ histograms → figures). See `README.md` for the pipeline.*
