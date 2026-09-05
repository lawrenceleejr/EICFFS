# Physics audit of the EICFFS framework

*Audit of the code as merged in PR #1, the fixes applied on this branch, what the
corrected simulation shows, and whether this study has been proposed for the EIC
before.  Reference: L. Lee, C. Bell, J. Lawless, C. Nash, E. Nibigira,
"Experimental impact of jet fragmentation reference frames at particle colliders",
Phys. Lett. B 866 (2025) 139561, [arXiv:2308.10951](https://arxiv.org/abs/2308.10951).*

---

## 1. Summary

The framework's premise is sound: in neutral-current DIS the colour-connected
system is the whole hadronic final state, with four-momentum P + q and invariant
mass W, so its rest frame (the γ*p or hadronic centre-of-mass frame) is the
"colour rest frame" of the reference paper, and W is the DIS analogue of m_Z in
the paper's ZZ → 4j example.  The original implementation, however, could not
have measured the effect:

| # | Severity | Finding | Status |
|---|----------|---------|--------|
| 1 | **Blocking** | Hard-process lepton looked up with `status == 23` in `pythia.event`; the shower copies it, so the record carries `−23` in ~96 % of events. Those events were silently dropped and the surviving 4 % (no lepton FSR/recoil) are a biased subset. | Fixed: read `pythia.process`. 100 % of events now pass. |
| 2 | **Major** | Lepton-beam ISR left on (Pythia default `PDF:lepton = on`). The stored q = k_beam − k′ then includes the ISR photon and is not the exchanged boson; W, Q², x, y are all mis-reconstructed for radiative events. | Fixed: `PDF:lepton = off` (standard for EIC studies); `--lepton-isr` restores it. |
| 3 | **Major** | The Breit frame and the γ*p CM frame were conflated in the README and docstrings ("Breit/γ*p frame"). They differ by a boost along the boson axis; the struck quark carries Q/2 in the Breit frame but W/2 in the γ*p frame. Only the latter is the colour rest frame. | Fixed in code and text; both boosts are now computed and stored per jet. |
| 4 | **Major** | Sign of the predicted effect stated backwards in `analyze_events.py` ("higher W → lower CM-frame momentum → fewer particles"). The colour rest frame moves forward with rapidity y_cm that *decreases* with W, so at fixed lab momentum a higher-W jet is *harder* in that frame and has *more* particles. | Fixed; confirmed by the simulation (Sec. 3). |
| 5 | **Major** | All jets were used, including proton-remnant fragments at forward rapidity, which have nothing to do with the struck quark's fragmentation. | Fixed: jets are flagged by Breit-frame hemisphere; the figures use current-hemisphere jets. |
| 6 | Moderate | FastJet was installed but never used by default: `run.sh` and `docker-compose.yml` never passed `--use-fastjet`, so the greedy cone fallback was the production algorithm. | Fixed: anti-kT via FastJet is the default (`--no-fastjet` to opt out). |
| 7 | Moderate | No `SpaceShower:dipoleRecoil = on`, the Pythia-recommended ISR recoil scheme for DIS; an arbitrary `SpaceShower:rapidityOrder = off` was set instead. | Fixed: dipole recoil on, rapidity ordering back to default. |
| 8 | Moderate | Uncertainty on ⟨n₉₀⟩ taken as √(mean/N), a Poisson assumption that does not hold for a fractional, interpolated observable. | Fixed: standard error of the mean from the per-bin sample. |
| 9 | Moderate | Only Q², W, x, y were stored per event, not q or k′, so no frame boost and no current-jet identification was possible downstream. | Fixed: q, k′, beam proton and struck parton four-vectors stored. |
| 10 | Minor | All final-state e, μ removed from the hadronic final state, rather than the scattered lepton and its QED FSR photons. Removes Dalitz electrons; a collinear FSR photon could have formed a spurious "jet". | Fixed: removal by ancestry of the hard-process lepton; neutrinos removed. |
| 12 | **Major** | ⟨n₉₀⟩ plotted against the jet's *total* lab momentum \|p\|_lab conflates momentum with angle: at fixed \|p\|, a higher-*W* jet is more central and therefore harder in *p*_T, and *p*_T is what sets how much radiation a fixed-*R* lab cone collects. Against *p*_T^lab the *W* slices differ by only a few percent (Sec. 3). | Fixed: lab-frame figures now use *p*_T^lab, and the frame test below replaces the \|p\|_lab fan as the primary result. |
| 11 | Minor | The 200k-event default sample has median W ≈ 8 GeV; events with W > 30 GeV and a current jet are ~1 % of it, so the high-W bins that carry the effect would have been empty. | Fixed: generation-time `--Wmin`, eight parallel seeds. |

Items already correct and verified: the n_x definition and its interpolation
match arXiv:2308.10951 Sec. 2 (vectorised version agrees with the scalar one to
1e-11); the Lorentz-boost formula; the Breit-frame construction (q + 2xP at
rest, boson energy exactly zero, |q| = Q); the DIS invariants (W reproduced by
the stored hadronic final state to 0.03 GeV on average, the residual being
neutrinos from heavy-flavour decays).

---

## 2. Details of the blocking finding

A 300-event test with the original generator settings:

```
electron status counts: {-12: 300, -21: 300, -23: 288, 23: 12}
```

In `pythia.event` the hard-process lepton has status −23 whenever the shower
makes a copy of it (QED FSR or recoil), i.e. in 288 of 300 events; the original
`extract_kinematics` returned `None` for those and the event loop `continue`d.
The reported "efficiency" would have been ~4 %.  Reading the same particle from
`pythia.process` (status +23 always) gives 100 % of events and, with
`PDF:lepton = off`, an exact q = k − k′.

---

## 3. What the corrected simulation shows

Sample: Pythia 8.317, e(10 GeV) p(100 GeV), NC DIS Q² > 1 GeV², W > 10 GeV,
2.4 M events (8 seeds × 300 k).  Two jet collections are built from the same
events:

* **Lab jets** — anti-kT R = 0.4 in the laboratory, p_T > 2 GeV, |η| < 3.5,
  Breit current hemisphere: 291 k jets.
* **Colour-frame jets** — every final-state particle boosted into the γ*p
  frame and clustered there with FastJet's e⁺e⁻ generalised-kT (p = −1) at an
  angular radius R = 0.4 rad, E > 1 GeV, current hemisphere: 4.05 M jets.
  An angular algorithm is required because in that frame the struck quark lies
  on the boson axis, where η–φ clustering is singular.

### The test: which variable makes the curves flat?

If fragmentation is set in the colour rest frame, then labelling jets by their
energy *in that frame* should make ⟨n₉₀⟩ independent of how hard the lab sees
them — one flat line per label, the DIS analogue of the flat ⟨n₉₀⟩ of the ZZ
jets in arXiv:2308.10951.  That is what happens, and only for the colour-frame
jets (`figures/flat_cmjets.pdf`):

| E_cm (GeV) | p_T^lab ≈ 0.8 | 1.7 | 3.1 | 5.2 | 9.2 | spread |
|---|---|---|---|---|---|---|
| 4–6 | 2.30 | 2.33 | 2.37 | 2.40 | 2.37 | 3.9 % |
| 6–9 | 2.79 | 2.82 | 2.89 | 2.93 | 2.92 | 5.3 % |
| 9–14 | 3.40 | 3.42 | 3.50 | 3.56 | 3.64 | 6.7 % |
| 14–22 | 4.17 | 4.12 | 4.15 | 4.24 | 4.24 | 2.8 % |

Across the six E_cm slices the residual variation over a factor of ten in lab
transverse momentum is 2.8 % to 12.1 %, the largest values belonging to the
lowest-energy slice where the jet is often a single hadron.  The same jets
plotted against E_cm with lab-p_T slices overlaid fall on one curve
(`universal_cm.pdf`).

The identical test on lab-clustered jets fails (`flat_labjets.pdf`).  Labelled
by the same colour-frame momentum, the slices collapse onto each other and rise
together with p_T^lab, spanning 41 % to 85 %:

| \|p\|_cm (GeV) | p_T^lab ≈ 2.2 | 3.6 | 5.7 | 9.0 | 11.4 | spread |
|---|---|---|---|---|---|---|
| 3.3–5 | 1.59 | 2.06 | 2.62 | 3.24 | — | 70 % |
| 7.5–11 | 1.56 | 2.01 | 2.62 | 3.26 | 3.58 | 81 % |
| 16–24 | 1.51 | 1.95 | 2.54 | 3.20 | 3.56 | 84 % |

A lab-frame jet's measured fragmentation is therefore controlled by its lab
transverse momentum and carries almost no memory of the colour rest frame,
while the intrinsic fragmentation of the same events is controlled by the
colour-frame energy and carries almost no memory of the lab.  This is the
frame-dependent fragmentation shift at the EIC stated as a measurement rather
than as an analogy.

### Where the effect lives: clustering, not counting

Computing n₉₀ from the *same* constituents but ordering them by lab momentum
instead of colour-frame momentum changes almost nothing (E_cm = 2.5–4 GeV:
1.87 → 1.80; E_cm = 9–14 GeV: 3.40 → 3.55).  The observable's ordering step is
essentially frame-stable.  What differs between the two jet collections is
*which particles end up in the jet*: a fixed lab cone gathers a boost-dependent
slice of the colour string, so the lab jet is a different object, not the same
object measured differently.

### The W dependence at fixed lab momentum

Against total \|p\|_lab the W slices fan out by 35–83 % (`ffs_fan.pdf`): at
\|p\|_lab = 7–10 GeV, ⟨n₉₀⟩ goes from 1.76 at W = 10–15 GeV to 3.09 at
W = 32–45 GeV.  Most of that is angle rather than frame.  At fixed p_T^lab the
same slices differ by a few percent (`pt_fan.pdf`: 1.56 against 1.61 at
p_T = 2–2.5 GeV).  A jet of fixed \|p\| in a higher-W event sits at smaller
rapidity and therefore larger p_T, and p_T determines how much radiation an
R = 0.4 lab cone collects.  The \|p\|_lab fan is a true statement about jets of
equal lab momentum and is the direct transcription of the paper's e⁺e⁻
comparison, but at the EIC it should be shown next to the p_T version so the
angular part is visible rather than hidden.  The colour-frame test above is the
cleaner claim.

### Supporting observations

* Within a W slice the current jet's colour-frame momentum is nearly fixed:
  \|p\|_cm/(W/2) has median 0.66 (16–84 %: 0.47–0.83).
* ⟨n₉₀⟩ of lab jets rises with the boost factor \|p\|_lab/\|p\|_cm from ≈1.5
  to ≈3 (`ffs_boost_factor.pdf`), with ⟨N_constituents⟩ following.
* The colour rest frame's rapidity relative to the lab runs from ≈3 at
  W = 10 GeV to ≈1.3 at W = 55 GeV (`boost_map.pdf`); the charged-hadron
  rapidity plateau in that frame grows like ln W² (`plateau.pdf`).
* Slicing colour-frame jets by W rather than by E_cm does *not* flatten them.
  W bounds the available energy; the jet's own colour-frame energy is what sets
  its fragmentation.

### Caveats to carry into a paper

* The angular radius used in the colour rest frame (0.4 rad) is a choice; the
  flatness holds at 0.8 and 1.0 rad too, with larger ⟨n₉₀⟩ throughout.
* Colour-frame jets are not directly measurable without the scattered lepton,
  which fixes the boost.  The EIC detectors provide it, but the resolution on
  that boost propagates into E_cm and has not been studied here.
* Pythia only; the paper compared Pythia, Vincia and Herwig.  A Herwig 7 or
  Sherpa cross-check of the fan would strengthen the claim.  The framework's
  Parquet interface makes this a generator swap.
* Q² > 1 GeV² includes a region where the "current jet" with p_T > 2 GeV comes
  from O(α_s) boson–gluon fusion and QCD Compton; the colour rest frame is
  still the γ*p frame (everything is one colour-singlet system), but the jet
  is not the LO struck quark.  Cutting at Q² > 5 GeV² or requiring a parton
  match (`dR_parton` is stored) isolates the Born-like sample.
* Backward jets (η < −1) at high W sit in the electron-endcap region; a
  detector-level study needs the ePIC tracking and calorimeter acceptance there.
* The jet radius sets the size of the cone effect; repeating with R = 0.8 or
  with Breit-frame (Centauro) clustering would separate "cone" from
  "fragmentation" contributions cleanly.
* Only e(10) × p(100) was simulated; 18 × 275 GeV extends W to ≈ 140 GeV and
  moves y_cm.

---

## 4. Has this been proposed for the EIC before?

Checked on 2026-09-04:

* **INSPIRE citations of arXiv:2308.10951:** 3 (a heavy-ion jet
  background-subtraction paper, a top-jet classification paper, and a
  quark/gluon tagging paper on CMS open data).  None mentions DIS, ep, the
  EIC, HERA or the Breit frame.
* **arXiv full-text and abstract searches** for "colour/color rest frame",
  "frame-dependent fragmentation" with jets, and EIC + jet + fragmentation +
  frame return only the reference paper and unrelated hits.  The EIC Yellow
  Report ([arXiv:2103.05419](https://arxiv.org/abs/2103.05419)), the EIC jet
  overview by Page, Chu and Aschenauer
  ([arXiv:1911.00657](https://arxiv.org/abs/1911.00657)) and the NC-DIS jet
  substructure study ([arXiv:2302.06941](https://arxiv.org/abs/2302.06941))
  discuss lab-frame versus Breit-frame *jet finding* and yields, but none
  studies a fragmentation observable at fixed lab momentum as a function of W,
  nor frames anything in terms of the colour rest frame.
* **Closest prior art is from HERA, and it is complementary rather than
  overlapping.**  ZEUS ([arXiv:0803.3878](https://arxiv.org/abs/0803.3878))
  and H1 ([hep-ex/9707005](https://arxiv.org/abs/hep-ex/9707005)) measured
  charged multiplicities in the Breit current hemisphere as a function of Q and
  in the γ*p current region as a function of W, and tested universality against
  e⁺e⁻ at √s = Q or W.  Those are *frame-corrected* measurements designed to
  remove the boost.  The present study asks the opposite question, the one
  posed by arXiv:2308.10951: what a fixed-momentum *lab-frame* jet looks like
  when its colour rest frame is boosted differently.  Breit-frame jet
  algorithms for the EIC (Centauro,
  [arXiv:2006.10751](https://arxiv.org/abs/2006.10751)) are likewise a way to
  undo the boost, not to measure its imprint.

Conclusion: no EIC (or HERA) proposal for a frame-dependent fragmentation shift
at fixed lab-frame jet momentum was found.  The idea appears to be new, and the
EIC adds something the e⁺e⁻ examples in the paper cannot: it varies the
colour-frame boost continuously over a factor of ≈5 within one dataset, and it
decouples the colour-frame energy (W) from the hard scale (Q).
