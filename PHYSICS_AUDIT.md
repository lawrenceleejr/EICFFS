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
2.4 M events (8 seeds × 300 k), anti-kT R = 0.4 in the lab, p_T > 2 GeV,
|η| < 3.5, jet in the Breit current hemisphere: 291 k jets.

**The effect is present and large.**  ⟨n₉₀⟩ at fixed lab-frame momentum rises
monotonically with W over most of the phase space (`figures/ffs_fan.pdf`,
`ffs_slopegraph.pdf`, `ffs_ratio.pdf`):

| \|p\|_lab (GeV) | W = 10–15 | W = 32–45 | change |
|---|---|---|---|
| 2–3   | 1.25 | 1.69 | +35 % |
| 4.5–7 | 1.58 | 2.44 | +54 % |
| 7–10  | 1.76 | 3.09 | +76 % |
| 10–15 | 1.99 | 3.57 | +79 % |
| 15–22 | 2.29 | 4.19 | +83 % |

These are shifts of the same order as the 50 % quoted in the reference paper
for 200 GeV jets, but for 5–20 GeV jets.

**Within a W slice the colour-frame momentum is nearly fixed.**  For current
jets |p|_cm / (W/2) has median 0.66 (16–84 %: 0.47–0.83), so a W slice is a
sample of jets with essentially one colour-frame momentum seen under different
boosts.  Plotted against the boost factor |p|_lab/|p|_cm
(`ffs_boost_factor.pdf`) each slice rises steeply from ≈1.5 when the jet recoils
against the frame's motion (backward jets, high y) to ≈3 when it is boosted
forward with it.  Plotted against |p|_cm itself (`ffs_collapse.pdf`, same
vertical scale as the fan) the slices agree with each other where they overlap
and are nearly flat, 1.8 → 2.2 from 1 to 20 GeV.  A fit of
log n₉₀ against log|p|_cm and log Q over all current jets gives exponents
0.03 and 0.33 respectively.

Interpretation: n₉₀ computed from the same constituents after boosting them
into the γ*p frame is indistinguishable from the lab value (correlation 0.993,
means 1.927 vs 1.922), so the ordering step of the observable is frame-stable.
The frame dependence enters through *which particles the fixed lab cone
contains*: a jet boosted forward is compressed into R = 0.4 and the cone
captures more of the string's current end; a jet recoiling against the frame's
motion opens up and the cone loses particles (`ncon_boost_factor.pdf` shows
⟨N_constituents⟩ following the same pattern).  This is the paper's point that a
jet is not a factorisable object, realised at the EIC in a regime where jets
have a median of two constituents and n₉₀ is a leading-hadron observable rather
than a shower-multiplicity one.

**W and Q both matter, and DIS separates them.**  In e⁺e⁻ the colour-frame
energy and the shower scale coincide; in DIS they are W and Q.  At fixed
|p|_lab = 4.5–10 GeV the (W, Q) table (`ffs_wq_table.pdf`) shows ⟨n₉₀⟩ rising
with Q at fixed W in every column, and rising then falling with W at fixed Q.
The turnover at W ≳ 40 GeV is the point where, at this beam energy, current jets
of fixed lab momentum are backward-going (y > 0.5) and de-boosted.  The
non-monotonic behaviour is visible in the full n₉₀ distribution as well
(`ffs_distribution.pdf`), not only in the mean.

**Where in the EIC phase space.**  `boost_map.pdf` gives the analytic rapidity
of the colour rest frame relative to the lab across the (x, Q²) plane with the
selected jets overlaid; at 10 × 100 GeV it runs from y_cm ≈ 1.3 at W = 55 GeV to
≈ 3 at W = 10 GeV, so lab jets of equal momentum are compared across frames
differing by a factor e^{1.7} ≈ 5 in longitudinal boost.  `plateau.pdf` shows
the charged-hadron rapidity plateau in the γ*p frame growing like ln W², the
string whose current end the lab jet samples.

### Caveats to carry into a paper

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
