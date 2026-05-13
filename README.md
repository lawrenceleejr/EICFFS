# EICFFS — EIC Frame-dependent Fragmentation Shift Study

A self-contained Docker framework for studying the **Frame-dependent Fragmentation
Shift (FFS) effect** at the Electron–Ion Collider (EIC), based on the paper:

**Reference:** [arXiv:2308.10951](https://arxiv.org/abs/2308.10951)  
*Phys.Lett.B 866, 2025, 139561 — Lee et al., University of Tennessee Knoxville*

---

## Physics Overview

### The FFS Effect

In QCD, jet fragmentation occurs in the **colour rest frame** — the rest frame
of all colour-connected particles — not the lab frame.  When a jet has fixed
lab-frame momentum |p|_lab, the boost between the lab frame and the colour rest
frame depends on the event kinematics.  As a result, two jets with the *same*
lab-frame momentum can have very different internal structure if they come from
different kinematic configurations (different colour rest frames).

This **Frame-dependent Fragmentation Shift (FFS)** is illustrated in
arXiv:2308.10951 using e⁺e⁻ → 3j (colour rest frame = lab frame) vs.
e⁺e⁻ → ZZ → 4j (colour rest frame = boosted Z rest frame).  For the ZZ
process, ⟨n₉₀⟩ is set by m_Z/2 ≈ 45 GeV *regardless* of the lab-frame
momentum, leading to 50% differences from the 3j case at high momentum.

### Primary observable: n₉₀

Following arXiv:2308.10951 Sec. 2, the main observable is **n_x** — the
fractional minimum number of jet constituents (ordered by decreasing 3-momentum
magnitude) needed to recover *x*% of the total jet momentum:

1. Sort constituents by decreasing |p|  
2. Compute cumulative momentum fraction  
3. Interpolate to find fractional n needed to reach the threshold  

We use **n₉₀** (threshold = 90%), which is IRC-safe under momentum ordering
and strongly discriminates between jet populations from different colour topologies.
The simple charged-particle count N_charged is also stored as a secondary observable.

### The EIC analogy

In NC-DIS (*ep → eX*) the relevant colour rest frame is the **photon-proton
CM frame** (Breit/γ\*p frame), characterised by the invariant mass *W*:

| Symbol | Definition | Meaning |
|--------|-----------|---------|
| Q²   | −q²       | Photon virtuality |
| *x*   | Q²/(2P·q) | Bjorken-x (struck parton momentum fraction) |
| *y*   | P·q / P·k | Inelasticity |
| *W*   | √((P+q)²) | Photon-proton CM energy = colour rest frame energy |

For **fixed lab-frame jet |p|**, the boost to the colour rest frame changes
with *W*.  Therefore the same lab-frame jet corresponds to a *different*
CM-frame momentum at different *W*, → the fragmentation (and hence n₉₀)
varies with *W* at fixed |p|_lab.  This is the FFS effect at the EIC.

---

## Repository Structure

```
EICFFS/
├── Dockerfile               # Container image definition
├── docker-compose.yml       # Multi-service orchestration
├── environment.yml          # Conda/mamba environment spec
├── run.sh                   # End-to-end pipeline script
│
├── generate_events.py       # Step 1 — Pythia8 NC-DIS event generation
├── analyze_events.py        # Step 2 — FFS analysis → ROOT histograms
├── make_plots.py            # Step 3 — Publication-quality matplotlib figures
│
├── utils/
│   ├── __init__.py
│   └── dis_kinematics.py    # 4-vector math, DIS invariants, frame boosts
│
├── data/                    # Generated events & histograms (gitignored)
└── plots/                   # Output figures (gitignored)
```

---

## Quick Start (Docker — recommended)

### Prerequisites
- [Docker](https://docs.docker.com/get-docker/) ≥ 24
- [Docker Compose](https://docs.docker.com/compose/) ≥ 2

### Full pipeline (generate → analyse → plot)

```bash
# Build the image (first time, ~5–10 min for package installation)
docker compose build

# Run the complete pipeline (creates data/ and plots/ on the host)
docker compose run eicffs

# Or run steps individually:
docker compose run generate   # generate 200k events → data/events.parquet
docker compose run analyze    # fill histograms     → data/histograms.root
docker compose run plot       # make figures        → plots/*.pdf
```

### Customise via environment variables

```bash
N_EVENTS=50000 docker compose run generate   # quick test with 50k events
```

---

## Quick Start (local Python environment)

### Install (conda/mamba recommended)

```bash
mamba env create -f environment.yml
conda activate eicffs
```

### Run

```bash
# Full pipeline
./run.sh

# Step by step
python generate_events.py --n-events 200000 --output data/events.parquet --quiet
python analyze_events.py  data/events.parquet --output data/histograms.root
python make_plots.py      data/histograms.root --outdir plots/
```

---

## Scripts

### `generate_events.py`

Generates NC-DIS events at EIC kinematics using **Pythia8** and saves
particle-level data (all final-state hadrons per event) to an **Apache
Parquet** file via `awkward-array`.

Key options:

| Flag | Default | Description |
|------|---------|-------------|
| `--n-events` | 200 000 | Events to save |
| `--electron-energy` | 10 GeV | Electron beam energy |
| `--proton-energy` | 100 GeV | Proton beam energy |
| `--Q2min` / `--Q2max` | 1 / 1000 GeV² | Phase-space cuts |
| `--output` | `data/events.parquet` | Output path |
| `--quiet` | off | Suppress Pythia8 output |

### `analyze_events.py`

Reads the Parquet event file, reconstructs DIS kinematics, finds jets
(anti-*k*_T, **R = 0.4** via **FastJet** or built-in fallback), and fills
the following histograms, written to a **ROOT** file via `uproot`:

| Histogram | Axes | Observable |
|-----------|------|-----------|
| `n90_3d` | W × \|p\|_lab × n₉₀ | **Primary FFS 3D histogram** |
| `sum_n90_vs_W`, `count_n90_vs_W` | W × \|p\|_lab | ⟨n₉₀⟩ profile (sum + count) |
| `mult_3d` | W × \|p\|_lab × N_charged | Secondary FFS histogram |
| `sum_N_vs_W`, `count_vs_W` | W × \|p\|_lab | ⟨N_charged⟩ profile (sum + count) |
| `Q2`, `x`, `y`, `W` | 1D | DIS kinematic distributions |
| `Q2_vs_W` | W × Q² | Kinematic plane |
| `jet_eta_pt` | η × pT | Jet landscape |
| `n_jets` | N_jets | Jet multiplicity per event |

### `make_plots.py`

Reads ROOT histograms and produces five publication-quality figures:

| File | Content |
|------|---------|
| `ffs_main.pdf` | **Primary result**: ⟨n₉₀⟩ vs *W* for each \|p\|_lab bin |
| `ffs_ratio.pdf` | ⟨n₉₀⟩ ratio to lowest-|p|_lab bin — FFS magnitude |
| `kinematics.pdf` | DIS kinematic plane (Q² vs W) + marginal distributions |
| `jet_landscape.pdf` | Jet η–pT heat-map + jet multiplicity per event |
| `ffs_heatmap.pdf` | 2D colour map of ⟨n₉₀⟩(*W*, \|p\|_lab) |

---

## Software Stack

| Package | Role |
|---------|------|
| [Pythia8](https://pythia.org) | Hard-process MC event generation |
| [FastJet](https://fastjet.fr) | Anti-*k*_T jet finding (R=0.4) |
| [uproot](https://uproot.readthedocs.io) | ROOT I/O (no ROOT dependency) |
| [awkward-array](https://awkward-array.org) | Ragged array handling |
| [hist](https://hist.readthedocs.io) | Histogram filling & manipulation |
| [matplotlib](https://matplotlib.org) | Plotting |
| [mplhep](https://mplhep.readthedocs.io) | HEP plot styling |
| [scipy](https://scipy.org) | Gaussian smoothing for heat maps |

