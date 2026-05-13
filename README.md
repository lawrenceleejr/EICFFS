# EICFFS — EIC Frame-dependent Fragmentation Shift Study

A self-contained Docker framework for studying the **Frame-dependent Fragmentation
Shift (FFS) effect** at the Electron–Ion Collider (EIC).  The FFS effect
describes how jet charged-particle multiplicity, observed at a fixed lab-frame
momentum, varies with the photon-proton invariant mass *W* — a direct
consequence of the Lorentz boost between the lab frame and the colour rest
frame changing as a function of *W*.

**Reference:** [arXiv:2308.10951](https://arxiv.org/abs/2308.10951)

---

## Physics Overview

In neutral-current DIS (*ep → eX*) the virtual photon carries 4-momentum
*q = k − k′*.  The key invariants are:

| Symbol | Definition | Meaning |
|--------|-----------|---------|
| Q²   | −q²       | Photon virtuality |
| *x*   | Q²/(2P·q) | Bjorken-x (struck parton momentum fraction) |
| *y*   | P·q / P·k | Inelasticity |
| *W*   | √((P+q)²) | Photon-proton CM energy |

For **fixed lab-frame jet momentum** |p|_lab, the boost from the lab to
the photon-proton colour rest frame (Breit / γ*p CM frame) depends on *W*.
Therefore the same lab-frame object probes a *different* part of the
fragmentation function at different *W* values — this is the FFS effect,
visible as a rise in ⟨*N*_charged⟩ per jet with increasing *W* at fixed |p|_lab.

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
(anti-*k*_T, *R* = 0.8 via **FastJet** or built-in fallback), and fills
the following histograms, written to a **ROOT** file via `uproot`:

| Histogram | Axes | Observable |
|-----------|------|-----------|
| `mult_3d` | W × \|p\|_lab × N_charged | 3D FFS histogram |
| `mean_mult_vs_W` | W × \|p\|_lab (Mean storage) | ⟨N_charged⟩ profile |
| `Q2`, `x`, `y`, `W` | 1D | DIS kinematic distributions |
| `Q2_vs_W` | W × Q² | Kinematic plane |
| `jet_eta_pt` | η × pT | Jet landscape |
| `n_jets` | N_jets | Jet multiplicity per event |

### `make_plots.py`

Reads ROOT histograms and produces five publication-quality figures:

| File | Content |
|------|---------|
| `ffs_main.pdf` | **Primary result**: ⟨N_charged⟩ vs *W* for each \|p\|_lab bin |
| `ffs_ratio.pdf` | Ratio to lowest-|p|_lab bin — FFS magnitude |
| `kinematics.pdf` | DIS kinematic plane (Q² vs W) + marginal distributions |
| `jet_landscape.pdf` | Jet η–pT heat-map + jet multiplicity per event |
| `ffs_heatmap.pdf` | 2D colour map of ⟨N_charged⟩(*W*, \|p\|_lab) |

---

## Software Stack

| Package | Role |
|---------|------|
| [Pythia8](https://pythia.org) | Hard-process MC event generation |
| [FastJet](https://fastjet.fr) | Anti-*k*_T jet finding |
| [uproot](https://uproot.readthedocs.io) | ROOT I/O (no ROOT dependency) |
| [awkward-array](https://awkward-array.org) | Ragged array handling |
| [hist](https://hist.readthedocs.io) | Histogram filling & manipulation |
| [matplotlib](https://matplotlib.org) | Plotting |
| [mplhep](https://mplhep.readthedocs.io) | HEP plot styling |
| [scipy](https://scipy.org) | Gaussian smoothing for heat maps |
