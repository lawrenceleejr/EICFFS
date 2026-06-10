#!/usr/bin/env python3
"""
EIC DIS Event Generator for the FFS (Frame-dependent Fragmentation Shift) study.

Generates Pythia8 neutral-current DIS events at EIC kinematics and saves
particle-level data to a Parquet file for downstream analysis.

Implements the simulation campaign of ANALYSIS_DESIGN.md Sec. 4:
  * three EIC beam configurations (5x41, 10x100, 18x275),
  * hadronization variations (Lund string parameter envelope),
  * recommended Pythia 8.3 DIS settings (dipole recoil, no lepton PDF),
  * scattered-electron four-vector and generator cross section stored
    for the detector-smearing and luminosity-projection stages.

Usage
-----
    python generate_events.py --config 10x100 --n-events 300000 --quiet
    python generate_events.py --config 18x275 --variation lund-soft

Reference: arXiv:2308.10951
"""

import argparse
import json
import sys
import os
import numpy as np
import awkward as ak

sys.path.insert(0, os.path.dirname(__file__))
from utils.dis_kinematics import DISKinematics


# ---------------------------------------------------------------------------
# EIC beam configurations (ePIC reference points)
# ---------------------------------------------------------------------------

BEAM_CONFIGS = {
    "5x41":   (5.0, 41.0),
    "10x100": (10.0, 100.0),
    "18x275": (18.0, 275.0),
}

# Hadronization variations: Lund string parameter envelope
# (proxy for a multi-generator envelope; see ANALYSIS_DESIGN.md Sec. 4.2)
VARIATIONS = {
    "baseline": [],
    "lund-soft": [          # more, softer hadrons
        "StringZ:aLund = 0.95",
        "StringPT:sigma = 0.37",
    ],
    "lund-hard": [          # fewer, harder hadrons
        "StringZ:aLund = 0.45",
        "StringPT:sigma = 0.30",
    ],
}


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Generate EIC DIS events with Pythia8 for FFS effect study",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--config", type=str, default="10x100",
                   choices=sorted(BEAM_CONFIGS),
                   help="EIC beam configuration (Ee x Ep, GeV)")
    p.add_argument("--variation", type=str, default="baseline",
                   choices=sorted(VARIATIONS),
                   help="Hadronization variation")
    p.add_argument("--n-events", type=int, default=200_000,
                   help="Number of events to save")
    p.add_argument("--Q2min", type=float, default=25.0,
                   help="Minimum Q² (GeV²)")
    p.add_argument("--Q2max", type=float, default=1000.0,
                   help="Maximum Q² (GeV²)")
    p.add_argument("--mpi", action="store_true",
                   help="Enable multiparton interactions (contamination check)")
    p.add_argument("--output", type=str, default=None,
                   help="Output Parquet path "
                        "(default data/events_<config>_<variation>.parquet)")
    p.add_argument("--seed", type=int, default=42,
                   help="Pythia8 random seed")
    p.add_argument("--quiet", action="store_true",
                   help="Suppress Pythia8 banner and progress output")
    return p.parse_args()


# ---------------------------------------------------------------------------
# Pythia8 initialisation
# ---------------------------------------------------------------------------

def import_pythia8():
    """Import Pythia8 python bindings (conda 'pythia8' or PyPI 'pythia8mc')."""
    try:
        import pythia8
        return pythia8
    except ImportError:
        pass
    try:
        import pythia8mc
        return pythia8mc
    except ImportError:
        sys.exit(
            "pythia8 is not installed.\n"
            "Install with:  pip install pythia8mc\n"
            "          or:  conda install -c conda-forge pythia8\n"
        )


def init_pythia(args):
    """Initialise Pythia8 for neutral-current DIS at EIC kinematics."""
    pythia8 = import_pythia8()
    pythia = pythia8.Pythia("", False) if args.quiet else pythia8.Pythia()

    def cfg(s):
        pythia.readString(s)

    if args.quiet:
        cfg("Print:quiet = on")
    cfg("Next:numberCount = 50000")

    # ── EIC beam layout ─────────────────────────────────────────────────────
    e_e, e_p = BEAM_CONFIGS[args.config]
    cfg("Beams:frameType = 2")            # fixed beam energies
    cfg("Beams:idA = 2212")               # proton, +z
    cfg("Beams:idB = 11")                 # electron, -z
    cfg(f"Beams:eA = {e_p:.4f}")
    cfg(f"Beams:eB = {e_e:.4f}")

    # ── Neutral-current DIS, recommended Pythia 8.3 settings ────────────────
    cfg("WeakBosonExchange:ff2ff(t:gmZ) = on")
    cfg("SpaceShower:dipoleRecoil = on")  # recommended for DIS (Pythia main343)
    cfg("PDF:lepton = off")               # monochromatic lepton beam
    cfg("TimeShower:QEDshowerByL = off")  # no FSR off the scattered lepton:
                                          # true kinematics from e' (radiative
                                          # corrections deferred, see design)

    # ── Phase-space cuts ────────────────────────────────────────────────────
    cfg(f"PhaseSpace:Q2min = {args.Q2min}")
    cfg(f"PhaseSpace:Q2max = {args.Q2max}")

    # ── PDF ─────────────────────────────────────────────────────────────────
    cfg("PDF:pSet = 13")                  # NNPDF2.3 QCD+QED LO

    # ── Shower / MPI ────────────────────────────────────────────────────────
    cfg(f"PartonLevel:MPI = {'on' if args.mpi else 'off'}")

    # ── Hadronization variation ─────────────────────────────────────────────
    for s in VARIATIONS[args.variation]:
        cfg(s)
    cfg("HadronLevel:all = on")

    # ── Random seed ─────────────────────────────────────────────────────────
    cfg("Random:setSeed = on")
    cfg(f"Random:seed = {args.seed}")

    if not pythia.init():
        sys.exit("Pythia8 initialisation failed.")
    return pythia


# ---------------------------------------------------------------------------
# DIS kinematics reconstruction from Pythia8 event record
# ---------------------------------------------------------------------------

def extract_kinematics(event, e_e, e_p):
    """
    Reconstruct true DIS kinematics from the Pythia8 event record.

    The scattered lepton is the final-state lepton descending from the
    hard process (status 23 hard-process copy traced to its final copy).
    Beam four-vectors are taken as nominal (monochromatic beams).
    """
    # Hard-process outgoing lepton (|status| = 23; the record entry may be
    # an intermediate copy), then walk to its final copy via iBotCopyId.
    k_out = None
    for j in range(event.size()):
        p = event[j]
        if abs(p.status()) == 23 and abs(p.id()) == 11:
            q = event[p.iBotCopyId()]
            k_out = (q.px(), q.py(), q.pz(), q.e())
            break
    if k_out is None:
        return None, None

    m_e = 0.000511
    m_p = 0.938272
    k_in = (0.0, 0.0, -np.sqrt(max(e_e**2 - m_e**2, 0.0)), e_e)
    P_in = (0.0, 0.0,  np.sqrt(max(e_p**2 - m_p**2, 0.0)), e_p)
    return DISKinematics(k_in, k_out, P_in), k_out


# ---------------------------------------------------------------------------
# Main generation loop
# ---------------------------------------------------------------------------

def generate_and_save(args):
    """Run the event loop and write output to Parquet."""
    pythia = init_pythia(args)
    e_e, e_p = BEAM_CONFIGS[args.config]

    if args.output is None:
        tag = f"{args.config}_{args.variation}" + ("_mpi" if args.mpi else "")
        args.output = f"data/events_{tag}.parquet"
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)

    # Event-level scalars
    ev_Q2, ev_W, ev_x, ev_y = [], [], [], []
    ev_kout = []                       # scattered electron lab 4-vector

    # Particle-level ragged arrays (one sub-list per event)
    par_px, par_py, par_pz, par_e, par_charge = [], [], [], [], []

    n_saved = 0
    n_tried = 0

    while n_saved < args.n_events:
        if not pythia.next():
            continue
        n_tried += 1

        kin, k_out = extract_kinematics(pythia.event, e_e, e_p)
        if kin is None or not kin.valid:
            continue

        # Collect final-state hadrons (exclude all leptons/neutrinos,
        # hence also the scattered electron)
        fpx, fpy, fpz, fe, fch = [], [], [], [], []
        for p in pythia.event:
            if not p.isFinal():
                continue
            pid = abs(p.id())
            if pid in {11, 12, 13, 14, 15, 16}:
                continue
            fpx.append(p.px())
            fpy.append(p.py())
            fpz.append(p.pz())
            fe.append(p.e())
            fch.append(int(round(p.charge())))

        ev_Q2.append(kin.Q2)
        ev_W.append(kin.W)
        ev_x.append(kin.x)
        ev_y.append(kin.y)
        ev_kout.append(list(k_out))

        par_px.append(fpx)
        par_py.append(fpy)
        par_pz.append(fpz)
        par_e.append(fe)
        par_charge.append(fch)

        n_saved += 1
        if n_saved % 50_000 == 0 and not args.quiet:
            print(f"  Saved {n_saved:>7d} / {args.n_events} events "
                  f"  (efficiency {n_saved/n_tried*100:.1f}%)", flush=True)

    pythia.stat()
    info = pythia.infoPython()
    sigma_gen_mb = info.sigmaGen()      # mb

    # ── Build awkward array ────────────────────────────────────────────────
    print("\nBuilding awkward array …", flush=True)
    kout = np.array(ev_kout, dtype=np.float32)
    array = ak.Array({
        "Q2": np.array(ev_Q2, dtype=np.float32),
        "W":  np.array(ev_W,  dtype=np.float32),
        "x":  np.array(ev_x,  dtype=np.float32),
        "y":  np.array(ev_y,  dtype=np.float32),
        "e_px": np.ascontiguousarray(kout[:, 0]),
        "e_py": np.ascontiguousarray(kout[:, 1]),
        "e_pz": np.ascontiguousarray(kout[:, 2]),
        "e_e":  np.ascontiguousarray(kout[:, 3]),
        "particles": ak.zip({
            "px":     ak.Array(par_px),
            "py":     ak.Array(par_py),
            "pz":     ak.Array(par_pz),
            "e":      ak.Array(par_e),
            "charge": ak.Array(par_charge),
        }),
    })

    meta = {
        "config": args.config,
        "variation": args.variation,
        "e_e": e_e, "e_p": e_p,
        "Q2min": args.Q2min, "Q2max": args.Q2max,
        "mpi": args.mpi,
        "n_events": n_saved,
        "n_tried": n_tried,
        "sigma_gen_mb": sigma_gen_mb,
        "seed": args.seed,
    }

    print(f"Writing {n_saved} events to {args.output} …", flush=True)
    ak.to_parquet(array, args.output)
    with open(args.output.replace(".parquet", ".json"), "w") as f:
        json.dump(meta, f, indent=2)
    print(f"Done.  sigma_gen = {sigma_gen_mb*1e12:.4g} fb   "
          f"({n_tried} trials, efficiency {n_saved/n_tried*100:.1f}%)")
    return array


if __name__ == "__main__":
    generate_and_save(parse_args())
