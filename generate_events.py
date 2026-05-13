#!/usr/bin/env python3
"""
EIC DIS Event Generator for the FFS (Frame-dependent Fragmentation Shift) study.

Generates Pythia8 neutral-current DIS events at EIC kinematics and saves
particle-level data to a Parquet file for downstream analysis.

Default EIC configuration: 10 GeV electrons on 100 GeV protons (√s ≈ 63 GeV).

Usage
-----
    python generate_events.py [options]
    python generate_events.py --n-events 500000 --output data/events.parquet --quiet

Reference: arXiv:2308.10951
"""

import argparse
import sys
import os
import numpy as np
import awkward as ak

sys.path.insert(0, os.path.dirname(__file__))
from utils.dis_kinematics import DISKinematics


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Generate EIC DIS events with Pythia8 for FFS effect study",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--n-events", type=int, default=200_000,
                   help="Number of events to save")
    p.add_argument("--electron-energy", type=float, default=10.0,
                   help="Electron beam energy (GeV)")
    p.add_argument("--proton-energy", type=float, default=100.0,
                   help="Proton beam energy (GeV)")
    p.add_argument("--Q2min", type=float, default=1.0,
                   help="Minimum Q² (GeV²)")
    p.add_argument("--Q2max", type=float, default=1000.0,
                   help="Maximum Q² (GeV²)")
    p.add_argument("--output", type=str, default="data/events.parquet",
                   help="Output Parquet file path")
    p.add_argument("--seed", type=int, default=42,
                   help="Pythia8 random seed")
    p.add_argument("--quiet", action="store_true",
                   help="Suppress Pythia8 banner and progress output")
    return p.parse_args()


# ---------------------------------------------------------------------------
# Pythia8 initialisation
# ---------------------------------------------------------------------------

def init_pythia(args):
    """Initialise Pythia8 for neutral-current DIS at EIC kinematics."""
    try:
        import pythia8
    except ImportError:
        sys.exit(
            "pythia8 is not installed.\n"
            "Install it with:  conda install -c conda-forge pythia8\n"
        )

    pythia = pythia8.Pythia()

    def cfg(s):
        pythia.readString(s)

    # ── Output verbosity ────────────────────────────────────────────────────
    if args.quiet:
        cfg("Print:quiet = on")
    cfg("Next:numberCount = 50000")

    # ── EIC beam layout ─────────────────────────────────────────────────────
    # idA (proton) travels in the +z direction
    # idB (electron) travels in the −z direction
    cfg("Beams:frameType = 2")                  # fixed beam energies
    cfg("Beams:idA = 2212")                     # proton
    cfg("Beams:idB = 11")                       # electron
    cfg(f"Beams:eA = {args.proton_energy:.4f}") # proton energy (GeV)
    cfg(f"Beams:eB = {args.electron_energy:.4f}") # electron energy (GeV)

    # ── Neutral-current DIS ─────────────────────────────────────────────────
    cfg("WeakBosonExchange:ff2ff(t:gmZ) = on")

    # ── Phase-space cuts ────────────────────────────────────────────────────
    cfg(f"PhaseSpace:Q2min = {args.Q2min}")
    cfg(f"PhaseSpace:Q2max = {args.Q2max}")

    # ── PDF ─────────────────────────────────────────────────────────────────
    # PDF set 13 = NNPDF2.3 QCD+QED LO in Pythia 8.2+.
    # Verify with your Pythia8 version; alternatively use PDF:pSet = "LHAPDF6:NNPDF23_lo_as_0130_qed".
    cfg("PDF:pSet = 13")

    # ── Shower / MPI ────────────────────────────────────────────────────────
    cfg("SpaceShower:rapidityOrder = off")
    cfg("PartonLevel:MPI = off")   # pure perturbative fragmentation study

    # ── Hadronisation ───────────────────────────────────────────────────────
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

def extract_kinematics(event):
    """
    Reconstruct DIS kinematics from the Pythia8 event record.

    Incoming beams have status −12.
    The hard-process outgoing lepton has status 23.

    Returns a DISKinematics object (may have .valid == False on failure).
    """
    k_in = P_in = k_out = None

    for j in range(event.size()):
        p = event[j]
        pid = p.id()
        st = p.status()

        if st == -12:                          # incoming beam particle
            if abs(pid) in {11, 13, 15}:      # lepton flavour
                k_in = (p.px(), p.py(), p.pz(), p.e())
            elif pid == 2212:                  # proton
                P_in = (p.px(), p.py(), p.pz(), p.e())

        elif st == 23 and abs(pid) in {11, 13, 15}:  # hard-process outgoing lepton
            k_out = (p.px(), p.py(), p.pz(), p.e())

    if k_in is None or P_in is None or k_out is None:
        return None

    return DISKinematics(k_in, k_out, P_in)


# ---------------------------------------------------------------------------
# Main generation loop
# ---------------------------------------------------------------------------

def generate_and_save(args):
    """Run the event loop and write output to Parquet."""
    pythia = init_pythia(args)
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)

    # Event-level scalars
    ev_Q2, ev_W, ev_x, ev_y = [], [], [], []

    # Particle-level ragged arrays (one sub-list per event)
    par_px, par_py, par_pz, par_e = [], [], [], []
    par_pdg, par_charge = [], []

    n_saved = 0
    n_tried = 0

    while n_saved < args.n_events:
        if not pythia.next():
            continue
        n_tried += 1

        kin = extract_kinematics(pythia.event)
        if kin is None or not kin.valid:
            continue

        # Collect final-state hadrons (exclude the scattered lepton)
        fpx, fpy, fpz, fe = [], [], [], []
        fpdg, fcharge = [], []

        for j in range(pythia.event.size()):
            p = pythia.event[j]
            if not p.isFinal():
                continue
            pid = abs(p.id())
            if pid in {11, 12, 13, 14, 15, 16}:  # skip all leptons/neutrinos
                continue

            fpx.append(p.px())
            fpy.append(p.py())
            fpz.append(p.pz())
            fe.append(p.e())
            fpdg.append(p.id())
            fcharge.append(int(round(p.charge())))

        ev_Q2.append(kin.Q2)
        ev_W.append(kin.W)
        ev_x.append(kin.x)
        ev_y.append(kin.y)

        par_px.append(fpx)
        par_py.append(fpy)
        par_pz.append(fpz)
        par_e.append(fe)
        par_pdg.append(fpdg)
        par_charge.append(fcharge)

        n_saved += 1
        if n_saved % 50_000 == 0 and not args.quiet:
            print(f"  Saved {n_saved:>7d} / {args.n_events} events "
                  f"  (efficiency {n_saved/n_tried*100:.1f}%)")

    pythia.stat()

    # ── Build awkward array ────────────────────────────────────────────────
    print(f"\nBuilding awkward array …", flush=True)
    array = ak.Array({
        "Q2":      np.array(ev_Q2,  dtype=np.float32),
        "W":       np.array(ev_W,   dtype=np.float32),
        "x":       np.array(ev_x,   dtype=np.float32),
        "y":       np.array(ev_y,   dtype=np.float32),
        "particles": {
            "px":     ak.Array(par_px),
            "py":     ak.Array(par_py),
            "pz":     ak.Array(par_pz),
            "e":      ak.Array(par_e),
            "pdg":    ak.Array(par_pdg),
            "charge": ak.Array(par_charge),
        },
    })

    print(f"Writing {n_saved} events to {args.output} …", flush=True)
    ak.to_parquet(array, args.output)
    print(f"Done.  ({n_tried} trials, efficiency {n_saved/n_tried*100:.1f}%)")
    return array


if __name__ == "__main__":
    args = parse_args()
    generate_and_save(args)
