#!/usr/bin/env python3
"""
EIC DIS Event Generator for the FFS (Frame-dependent Fragmentation Shift) study.

Generates Pythia8 neutral-current DIS events at EIC kinematics and saves
particle-level data plus the hard-process four-vectors to a Parquet file for
downstream analysis.

Default EIC configuration: 10 GeV electrons on 100 GeV protons (√s ≈ 63 GeV).

Usage
-----
    python generate_events.py [options]
    python generate_events.py --n-events 250000 --seed 1 --output data/events_1.parquet --quiet

Reference: arXiv:2308.10951

Physics notes (see PHYSICS_AUDIT.md)
------------------------------------
* Lepton-beam ISR is switched off (PDF:lepton = off), the standard choice for
  EIC studies.  With it on, q = k_beam − k' is *not* the exchanged boson.
* The hard-process lepton is taken from ``pythia.process``, where its status
  is always +23.  In ``pythia.event`` the same particle usually carries
  status −23 because the shower makes a copy of it; testing ``status == 23``
  there silently drops ~96 % of events and keeps a biased remnant.
* ``SpaceShower:dipoleRecoil = on`` is the Pythia-recommended recoil scheme
  for DIS.
* The scattered lepton and every particle descending from it (QED FSR
  photons) are removed from the stored hadronic final state.  Neutrinos are
  removed.  Everything else — charged and neutral hadrons, photons — is kept,
  matching the particle-level treatment in the reference paper.
* Per event we store q, k', the struck outgoing parton and the beam proton so
  the analysis can boost into the γ*p (hadronic CM) and Breit frames and
  identify the current jet.
"""

import argparse
import os
import sys
import time

import numpy as np
import awkward as ak

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
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
    p.add_argument("--Wmin", type=float, default=0.0,
                   help="Keep only events with W above this value (GeV); applied "
                        "before the final-state loop so it is cheap. Zero = keep all.")
    p.add_argument("--output", type=str, default="data/events.parquet",
                   help="Output Parquet file path")
    p.add_argument("--seed", type=int, default=42,
                   help="Pythia8 random seed")
    p.add_argument("--lepton-isr", action="store_true",
                   help="Keep Pythia's lepton-beam ISR on (PDF:lepton = on). "
                        "Off by default, as is standard for EIC studies.")
    p.add_argument("--quiet", action="store_true",
                   help="Suppress Pythia8 banner and progress output")
    return p.parse_args()


# ---------------------------------------------------------------------------
# Pythia8 import (conda package is ``pythia8``, PyPI wheel is ``pythia8mc``)
# ---------------------------------------------------------------------------

def import_pythia():
    try:
        import pythia8
        return pythia8
    except ImportError:
        pass
    try:
        import pythia8mc as pythia8
        return pythia8
    except ImportError:
        sys.exit(
            "pythia8 is not installed.\n"
            "Install it with:  conda install -c conda-forge pythia8\n"
            "            or:  pip install pythia8mc\n"
        )


# ---------------------------------------------------------------------------
# Pythia8 initialisation
# ---------------------------------------------------------------------------

def init_pythia(args):
    """Initialise Pythia8 for neutral-current DIS at EIC kinematics."""
    pythia8 = import_pythia()
    pythia = pythia8.Pythia("", not args.quiet)

    def cfg(s):
        pythia.readString(s)

    # ── Output verbosity ────────────────────────────────────────────────────
    if args.quiet:
        cfg("Print:quiet = on")
    cfg("Next:numberCount = 0")
    cfg("Next:numberShowEvent = 0")
    cfg("Next:numberShowProcess = 0")
    cfg("Next:numberShowInfo = 0")

    # ── EIC beam layout ─────────────────────────────────────────────────────
    # idA (proton) travels in the +z direction, idB (electron) in −z.
    cfg("Beams:frameType = 2")
    cfg("Beams:idA = 2212")
    cfg("Beams:idB = 11")
    cfg(f"Beams:eA = {args.proton_energy:.4f}")
    cfg(f"Beams:eB = {args.electron_energy:.4f}")

    # ── Neutral-current DIS ─────────────────────────────────────────────────
    cfg("WeakBosonExchange:ff2ff(t:gmZ) = on")

    # ── Phase-space cuts ────────────────────────────────────────────────────
    cfg(f"PhaseSpace:Q2min = {args.Q2min}")
    cfg(f"PhaseSpace:Q2max = {args.Q2max}")

    # ── Lepton beam: no ISR so that q = k − k' is the exchanged boson ──────
    cfg(f"PDF:lepton = {'on' if args.lepton_isr else 'off'}")

    # ── Shower / MPI ────────────────────────────────────────────────────────
    # Dipole recoil is the Pythia-recommended ISR recoil scheme for DIS.
    cfg("SpaceShower:dipoleRecoil = on")
    cfg("PartonLevel:MPI = off")          # no MPI in the direct DIS process

    # ── Hadronisation ───────────────────────────────────────────────────────
    cfg("HadronLevel:all = on")

    # ── Random seed ─────────────────────────────────────────────────────────
    cfg("Random:setSeed = on")
    cfg(f"Random:seed = {args.seed}")

    if not pythia.init():
        sys.exit("Pythia8 initialisation failed.")
    return pythia


# ---------------------------------------------------------------------------
# Hard-process kinematics from the Pythia8 process record
# ---------------------------------------------------------------------------

def _p4(p):
    return (p.px(), p.py(), p.pz(), p.e())


def extract_hard_process(pythia):
    """
    Read the hard 2→2 scattering from ``pythia.process``.

    In the process record the incoming beams have status −12, the particles
    entering the hard scattering have status 21 and those leaving it have
    status 23 (always positive here, unlike in ``pythia.event``).

    Returns (k_in, k_out, P_in, parton_out, parton_id) or None.
    """
    proc = pythia.process
    k_in = k_out = P_in = q_out = None
    q_id = 0
    for j in range(proc.size()):
        p = proc[j]
        pid = p.id()
        st = p.status()
        if st == -12:
            if abs(pid) in (11, 13, 15):
                k_in = _p4(p)
            elif pid == 2212:
                P_in = _p4(p)
        elif st == 21 and abs(pid) in (11, 13, 15):
            # Lepton entering the hard scattering (equals the beam lepton
            # when lepton ISR is off; differs when it is on).
            k_in = _p4(p)
        elif st == 23:
            if abs(pid) in (11, 13, 15):
                k_out = _p4(p)
            elif abs(pid) <= 21:
                q_out = _p4(p)
                q_id = pid
    if k_in is None or k_out is None or P_in is None or q_out is None:
        return None
    return k_in, k_out, P_in, q_out, q_id


def scattered_lepton_index(event):
    """Index in ``pythia.event`` of the hard-process scattered lepton."""
    for j in range(event.size()):
        p = event[j]
        if abs(p.status()) == 23 and abs(p.id()) in (11, 13, 15):
            return j
    return -1


# ---------------------------------------------------------------------------
# Main generation loop
# ---------------------------------------------------------------------------

def generate_and_save(args):
    """Run the event loop and write output to Parquet."""
    pythia = init_pythia(args)
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)

    # Event-level scalars and four-vectors
    ev_Q2, ev_W, ev_x, ev_y = [], [], [], []
    ev_q, ev_kout, ev_Pin, ev_parton = [], [], [], []
    ev_parton_id = []

    # Particle-level ragged arrays (one sub-list per event)
    par_px, par_py, par_pz, par_e = [], [], [], []
    par_pdg, par_charge = [], []

    n_saved = 0
    n_tried = 0
    n_bad_kin = 0
    n_below_W = 0
    t0 = time.time()

    while n_saved < args.n_events:
        if not pythia.next():
            continue
        n_tried += 1

        hard = extract_hard_process(pythia)
        if hard is None:
            n_bad_kin += 1
            continue
        k_in, k_out, P_in, parton, parton_id = hard

        kin = DISKinematics(k_in, k_out, P_in)
        if not kin.valid:
            n_bad_kin += 1
            continue
        if kin.W < args.Wmin:
            n_below_W += 1
            continue

        event = pythia.event
        i_lep = scattered_lepton_index(event)

        # Collect the hadronic final state: every final-state particle that
        # is neither the scattered lepton, one of its descendants (QED FSR
        # photons), nor a neutrino.
        fpx, fpy, fpz, fe = [], [], [], []
        fpdg, fcharge = [], []
        for j in range(event.size()):
            p = event[j]
            if not p.isFinal():
                continue
            pid = abs(p.id())
            if pid in (12, 14, 16):
                continue
            if i_lep >= 0 and (j == i_lep or p.isAncestor(i_lep)):
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
        ev_q.append(tuple(kin.q))
        ev_kout.append(k_out)
        ev_Pin.append(P_in)
        ev_parton.append(parton)
        ev_parton_id.append(parton_id)

        par_px.append(fpx)
        par_py.append(fpy)
        par_pz.append(fpz)
        par_e.append(fe)
        par_pdg.append(fpdg)
        par_charge.append(fcharge)

        n_saved += 1
        if n_saved % 25_000 == 0 and not args.quiet:
            rate = n_saved / (time.time() - t0)
            print(f"  Saved {n_saved:>8d} / {args.n_events}   "
                  f"({rate:.0f} ev/s, eta {(args.n_events - n_saved) / rate:.0f} s)",
                  flush=True)

    if not args.quiet:
        pythia.stat()

    def vec4(lst):
        a = np.asarray(lst, dtype=np.float32)
        return {"px": np.ascontiguousarray(a[:, 0]), "py": np.ascontiguousarray(a[:, 1]),
                "pz": np.ascontiguousarray(a[:, 2]), "e": np.ascontiguousarray(a[:, 3])}

    array = ak.Array({
        "Q2": np.array(ev_Q2, dtype=np.float32),
        "W":  np.array(ev_W,  dtype=np.float32),
        "x":  np.array(ev_x,  dtype=np.float32),
        "y":  np.array(ev_y,  dtype=np.float32),
        "q":       vec4(ev_q),        # exchanged boson  k − k'
        "lepton":  vec4(ev_kout),     # scattered lepton (hard process)
        "proton":  vec4(ev_Pin),      # beam proton
        "parton":  vec4(ev_parton),   # struck outgoing parton (hard process)
        "parton_id": np.array(ev_parton_id, dtype=np.int32),
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
    print(f"Done.  {n_tried} generated, {n_bad_kin} failed kinematics, "
          f"{n_below_W} below Wmin, {n_saved} saved in {time.time() - t0:.0f} s.")
    return array


if __name__ == "__main__":
    generate_and_save(parse_args())
