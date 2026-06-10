#!/usr/bin/env python3
"""
Convert Herwig HepMC2 ASCII output to the EICFFS event Parquet schema
(same columns as generate_events.py), so analyze_jets.py runs unchanged.

DIS kinematics are reconstructed from the beams and the scattered electron
(highest-energy final-state electron), mirroring the Pythia treatment.

Usage:
    python convert_hepmc.py run.hepmc --config 10x100 --output data/events_10x100_herwig.parquet \
        --sigma-fb 1.1e7
"""

import argparse
import json
import os
import sys

import numpy as np
import awkward as ak

sys.path.insert(0, os.path.dirname(__file__))
from generate_events import BEAM_CONFIGS
from utils.dis_kinematics import DISKinematics

M_E = 0.000511
M_P = 0.938272


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("input", help="HepMC2 ASCII file from Herwig")
    p.add_argument("--config", required=True, choices=sorted(BEAM_CONFIGS))
    p.add_argument("--output", required=True)
    p.add_argument("--sigma-fb", type=float, required=True,
                   help="Generated cross section in fb (from Herwig .out)")
    p.add_argument("--max-events", type=int, default=None)
    return p.parse_args()


def main():
    args = parse_args()
    import pyhepmc

    e_e, e_p = BEAM_CONFIGS[args.config]
    k_in = np.array([0.0, 0.0, -np.sqrt(e_e**2 - M_E**2), e_e])

    ev_Q2, ev_W, ev_x, ev_y, ev_kout = [], [], [], [], []
    par_px, par_py, par_pz, par_e, par_charge = [], [], [], [], []

    # PDG charge lookup (cached; covers any hadron Herwig declares stable)
    from functools import lru_cache

    @lru_cache(maxsize=None)
    def charge_of(pid):
        try:
            from particle import Particle
            q = Particle.from_pdgid(pid).charge
            return int(round(q)) if q is not None else 0
        except Exception:
            return 0

    n_read = 0
    n_saved = 0
    with pyhepmc.open(args.input) as f:
        for evt in f:
            n_read += 1
            if args.max_events and n_saved >= args.max_events:
                break

            k_out = None
            best_E = -1.0
            hadrons = []
            for prt in evt.particles:
                if prt.status != 1:
                    continue
                pid = prt.pid
                m = prt.momentum
                if pid == 11 and m.e > best_E:        # scattered electron
                    best_E = m.e
                    k_out = (m.px, m.py, m.pz, m.e)
                    continue
                if abs(pid) in {11, 12, 13, 14, 15, 16}:
                    continue
                hadrons.append((m.px, m.py, m.pz, m.e, charge_of(pid)))
            if k_out is None:
                continue

            P_in = np.array([0.0, 0.0, np.sqrt(e_p**2 - M_P**2), e_p])
            kin = DISKinematics(k_in, np.array(k_out), P_in)
            if not kin.valid:
                continue

            ev_Q2.append(kin.Q2)
            ev_W.append(kin.W)
            ev_x.append(kin.x)
            ev_y.append(kin.y)
            ev_kout.append(list(k_out))
            par_px.append([h[0] for h in hadrons])
            par_py.append([h[1] for h in hadrons])
            par_pz.append([h[2] for h in hadrons])
            par_e.append([h[3] for h in hadrons])
            par_charge.append([h[4] for h in hadrons])
            n_saved += 1
            if n_saved % 50_000 == 0:
                print(f"  {n_saved:,} events", flush=True)

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
            "px": ak.Array(par_px), "py": ak.Array(par_py),
            "pz": ak.Array(par_pz), "e": ak.Array(par_e),
            "charge": ak.Array(par_charge),
        }),
    })

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    ak.to_parquet(array, args.output)
    meta = {
        "config": args.config,
        "variation": "herwig",
        "generator": "herwig7",
        "e_e": e_e, "e_p": e_p,
        "Q2min": 25.0, "Q2max": 1000.0,
        "mpi": False,
        "n_events": n_saved,
        "n_tried": n_read,
        "sigma_gen_mb": args.sigma_fb / 1e12,
        "seed": 0,
    }
    with open(args.output.replace(".parquet", ".json"), "w") as f:
        json.dump(meta, f, indent=2)
    print(f"Wrote {n_saved:,} events ({n_read:,} read) -> {args.output}")


if __name__ == "__main__":
    main()
