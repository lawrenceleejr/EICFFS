#!/usr/bin/env python3
"""
FFS jet analysis: event Parquet  ->  per-jet table (Parquet).

Implements ANALYSIS_DESIGN.md Sec. 5-6:

  * boost all final-state particles to the gamma*p CM frame (the color
    rest frame of the struck-quark--remnant string);
  * cluster with the e+e- generalised-kT algorithm (p = -1, R = 1.0) in
    that frame -- the current jet lies along the photon axis, where
    beam-axis pT algorithms are degenerate;
  * select the leading-energy jet in the current hemisphere of the
    Breit frame (p_z^Breit < 0), removing target fragmentation;
  * compute the frame-resolved observables: n90/n75/n95 from constituent
    |p| in the lab frame and in the CM frame, charged multiplicity, pTD;
  * optional reco-level pass (--smear): charged-track-only jets with
    parametric ePIC-like smearing, and DIS kinematics from the smeared
    scattered electron (electron method).

Output: one row per selected jet with event kinematics attached.

Usage
-----
    python analyze_jets.py data/events_10x100_baseline.parquet
    python analyze_jets.py data/events_10x100_baseline.parquet --smear
"""

import argparse
import json
import os
import sys

import numpy as np
import awkward as ak

sys.path.insert(0, os.path.dirname(__file__))
from utils.dis_kinematics import four_dot, boost_to_breit_frame
from utils.smearing import smear_tracks, smear_electron

M_E = 0.000511
M_P = 0.938272

JET_R = 1.0          # ee_genkt radius in the gamma*p CM frame
PCM_MIN = 2.0        # GeV, minimum CM-frame jet momentum
PLAB_MIN = 3.0       # GeV, loose lab-frame momentum window (binned later)
PLAB_MAX = 80.0
ETA_LAB_MAX = 3.5    # EIC acceptance

Q2_MIN, Q2_MAX = 25.0, 1000.0
Y_MIN, Y_MAX = 0.05, 0.95
W_MIN = 4.0


# ---------------------------------------------------------------------------
# Vectorised Lorentz boost
# ---------------------------------------------------------------------------

def boost_matrix_to_rest(p_sys):
    """
    Return a function boosting (N,4) arrays of 4-vectors (px,py,pz,E)
    into the rest frame of p_sys.
    """
    p_sys = np.asarray(p_sys, dtype=float)
    E = p_sys[3]
    m = np.sqrt(max(E**2 - p_sys[0]**2 - p_sys[1]**2 - p_sys[2]**2, 1e-12))
    beta = p_sys[:3] / E
    gamma = E / m

    def _boost(v):
        v = np.atleast_2d(np.asarray(v, dtype=float))
        bp = v[:, :3] @ beta
        factor = gamma * (gamma * bp / (gamma + 1.0) - v[:, 3])
        out = np.empty_like(v)
        out[:, :3] = v[:, :3] + factor[:, None] * beta
        out[:, 3] = gamma * (v[:, 3] - bp)
        return out

    return _boost


# ---------------------------------------------------------------------------
# n_x observable (paper's primary FFS quantity), arXiv:2308.10951 Sec. 2
# ---------------------------------------------------------------------------

def compute_n_x(pmags, thresholds):
    """
    Fractional number of momentum-ordered constituents carrying each
    momentum-fraction threshold.  Returns one value per threshold.
    """
    pmags = np.sort(np.asarray(pmags, dtype=float))[::-1]
    total = pmags.sum()
    if len(pmags) == 0 or total <= 0:
        return [np.nan] * len(thresholds)
    cumfrac = np.cumsum(pmags) / total
    out = []
    for thr in thresholds:
        idx = int(np.searchsorted(cumfrac, thr, side="left"))
        if idx >= len(cumfrac):
            out.append(float(len(cumfrac)))
        elif idx == 0:
            out.append(float(thr / cumfrac[0]))
        else:
            frac = (thr - cumfrac[idx - 1]) / (cumfrac[idx] - cumfrac[idx - 1])
            out.append(float(idx + frac))
    return out


# ---------------------------------------------------------------------------
# Electron-method DIS kinematics
# ---------------------------------------------------------------------------

def electron_method(k_in, k_out, P_in):
    """Return (Q2, W, x, y, q) from beam and scattered-lepton 4-vectors."""
    q = k_in - k_out
    Q2 = -four_dot(q, q)
    Pq = four_dot(P_in, q)
    Pk = four_dot(P_in, k_in)
    if Q2 <= 0 or Pq <= 0 or Pk <= 0:
        return None
    y = Pq / Pk
    x = Q2 / (2.0 * Pq)
    W2 = M_P**2 + 2.0 * Pq - Q2
    if W2 <= 0:
        return None
    return Q2, float(np.sqrt(W2)), x, y, q


# ---------------------------------------------------------------------------
# Main loop
# ---------------------------------------------------------------------------

def analyze(args):
    import fastjet as fj

    meta = json.load(open(args.input.replace(".parquet", ".json")))
    e_e, e_p = meta["e_e"], meta["e_p"]
    k_in = np.array([0.0, 0.0, -np.sqrt(e_e**2 - M_E**2), e_e])
    P_in = np.array([0.0, 0.0, np.sqrt(e_p**2 - M_P**2), e_p])

    print(f"Reading {args.input}  "
          f"[{meta['config']}, {meta['variation']}"
          f"{', MPI' if meta.get('mpi') else ''}, "
          f"{'reco-smeared' if args.smear else 'particle level'}]", flush=True)
    events = ak.from_parquet(args.input)
    n_events = len(events)

    # Flatten particle lists once (per-event numpy views are then cheap)
    counts = np.asarray(ak.num(events.particles.px))
    offs = np.zeros(n_events + 1, dtype=np.int64)
    np.cumsum(counts, out=offs[1:])
    f_px = np.asarray(ak.flatten(events.particles.px), dtype=float)
    f_py = np.asarray(ak.flatten(events.particles.py), dtype=float)
    f_pz = np.asarray(ak.flatten(events.particles.pz), dtype=float)
    f_e = np.asarray(ak.flatten(events.particles.e), dtype=float)
    f_ch = np.asarray(ak.flatten(events.particles.charge), dtype=int)
    ev_Q2 = np.asarray(events.Q2, dtype=float)
    ev_W = np.asarray(events.W, dtype=float)
    ev_x = np.asarray(events.x, dtype=float)
    ev_y = np.asarray(events.y, dtype=float)
    ev_ke = np.stack([np.asarray(events[f], dtype=float)
                      for f in ("e_px", "e_py", "e_pz", "e_e")], axis=1)
    del events

    jdef = fj.JetDefinition(fj.ee_genkt_algorithm, JET_R, -1.0)
    rng = np.random.default_rng(args.seed)

    cols = {k: [] for k in [
        "Q2", "W", "x", "y",
        "plab", "ptlab", "etalab", "pcm", "ecm",
        "n90lab", "n75lab", "n95lab", "n90cm", "n75cm", "n95cm",
        "nconst", "nch", "ptd", "gboost",
    ]}
    n_sel_event = 0
    n_jet = 0

    for i in range(n_events):
        sl = slice(offs[i], offs[i + 1])
        px, py, pz, e = f_px[sl], f_py[sl], f_pz[sl], f_e[sl]
        ch = f_ch[sl]
        k_out = ev_ke[i]

        if args.smear:
            px, py, pz, e, ch = smear_tracks(px, py, pz, e, ch, rng)
            k_out = smear_electron(k_out, rng)
            kin = electron_method(k_in, k_out, P_in)
            if kin is None:
                continue
            Q2, W, x, y, q = kin
        else:
            Q2, W, x, y = ev_Q2[i], ev_W[i], ev_x[i], ev_y[i]
            q = k_in - k_out

        # Event selection (reconstructed kinematics in the smeared pass)
        if not (Q2_MIN <= Q2 <= Q2_MAX and Y_MIN < y < Y_MAX and W > W_MIN):
            continue
        if len(px) < 1:
            continue
        n_sel_event += 1

        # Boost particles to the gamma*p CM frame
        boost = boost_matrix_to_rest(P_in + q)
        v_lab = np.stack([px, py, pz, e], axis=1)
        v_cm = boost(v_lab)

        # Cluster with ee_genkt (anti-kT-like, p=-1) in the CM frame
        pjs = []
        for j in range(len(v_cm)):
            pj = fj.PseudoJet(v_cm[j, 0], v_cm[j, 1], v_cm[j, 2], v_cm[j, 3])
            pj.set_user_index(j)
            pjs.append(pj)
        cs = fj.ClusterSequence(pjs, jdef)
        jets = sorted(cs.inclusive_jets(0.0), key=lambda jj: -jj.e())
        if not jets:
            continue

        # Leading-energy jet in the current hemisphere of the Breit frame
        jet = None
        for cand in jets:
            idx = [c.user_index() for c in cand.constituents()]
            jlab = v_lab[idx].sum(axis=0)
            jbreit = boost_to_breit_frame(P_in, q, jlab)
            if jbreit[2] < 0.0:          # current hemisphere
                jet = cand
                jet_idx = idx
                jet_lab = jlab
                break
        if jet is None:
            continue

        pcm = np.sqrt(jet.px()**2 + jet.py()**2 + jet.pz()**2)
        plab = np.sqrt(jet_lab[0]**2 + jet_lab[1]**2 + jet_lab[2]**2)
        ptlab = np.sqrt(jet_lab[0]**2 + jet_lab[1]**2)
        if plab < 1e-9 or pcm < PCM_MIN:
            continue
        etalab = np.arctanh(np.clip(jet_lab[2] / plab, -1 + 1e-9, 1 - 1e-9))
        if not (PLAB_MIN <= plab <= PLAB_MAX and abs(etalab) < ETA_LAB_MAX):
            continue

        # Frame-resolved constituent observables
        c_lab = v_lab[jet_idx]
        c_cm = v_cm[jet_idx]
        pm_lab = np.sqrt((c_lab[:, :3]**2).sum(axis=1))
        pm_cm = np.sqrt((c_cm[:, :3]**2).sum(axis=1))
        n75l, n90l, n95l = compute_n_x(pm_lab, (0.75, 0.90, 0.95))
        n75c, n90c, n95c = compute_n_x(pm_cm, (0.75, 0.90, 0.95))
        pt_lab = np.sqrt((c_lab[:, :2]**2).sum(axis=1))
        sum_pt = pt_lab.sum()
        ptd = float(np.sqrt((pt_lab**2).sum()) / sum_pt) if sum_pt > 0 else np.nan

        gboost = (P_in + q)[3] / W      # Lorentz gamma of the gamma*p frame
        row = (Q2, W, x, y, plab, ptlab, etalab, pcm, jet.e(),
               n90l, n75l, n95l, n90c, n75c, n95c,
               len(jet_idx), int(np.sum(np.abs(ch[jet_idx]) > 0)), ptd, gboost)
        for k, val in zip(cols, row):
            cols[k].append(val)
        n_jet += 1

        if (i + 1) % 100_000 == 0:
            print(f"  {i+1:>8,}/{n_events:,} events, {n_jet:,} jets", flush=True)

    # ── Write jet table ────────────────────────────────────────────────────
    out = args.output
    if out is None:
        base = os.path.basename(args.input).replace("events_", "jets_")
        base = base.replace(".parquet",
                            "_reco.parquet" if args.smear else "_truth.parquet")
        out = os.path.join(os.path.dirname(args.input), base)
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)

    table = ak.Array({k: np.array(v, dtype=np.float32) for k, v in cols.items()})
    ak.to_parquet(table, out)

    meta_out = dict(meta)
    meta_out.update({
        "level": "reco" if args.smear else "truth",
        "n_events_input": n_events,
        "n_events_selected": n_sel_event,
        "n_jets": n_jet,
        "jet_R": JET_R,
        "jet_algorithm": "ee_genkt(p=-1) in gamma*p CM frame",
    })
    with open(out.replace(".parquet", ".json"), "w") as f:
        json.dump(meta_out, f, indent=2)

    print(f"Wrote {n_jet:,} jets ({n_sel_event:,} selected events) -> {out}")
    return out


def parse_args():
    p = argparse.ArgumentParser(
        description="FFS jet analysis: events -> per-jet table",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("input", help="Input event Parquet from generate_events.py")
    p.add_argument("--output", default=None, help="Output jet-table path")
    p.add_argument("--smear", action="store_true",
                   help="Reco-level pass: parametric smearing + track jets")
    p.add_argument("--seed", type=int, default=7, help="Smearing RNG seed")
    return p.parse_args()


if __name__ == "__main__":
    analyze(parse_args())
