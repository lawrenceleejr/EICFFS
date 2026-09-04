#!/usr/bin/env python3
"""
FFS Effect Analysis Script
===========================
Reads particle-level EIC DIS events (Parquet) produced by generate_events.py,
finds anti-kT R = 0.4 jets in the laboratory frame, and writes

  * a flat per-jet TTree ``jets`` with lab-, γ*p-CM- and Breit-frame
    quantities and the n₉₀ observable, and
  * a per-event TTree ``events``, and
  * the legacy histograms used by make_plots.py (filled with current jets),

to a ROOT file via uproot.

Primary observable
------------------
n₉₀ — the *fractional* minimum number of jet constituents (ordered by
decreasing 3-momentum magnitude) needed to account for 90 % of the jet's
total scalar momentum, with per-jet linear interpolation, exactly as in
arXiv:2308.10951 Sec. 2.  It is computed from the lab-frame constituent
momenta (n90), from the same constituents boosted into the γ*p frame
(n90_hcm), and from charged constituents only (n90_ch).

Frames
------
The colour-connected system in DIS is the whole hadronic final state, with
four-momentum P + q and invariant mass W.  Its rest frame — the γ*p or
hadronic centre-of-mass (HCM) frame — is therefore the colour rest frame of
arXiv:2308.10951.  It is *not* the Breit frame: the two are related by a
boost along the boson axis, and the struck quark carries W/2 in the HCM
frame but Q/2 in the Breit frame.  The Breit frame is used only to define
the current hemisphere (particles moving along the boson direction).

Usage
-----
    python analyze_events.py data/events_*.parquet --output data/analysis.root
    python analyze_events.py data/events.parquet --no-fastjet     # cone fallback

Reference: arXiv:2308.10951  (Phys.Lett.B 866, 2025, 139561)
"""

import argparse
import glob
import os
import sys
import time

import numpy as np
import awkward as ak
import uproot
import hist
import vector

vector.register_awkward()

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils.dis_kinematics import (hcm_boost_matrix, breit_boost_matrix,
                                  apply_boost, rapidity)


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

# Legacy histogram binning (kept for make_plots.py)
W_BINS = np.array([5.0, 10.0, 20.0, 30.0, 40.0, 55.0], dtype=float)
P_LAB_BINS = np.array([2.0, 5.0, 10.0, 20.0], dtype=float)
NMAX = 60
N90_MAX = 25.0
N90_BINS = 50

# Jet finding
JET_R = 0.4          # anti-kT radius, as in arXiv:2308.10951
JET_ETA_MAX = 3.5    # EIC central-detector acceptance (lab)
JET_PT_MIN = 2.0     # GeV, lab frame
CHUNK = 20_000       # events per processing chunk


# ---------------------------------------------------------------------------
# n_x observable (arXiv:2308.10951 Sec. 2)
# ---------------------------------------------------------------------------

def compute_n_x(const_pmags, threshold=0.90):
    """
    Scalar reference implementation of n_x for one jet.

    1. Sort constituents by decreasing |p|.
    2. Cumulatively sum, normalised to the jet total.
    3. Linearly interpolate to the threshold.

    >>> round(compute_n_x([10., 5., 3., 1.]), 3)      # cumfrac .526 .789 .947 1
    2.7
    """
    pmags = np.sort(np.asarray(const_pmags, dtype=float))[::-1]
    if len(pmags) == 0 or pmags.sum() <= 0:
        return np.nan
    cumfrac = np.cumsum(pmags) / pmags.sum()
    idx = int(np.searchsorted(cumfrac, threshold, side="left"))
    if idx >= len(cumfrac):
        return float(len(cumfrac))
    prev = cumfrac[idx - 1] if idx > 0 else 0.0
    return float(idx + (threshold - prev) / (cumfrac[idx] - prev))


def n_x_segments(pmag, seg_id, n_seg, threshold=0.90):
    """
    Vectorised n_x for many jets at once.

    Parameters
    ----------
    pmag   : (K,) constituent |p| values, any order
    seg_id : (K,) jet index of each constituent, in [0, n_seg)
    n_seg  : number of jets

    Returns (n_seg,) array of n_x, NaN for jets without constituents.
    """
    out = np.full(n_seg, np.nan)
    if len(pmag) == 0:
        return out
    order = np.lexsort((-pmag, seg_id))
    p = pmag[order]
    s = seg_id[order]
    counts = np.bincount(s, minlength=n_seg)
    starts = np.concatenate([[0], np.cumsum(counts)[:-1]])
    nonempty = counts > 0
    total = np.bincount(s, weights=p, minlength=n_seg)

    cum = np.cumsum(p)
    base = cum[starts[nonempty]] - p[starts[nonempty]]
    base_full = np.zeros(n_seg)
    base_full[nonempty] = base
    cumfrac = (cum - base_full[s]) / np.where(total[s] > 0, total[s], 1.0)

    pos = np.arange(len(p)) - starts[s]
    big = len(p) + 1
    pos_hit = np.where(cumfrac >= threshold - 1e-12, pos, big)
    first = np.minimum.reduceat(pos_hit, starts[nonempty])
    first = np.minimum(first, counts[nonempty] - 1)          # numerical guard

    st = starts[nonempty]
    cf_idx = cumfrac[st + first]
    cf_prev = np.where(first > 0, cumfrac[np.maximum(st + first - 1, 0)], 0.0)
    denom = np.where(cf_idx - cf_prev > 0, cf_idx - cf_prev, 1.0)
    out[nonempty] = first + (threshold - cf_prev) / denom
    return out


# ---------------------------------------------------------------------------
# Jet finding
# ---------------------------------------------------------------------------

def cluster_fastjet(parts, ptmin=JET_PT_MIN, R=JET_R):
    """Anti-kT clustering of a jagged Momentum4D array; returns (jets, constituent_index)."""
    import fastjet as fj
    jetdef = fj.JetDefinition(fj.antikt_algorithm, R)
    cs = fj.ClusterSequence(parts, jetdef)
    jets = cs.inclusive_jets(min_pt=ptmin)
    cidx = cs.constituent_index(min_pt=ptmin)
    return jets, cidx


def cluster_cone(parts, ptmin=JET_PT_MIN, R=JET_R):
    """
    Greedy fixed-cone fallback (debugging / CI only).  Returns the same
    structure as cluster_fastjet.
    """
    jets_out, cidx_out = [], []
    px_all, py_all, pz_all, e_all = (ak.to_list(parts.px), ak.to_list(parts.py),
                                     ak.to_list(parts.pz), ak.to_list(parts.E))
    for px, py, pz, e in zip(px_all, py_all, pz_all, e_all):
        px, py, pz, e = map(lambda v: np.asarray(v, float), (px, py, pz, e))
        n = len(px)
        pt = np.hypot(px, py)
        p = np.sqrt(pt**2 + pz**2)
        eta = np.arctanh(np.clip(pz / np.maximum(p, 1e-9), -1 + 1e-7, 1 - 1e-7))
        phi = np.arctan2(py, px)
        used = np.zeros(n, bool)
        ev_jets, ev_cidx = [], []
        for seed in np.argsort(-pt):
            if used[seed] or pt[seed] < 0.5 * ptmin:
                continue
            dphi = (phi - phi[seed] + np.pi) % (2 * np.pi) - np.pi
            members = np.where(~used & (np.hypot(eta - eta[seed], dphi) < R))[0]
            used[members] = True
            jpx, jpy, jpz, je = (px[members].sum(), py[members].sum(),
                                 pz[members].sum(), e[members].sum())
            if np.hypot(jpx, jpy) < ptmin:
                continue
            ev_jets.append({"px": jpx, "py": jpy, "pz": jpz, "E": je})
            ev_cidx.append(members.tolist())
        jets_out.append(ev_jets)
        cidx_out.append(ev_cidx)
    jets = ak.Array(jets_out, with_name="Momentum4D") if any(jets_out) else \
        ak.zip({"px": ak.Array([[]] * len(jets_out)), "py": ak.Array([[]] * len(jets_out)),
                "pz": ak.Array([[]] * len(jets_out)), "E": ak.Array([[]] * len(jets_out))},
               with_name="Momentum4D")
    return jets, ak.Array(cidx_out)


# ---------------------------------------------------------------------------
# Legacy histograms (for make_plots.py)
# ---------------------------------------------------------------------------

def make_histograms():
    W_ax = hist.axis.Variable(W_BINS, name="W", label=r"$W$ [GeV]")
    PL_ax = hist.axis.Variable(P_LAB_BINS, name="plab", label=r"$|p|_{\rm lab}$ [GeV]")
    N_ax = hist.axis.Regular(NMAX, 0, NMAX, name="N", label=r"$N_{\rm charged}$")
    N90_ax = hist.axis.Regular(N90_BINS, 0, N90_MAX, name="n90", label=r"$n_{90}$")
    Q2_ax = hist.axis.Regular(50, 1, 1000, name="Q2", label=r"$Q^2$ [GeV$^2$]",
                              transform=hist.axis.transform.log)
    x_ax = hist.axis.Regular(50, 1e-4, 1, name="x", label=r"Bjorken $x$",
                             transform=hist.axis.transform.log)
    y_ax = hist.axis.Regular(50, 0.01, 1, name="y", label=r"Inelasticity $y$")
    W_fine = hist.axis.Regular(50, 5, 55, name="W_fine", label=r"$W$ [GeV]")
    D = hist.storage.Double()
    return {
        "n90_3d": hist.Hist(W_ax, PL_ax, N90_ax, storage=D),
        "sum_n90_vs_W": hist.Hist(W_fine, PL_ax, storage=D),
        "count_n90_vs_W": hist.Hist(W_fine, PL_ax, storage=D),
        "mult_3d": hist.Hist(W_ax, PL_ax, N_ax, storage=D),
        "sum_N_vs_W": hist.Hist(W_fine, PL_ax, storage=D),
        "count_vs_W": hist.Hist(W_fine, PL_ax, storage=D),
        "Q2": hist.Hist(Q2_ax, storage=D),
        "x": hist.Hist(x_ax, storage=D),
        "y": hist.Hist(y_ax, storage=D),
        "W": hist.Hist(W_fine, storage=D),
        "Q2_vs_W": hist.Hist(
            hist.axis.Regular(50, 5, 55, name="W2", label=r"$W$ [GeV]"),
            hist.axis.Regular(50, 1, 1000, name="Q22", label=r"$Q^2$ [GeV$^2$]",
                              transform=hist.axis.transform.log), storage=D),
        "jet_eta_pt": hist.Hist(
            hist.axis.Regular(60, -4, 4, name="eta", label=r"Jet $\eta_{\rm lab}$"),
            hist.axis.Regular(50, 0, 30, name="pt", label=r"Jet $p_T$ [GeV]"), storage=D),
        "n_jets": hist.Hist(
            hist.axis.Regular(20, 0, 20, name="nj", label=r"$N_{\rm jets}$ per event"),
            storage=D),
    }


def fill_legacy(hists, ev, jets):
    """Fill the legacy histograms from flat event/jet dictionaries."""
    W, Q2, x, y = ev["W"], ev["Q2"], ev["x"], ev["y"]
    inW = (W >= 5.0) & (W <= 55.0)
    hists["W"].fill(W_fine=W[inW])
    hists["Q2"].fill(Q2=Q2)
    hists["x"].fill(x=x)
    hists["y"].fill(y=y)
    hists["Q2_vs_W"].fill(W2=np.clip(W, 5.01, 54.99), Q22=np.clip(Q2, 1.01, 999.9))
    hists["n_jets"].fill(nj=np.minimum(ev["n_jets"], 19))

    hists["jet_eta_pt"].fill(eta=np.clip(jets["eta"], -3.99, 3.99),
                             pt=np.clip(jets["pt"], 0.01, 29.99))
    sel = (jets["current"] & (jets["W"] >= 5.0) & (jets["W"] <= 55.0)
           & (jets["plab"] >= P_LAB_BINS[0]) & (jets["plab"] < P_LAB_BINS[-1]))
    jW, jp, n90, nch = jets["W"][sel], jets["plab"][sel], jets["n90"][sel], jets["n_charged"][sel]
    ok = np.isfinite(n90)
    hists["n90_3d"].fill(W=jW[ok], plab=jp[ok],
                         n90=np.minimum(n90[ok], N90_MAX - N90_MAX / N90_BINS))
    hists["sum_n90_vs_W"].fill(W_fine=jW[ok], plab=jp[ok], weight=n90[ok])
    hists["count_n90_vs_W"].fill(W_fine=jW[ok], plab=jp[ok])
    hists["mult_3d"].fill(W=jW, plab=jp, N=np.minimum(nch, NMAX - 1))
    hists["sum_N_vs_W"].fill(W_fine=jW, plab=jp, weight=nch.astype(float))
    hists["count_vs_W"].fill(W_fine=jW, plab=jp)


# ---------------------------------------------------------------------------
# Per-chunk analysis
# ---------------------------------------------------------------------------

def _vec4(rec):
    return np.stack([ak.to_numpy(rec.px), ak.to_numpy(rec.py),
                     ak.to_numpy(rec.pz), ak.to_numpy(rec.e)], axis=1).astype(float)


def analyze_chunk(events, use_fastjet=True):
    """Return (event_dict, jet_dict) of flat numpy arrays for one chunk."""
    n_ev = len(events)
    W = ak.to_numpy(events.W).astype(float)
    Q2 = ak.to_numpy(events.Q2).astype(float)
    x = ak.to_numpy(events.x).astype(float)
    y = ak.to_numpy(events.y).astype(float)
    P_in = _vec4(events.proton)
    q = _vec4(events.q)
    k_out = _vec4(events.lepton)
    parton = _vec4(events.parton)
    parton_id = ak.to_numpy(events.parton_id)

    # Frames --------------------------------------------------------------
    L_hcm = hcm_boost_matrix(P_in, q)
    L_breit = breit_boost_matrix(P_in, q)
    y_hcm = rapidity(P_in + q)                      # rapidity of the colour rest frame in the lab
    gamma_hcm = (P_in + q)[:, 3] / W                # Lorentz factor of that frame
    q_hcm = apply_boost(L_hcm, q)
    qhat_hcm = q_hcm[:, :3] / np.linalg.norm(q_hcm[:, :3], axis=1)[:, None]
    q_breit = apply_boost(L_breit, q)
    qhat_breit = q_breit[:, :3] / np.linalg.norm(q_breit[:, :3], axis=1)[:, None]

    # Jets ----------------------------------------------------------------
    parts = ak.zip({"px": events.particles.px, "py": events.particles.py,
                    "pz": events.particles.pz, "E": events.particles.e},
                   with_name="Momentum4D")
    charge = events.particles.charge
    if use_fastjet:
        jets, cidx = cluster_fastjet(parts)
    else:
        jets, cidx = cluster_cone(parts)

    keep = abs(jets.eta) < JET_ETA_MAX
    jets, cidx = jets[keep], cidx[keep]
    order = ak.argsort(jets.pt, ascending=False)
    jets, cidx = jets[order], cidx[order]
    n_jets = ak.to_numpy(ak.num(jets))
    ev_of_jet = np.repeat(np.arange(n_ev), n_jets)
    rank = ak.to_numpy(ak.flatten(ak.local_index(jets)))

    J = np.stack([ak.to_numpy(ak.flatten(jets.px)), ak.to_numpy(ak.flatten(jets.py)),
                  ak.to_numpy(ak.flatten(jets.pz)), ak.to_numpy(ak.flatten(jets.E))], axis=1)
    n_j = len(J)
    plab = np.linalg.norm(J[:, :3], axis=1)
    pt = np.hypot(J[:, 0], J[:, 1])
    eta = np.arcsinh(J[:, 2] / np.maximum(pt, 1e-9))
    phi = np.arctan2(J[:, 1], J[:, 0])
    m2 = J[:, 3]**2 - plab**2
    mass = np.sqrt(np.maximum(m2, 0.0))

    J_hcm = np.einsum("nij,nj->ni", L_hcm[ev_of_jet], J)
    p_hcm = np.linalg.norm(J_hcm[:, :3], axis=1)
    cos_hcm = np.einsum("ni,ni->n", J_hcm[:, :3], qhat_hcm[ev_of_jet]) / np.maximum(p_hcm, 1e-9)
    J_breit = np.einsum("nij,nj->ni", L_breit[ev_of_jet], J)
    p_breit = np.linalg.norm(J_breit[:, :3], axis=1)
    pz_breit = np.einsum("ni,ni->n", J_breit[:, :3], qhat_breit[ev_of_jet])
    current = pz_breit > 0                           # current hemisphere of the Breit frame

    # Struck parton / lepton matching (lab) ------------------------------
    p_pt = np.hypot(parton[:, 0], parton[:, 1])
    p_eta = np.arcsinh(parton[:, 2] / np.maximum(p_pt, 1e-9))
    p_phi = np.arctan2(parton[:, 1], parton[:, 0])
    dphi = (phi - p_phi[ev_of_jet] + np.pi) % (2 * np.pi) - np.pi
    dR_parton = np.hypot(eta - p_eta[ev_of_jet], dphi)
    l_phi = np.arctan2(k_out[:, 1], k_out[:, 0])
    dphi_lepton = np.abs((phi - l_phi[ev_of_jet] + np.pi) % (2 * np.pi) - np.pi)

    # Constituents ----------------------------------------------------------
    # Gather constituents: cidx is [ev][jet][c]; index the [ev][c] arrays
    # with the per-event flattened indices and restore the jet structure.
    flat_idx = ak.flatten(cidx, axis=2)
    per_jet = ak.flatten(ak.num(cidx, axis=2))
    consts = ak.unflatten(parts[flat_idx], per_jet, axis=1)    # [ev][jet][c]
    cch = ak.unflatten(charge[flat_idx], per_jet, axis=1)
    c_flat = ak.flatten(consts, axis=1)               # [jet][c]
    cch_flat = ak.flatten(cch, axis=1)
    n_const = ak.to_numpy(ak.num(c_flat))
    n_charged = ak.to_numpy(ak.sum(cch_flat != 0, axis=1))
    jet_of_c = np.repeat(np.arange(n_j), n_const)
    C = np.stack([ak.to_numpy(ak.flatten(c_flat.px)), ak.to_numpy(ak.flatten(c_flat.py)),
                  ak.to_numpy(ak.flatten(c_flat.pz)), ak.to_numpy(ak.flatten(c_flat.E))], axis=1)
    c_charged = ak.to_numpy(ak.flatten(cch_flat)) != 0
    c_p = np.linalg.norm(C[:, :3], axis=1)
    C_hcm = np.einsum("nij,nj->ni", L_hcm[ev_of_jet[jet_of_c]], C)
    c_p_hcm = np.linalg.norm(C_hcm[:, :3], axis=1)

    n90 = n_x_segments(c_p, jet_of_c, n_j)
    n90_hcm = n_x_segments(c_p_hcm, jet_of_c, n_j)
    n90_ch = n_x_segments(c_p[c_charged], jet_of_c[c_charged], n_j)
    # Momentum fraction of the leading constituent (a classic FF-like quantity)
    z_lead = np.zeros(n_j)
    if len(c_p):
        np.maximum.at(z_lead, jet_of_c, c_p)
        z_lead /= np.bincount(jet_of_c, weights=c_p, minlength=n_j).clip(1e-9)

    n_current = np.bincount(ev_of_jet[current], minlength=n_ev)

    ev_out = {
        "W": W, "Q2": Q2, "x": x, "y": y, "y_hcm": y_hcm, "gamma_hcm": gamma_hcm,
        "n_jets": n_jets.astype(np.int32), "n_current": n_current.astype(np.int32),
        "parton_id": parton_id.astype(np.int32),
    }
    jet_out = {
        "W": W[ev_of_jet], "Q2": Q2[ev_of_jet], "x": x[ev_of_jet], "y": y[ev_of_jet],
        "y_hcm": y_hcm[ev_of_jet], "parton_id": parton_id[ev_of_jet].astype(np.int32),
        "rank": rank.astype(np.int32), "plab": plab, "pt": pt, "eta": eta, "phi": phi,
        "e": J[:, 3], "mass": mass,
        "p_hcm": p_hcm, "e_hcm": J_hcm[:, 3], "cos_hcm": cos_hcm,
        "p_breit": p_breit, "pz_breit": pz_breit, "current": current,
        "n_const": n_const.astype(np.int32), "n_charged": n_charged.astype(np.int32),
        "n90": n90, "n90_hcm": n90_hcm, "n90_ch": n90_ch, "z_lead": z_lead,
        "dR_parton": dR_parton, "dphi_lepton": dphi_lepton,
    }
    return ev_out, jet_out


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def analyze(args):
    files = []
    for pattern in args.inputs:
        files.extend(sorted(glob.glob(pattern)) or [pattern])
    hists = make_histograms()
    ev_parts, jet_parts = [], []
    t0 = time.time()
    n_done = 0

    for fname in files:
        events = ak.from_parquet(fname)
        n = len(events)
        print(f"{fname}: {n:,} events", flush=True)
        for start in range(0, n, CHUNK):
            chunk = events[start:start + CHUNK]
            ev_out, jet_out = analyze_chunk(chunk, use_fastjet=not args.no_fastjet)
            fill_legacy(hists, ev_out, jet_out)
            ev_parts.append(ev_out)
            jet_parts.append(jet_out)
            n_done += len(chunk)
            rate = n_done / (time.time() - t0)
            print(f"  {n_done:>9,} events   {len(jet_out['W']):>7,} jets in chunk   "
                  f"({rate:.0f} ev/s)", flush=True)

    ev_all = {k: np.concatenate([p[k] for p in ev_parts]) for k in ev_parts[0]}
    jet_all = {k: np.concatenate([p[k] for p in jet_parts]) for k in jet_parts[0]}

    n_ev = len(ev_all["W"])
    n_jet = len(jet_all["W"])
    n_cur = int(jet_all["current"].sum())
    print(f"\n{n_ev:,} events, {n_jet:,} jets (|η|<{JET_ETA_MAX}, pT>{JET_PT_MIN} GeV), "
          f"{n_cur:,} in the Breit current hemisphere")
    print(f"Events with ≥1 current jet: {100 * (ev_all['n_current'] > 0).mean():.1f}%")

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    print(f"Writing {args.output} …", flush=True)
    with uproot.recreate(args.output) as f:
        f["jets"] = jet_all
        f["events"] = ev_all
        for name, h in hists.items():
            f[name] = h
    print(f"Done in {time.time() - t0:.0f} s.")


def parse_args():
    p = argparse.ArgumentParser(
        description="Analyze EIC DIS events for the FFS effect",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("inputs", nargs="+",
                   help="Input Parquet file(s) or glob pattern(s) from generate_events.py")
    p.add_argument("--output", type=str, default="data/analysis.root",
                   help="Output ROOT file")
    p.add_argument("--no-fastjet", action="store_true",
                   help="Use the simple cone fallback instead of FastJet anti-kT")
    return p.parse_args()


if __name__ == "__main__":
    analyze(parse_args())
