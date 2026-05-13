#!/usr/bin/env python3
"""
FFS Effect Analysis Script
===========================
Reads particle-level EIC DIS events (Parquet) produced by generate_events.py,
reconstructs jet observables, and measures the Frame-dependent Fragmentation
Shift (FFS) effect.

Primary observable
------------------
n₉₀ — the *fractional* minimum number of jet constituents (ordered by
decreasing 3-momentum magnitude) needed to account for 90% of the jet's
total momentum.  A per-jet linear interpolation is used so that n₉₀ is
a continuous quantity (e.g. 3.7 particles).  This exactly mirrors the n_x
observable defined in arXiv:2308.10951 Sec. 2.

The FFS prediction: at fixed lab-frame jet |p|, the colour rest frame
(the γ*p CM frame, characterised by W) is more boosted relative to the
lab at higher W.  Therefore the same lab-frame jet corresponds to a lower
CM-frame momentum at higher W → fewer particles needed to carry 90% of
the jet energy → ⟨n₉₀⟩ varies with W at fixed |p_lab|.

Secondary observable
--------------------
N_charged — simple count of charged jet constituents (IRC unsafe, but
widely used in practice as discussed in the paper).

Histograms are written to a ROOT file via uproot.

Usage
-----
    python analyze_events.py data/events.parquet [--output data/histograms.root]
    python analyze_events.py data/events.parquet --use-fastjet

Reference: arXiv:2308.10951  (Phys.Lett.B 866, 2025, 139561)
"""

import argparse
import sys
import os
import warnings

import numpy as np
import awkward as ak
import uproot
import hist

sys.path.insert(0, os.path.dirname(__file__))
from utils.dis_kinematics import DISKinematics, four_dot


# ---------------------------------------------------------------------------
# Configuration – bin edges
# ---------------------------------------------------------------------------

# W bins (GeV): study the FFS effect across the accessible EIC range
W_BINS = np.array([5.0, 10.0, 20.0, 30.0, 40.0, 55.0], dtype=float)

# Lab-frame jet momentum |p| bins (GeV)
P_LAB_BINS = np.array([2.0, 5.0, 10.0, 20.0], dtype=float)

# Charged-particle multiplicity axis (for 2D histograms)
NMAX = 60
N_BINS = np.arange(0, NMAX + 1, dtype=float)

# n_90 observable axis — fractional number of particles for 90% of jet |p|
# Values typically 1–10 for R=0.4 anti-kt DIS jets at EIC energies
N90_MAX  = 25.0
N90_BINS = 50

# Jet finding parameters
JET_R = 0.4          # anti-kt cone radius (matches paper: arXiv:2308.10951)
JET_ETA_MAX = 3.5    # EIC detector acceptance
JET_PT_MIN = 2.0     # GeV



# ---------------------------------------------------------------------------
# n_x observable (paper's primary FFS quantity)
# ---------------------------------------------------------------------------

def compute_n_x(const_pmags, threshold=0.90):
    """
    Compute fractional n_x: the minimum (fractional) number of jet
    constituents needed to account for *threshold* of the jet's total
    3-momentum magnitude.

    Algorithm (per arXiv:2308.10951 Sec. 2):
      1. Sort constituents by decreasing |p|.
      2. Cumulatively sum |p| values normalised to the jet total.
      3. Use linear interpolation to obtain a non-integer result.

    Parameters
    ----------
    const_pmags : array-like
        3-momentum magnitudes (|p|) of all jet constituents.
    threshold : float
        Fraction of total momentum to recover (default 0.90 for n₉₀).

    Returns
    -------
    float
        Fractional multiplicity in [0, N_constituents].
        Returns nan for empty or zero-momentum jets.

    Examples
    --------
    >>> compute_n_x([10., 5., 3., 1.])   # threshold=0.90
    # cumfrac = [0.526, 0.789, 0.947, 1.0] → need ~2.6 particles
    """
    pmags = np.asarray(const_pmags, dtype=float)
    if len(pmags) == 0:
        return np.nan
    pmags = np.sort(pmags)[::-1]          # sort descending
    total = pmags.sum()
    if total <= 0.0:
        return np.nan

    cumfrac = np.cumsum(pmags) / total

    # Index of first bin where cumulative fraction reaches threshold
    idx = int(np.searchsorted(cumfrac, threshold, side="left"))
    if idx >= len(cumfrac):
        return float(len(cumfrac))
    if idx == 0:
        # First constituent alone exceeds threshold
        return float(threshold / cumfrac[0])
    # Interpolate within the idx-th constituent
    frac = (threshold - cumfrac[idx - 1]) / (cumfrac[idx] - cumfrac[idx - 1])
    return float(idx + frac)


# ---------------------------------------------------------------------------
# Jet finding (anti-kt via fastjet, with lightweight fallback)
# ---------------------------------------------------------------------------

def _fastjet_jets(px_arr, py_arr, pz_arr, e_arr, charge_arr, R, ptmin, etamax):
    """Run anti-kt jet finding with FastJet."""
    import fastjet as fj

    pjets = []
    for i in range(len(px_arr)):
        pj = fj.PseudoJet(float(px_arr[i]), float(py_arr[i]),
                           float(pz_arr[i]), float(e_arr[i]))
        pj.set_user_index(i)
        pjets.append(pj)

    jdef = fj.JetDefinition(fj.antikt_algorithm, R)
    cs = fj.ClusterSequence(pjets, jdef)
    raw_jets = fj.sorted_by_pt(cs.inclusive_jets(ptmin))

    jets = []
    for jet in raw_jets:
        if abs(jet.eta()) > etamax:
            continue
        constituents = jet.constituents()
        idxs = [c.user_index() for c in constituents]
        n_ch = sum(1 for i in idxs if abs(charge_arr[i]) > 0)
        pmag = np.sqrt(jet.px()**2 + jet.py()**2 + jet.pz()**2)
        # Compute n_90: fractional particle count for 90% of jet |p|
        const_pmags = [
            np.sqrt(float(px_arr[i])**2 + float(py_arr[i])**2
                    + float(pz_arr[i])**2)
            for i in idxs
        ]
        n90 = compute_n_x(const_pmags, threshold=0.90)
        jets.append({
            "px": jet.px(), "py": jet.py(), "pz": jet.pz(), "e": jet.e(),
            "pt": jet.pt(), "eta": jet.eta(), "pmag": pmag,
            "n_const": len(constituents), "n_charged": n_ch, "n90": n90,
        })
    return jets


def _simple_jets(px_arr, py_arr, pz_arr, e_arr, charge_arr, ptmin, etamax):
    """
    Lightweight jet proxy when FastJet is unavailable.

    Groups particles into cones of radius R (≈ 0.4) by iterative nearest-
    neighbour clustering in (η, φ) space.  Intended for debugging / CI only;
    production runs should use FastJet.
    """
    n = len(px_arr)
    if n == 0:
        return []

    pt = np.sqrt(px_arr**2 + py_arr**2)
    p = np.sqrt(px_arr**2 + py_arr**2 + pz_arr**2)

    # Guard against zero-pt particles
    eps = 1e-9
    pt_safe = np.where(pt > eps, pt, eps)
    p_safe  = np.where(p  > eps, p,  eps)

    eta = np.arctanh(np.clip(pz_arr / p_safe, -1 + 1e-7, 1 - 1e-7))
    phi = np.arctan2(py_arr, px_arr)

    used = np.zeros(n, dtype=bool)
    jets = []

    # Process particles in descending pT order (greedy cone)
    order = np.argsort(-pt)
    for seed in order:
        if used[seed] or pt[seed] < ptmin or abs(eta[seed]) > etamax:
            continue

        deta = eta - eta[seed]
        dphi = phi - phi[seed]
        dphi = (dphi + np.pi) % (2 * np.pi) - np.pi
        dR   = np.sqrt(deta**2 + dphi**2)

        members = np.where((~used) & (dR < JET_R))[0]
        if len(members) == 0:
            continue
        used[members] = True

        jpx = float(np.sum(px_arr[members]))
        jpy = float(np.sum(py_arr[members]))
        jpz = float(np.sum(pz_arr[members]))
        je  = float(np.sum(e_arr[members]))
        jpt = float(np.sqrt(jpx**2 + jpy**2))
        jpmag = float(np.sqrt(jpx**2 + jpy**2 + jpz**2))
        jeta = float(np.arctanh(np.clip(jpz / max(jpmag, eps), -1 + 1e-7, 1 - 1e-7)))
        if jpt < ptmin or abs(jeta) > etamax:
            continue

        n_ch = int(np.sum(np.abs(charge_arr[members]) > 0))
        const_pmags = np.sqrt(px_arr[members]**2 + py_arr[members]**2
                              + pz_arr[members]**2)
        n90 = compute_n_x(const_pmags, threshold=0.90)
        jets.append({
            "px": jpx, "py": jpy, "pz": jpz, "e": je,
            "pt": jpt, "eta": jeta, "pmag": jpmag,
            "n_const": len(members), "n_charged": n_ch, "n90": n90,
        })
    return jets


def find_jets(px_arr, py_arr, pz_arr, e_arr, charge_arr,
              R=JET_R, ptmin=JET_PT_MIN, etamax=JET_ETA_MAX,
              use_fastjet=True):
    """Find jets; fall back to simple algorithm if FastJet is unavailable."""
    px_arr = np.asarray(px_arr, dtype=float)
    py_arr = np.asarray(py_arr, dtype=float)
    pz_arr = np.asarray(pz_arr, dtype=float)
    e_arr  = np.asarray(e_arr,  dtype=float)
    charge_arr = np.asarray(charge_arr, dtype=int)

    if use_fastjet:
        try:
            return _fastjet_jets(px_arr, py_arr, pz_arr, e_arr, charge_arr,
                                 R, ptmin, etamax)
        except ImportError:
            warnings.warn(
                "fastjet not available; using simplified cone algorithm. "
                "Install with: conda install -c conda-forge fastjet",
                stacklevel=2,
            )
    return _simple_jets(px_arr, py_arr, pz_arr, e_arr, charge_arr,
                        ptmin, etamax)


# ---------------------------------------------------------------------------
# Histogram definitions (using hist / boost-histogram)
# ---------------------------------------------------------------------------

def make_histograms():
    """Create and return all analysis histograms."""
    W_ax  = hist.axis.Variable(W_BINS,    name="W",    label=r"$W$ [GeV]")
    PL_ax = hist.axis.Variable(P_LAB_BINS, name="plab", label=r"$|p|_{\rm lab}$ [GeV]")
    N_ax  = hist.axis.Regular(NMAX, 0, NMAX, name="N",  label=r"$N_{\rm charged}$")
    N90_ax = hist.axis.Regular(N90_BINS, 0, N90_MAX, name="n90", label=r"$n_{90}$")

    Q2_ax   = hist.axis.Regular(50, 1,  1000, name="Q2",  label=r"$Q^2$ [GeV$^2$]",
                                transform=hist.axis.transform.log)
    x_ax    = hist.axis.Regular(50, 1e-4, 1, name="x",   label=r"Bjorken $x$",
                                transform=hist.axis.transform.log)
    y_ax    = hist.axis.Regular(50, 0.01, 1, name="y",   label=r"Inelasticity $y$")
    W_fine  = hist.axis.Regular(50, 5, 55,  name="W_fine", label=r"$W$ [GeV]")

    hists = {
        # ── Primary FFS observable: n₉₀ (arXiv:2308.10951, primary figure) ──
        # 3D: (W bin) × (|p_lab| bin) × (n₉₀)
        "n90_3d": hist.Hist(W_ax, PL_ax, N90_ax, storage=hist.storage.Double()),

        # Profile: ⟨n₉₀⟩ as function of W (fine binning), per p_lab bin
        # Stored as two histograms: sum_n90 and count (mean = sum / count)
        "sum_n90_vs_W":   hist.Hist(W_fine, PL_ax, storage=hist.storage.Double()),
        "count_n90_vs_W": hist.Hist(W_fine, PL_ax, storage=hist.storage.Double()),

        # ── Secondary FFS observable: N_charged (IRC unsafe, but common) ───
        # 3D: (W bin) × (|p_lab| bin) × (N_charged)
        "mult_3d": hist.Hist(W_ax, PL_ax, N_ax, storage=hist.storage.Double()),

        # Profile: ⟨N_charged⟩ as function of W (fine binning), per p_lab bin
        "sum_N_vs_W": hist.Hist(W_fine, PL_ax, storage=hist.storage.Double()),
        "count_vs_W": hist.Hist(W_fine, PL_ax, storage=hist.storage.Double()),

        # ── DIS kinematics cross-checks ─────────────────────────────────────
        "Q2":  hist.Hist(Q2_ax,  storage=hist.storage.Double()),
        "x":   hist.Hist(x_ax,   storage=hist.storage.Double()),
        "y":   hist.Hist(y_ax,   storage=hist.storage.Double()),
        "W":   hist.Hist(W_fine, storage=hist.storage.Double()),

        # 2D: Q² vs W (kinematic plane)
        "Q2_vs_W": hist.Hist(
            hist.axis.Regular(50, 5, 55,   name="W2",  label=r"$W$ [GeV]"),
            hist.axis.Regular(50, 1, 1000, name="Q22", label=r"$Q^2$ [GeV$^2$]",
                              transform=hist.axis.transform.log),
            storage=hist.storage.Double(),
        ),

        # 2D: jet η vs pT (lab frame)
        "jet_eta_pt": hist.Hist(
            hist.axis.Regular(60, -4, 4,    name="eta", label=r"Jet $\eta_{\rm lab}$"),
            hist.axis.Regular(50,  0,  30,  name="pt",  label=r"Jet $p_T$ [GeV]"),
            storage=hist.storage.Double(),
        ),

        # 1D: jet multiplicity per event (all accepted jets)
        "n_jets": hist.Hist(
            hist.axis.Regular(20, 0, 20, name="nj", label=r"$N_{\rm jets}$ per event"),
            storage=hist.storage.Double(),
        ),
    }
    return hists


# ---------------------------------------------------------------------------
# Main analysis loop
# ---------------------------------------------------------------------------

def analyze(args):
    print(f"Reading events from  {args.input}", flush=True)
    events = ak.from_parquet(args.input)
    n_events = len(events)
    print(f"Loaded {n_events:,} events", flush=True)

    hists = make_histograms()

    # Counters for diagnostics
    n_with_jets = 0
    n_jet_total = 0

    for i in range(n_events):
        ev = events[i]

        W  = float(ev["W"])
        Q2 = float(ev["Q2"])
        x  = float(ev["x"])
        y  = float(ev["y"])

        # Fill DIS kinematics histograms
        if 5.0 <= W <= 55.0:
            hists["W"].fill(W_fine=W)
        hists["Q2"].fill(Q2=Q2)
        hists["x"].fill(x=x)
        hists["y"].fill(y=y)
        hists["Q2_vs_W"].fill(W2=np.clip(W, 5.01, 54.99),
                               Q22=np.clip(Q2, 1.01, 999.9))

        px  = np.asarray(ev["particles"]["px"])
        py  = np.asarray(ev["particles"]["py"])
        pz  = np.asarray(ev["particles"]["pz"])
        e   = np.asarray(ev["particles"]["e"])
        charge = np.asarray(ev["particles"]["charge"], dtype=int)

        jets = find_jets(px, py, pz, e, charge,
                         use_fastjet=args.use_fastjet)

        # ── Jet-level observables ──────────────────────────────────────────
        n_accepted = 0
        for jet in jets:
            jpt   = jet["pt"]
            jeta  = jet["eta"]
            jpmag = jet["pmag"]
            n_ch  = jet["n_charged"]
            n90   = jet["n90"]

            hists["jet_eta_pt"].fill(eta=np.clip(jeta, -3.99, 3.99),
                                     pt=np.clip(jpt, 0.01, 29.99))

            # FFS observables: fill when in the analysis kinematic range
            if 5.0 <= W <= 55.0 and P_LAB_BINS[0] <= jpmag < P_LAB_BINS[-1]:
                # ── Primary observable: n₉₀ ───────────────────────────────
                if np.isfinite(n90):
                    n90_clipped = min(float(n90), N90_MAX - N90_MAX / N90_BINS)
                    hists["n90_3d"].fill(W=W, plab=jpmag, n90=n90_clipped)
                    hists["sum_n90_vs_W"].fill(
                        W_fine=W, plab=jpmag, weight=float(n90))
                    hists["count_n90_vs_W"].fill(W_fine=W, plab=jpmag)

                # ── Secondary observable: N_charged ───────────────────────
                n_ch_clipped = min(n_ch, NMAX - 1)
                hists["mult_3d"].fill(W=W, plab=jpmag, N=n_ch_clipped)
                hists["sum_N_vs_W"].fill(
                    W_fine=W, plab=jpmag, weight=float(n_ch))
                hists["count_vs_W"].fill(W_fine=W, plab=jpmag)

            n_accepted += 1

        hists["n_jets"].fill(nj=min(n_accepted, 19))
        if n_accepted > 0:
            n_with_jets += 1
            n_jet_total += n_accepted

        if (i + 1) % 50_000 == 0:
            print(f"  Processed {i+1:>8,} / {n_events:,} events  "
                  f"[jets/ev = {n_jet_total/(i+1):.2f}]", flush=True)

    print(f"\nEvents with ≥1 jet: {n_with_jets}/{n_events} "
          f"({100*n_with_jets/n_events:.1f}%)")
    print(f"Mean jets per event: {n_jet_total/n_events:.2f}")

    # ── Write ROOT file ────────────────────────────────────────────────────
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    print(f"\nWriting histograms to {args.output} …", flush=True)

    with uproot.recreate(args.output) as f:
        for name, h in hists.items():
            f[name] = h

    print("Done.")
    return hists


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Analyze EIC DIS events for the FFS effect",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("input",  type=str,
                   help="Input Parquet file from generate_events.py")
    p.add_argument("--output", type=str, default="data/histograms.root",
                   help="Output ROOT file path")
    p.add_argument("--use-fastjet", action="store_true",
                   help="Use FastJet for jet finding (strongly recommended)")
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    analyze(args)
