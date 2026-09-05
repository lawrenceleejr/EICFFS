#!/usr/bin/env python3
"""
FFS figures — Tufte edition
===========================
Reads the per-jet TTree written by analyze_events.py and draws one figure per
PDF, following Tufte: no gridlines, no boxes, no legends where a direct label
will do, range frames, muted ink, small multiples as separate files.

Figures (all in --outdir, PDF + PNG preview)
--------------------------------------------
ffs_fan              <n90> vs |p|_lab, one line per W slice — the EIC FFS effect
ffs_collapse         the same jets vs their momentum in the colour rest frame (γ*p)
ffs_ratio            <n90>(W) / <n90>(W_ref) vs |p|_lab — size of the shift
ffs_slopegraph       low-W → high-W shift of <n90> per |p|_lab bin
ffs_distribution     n90 quantile bands vs W at fixed |p|_lab
ffs_wq_table         <n90> over the (W, Q) plane at fixed |p|_lab — which scale does the jet know?
nch_fan              <N_charged> vs |p|_lab per W slice (secondary observable)
boost_map            rapidity of the colour rest frame across the (x, Q²) plane, with the jets
plateau              charged-hadron rapidity density in the γ*p frame per W slice (needs --events)

Usage
-----
    python make_figures.py data/analysis.root --outdir figures/
    python make_figures.py data/analysis*.root --events data/events_1.parquet --outdir figures/
"""

import argparse
import glob
import os
import sys

import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
plt.style.use(os.path.join(HERE, "utils", "tufte.mplstyle"))

# ---------------------------------------------------------------------------
# Selection and binning
# ---------------------------------------------------------------------------

E_LEPTON, E_PROTON = 10.0, 100.0
SQRT_S = np.sqrt(4 * E_LEPTON * E_PROTON)

W_SLICES = [(10, 15), (15, 22), (22, 32), (32, 45), (45, 60)]
P_EDGES = np.array([2.0, 3.0, 4.5, 7.0, 10.0, 15.0, 22.0, 33.0])
P_WINDOW = (4.5, 10.0)          # |p|_lab window used for the fixed-momentum figures
MIN_ENTRIES = 40

# Single-hue sequential ramp for the ordered variable W (light = low W)
W_COLORS = ["#9ecae1", "#6baed6", "#3182bd", "#08519c", "#08306b"]
# Sequential ramp for the ordered variable E_cm (six steps, light = low energy)
E_COLORS = ["#c6dbef", "#9ecae1", "#6baed6", "#3182bd", "#08519c", "#08306b"]
E_SLICES = [(1.5, 2.5), (2.5, 4), (4, 6), (6, 9), (9, 14), (14, 22)]
PCM_SLICES = [(2, 3.3), (3.3, 5), (5, 7.5), (7.5, 11), (11, 16), (16, 24)]
PT_EDGES_CM = np.array([0.6, 1.0, 1.5, 2.0, 2.7, 3.5, 4.5, 6.0, 8.0, 11.0])
PT_EDGES_LAB = np.array([2.0, 2.5, 3.2, 4.0, 5.0, 6.5, 8.0, 10.0, 13.0])


def e_label(lo, hi):
    return rf"$E_{{\rm cm}}$ = {lo:g}$-${hi:g} GeV"
INK = "#1f1f1f"
MUTED = "#8a8a8a"
FAINT = "#d9d9d9"
ACCENT = "#c44e52"


def w_label(lo, hi):
    return rf"$W$ = {lo}$-${hi} GeV"


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------

def load_jets(paths):
    files = []
    for p in paths:
        files.extend(sorted(glob.glob(p)) or [p])
    cols = ["W", "Q2", "x", "y", "plab", "pt", "eta", "p_hcm", "current",
            "n_const", "n_charged", "n90", "n90_hcm", "n90_ch", "z_lead", "y_hcm"]
    parts = [uproot.open(f)["jets"].arrays(cols, library="np") for f in files]
    jets = {k: np.concatenate([p[k] for p in parts]) for k in cols}
    sel = jets["current"] & np.isfinite(jets["n90"])
    return {k: v[sel] for k, v in jets.items()}


def load_cmjets(paths):
    """Jets clustered in the colour rest frame (tree ``cmjets``)."""
    files = []
    for p in paths:
        files.extend(sorted(glob.glob(p)) or [p])
    cols = ["W", "Q2", "x", "y", "e_cm", "p_cm", "cos_cm", "pt_lab", "p_lab",
            "eta_lab", "n_const", "n_charged", "n90_cm", "n90_labmom"]
    parts = [uproot.open(f)["cmjets"].arrays(cols, library="np") for f in files]
    cj = {k: np.concatenate([p[k] for p in parts]) for k in cols}
    ok = np.isfinite(cj["n90_cm"])
    return {k: v[ok] for k, v in cj.items()}


def profile(xvals, yvals, edges, min_entries=MIN_ENTRIES):
    """Mean, standard error and bin centres (geometric) of y in bins of x."""
    idx = np.digitize(xvals, edges) - 1
    n = len(edges) - 1
    mean = np.full(n, np.nan)
    err = np.full(n, np.nan)
    for i in range(n):
        m = idx == i
        if m.sum() >= min_entries:
            v = yvals[m]
            mean[i] = v.mean()
            err[i] = v.std(ddof=1) / np.sqrt(len(v))
    centres = np.sqrt(edges[:-1] * edges[1:])
    return centres, mean, err


def range_frame(ax, x=None, y=None):
    """Tufte range frame: spines span the data range only."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if x is not None and np.isfinite(x).any():
        ax.spines["bottom"].set_bounds(np.nanmin(x), np.nanmax(x))
    if y is not None and np.isfinite(y).any():
        ax.spines["left"].set_bounds(np.nanmin(y), np.nanmax(y))


def label_line_end(ax, x, y, text, color, dy=0, fontsize=8):
    """Direct label at the right end of a line; dy is a vertical offset in points."""
    ok = np.isfinite(y)
    if not ok.any():
        return
    xe, ye = x[ok][-1], y[ok][-1]
    ax.annotate(text, (xe, ye), xytext=(6, dy), textcoords="offset points",
                va="center", ha="left", fontsize=fontsize, color=color)


class EndLabels:
    """Collect direct labels for line ends and draw them without collisions."""

    def __init__(self, ax, min_sep_pt=9.0, fontsize=8):
        self.ax, self.items, self.min_sep, self.fs = ax, [], min_sep_pt, fontsize

    def add(self, x, y, text, color):
        ok = np.isfinite(np.asarray(y))
        if ok.any():
            self.items.append((np.asarray(x)[ok][-1], np.asarray(y)[ok][-1], text, color))

    def draw(self):
        if not self.items:
            return
        fig = self.ax.figure
        fig.canvas.draw()                     # axes transforms are final now
        pts = np.array([self.ax.transData.transform((x, y)) for x, y, _, _ in self.items])
        disp_y = pts[:, 1] * 72.0 / fig.dpi   # points
        order = np.argsort(disp_y)
        target = disp_y[order].copy()
        for i in range(1, len(target)):       # push upwards to keep min separation
            target[i] = max(target[i], target[i - 1] + self.min_sep)
        shift = 0.5 * (target[-1] - disp_y[order][-1])  # recentre the stack
        target -= shift * (len(target) > 1)
        for k, idx in enumerate(order):
            x, y, text, color = self.items[idx]
            dy = target[k] - disp_y[idx]
            self.ax.annotate(text, (x, y), xytext=(6, dy), textcoords="offset points",
                             va="center", ha="left", fontsize=self.fs, color=color)


def caption(ax, text, fontsize=7.5):
    """Small explanatory note under the axes — the figure explains itself."""
    ax.text(0.0, -0.22, text, transform=ax.transAxes, fontsize=fontsize,
            color=MUTED, va="top", ha="left", wrap=True)


def save(fig, outdir, name):
    os.makedirs(outdir, exist_ok=True)
    pdf = os.path.join(outdir, name + ".pdf")
    fig.savefig(pdf, format="pdf", bbox_inches="tight")
    fig.savefig(os.path.join(outdir, name + ".png"), dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f"  → {pdf}")


def eic_note_cm():
    return (rf"Pythia 8 NC DIS, {E_LEPTON:.0f}$\times${E_PROTON:.0f} GeV, angular anti-$k_T$ "
            rf"$R$ = 0.4 in the $\gamma^*p$ frame, $E_{{\rm cm}} > 1$ GeV, current hemisphere")


def eic_note():
    return (rf"Pythia 8 NC DIS, {E_LEPTON:.0f}$\times${E_PROTON:.0f} GeV, anti-$k_T$ $R$ = 0.4, "
            rf"$p_T^{{\rm lab}} > 2$ GeV, $|\eta| < 3.5$, Breit current hemisphere")


# ---------------------------------------------------------------------------
# Figure 1: the fan  — <n90> vs |p|_lab per W slice
# ---------------------------------------------------------------------------

def fig_fan(jets, outdir, obs="n90", name="ffs_fan",
            ylabel=r"$\langle n_{90}\rangle$"):
    fig, ax = plt.subplots(figsize=(4.6, 3.4))
    allx, ally = [], []
    labels = EndLabels(ax)
    for (lo, hi), col in zip(W_SLICES, W_COLORS):
        m = (jets["W"] >= lo) & (jets["W"] < hi)
        xc, mu, se = profile(jets["plab"][m], jets[obs][m], P_EDGES)
        ok = np.isfinite(mu)
        ax.errorbar(xc[ok], mu[ok], yerr=se[ok], color=col, lw=1.2, elinewidth=0.6,
                    capsize=0, marker="o", ms=3, mec="white", mew=0.4)
        labels.add(xc[ok], mu[ok], w_label(lo, hi), col)
        allx.append(xc[ok]); ally.append(mu[ok])
    allx, ally = np.concatenate(allx), np.concatenate(ally)
    ax.set_xscale("log")
    ax.set_xticks([2, 3, 5, 10, 20, 30])
    ax.set_xticklabels(["2", "3", "5", "10", "20", "30"])
    ax.minorticks_off()
    ax.set_xlabel(r"jet momentum in the lab, $|\vec p|_{\rm lab}$  [GeV]")
    ax.set_ylabel(ylabel)
    ax.set_xlim(1.8, 60)
    range_frame(ax, allx, ally)
    labels.draw()
    ylim = ax.get_ylim()
    caption(ax, "Same lab-frame jet momentum, different colour rest frames.  "
                "The colour-connected system in DIS has invariant mass $W$ and moves "
                "with rapidity $y_{\\rm cm}$ that falls with $W$; at fixed "
                "$|\\vec p|_{\\rm lab}$ a jet from a higher-$W$ event is harder in that "
                "frame and fragments into more particles.\n" + eic_note())
    save(fig, outdir, name)
    return ylim


# ---------------------------------------------------------------------------
# Figure 2: the collapse — same jets vs momentum in the colour rest frame
# ---------------------------------------------------------------------------

def fig_collapse(jets, outdir, obs="n90", name="ffs_collapse",
                 ylabel=r"$\langle n_{90}\rangle$", ylim=None):
    fig, ax = plt.subplots(figsize=(4.6, 3.4))
    edges = np.array([1.0, 1.5, 2.2, 3.3, 5.0, 7.5, 11.0, 16.0, 24.0])
    allx, ally = [], []
    labels = EndLabels(ax)
    for (lo, hi), col in zip(W_SLICES, W_COLORS):
        m = (jets["W"] >= lo) & (jets["W"] < hi)
        xc, mu, se = profile(jets["p_hcm"][m], jets[obs][m], edges, min_entries=200)
        ok = np.isfinite(mu)
        ax.errorbar(xc[ok], mu[ok], yerr=se[ok], color=col, lw=1.2, elinewidth=0.6,
                    capsize=0, marker="o", ms=3, mec="white", mew=0.4)
        labels.add(xc[ok], mu[ok], w_label(lo, hi), col)
        allx.append(xc[ok]); ally.append(mu[ok])
    allx, ally = np.concatenate(allx), np.concatenate(ally)
    ax.set_xscale("log")
    ax.set_xticks([1, 2, 3, 5, 10, 20])
    ax.set_xticklabels(["1", "2", "3", "5", "10", "20"])
    ax.minorticks_off()
    ax.set_xlabel(r"jet momentum in the colour rest frame ($\gamma^*p$ CM), $|\vec p|_{\rm cm}$  [GeV]")
    ax.set_ylabel(ylabel)
    ax.set_xlim(0.9, 40)
    if ylim is not None:
        ax.set_ylim(*ylim)
    range_frame(ax, allx, ally)
    labels.draw()
    caption(ax, "The same jets, now placed at their momentum in the frame where the "
                "colour-connected system is at rest, on the same vertical scale as the "
                "lab-frame figure.  The $W$ slices agree where they overlap, and the "
                "dependence on colour-frame momentum is weak: most of the lab-frame "
                "spread is the boost, not the jet.")
    save(fig, outdir, name)


# ---------------------------------------------------------------------------
# The flatness test: which frame organises fragmentation?
# ---------------------------------------------------------------------------

def _flat_panel(ax, slice_var, x, y, slices, colors, edges, label_fn,
                label_extremes_only=False):
    """Profile y against x within slices of slice_var; draw and measure flatness."""
    labels = EndLabels(ax)
    drawn = []
    allx, ally, spreads = [], [], {}
    for (lo, hi), col in zip(slices, colors):
        m = (slice_var >= lo) & (slice_var < hi)
        xc, mu, se = profile(x[m], y[m], edges, min_entries=150)
        ok = np.isfinite(mu)
        if ok.sum() < 3:
            continue
        ax.errorbar(xc[ok], mu[ok], yerr=se[ok], color=col, lw=1.2, elinewidth=0.6,
                    capsize=0, marker="o", ms=3, mec="white", mew=0.4)
        drawn.append((xc[ok], mu[ok], label_fn(lo, hi), col))
        spreads[(lo, hi)] = 100 * (mu[ok].max() - mu[ok].min()) / mu[ok].mean()
        allx.append(xc[ok])
        ally.append(mu[ok])
    for k, (xx, yy, txt, col) in enumerate(drawn):
        if label_extremes_only and 0 < k < len(drawn) - 1:
            continue
        labels.add(xx, yy, txt, col)
    return labels, spreads, np.concatenate(allx), np.concatenate(ally)


def fig_flat_cm(cj, outdir):
    """n90 of colour-frame jets against lab pT, sliced in colour-frame energy."""
    fig, ax = plt.subplots(figsize=(4.6, 3.4))
    labels, spreads, allx, ally = _flat_panel(
        ax, cj["e_cm"], cj["pt_lab"], cj["n90_cm"], E_SLICES, E_COLORS,
        PT_EDGES_CM, e_label)
    ax.set_xscale("log")
    ax.set_xticks([0.7, 1, 2, 3, 5, 8])
    ax.set_xticklabels(["0.7", "1", "2", "3", "5", "8"])
    ax.minorticks_off()
    ax.set_xlabel(r"jet transverse momentum in the lab, $p_T^{\rm lab}$  [GeV]")
    ax.set_ylabel(r"$\langle n_{90}\rangle$ in the colour rest frame")
    ax.set_xlim(0.55, 24)
    range_frame(ax, allx, ally)
    labels.draw()
    worst = max(spreads.values())
    caption(ax, "Jets clustered in the colour rest frame with an angular algorithm, then "
                "labelled by their energy in that frame.  Every line is flat: across a factor "
                f"of ten in lab transverse momentum no slice varies by more than {worst:.0f}%.  "
                "This is the frame that organises fragmentation.\n" + eic_note_cm())
    save(fig, outdir, "flat_cmjets")
    return ax.get_ylim(), spreads


def fig_flat_lab(jets, outdir, ylim=None):
    """The same test for lab-clustered jets, sliced in colour-frame momentum."""
    fig, ax = plt.subplots(figsize=(4.6, 3.4))
    labels, spreads, allx, ally = _flat_panel(
        ax, jets["p_hcm"], jets["pt"], jets["n90"], PCM_SLICES, E_COLORS,
        PT_EDGES_LAB, lambda lo, hi: rf"$|\vec p|_{{\rm cm}}$ = {lo:g}$-${hi:g} GeV",
        label_extremes_only=True)
    ax.set_xscale("log")
    ax.set_xticks([2, 3, 5, 8, 12])
    ax.set_xticklabels(["2", "3", "5", "8", "12"])
    ax.minorticks_off()
    ax.set_xlabel(r"jet transverse momentum in the lab, $p_T^{\rm lab}$  [GeV]")
    ax.set_ylabel(r"$\langle n_{90}\rangle$ of lab-frame jets")
    ax.set_xlim(1.8, 34)
    if ylim is not None:
        ax.set_ylim(*ylim)
    range_frame(ax, allx, ally)
    labels.draw()
    ax.annotate("four more slices between,\nindistinguishable", (6.0, 1.85),
                fontsize=7, color=MUTED, ha="left", va="top")
    worst = max(spreads.values())
    caption(ax, "The identical test for jets clustered in the lab with anti-$k_T$ $R$ = 0.4, "
                "carrying the same colour-frame labels and drawn on the same vertical scale.  "
                f"The slices collapse onto one rising curve and span {worst:.0f}%: a lab-frame "
                "jet definition measures the lab cone, not the fragmenting system.\n" + eic_note())
    save(fig, outdir, "flat_labjets")
    return spreads


def fig_universal_cm(cj, outdir):
    """n90 against colour-frame energy, sliced in lab pT: a single curve."""
    fig, ax = plt.subplots(figsize=(4.6, 3.4))
    edges = np.array([1.5, 2.2, 3.2, 4.5, 6.5, 9.0, 13.0, 18.0, 25.0])
    slices = [(0.6, 1.2), (1.2, 2.0), (2.0, 3.2), (3.2, 5.0), (5.0, 9.0)]
    labels, spreads, allx, ally = _flat_panel(
        ax, cj["pt_lab"], cj["e_cm"], cj["n90_cm"], slices, W_COLORS, edges,
        lambda lo, hi: rf"$p_T^{{\rm lab}}$ = {lo:g}$-${hi:g} GeV",
        label_extremes_only=True)
    ax.set_xscale("log")
    ax.set_xticks([2, 3, 5, 8, 12, 20])
    ax.set_xticklabels(["2", "3", "5", "8", "12", "20"])
    ax.minorticks_off()
    ax.set_xlabel(r"jet energy in the colour rest frame, $E_{\rm cm}$  [GeV]")
    ax.set_ylabel(r"$\langle n_{90}\rangle$ in the colour rest frame")
    ax.set_xlim(1.4, 50)
    range_frame(ax, allx, ally)
    labels.draw()
    ax.annotate("three more slices between,\nall within a few percent", (5.0, 2.15),
                fontsize=7, color=MUTED, ha="left", va="top")
    caption(ax, "The converse view of the same jets: against colour-frame energy the lab-momentum "
                "slices fall on one curve, a factor of ten in lab momentum collapsed onto a single "
                "line.  Measured in the right frame, one variable carries the whole dependence.")
    save(fig, outdir, "universal_cm")


def fig_pt_fan(jets, outdir):
    """Lab jets against lab pT, sliced in W: the fan nearly closes."""
    fig, ax = plt.subplots(figsize=(4.6, 3.2))
    labels, spreads, allx, ally = _flat_panel(
        ax, jets["W"], jets["pt"], jets["n90"], W_SLICES, W_COLORS,
        PT_EDGES_LAB, w_label)
    ax.set_xscale("log")
    ax.set_xticks([2, 3, 5, 8, 12])
    ax.set_xticklabels(["2", "3", "5", "8", "12"])
    ax.minorticks_off()
    ax.set_xlabel(r"$p_T^{\rm lab}$  [GeV]")
    ax.set_ylabel(r"$\langle n_{90}\rangle$")
    ax.set_xlim(1.8, 30)
    range_frame(ax, allx, ally)
    labels.draw()
    caption(ax, "Why lab transverse momentum is the honest abscissa for a lab jet: at fixed "
                "$p_T^{\\rm lab}$ the $W$ slices differ by a few percent.  The wide fan against "
                "total $|\\vec p|_{\\rm lab}$ is largely a change of jet angle, since a "
                "fixed-$|\\vec p|$ jet at higher $W$ is more central and therefore harder in "
                "$p_T$.")
    save(fig, outdir, "pt_fan")


# ---------------------------------------------------------------------------
# Radius scan: how much of the current system does a lab cone hold?
# ---------------------------------------------------------------------------

R_TREES = [("jets_R0p4", 0.4), ("jets_R0p8", 0.8), ("jets_R1p2", 1.2),
           ("jets_R1p6", 1.6), ("jets_R2p4", 2.4)]
R_COLORS = ["#c6dbef", "#9ecae1", "#6baed6", "#3182bd", "#08519c"]
Q_WINDOW = (5.0, 7.5)          # GeV, the Q band used for the fixed-Q figures
QUOTE_PT = (1.9, 7.3)          # GeV, pT range both jet definitions populate
PT_EDGES_Q = np.array([1.0, 1.4, 1.9, 2.5, 3.3, 4.3, 5.6, 7.3, 9.5])


def load_radius_trees(paths):
    """The lab radius scan and the hemisphere, keyed by tree name."""
    files = []
    for p in paths:
        files.extend(sorted(glob.glob(p)) or [p])
    out = {}
    for name, _ in R_TREES:
        cols = ["W", "Q2", "pt", "plab", "e_hcm", "current", "captured", "lead",
                "n_const", "n90"]
        parts = [uproot.open(f)[name].arrays(cols, library="np") for f in files]
        d = {k: np.concatenate([p[k] for p in parts]) for k in cols}
        sel = d["lead"] & np.isfinite(d["n90"])
        out[name] = {k: v[sel] for k, v in d.items()}
    cols = ["W", "Q2", "pt", "plab", "e_hcm", "n_const", "n90", "n90_cm"]
    parts = [uproot.open(f)["hemisphere"].arrays(cols, library="np") for f in files]
    d = {k: np.concatenate([p[k] for p in parts]) for k in cols}
    sel = np.isfinite(d["n90"]) & (d["pt"] > 1.0)
    out["hemisphere"] = {k: v[sel] for k, v in d.items()}
    return out


def fig_capture(trees, outdir):
    """Fraction of the current hemisphere held by the leading lab jet."""
    fig, ax = plt.subplots(figsize=(4.6, 3.2))
    labels = EndLabels(ax)
    allx, ally = [], []
    qlo, qhi = Q_WINDOW
    for (name, R), col in zip(R_TREES, R_COLORS):
        d = trees[name]
        Q = np.sqrt(d["Q2"])
        m = (Q >= qlo) & (Q < qhi)
        xc, mu, se = profile(d["pt"][m], d["captured"][m], PT_EDGES_Q, min_entries=150)
        ok = np.isfinite(mu)
        if ok.sum() < 3:
            continue
        ax.errorbar(xc[ok], mu[ok], yerr=se[ok], color=col, lw=1.2, elinewidth=0.6,
                    capsize=0, marker="o", ms=3, mec="white", mew=0.4)
        labels.add(xc[ok], mu[ok], rf"$R$ = {R}", col)
        allx.append(xc[ok]); ally.append(mu[ok])
    ax.axhline(1.0, color=FAINT, lw=0.8, zorder=0)
    ax.annotate("the whole current system", (1.05, 1.02), fontsize=7, color=MUTED, va="bottom")
    allx, ally = np.concatenate(allx), np.concatenate(ally)
    ax.set_xscale("log")
    ax.set_xticks([1, 2, 3, 5, 8]); ax.set_xticklabels(["1", "2", "3", "5", "8"])
    ax.minorticks_off()
    ax.set_xlabel(r"$p_T^{\rm lab}$  [GeV]")
    ax.set_ylabel("fraction of the current system\nheld by the jet")
    ax.set_xlim(0.9, 16)
    range_frame(ax, allx, ally)
    labels.draw()
    caption(ax, rf"Leading current jet, $Q$ = {qlo:g}$-${qhi:g} GeV.  A fixed lab cone does not "
                "hold a fixed piece of the shower: at $R$ = 0.4 it slides from about half the "
                "current system to all of it as the jet gets harder in the lab.  Above unity the "
                "cone is also sweeping in the target side.  Only near $R$ = 2.4 does the "
                "fraction stop depending on $p_T$.")
    save(fig, outdir, "capture_fraction")


def _fixed_q_profiles(d, edges, slices=E_SLICES, quote=QUOTE_PT):
    """Profiles of n90 vs lab pT within a fixed Q window, one per E_cm slice."""
    qlo, qhi = Q_WINDOW
    Q = np.sqrt(d["Q2"])
    inq = (Q >= qlo) & (Q < qhi)
    out, spreads, slopes = [], [], []
    for (lo, hi), col in zip(slices, E_COLORS):
        m = inq & (d["e_hcm"] >= lo) & (d["e_hcm"] < hi)
        xc, mu, se = profile(d["pt"][m], d["n90"][m], edges, min_entries=150)
        ok = np.isfinite(mu)
        if ok.sum() < 3:
            continue
        out.append((xc[ok], mu[ok], se[ok], e_label(lo, hi), col))
        q = ok & (xc >= quote[0]) & (xc <= quote[1])
        if q.sum() >= 3:
            spreads.append(100 * (mu[q].max() - mu[q].min()) / mu[q].mean())
            slopes.append(np.polyfit(np.log(xc[q]), np.log(mu[q]), 1)[0])
    return out, spreads, slopes


def _draw_fixed_q(profiles, ylabel, span, name, outdir, caption_text):
    fig, ax = plt.subplots(figsize=(4.6, 3.4))
    labels = EndLabels(ax)
    allv = []
    for xc, mu, se, txt, col in profiles:
        ax.errorbar(xc, mu, yerr=se, color=col, lw=1.2, elinewidth=0.6, capsize=0,
                    marker="o", ms=3, mec="white", mew=0.4)
        labels.add(xc, mu, txt, col)
        allv.append(mu)
    allv = np.concatenate(allv)
    gm = np.sqrt(allv.min() * allv.max())
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim(gm / span, gm * span)                 # identical log span in both figures
    ax.set_xticks([1, 2, 3, 5, 8]); ax.set_xticklabels(["1", "2", "3", "5", "8"])
    yt = [t for t in (1.5, 2, 2.5, 3, 4, 5, 6) if gm / span <= t <= gm * span]
    ax.set_yticks(yt); ax.set_yticklabels([f"{t:g}" for t in yt])
    ax.minorticks_off()
    ax.set_xlabel(r"$p_T^{\rm lab}$  [GeV]")
    ax.set_ylabel(ylabel)
    ax.set_xlim(0.9, 18)
    range_frame(ax, np.concatenate([p[0] for p in profiles]), allv)
    labels.draw()
    caption(ax, caption_text)
    save(fig, outdir, name)


def fig_fixed_q(trees, outdir):
    """At fixed Q: a small cone stays sloped, the whole current system goes flat."""
    qlo, qhi = Q_WINDOW
    p_r04, sp_r04, sl_r04 = _fixed_q_profiles(trees["jets_R0p4"], PT_EDGES_Q)
    p_hem, sp_hem, sl_hem = _fixed_q_profiles(trees["hemisphere"], PT_EDGES_Q)
    # one log span for both panels, so equal visual slope means equal fractional change
    span = 1.02 * max(np.sqrt(np.concatenate([p[1] for p in ps]).max()
                              / np.concatenate([p[1] for p in ps]).min())
                      for ps in (p_r04, p_hem))
    med_r, med_h = float(np.median(sp_r04)), float(np.median(sp_hem))
    _draw_fixed_q(p_r04, r"$\langle n_{90}\rangle$,  lab jets $R$ = 0.4", span,
                  "fixed_q_R04", outdir,
                  rf"Leading lab jet, $R$ = 0.4, with $Q$ = {qlo:g}$-${qhi:g} GeV held fixed as well "
                  rf"as the colour-frame energy.  The lines still climb by {med_r:.0f}% over "
                  rf"$p_T^{{\rm lab}}$ = {QUOTE_PT[0]:g}$-${QUOTE_PT[1]:g} GeV: the cone keeps a "
                  "$p_T$-dependent fraction of the shower.  Both axes are logarithmic and this "
                  "figure and the next share a vertical span, so equal visual slope means equal "
                  "fractional change.")
    _draw_fixed_q(p_hem, r"$\langle n_{90}\rangle$,  whole current hemisphere", span,
                  "fixed_q_hemisphere", outdir,
                  "The same events with the whole Breit current hemisphere taken as the jet, on the "
                  f"same logarithmic span.  Nothing is left outside the cone and the lines flatten "
                  f"to {med_h:.0f}%.  Lab-frame clustering does show the shower pattern, but only "
                  "once it stops cutting the shower and $Q$ is held fixed.")
    return med_r, med_h


def fig_slope_vs_radius(trees, outdir):
    """Residual dependence on lab pT, at fixed Q and colour-frame energy, against jet radius."""
    fig, ax = plt.subplots(figsize=(4.6, 3.2))
    xs, ys, es = [], [], []
    for name, R in R_TREES:
        _, _, slopes = _fixed_q_profiles(trees[name], PT_EDGES_Q)
        if not slopes:
            continue
        xs.append(R); ys.append(np.median(slopes))
        es.append(np.std(slopes, ddof=1) / np.sqrt(len(slopes)) if len(slopes) > 1 else 0.0)
    _, _, sl_h = _fixed_q_profiles(trees["hemisphere"], PT_EDGES_Q)
    y_h = float(np.median(sl_h))
    X_HEMI = 3.6                                        # plotted beyond the radius scan
    ax.axhline(0.0, color=FAINT, lw=0.8, zorder=0)
    ax.errorbar(xs, ys, yerr=es, color=INK, lw=1.2, elinewidth=0.6, capsize=0,
                marker="o", ms=4, mec="white", mew=0.5)
    ax.plot([X_HEMI], [y_h], marker="D", ms=5, color=ACCENT, mec="white", mew=0.5)
    ax.annotate("whole current\nhemisphere", (X_HEMI, y_h), xytext=(0, -14),
                textcoords="offset points", fontsize=7.5, color=ACCENT, ha="center", va="top")
    ax.annotate("flat", (0.42, 0.012), fontsize=7.5, color=MUTED, va="bottom")
    ax.annotate(r"anti-$k_T$ radius in the lab", (1.2, ys[-1] + 0.03), fontsize=7.5,
                color=MUTED, ha="center")
    ax.set_xlabel(r"jet radius $R$   (rightmost point: no cone at all)")
    ax.set_ylabel(r"residual slope  $\mathrm{d}\ln\langle n_{90}\rangle\,/\,\mathrm{d}\ln p_T^{\rm lab}$")
    ax.set_xlim(0.2, 4.1)
    ax.set_xticks([0.4, 0.8, 1.2, 1.6, 2.4, X_HEMI])
    ax.set_xticklabels(["0.4", "0.8", "1.2", "1.6", "2.4", "all"])
    range_frame(ax, np.array(xs), np.array(ys + [y_h, 0.0]))
    qlo, qhi = Q_WINDOW
    caption(ax, rf"With $Q$ = {qlo:g}$-${qhi:g} GeV and the colour-frame energy both held fixed, "
                "how much dependence on lab transverse momentum survives.  Zero is a flat curve.  "
                "Widening the cone barely helps until it is wide enough to hold the whole current "
                "system, which in the lab means no cone at all.")
    save(fig, outdir, "slope_vs_radius")
    return dict(zip([r for _, r in R_TREES], ys)), y_h


# ---------------------------------------------------------------------------
# Hemisphere against full lab momentum: sliced in E_cm, and inclusive
# ---------------------------------------------------------------------------

P_HEMI = np.array([1.0, 1.6, 2.5, 4.0, 6.3, 10.0, 16.0, 25.0, 40.0, 63.0])
E_SLICES_H = [(2.5, 4), (4, 6), (6, 9), (9, 14), (14, 22)]
E_COLORS_H = ["#9ecae1", "#6baed6", "#3182bd", "#08519c", "#08306b"]


def _hemi_lines(ax, d, mask, labels, allx, ally, annotate_q_for=(4, 6), min_entries=200):
    """E_cm-sliced <n90> vs |p|_lab for hemispheres passing mask; returns spreads."""
    spreads = []
    for (lo, hi), col in zip(E_SLICES_H, E_COLORS_H):
        m = mask & (d["e_hcm"] >= lo) & (d["e_hcm"] < hi)
        xc, mu, se = profile(d["plab"][m], d["n90"][m], P_HEMI, min_entries=min_entries)
        ok = np.isfinite(mu)
        if ok.sum() < 3:
            continue
        ax.errorbar(xc[ok], mu[ok], yerr=se[ok], color=col, lw=1.2, elinewidth=0.6,
                    capsize=0, marker="o", ms=3, mec="white", mew=0.4)
        labels.add(xc[ok], mu[ok], e_label(lo, hi), col)
        spreads.append(100 * (mu[ok].max() - mu[ok].min()) / mu[ok].mean())
        allx.append(xc[ok]); ally.append(mu[ok])
        if annotate_q_for and (lo, hi) == annotate_q_for:
            idx = np.digitize(d["plab"][m], P_HEMI) - 1
            for k in np.where(ok)[0]:
                qmed = np.median(np.sqrt(d["Q2"][m][idx == k]))
                ax.annotate(rf"$Q\approx${qmed:.1f}", (xc[k], mu[k]), xytext=(0, -11),
                            textcoords="offset points", ha="center", fontsize=6, color=MUTED)
    return spreads


def fig_hemisphere_vs_p(trees, outdir):
    """Inclusive slope, and E_cm slices that do not flatten because |p|_lab is Q."""
    d = trees["hemisphere"]
    fig, ax = plt.subplots(figsize=(4.6, 3.6))
    labels = EndLabels(ax)
    allx, ally = [], []
    every = np.ones(len(d["W"]), dtype=bool)
    spreads = _hemi_lines(ax, d, every, labels, allx, ally)
    xc, mu, se = profile(d["plab"], d["n90"], P_HEMI, min_entries=200)
    ok = np.isfinite(mu)
    ax.errorbar(xc[ok], mu[ok], yerr=se[ok], color=INK, lw=1.8, elinewidth=0.6, capsize=0,
                marker="o", ms=3.5, mec="white", mew=0.4, zorder=5)
    labels.add(xc[ok], mu[ok], "all hemispheres", INK)
    allx.append(xc[ok]); ally.append(mu[ok])
    ax.set_xscale("log")
    ax.set_xticks([1, 2, 5, 10, 20, 50]); ax.set_xticklabels(["1", "2", "5", "10", "20", "50"])
    ax.minorticks_off()
    ax.set_xlabel(r"full lab momentum of the current hemisphere, $|\vec p|_{\rm lab}$  [GeV]")
    ax.set_ylabel(r"$\langle n_{90}\rangle$")
    ax.set_xlim(0.9, 130)
    range_frame(ax, np.concatenate(allx), np.concatenate(ally))
    labels.draw()
    caption(ax, "Whole Breit current hemisphere, every quantity measured in the lab.  The inclusive "
                "curve climbs steeply, and slicing in colour-frame energy does not flatten it "
                f"(median spread {np.median(spreads):.0f}%).  The small grey numbers give the median "
                "$Q$ in each bin of one slice: at fixed $E_{\\rm cm}$ the lab momentum of the current "
                "system is set by $Q$, so this axis is a $Q$ scan in disguise.")
    save(fig, outdir, "hemisphere_vs_p")
    return float(np.median(spreads))


def fig_hemisphere_vs_p_fixed_q(trees, outdir, q_window=(5.0, 7.5)):
    d = trees["hemisphere"]
    Q = np.sqrt(d["Q2"])
    qlo, qhi = q_window
    fig, ax = plt.subplots(figsize=(4.6, 3.6))
    labels = EndLabels(ax)
    allx, ally = [], []
    spreads = _hemi_lines(ax, d, (Q >= qlo) & (Q < qhi), labels, allx, ally,
                          annotate_q_for=None, min_entries=150)
    ax.set_xscale("log")
    ax.set_xticks([1, 2, 5, 10, 20, 50]); ax.set_xticklabels(["1", "2", "5", "10", "20", "50"])
    ax.minorticks_off()
    ax.set_xlabel(r"$|\vec p|_{\rm lab}$ of the current hemisphere  [GeV]")
    ax.set_ylabel(r"$\langle n_{90}\rangle$")
    ax.set_xlim(0.9, 130)
    range_frame(ax, np.concatenate(allx), np.concatenate(ally))
    labels.draw()
    caption(ax, rf"The same object with $Q$ = {qlo:g}$-${qhi:g} GeV held fixed.  The slices with most "
                "of the current system's energy flatten (a few percent); the low-energy slices keep a "
                "residual rise.  Within a fixed $(E_{\\rm cm}, Q)$ cell the lab momentum has little room "
                "left to vary, and what varies is the share of $W$ the hemisphere carries.")
    save(fig, outdir, "hemisphere_vs_p_fixed_q")
    return float(np.median(spreads))


def fig_hemisphere_p_vs_q(trees, outdir):
    """Why: at fixed W the boost of the current system is a function of Q."""
    d = trees["hemisphere"]
    Q = np.sqrt(d["Q2"])
    fig, ax = plt.subplots(figsize=(4.6, 3.4))
    labels = EndLabels(ax)
    q_edges = np.array([1.0, 1.5, 2.2, 3.3, 5.0, 7.5, 11.0, 17.0])
    allx, ally = [], []
    for (lo, hi), col in zip(W_SLICES, W_COLORS):
        m = (d["W"] >= lo) & (d["W"] < hi)
        idx = np.digitize(Q[m], q_edges) - 1
        med = np.array([np.median(d["plab"][m][idx == k]) if (idx == k).sum() >= 300 else np.nan
                        for k in range(len(q_edges) - 1)])
        lo16 = np.array([np.percentile(d["plab"][m][idx == k], 16) if (idx == k).sum() >= 300 else np.nan
                         for k in range(len(q_edges) - 1)])
        hi84 = np.array([np.percentile(d["plab"][m][idx == k], 84) if (idx == k).sum() >= 300 else np.nan
                         for k in range(len(q_edges) - 1)])
        xc = np.sqrt(q_edges[:-1] * q_edges[1:])
        ok = np.isfinite(med)
        if ok.sum() < 2:
            continue
        ax.fill_between(xc[ok], lo16[ok], hi84[ok], color=col, alpha=0.18, lw=0)
        ax.plot(xc[ok], med[ok], color=col, lw=1.2, marker="o", ms=3, mec="white", mew=0.4)
        labels.add(xc[ok], med[ok], w_label(lo, hi), col)
        allx.append(xc[ok]); ally.append(med[ok])
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xticks([1, 2, 3, 5, 10, 15]); ax.set_xticklabels(["1", "2", "3", "5", "10", "15"])
    ax.set_yticks([1, 2, 5, 10, 20, 50]); ax.set_yticklabels(["1", "2", "5", "10", "20", "50"])
    ax.minorticks_off()
    ax.set_xlabel(r"$Q$  [GeV]")
    ax.set_ylabel(r"$|\vec p|_{\rm lab}$ of the current hemisphere  [GeV]")
    ax.set_xlim(0.9, 30)
    range_frame(ax, np.concatenate(allx), np.concatenate(ally))
    labels.draw()
    caption(ax, "Median and central 68% of the current hemisphere's lab momentum against $Q$ in "
                "slices of $W$.  DIS has two kinematic degrees of freedom; once $W$ (hence the "
                "colour-frame energy) is fixed, the boost of the current system into the lab is a "
                "function of $Q$ alone.  Unlike $e^+e^-\\to ZZ$, there is no independent boost knob.")
    save(fig, outdir, "hemisphere_p_vs_q")


# ---------------------------------------------------------------------------
# Beam-energy test: same (W, Q), three lab frames
# ---------------------------------------------------------------------------

BEAM_CELLS_W = [(10, 15), (15, 22), (22, 28)]
BEAM_CELLS_Q = [(2.2, 3.3), (3.3, 5.0), (5.0, 7.5)]
Q_MARKERS = {0: "o", 1: "s", 2: "^"}
BEAM_MIN = 400


def _load_tree(path, tree, cols, sel_fn):
    d = uproot.open(path)[tree].arrays(cols, library="np")
    m = sel_fn(d)
    return {k: v[m] for k, v in d.items()}


def _beam_cells(beams, xkey="plab", ecm_cells=False):
    """
    For each (W or E_cm, Q) cell and each beam: median lab momentum, <n90>, sem.
    Returns list of (cell label, W index, Q index, [(x, y, e, beam label), ...]).
    """
    rows = []
    kcells = E_SLICES if ecm_cells else BEAM_CELLS_W
    kvar = "e_hcm" if ecm_cells else "W"
    for iw, (wlo, whi) in enumerate(kcells):
        for iq, (qlo, qhi) in enumerate(BEAM_CELLS_Q):
            pts = []
            for label, d in beams:
                Q = np.sqrt(d["Q2"])
                m = (d[kvar] >= wlo) & (d[kvar] < whi) & (Q >= qlo) & (Q < qhi)
                if m.sum() < BEAM_MIN:
                    pts = []
                    break
                v = d["n90"][m]
                pts.append((np.median(d[xkey][m]), v.mean(), v.std(ddof=1) / np.sqrt(len(v)), label))
            if len(pts) == len(beams):
                name = (rf"$E_{{\rm cm}}$ {wlo:g}$-${whi:g}, $Q$ {qlo:g}$-${qhi:g}" if ecm_cells
                        else rf"$W$ {wlo:g}$-${whi:g}, $Q$ {qlo:g}$-${qhi:g}")
                rows.append((name, iw, iq, pts))
    return rows


def _draw_beam_rows(ax, rows, colors):
    labels = EndLabels(ax, min_sep_pt=8.5, fontsize=7)
    allx, ally, flat = [], [], []
    for name, iw, iq, pts in rows:
        xs = np.array([p[0] for p in pts]); ys = np.array([p[1] for p in pts]); es = np.array([p[2] for p in pts])
        col = colors[iw]
        ax.errorbar(xs, ys, yerr=es, color=col, lw=1.1, elinewidth=0.6, capsize=0,
                    marker=Q_MARKERS[iq], ms=3.5, mec="white", mew=0.4)
        labels.add(xs, ys, name, col)
        allx.append(xs); ally.append(ys)
        flat.append(100 * (ys.max() - ys.min()) / ys.mean())
    return labels, np.concatenate(allx), np.concatenate(ally), flat


def fig_beam_energy(beam_paths, outdir, inclusive_path=None):
    """
    beam_paths: list of (label, path) in increasing sqrt(s).  Three figures:
    hemisphere at fixed (W, Q); leading R = 0.4 lab jet at fixed (W, Q);
    gamma*p-frame jets at fixed (E_cm, Q).
    """
    results = {}

    # -- whole current hemisphere ------------------------------------------
    beams = [(lab, _load_tree(p, "hemisphere", ["W", "Q2", "plab", "n90"],
                              lambda d: np.isfinite(d["n90"]) & (d["plab"] > 1.0)))
             for lab, p in beam_paths]
    rows = _beam_cells(beams)
    fig, ax = plt.subplots(figsize=(4.8, 3.8))
    if inclusive_path:
        d = beams[[l for l, _ in beam_paths].index(inclusive_path)][1]
        xc, mu, se = profile(d["plab"], d["n90"], P_HEMI, min_entries=200)
        ok = np.isfinite(mu)
        ax.plot(xc[ok], mu[ok], color=FAINT, lw=2.4, zorder=0)
        ax.annotate(f"all hemispheres, {inclusive_path}", (xc[ok][-1], mu[ok][-1]),
                    xytext=(6, 0), textcoords="offset points", fontsize=7, color=MUTED, va="center")
    labels, allx, ally, flat = _draw_beam_rows(ax, rows, ["#6baed6", "#3182bd", "#08306b"])
    # beam labels along one row
    if rows:
        _, _, _, pts = rows[len(rows) // 2]
        for x, y, _, lab in pts:
            ax.annotate(lab, (x, y), xytext=(0, -12), textcoords="offset points",
                        fontsize=6.5, color=MUTED, ha="center")
    ax.set_xscale("log")
    ax.set_xticks([1, 2, 5, 10, 20, 50]); ax.set_xticklabels(["1", "2", "5", "10", "20", "50"])
    ax.minorticks_off()
    ax.set_xlabel(r"$|\vec p|_{\rm lab}$ of the current hemisphere  [GeV]")
    ax.set_ylabel(r"$\langle n_{90}\rangle$")
    ax.set_xlim(0.9, 200)
    range_frame(ax, allx, ally)
    labels.draw()
    med = float(np.median(flat))
    caption(ax, "Each line is one cell of fixed $(W, Q)$, so fixed colour-frame kinematics, measured "
                "in three EIC beam configurations.  The beam energy changes only the lab frame; the "
                "current hemisphere's lab momentum moves by a large factor along each line and "
                f"$\\langle n_{{90}}\\rangle$ does not: median variation {med:.0f}%.  The grey band is "
                "the inclusive curve at one beam energy.  Marker shape encodes the $Q$ bin.")
    save(fig, outdir, "beam_energy_hemisphere")
    results["hemisphere"] = flat

    # -- leading R = 0.4 lab jet: the cone breaks the invariance -----------
    beams_j = [(lab, _load_tree(p, "jets_R0p4", ["W", "Q2", "plab", "n90", "lead", "current"],
                                lambda d: d["lead"] & d["current"] & np.isfinite(d["n90"])))
               for lab, p in beam_paths]
    rows = _beam_cells(beams_j)
    fig, ax = plt.subplots(figsize=(4.8, 3.8))
    labels, allx, ally, flat = _draw_beam_rows(ax, rows, ["#6baed6", "#3182bd", "#08306b"])
    ax.set_xscale("log")
    ax.set_xticks([1, 2, 5, 10, 20, 50]); ax.set_xticklabels(["1", "2", "5", "10", "20", "50"])
    ax.minorticks_off()
    ax.set_xlabel(r"$|\vec p|_{\rm lab}$ of the leading lab jet, $R$ = 0.4  [GeV]")
    ax.set_ylabel(r"$\langle n_{90}\rangle$")
    ax.set_xlim(0.9, 200)
    range_frame(ax, allx, ally)
    labels.draw()
    med_j = float(np.median(flat))
    caption(ax, "The same cells and beam configurations, now for the leading anti-$k_T$ $R$ = 0.4 jet "
                "clustered in the lab.  Fixed colour-frame kinematics no longer gives a fixed answer "
                f"(median variation {med_j:.0f}%): the boost changes which particles the cone holds.")
    save(fig, outdir, "beam_energy_labjet")
    results["labjet"] = flat

    # -- gamma*p-frame jets at fixed (E_cm, Q) ------------------------------
    beams_c = [(lab, _load_tree(p, "cmjets", ["e_hcm", "Q2", "plab", "n90"],
                                lambda d: np.isfinite(d["n90"])))
               for lab, p in beam_paths]
    for _, d in beams_c:
        pass
    rows = _beam_cells(beams_c, ecm_cells=True)
    if rows:
        fig, ax = plt.subplots(figsize=(4.8, 3.8))
        labels, allx, ally, flat = _draw_beam_rows(ax, rows, E_COLORS)
        ax.set_xscale("log")
        ax.set_xticks([1, 2, 5, 10, 20, 50]); ax.set_xticklabels(["1", "2", "5", "10", "20", "50"])
        ax.minorticks_off()
        ax.set_xlabel(r"$|\vec p|_{\rm lab}$ of the $\gamma^*p$-frame jet  [GeV]")
        ax.set_ylabel(r"$\langle n_{90}\rangle$")
        ax.set_xlim(0.9, 200)
        range_frame(ax, allx, ally)
        labels.draw()
        med_c = float(np.median(flat))
        caption(ax, "Jets clustered in the colour rest frame, cells of fixed jet energy in that frame "
                    f"and fixed $Q$, across the three beam configurations: median variation {med_c:.0f}%.")
        save(fig, outdir, "beam_energy_cmjet")
        results["cmjet"] = flat
    return results


# ---------------------------------------------------------------------------
# Figure 2b: the boost factor — same colour-frame jet, different lab boosts
# ---------------------------------------------------------------------------

def fig_boost_factor(jets, outdir, obs="n90", name="ffs_boost_factor",
                     ylabel=r"$\langle n_{90}\rangle$", ylim=None):
    fig, ax = plt.subplots(figsize=(4.6, 3.4))
    kappa = jets["plab"] / jets["p_hcm"]
    edges = np.array([0.15, 0.25, 0.4, 0.6, 0.85, 1.2, 1.7, 2.5, 3.5, 5.0, 7.0])
    allx, ally = [], []
    labels = EndLabels(ax)
    for (lo, hi), col in zip(W_SLICES, W_COLORS):
        m = (jets["W"] >= lo) & (jets["W"] < hi)
        xc, mu, se = profile(kappa[m], jets[obs][m], edges, min_entries=200)
        ok = np.isfinite(mu)
        ax.errorbar(xc[ok], mu[ok], yerr=se[ok], color=col, lw=1.2, elinewidth=0.6,
                    capsize=0, marker="o", ms=3, mec="white", mew=0.4)
        labels.add(xc[ok], mu[ok], w_label(lo, hi), col)
        allx.append(xc[ok]); ally.append(mu[ok])
    allx, ally = np.concatenate(allx), np.concatenate(ally)
    ax.axvline(1.0, color=FAINT, lw=0.6, zorder=0)
    ax.set_xscale("log")
    ax.set_xticks([0.2, 0.3, 0.5, 1, 2, 3, 5])
    ax.set_xticklabels(["0.2", "0.3", "0.5", "1", "2", "3", "5"])
    ax.minorticks_off()
    ax.set_xlabel(r"boost factor seen by the jet, $|\vec p|_{\rm lab}\,/\,|\vec p|_{\rm cm}$")
    ax.set_ylabel(ylabel)
    ax.set_xlim(0.14, 14)
    if ylim is not None:
        ax.set_ylim(*ylim)
    range_frame(ax, allx, ally)
    labels.draw()
    caption(ax, "Within a $W$ slice the current jet carries a nearly fixed momentum in the "
                "colour rest frame ($\\approx W/3$ here), so moving along a line changes only "
                "how the lab frame sees that same jet: compressed into the $R$ = 0.4 cone "
                "when boosted forward (right), opened up and partly lost when it recoils "
                "against the boost (left).")
    save(fig, outdir, name)


# ---------------------------------------------------------------------------
# Figure 3: ratio to the reference W slice — size of the shift
# ---------------------------------------------------------------------------

def fig_ratio(jets, outdir, ref=0):
    fig, ax = plt.subplots(figsize=(4.6, 3.0))
    lo_r, hi_r = W_SLICES[ref]
    m_r = (jets["W"] >= lo_r) & (jets["W"] < hi_r)
    xc, mu_r, se_r = profile(jets["plab"][m_r], jets["n90"][m_r], P_EDGES)
    ax.axhline(1.0, color=FAINT, lw=0.8, zorder=0)
    labels = EndLabels(ax)
    labels.add(xc[np.isfinite(mu_r)], np.ones(np.isfinite(mu_r).sum()), w_label(lo_r, hi_r) + " (reference)", MUTED)
    allx, ally = [xc[np.isfinite(mu_r)]], [np.ones(np.isfinite(mu_r).sum())]
    for (lo, hi), col in list(zip(W_SLICES, W_COLORS))[1:]:
        m = (jets["W"] >= lo) & (jets["W"] < hi)
        _, mu, se = profile(jets["plab"][m], jets["n90"][m], P_EDGES)
        r = mu / mu_r
        re = r * np.sqrt((se / mu)**2 + (se_r / mu_r)**2)
        ok = np.isfinite(r)
        ax.errorbar(xc[ok], r[ok], yerr=re[ok], color=col, lw=1.2, elinewidth=0.6,
                    capsize=0, marker="o", ms=3, mec="white", mew=0.4)
        labels.add(xc[ok], r[ok], w_label(lo, hi), col)
        allx.append(xc[ok]); ally.append(r[ok])
    allx, ally = np.concatenate(allx), np.concatenate(ally)
    ax.set_xscale("log")
    ax.set_xticks([2, 3, 5, 10, 20, 30])
    ax.set_xticklabels(["2", "3", "5", "10", "20", "30"])
    ax.minorticks_off()
    ax.set_xlim(1.8, 60)
    ax.set_xlabel(r"$|\vec p|_{\rm lab}$  [GeV]")
    ax.set_ylabel(r"$\langle n_{90}\rangle\,/\,\langle n_{90}\rangle_{\rm ref}$")
    range_frame(ax, allx, ally)
    labels.draw()
    caption(ax, "Fractional frame-dependent fragmentation shift at the EIC: the ratio of "
                "$\\langle n_{90}\\rangle$ at fixed lab momentum between $W$ slices.")
    save(fig, outdir, "ffs_ratio")


# ---------------------------------------------------------------------------
# Figure 4: slopegraph — low-W to high-W shift per |p|_lab bin
# ---------------------------------------------------------------------------

def fig_slopegraph(jets, outdir, left=0, right=3):
    fig, ax = plt.subplots(figsize=(3.2, 3.6))
    (lo_l, hi_l), (lo_r, hi_r) = W_SLICES[left], W_SLICES[right]
    m_l = (jets["W"] >= lo_l) & (jets["W"] < hi_l)
    m_r = (jets["W"] >= lo_r) & (jets["W"] < hi_r)
    _, mu_l, se_l = profile(jets["plab"][m_l], jets["n90"][m_l], P_EDGES)
    _, mu_r, se_r = profile(jets["plab"][m_r], jets["n90"][m_r], P_EDGES)
    ys = []
    for i in range(len(P_EDGES) - 1):
        if not (np.isfinite(mu_l[i]) and np.isfinite(mu_r[i])):
            continue
        shade = 0.25 + 0.75 * i / (len(P_EDGES) - 2)
        col = mcolors.to_hex((1 - shade) * np.array([1, 1, 1]) + shade * np.array([0.12, 0.12, 0.12]))
        ax.plot([0, 1], [mu_l[i], mu_r[i]], color=col, lw=1.0)
        ax.plot([0, 1], [mu_l[i], mu_r[i]], "o", color=col, ms=3)
        ax.text(-0.04, mu_l[i], f"{mu_l[i]:.2f}", ha="right", va="center", fontsize=7.5, color=col)
        pct = 100 * (mu_r[i] / mu_l[i] - 1)
        ax.text(1.04, mu_r[i],
                rf"{mu_r[i]:.2f}   {P_EDGES[i]:g}$-${P_EDGES[i+1]:g} GeV  ({pct:+.0f}%)",
                ha="left", va="center", fontsize=7.5, color=col)
        ys += [mu_l[i], mu_r[i]]
    ax.set_xlim(-0.35, 2.3)
    ax.set_xticks([0, 1])
    ax.set_xticklabels([w_label(lo_l, hi_l), w_label(lo_r, hi_r)], fontsize=8)
    ax.tick_params(axis="x", length=0)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.set_yticks([])
    ax.set_title(r"$\langle n_{90}\rangle$ per $|\vec p|_{\rm lab}$ bin", fontsize=9, loc="left")
    caption(ax, "Each line is one lab-momentum bin; the right-hand label gives the bin and "
                "the relative change from low to high $W$ at fixed lab momentum.", fontsize=7)
    save(fig, outdir, "ffs_slopegraph")


# ---------------------------------------------------------------------------
# Figure 5: distribution bands — n90 quantiles vs W at fixed |p|_lab
# ---------------------------------------------------------------------------

def fig_distribution(jets, outdir):
    fig, ax = plt.subplots(figsize=(4.6, 3.2))
    plo, phi = P_WINDOW
    m = (jets["plab"] >= plo) & (jets["plab"] < phi)
    W = jets["W"][m]; v = jets["n90"][m]
    edges = np.array([10, 13, 16, 20, 24, 29, 35, 42, 50, 60], float)
    centres = 0.5 * (edges[:-1] + edges[1:])
    qs = {q: np.full(len(centres), np.nan) for q in (10, 25, 50, 75, 90)}
    mean = np.full(len(centres), np.nan)
    idx = np.digitize(W, edges) - 1
    for i in range(len(centres)):
        s = v[idx == i]
        if len(s) >= MIN_ENTRIES:
            for q in qs:
                qs[q][i] = np.percentile(s, q)
            mean[i] = s.mean()
    ax.fill_between(centres, qs[10], qs[90], color="#e6e6e6", lw=0)
    ax.fill_between(centres, qs[25], qs[75], color="#c4c4c4", lw=0)
    ax.plot(centres, qs[50], color=INK, lw=1.2)
    ax.plot(centres, mean, color=ACCENT, lw=1.0, ls=(0, (3, 2)))
    for q, txt, dy in ((90, "90th percentile", 0), (75, "75th", 0), (50, "median", -4),
                       (25, "25th", 3), (10, "10th", -4)):
        label_line_end(ax, centres, qs[q], txt, INK if q == 50 else MUTED, dy=dy, fontsize=7.5)
    label_line_end(ax, centres, mean, "mean", ACCENT, dy=4, fontsize=7.5)
    ax.set_xlabel(r"$W$  [GeV]")
    ax.set_ylabel(r"$n_{90}$")
    ax.set_xlim(10, 62)
    ax.set_xticks([10, 20, 30, 40, 50])
    range_frame(ax, centres, np.concatenate([qs[10], qs[90]]))
    caption(ax, rf"Jets with $|\vec p|_{{\rm lab}}$ = {plo:g}$-${phi:g} GeV.  The whole distribution "
                "moves, not only its mean: bands are the central 50% and 80% of jets.")
    save(fig, outdir, "ffs_distribution")


# ---------------------------------------------------------------------------
# Figure 6: (W, Q) table-graphic — which scale does the jet know about?
# ---------------------------------------------------------------------------

def fig_wq_table(jets, outdir):
    fig, ax = plt.subplots(figsize=(4.6, 3.6))
    plo, phi = P_WINDOW
    m = (jets["plab"] >= plo) & (jets["plab"] < phi)
    W = jets["W"][m]; Q = np.sqrt(jets["Q2"][m]); v = jets["n90"][m]
    W_e = np.array([10, 15, 22, 32, 45, 60], float)
    Q_e = np.array([1.0, 1.5, 2.2, 3.3, 5.0, 7.5, 11.0, 17.0])
    iw = np.digitize(W, W_e) - 1
    iq = np.digitize(Q, Q_e) - 1
    tab = np.full((len(Q_e) - 1, len(W_e) - 1), np.nan)
    cnt = np.zeros_like(tab)
    for a in range(tab.shape[0]):
        for b in range(tab.shape[1]):
            s = v[(iq == a) & (iw == b)]
            cnt[a, b] = len(s)
            if len(s) >= MIN_ENTRIES:
                tab[a, b] = s.mean()
    vmin, vmax = np.nanmin(tab), np.nanmax(tab)
    cmap = mcolors.LinearSegmentedColormap.from_list("blues", ["#f7fbff", "#08306b"])
    for a in range(tab.shape[0]):
        for b in range(tab.shape[1]):
            if np.isfinite(tab[a, b]):
                t = (tab[a, b] - vmin) / max(vmax - vmin, 1e-9)
                ax.add_patch(plt.Rectangle((b, a), 1, 1, color=cmap(0.15 + 0.7 * t), lw=0))
                ax.text(b + 0.5, a + 0.5, f"{tab[a, b]:.2f}", ha="center", va="center",
                        fontsize=7.5, color="white" if t > 0.55 else INK)
    ax.set_xlim(0, tab.shape[1]); ax.set_ylim(0, tab.shape[0])
    ax.set_xticks(np.arange(tab.shape[1] + 1)); ax.set_xticklabels([f"{w:g}" for w in W_e], fontsize=8)
    ax.set_yticks(np.arange(tab.shape[0] + 1)); ax.set_yticklabels([f"{q:g}" for q in Q_e], fontsize=8)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_xlabel(r"$W$  [GeV]   (colour rest frame energy)")
    ax.set_ylabel(r"$Q$  [GeV]   (hard scale)")
    ax.set_title(rf"$\langle n_{{90}}\rangle$ for jets with $|\vec p|_{{\rm lab}}$ = {plo:g}$-${phi:g} GeV",
                 fontsize=9, loc="left")
    caption(ax, "DIS separates the two scales that coincide in $e^+e^-$: the energy of the colour "
                "rest frame ($W$) and the virtuality that starts the shower ($Q$).  Read along a "
                "row for the $W$ dependence at fixed $Q$, down a column for the reverse.  "
                "Cells with fewer than 40 jets are blank.", fontsize=7)
    save(fig, outdir, "ffs_wq_table")


# ---------------------------------------------------------------------------
# Figure 7: the boost map — where in the EIC phase space is the colour frame boosted?
# ---------------------------------------------------------------------------

def fig_boost_map(jets, outdir):
    fig, ax = plt.subplots(figsize=(4.6, 3.8))
    s = SQRT_S**2
    x = np.logspace(-4, 0, 300)
    Q2 = np.logspace(0, 3, 300)
    X, Q2g = np.meshgrid(x, Q2)
    Y = Q2g / (X * s)                                # inelasticity (massless beams)
    W2 = Q2g * (1 - X) / X
    Yr = np.clip(Y, 1e-6, 1 - 1e-6)
    # Massless beams: k = (0,0,-E_e,E_e), P = (0,0,E_p,E_p).
    # From y = 1 - P.k'/P.k and Q2 = 2 k.k':  E'-pz' = 2E_e(1-y),  E'+pz' = Q2/(2E_e).
    Ep = E_LEPTON * (1 - Yr) + Q2g / (4 * E_LEPTON)
    pz_p = -E_LEPTON * (1 - Yr) + Q2g / (4 * E_LEPTON)
    Eq, pzq = E_LEPTON - Ep, -E_LEPTON - pz_p            # q = k - k'
    Etot, pztot = E_PROTON + Eq, E_PROTON + pzq           # P + q
    yh = 0.5 * np.log(np.clip((Etot + pztot) / np.clip(Etot - pztot, 1e-9, None), 1e-9, None))
    yh = np.where((Y < 1) & (W2 > 0), yh, np.nan)

    def along_row(field, level, q2_row):
        """(x, Q2) where `field` crosses `level` on the row nearest Q2 = q2_row."""
        r = np.argmin(abs(Q2 - q2_row))
        row = field[r]
        ok = np.isfinite(row)
        if ok.sum() < 2 or not (row[ok].min() <= level <= row[ok].max()):
            return None
        xi = np.interp(level, row[ok][::-1], x[ok][::-1]) if row[ok][0] > row[ok][-1] \
            else np.interp(level, row[ok], x[ok])
        return (xi, Q2[r])

    def along_col(field, level, x_col):
        c = np.argmin(abs(x - x_col))
        col = field[:, c]
        ok = np.isfinite(col)
        if ok.sum() < 2 or not (col[ok].min() <= level <= col[ok].max()):
            return None
        cc, qq = col[ok], Q2[ok]
        if cc[0] > cc[-1]:
            cc, qq = cc[::-1], qq[::-1]
        return (x[c], np.interp(level, cc, qq))

    y_levels = [1.5, 2, 2.5, 3, 3.5, 4]
    cs = ax.contour(X, Q2g, yh, levels=y_levels, colors=INK, linewidths=0.7)
    pos = [p for p in ([along_row(yh, 1.5, 1.6), along_row(yh, 2.0, 1.6)]
                       + [along_col(yh, lv, 0.55) for lv in (2.5, 3, 3.5, 4)]) if p is not None]
    ax.clabel(cs, fmt=lambda v: rf"$y_{{\rm cm}}$ = {v:g}", fontsize=7, inline=True,
              manual=pos)
    W_levels = [10, 15, 22, 32, 45]
    Wl = ax.contour(X, Q2g, np.sqrt(W2), levels=W_levels, colors=W_COLORS,
                    linewidths=0.9, linestyles="--")
    pos = [p for p in (along_row(np.sqrt(W2), lv, 400.0) for lv in W_levels) if p is not None]
    ax.clabel(Wl, fmt=lambda v: rf"$W$ = {v:g}", fontsize=7, inline=True, manual=pos)
    ax.contour(X, Q2g, Y, levels=[0.95], colors=[MUTED], linewidths=0.6)
    ax.text(3e-4, 1.6e2, r"$y$ = 0.95", fontsize=7, color=MUTED, rotation=38)
    ax.scatter(jets["x"], jets["Q2"], s=0.5, color=ACCENT, alpha=0.08, lw=0, rasterized=True)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(1e-4, 1); ax.set_ylim(1, 1e3)
    ax.set_xlabel(r"Bjorken $x$")
    ax.set_ylabel(r"$Q^2$  [GeV$^2$]")
    range_frame(ax, np.array([1e-4, 1]), np.array([1, 1e3]))
    caption(ax, "Solid: rapidity of the colour rest frame ($\\gamma^*p$ CM) relative to the lab, "
                "$y_{\\rm cm}$; a jet at fixed lab momentum is seen in that frame with momentum "
                "reduced by roughly $e^{-y_{\\rm cm}}$ along the axis.  Dashed: lines of constant $W$.  "
                "Red: the selected current jets.", fontsize=7)
    save(fig, outdir, "boost_map")


# ---------------------------------------------------------------------------
# Figure 8: rapidity plateau in the γ*p frame per W slice (needs event file)
# ---------------------------------------------------------------------------

def fig_plateau(event_file, outdir, max_events=200_000):
    import awkward as ak
    from utils.dis_kinematics import hcm_boost_matrix, apply_boost
    ev = ak.from_parquet(event_file)[:max_events]
    W = ak.to_numpy(ev.W)
    P = np.stack([ak.to_numpy(ev.proton[c]) for c in ("px", "py", "pz", "e")], 1).astype(float)
    q = np.stack([ak.to_numpy(ev.q[c]) for c in ("px", "py", "pz", "e")], 1).astype(float)
    L = hcm_boost_matrix(P, q)
    qh = apply_boost(L, q); qhat = qh[:, :3] / np.linalg.norm(qh[:, :3], axis=1)[:, None]
    ch = ev.particles.charge != 0
    n = ak.to_numpy(ak.num(ev.particles.px[ch]))
    evi = np.repeat(np.arange(len(W)), n)
    C = np.stack([ak.to_numpy(ak.flatten(ev.particles[c][ch])) for c in ("px", "py", "pz", "e")], 1).astype(float)
    Ch = np.einsum("nij,nj->ni", L[evi], C)
    pl = np.einsum("ni,ni->n", Ch[:, :3], qhat[evi])
    ystar = 0.5 * np.log(np.clip((Ch[:, 3] + pl) / np.clip(Ch[:, 3] - pl, 1e-9, None), 1e-12, None))

    fig, ax = plt.subplots(figsize=(4.6, 3.2))
    edges = np.linspace(-6, 6, 49)
    cent = 0.5 * (edges[1:] + edges[:-1])
    for (lo, hi), col in zip(W_SLICES, W_COLORS):
        m_ev = (W >= lo) & (W < hi)
        if m_ev.sum() < 200:
            continue
        h, _ = np.histogram(ystar[m_ev[evi]], bins=edges)
        dens = h / m_ev.sum() / np.diff(edges)
        ax.plot(cent, dens, color=col, lw=1.2)
        k = W_SLICES.index((lo, hi))
        level = (0.9 - 0.2 * k) * dens.max()                 # graduated heights, no overlap
        right = cent > 0
        edge = np.max(cent[right][dens[right] > level])       # falling edge on the current side
        i = np.argmin(abs(cent - edge))
        ax.text(cent[i] + 0.15, dens[i], w_label(lo, hi), fontsize=7, color=col, ha="left", va="center")
    ax.axvline(0, color=FAINT, lw=0.6, zorder=0)
    ax.text(0.1, ax.get_ylim()[1] * 0.97, r"current $\rightarrow$", fontsize=7, color=MUTED, va="top")
    ax.text(-0.1, ax.get_ylim()[1] * 0.97, r"$\leftarrow$ remnant", fontsize=7, color=MUTED, va="top", ha="right")
    ax.set_xlabel(r"charged-hadron rapidity along the boson axis in the $\gamma^*p$ frame, $y^*$")
    ax.set_ylabel(r"$\mathrm{d}N_{\rm ch}/\mathrm{d}y^*$ per event")
    range_frame(ax, cent, np.array([0, ax.get_ylim()[1]]))
    caption(ax, "The colour string spans the whole hadronic system: its rapidity length grows "
                "like $\\ln W^2$.  A lab-frame jet is a window onto the current end of this "
                "plateau, and how much of the string it sees depends on the boost.", fontsize=7)
    save(fig, outdir, "plateau")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(description="Tufte-style FFS figures",
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("inputs", nargs="+", help="ROOT file(s) from analyze_events.py")
    p.add_argument("--outdir", default="figures/")
    p.add_argument("--events", default=None,
                   help="Optional event Parquet file for the rapidity-plateau figure")
    p.add_argument("--beams", nargs="*", default=None, metavar="LABEL=FILE",
                   help="Analysis files for other beam configurations, e.g. "
                        "5x41=data/analysis_e5p41.root; the main input is labelled 10x100")
    return p.parse_args()


def main():
    args = parse_args()
    jets = load_jets(args.inputs)
    print(f"{len(jets['W']):,} current jets loaded")
    trees = load_radius_trees(args.inputs)
    fig_capture(trees, args.outdir)
    med_r04, med_hemi = fig_fixed_q(trees, args.outdir)
    slopes_R, slope_h = fig_slope_vs_radius(trees, args.outdir)
    if args.beams:
        beam_paths = [("10x100", args.inputs[0])]
        for item in args.beams:
            lab, path = item.split("=", 1)
            beam_paths.append((lab, path))
        beam_paths.sort(key=lambda t: float(t[0].split("x")[0]) * float(t[0].split("x")[1]))
        res = fig_beam_energy(beam_paths, args.outdir, inclusive_path="10x100")
        for k, v in res.items():
            print(f"  beam-energy test, {k}: variation per cell "
                  + ", ".join(f"{x:.0f}%" for x in v) + f"  (median {np.median(v):.0f}%)")
    sp_incl = fig_hemisphere_vs_p(trees, args.outdir)
    sp_fq = fig_hemisphere_vs_p_fixed_q(trees, args.outdir)
    fig_hemisphere_p_vs_q(trees, args.outdir)
    print(f"  hemisphere vs |p|_lab: E_cm-sliced spread {sp_incl:.0f}%, at fixed Q {sp_fq:.0f}%")
    print(f"  at fixed Q: R=0.4 spread {med_r04:.1f}%, hemisphere {med_hemi:.1f}%")
    print("  residual slope vs R: "
          + ", ".join(f"{r}:{v:+.2f}" for r, v in slopes_R.items())
          + f", hemisphere:{slope_h:+.2f}")
    cj = load_cmjets(args.inputs)
    print(f"{len(cj['e_cm']):,} colour-frame jets loaded")
    ylim_cm, sp_cm = fig_flat_cm(cj, args.outdir)
    sp_lab = fig_flat_lab(jets, args.outdir, ylim=ylim_cm)
    fig_universal_cm(cj, args.outdir)
    fig_pt_fan(jets, args.outdir)
    print(f"  flatness: colour-frame slices vary by "
          f"{min(sp_cm.values()):.1f}-{max(sp_cm.values()):.1f}%, "
          f"lab slices by {min(sp_lab.values()):.1f}-{max(sp_lab.values()):.1f}%")
    ylim = fig_fan(jets, args.outdir)
    fig_collapse(jets, args.outdir, ylim=ylim)
    fig_boost_factor(jets, args.outdir, ylim=ylim)
    fig_boost_factor(jets, args.outdir, obs="n_const", name="ncon_boost_factor",
                     ylabel=r"$\langle N_{\rm constituents}\rangle$")
    fig_ratio(jets, args.outdir)
    fig_slopegraph(jets, args.outdir)
    fig_distribution(jets, args.outdir)
    fig_wq_table(jets, args.outdir)
    fig_fan(jets, args.outdir, obs="n_charged", name="nch_fan",
            ylabel=r"$\langle N_{\rm charged}\rangle$")
    fig_boost_map(jets, args.outdir)
    if args.events:
        fig_plateau(args.events, args.outdir)


if __name__ == "__main__":
    main()
