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
    return p.parse_args()


def main():
    args = parse_args()
    jets = load_jets(args.inputs)
    print(f"{len(jets['W']):,} current jets loaded")
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
