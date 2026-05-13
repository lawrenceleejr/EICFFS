#!/usr/bin/env python3
"""
FFS Effect Plotting Script
==========================
Produces publication-quality figures from the histogram ROOT file generated
by analyze_events.py, illustrating the Frame-dependent Fragmentation Shift
(FFS) effect in EIC deep-inelastic scattering.

Primary observable: n₉₀ — the fractional minimum number of jet constituents
(sorted by decreasing |p|) needed to carry 90% of the jet's total momentum.
This is the quantity studied in arXiv:2308.10951.  At fixed lab-frame jet
|p|, the FFS effect predicts ⟨n₉₀⟩ varies with W because different W values
correspond to different boosts between the lab and the photon-proton colour
rest frame.

Figures produced
----------------
1. ffs_main.pdf  — Primary result: ⟨n₉₀⟩ vs W for fixed |p_lab| bins
2. ffs_ratio.pdf — ⟨n₉₀⟩ ratio to lowest-|p_lab| reference bin (FFS magnitude)
3. kinematics.pdf — DIS kinematic plane (Q² vs W, x, y distributions)
4. jet_landscape.pdf — Jet η vs pT heat-map + jet multiplicity per event
5. ffs_heatmap.pdf — 2D: ⟨n₉₀⟩(|p_lab|, W) as a colour map

Usage
-----
    python make_plots.py data/histograms.root [--outdir plots/]
    python make_plots.py data/histograms.root --outdir plots/ --format pdf

Reference: arXiv:2308.10951  (Phys.Lett.B 866, 2025, 139561)
"""

import argparse
import os
import sys
import warnings

import numpy as np
import uproot
import hist
import matplotlib
matplotlib.use("Agg")                         # non-interactive backend
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.ticker as mticker
from matplotlib.gridspec import GridSpec
from scipy.ndimage import gaussian_filter

# ── Try to import mplhep for HEP-style axes; degrade gracefully ───────────
try:
    import mplhep as hep
    HEP_STYLE = True
except ImportError:
    HEP_STYLE = False
    warnings.warn("mplhep not installed; using default matplotlib style.")


# ---------------------------------------------------------------------------
# Design system  ─  colours, typography, geometry
# ---------------------------------------------------------------------------

# EIC brand palette (custom)
EIC_BLUE   = "#003865"
EIC_TEAL   = "#00A3AD"
EIC_GOLD   = "#F5A800"
EIC_CORAL  = "#E05C2B"
EIC_LILAC  = "#7B5EA7"
EIC_GREEN  = "#5CA35C"

# Ordered palette for multi-curve plots
PALETTE = [EIC_TEAL, EIC_GOLD, EIC_CORAL, EIC_LILAC, EIC_GREEN, EIC_BLUE]

MARKERS = ["o", "s", "^", "D", "v", "P"]

# Typography
FONT_TITLE  = 18
FONT_LABEL  = 15
FONT_TICK   = 13
FONT_LEGEND = 13
FONT_ANNOT  = 11

# W and p_lab bin edges (must match analyze_events.py)
W_BINS    = np.array([5.0, 10.0, 20.0, 30.0, 40.0, 55.0])
P_LAB_BINS = np.array([2.0, 5.0, 10.0, 20.0])


# ---------------------------------------------------------------------------
# Global style
# ---------------------------------------------------------------------------

def apply_style():
    """Configure matplotlib rcParams for a clean, publication-ready look."""
    rc = {
        # Figure
        "figure.facecolor":      "white",
        "figure.dpi":            150,
        "figure.autolayout":     False,

        # Axes
        "axes.facecolor":        "white",
        "axes.edgecolor":        "#333333",
        "axes.linewidth":        1.2,
        "axes.spines.top":       False,
        "axes.spines.right":     False,
        "axes.grid":             True,
        "grid.color":            "#DDDDDD",
        "grid.linestyle":        "-",
        "grid.linewidth":        0.6,
        "grid.alpha":            0.8,

        # Ticks
        "xtick.direction":       "out",
        "ytick.direction":       "out",
        "xtick.major.width":     1.2,
        "ytick.major.width":     1.2,
        "xtick.minor.visible":   True,
        "ytick.minor.visible":   True,
        "xtick.labelsize":       FONT_TICK,
        "ytick.labelsize":       FONT_TICK,

        # Lines / markers
        "lines.linewidth":       2.0,
        "lines.markersize":      7,
        "errorbar.capsize":      3,

        # Font
        "font.family":           "sans-serif",
        "font.size":             FONT_LABEL,
        "axes.titlesize":        FONT_TITLE,
        "axes.labelsize":        FONT_LABEL,
        "legend.fontsize":       FONT_LEGEND,
        "legend.frameon":        True,
        "legend.framealpha":     0.92,
        "legend.edgecolor":      "#CCCCCC",

        # Text rendering
        "text.usetex":           False,
        "mathtext.fontset":      "stix",
    }
    plt.rcParams.update(rc)

    if HEP_STYLE:
        hep.style.use("CMS")          # clean HEP style as base
        plt.rcParams.update(rc)       # re-apply our overrides


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def eic_label(ax, text="EIC MC  |  Pythia8 NC-DIS  |  10×100 GeV",
              loc="upper right"):
    """Add a standard experiment-label annotation."""
    ax.annotate(
        text,
        xy=(0.98, 0.97) if "right" in loc else (0.02, 0.97),
        xycoords="axes fraction",
        ha="right" if "right" in loc else "left",
        va="top",
        fontsize=FONT_ANNOT,
        color="#555555",
        style="italic",
    )


def arxiv_label(ax, ref="arXiv:2308.10951"):
    """Watermark the FFS paper reference."""
    ax.annotate(
        f"FFS effect  [{ref}]",
        xy=(0.5, -0.12),
        xycoords="axes fraction",
        ha="center",
        va="top",
        fontsize=FONT_ANNOT - 1,
        color="#888888",
    )


def save_fig(fig, path):
    fig.savefig(path, bbox_inches="tight", dpi=200)
    plt.close(fig)
    print(f"  → {path}")


def load_hist(f, name):
    """Load a histogram from an open uproot file as a hist.Hist object."""
    return f[name].to_hist()


# ---------------------------------------------------------------------------
# Mean-observable extractors from sum / count histogram pairs
# ---------------------------------------------------------------------------

def _mean_from_sum_count(hists, sum_key, count_key, plab_bin_idx,
                         min_count=10):
    """
    Generic extractor: ⟨obs⟩ ± σ_mean as a function of W for a given
    |p_lab| bin, from a pair of (sum, count) histograms stored by
    analyze_events.py.

    Both histograms are expected to have axes [W_fine, plab].
    """
    h_sum   = hists[sum_key]
    h_count = hists[count_key]
    ip = plab_bin_idx

    sum_vals   = h_sum.values()[:, ip]
    count_vals = h_count.values()[:, ip]
    W_centers  = h_sum.axes[0].centers

    with np.errstate(invalid="ignore", divide="ignore"):
        means = np.where(count_vals > 0, sum_vals / count_vals, np.nan)
        # Poisson approximation for σ_mean
        sigma_mean = np.where(
            count_vals > 1,
            np.sqrt(np.abs(means) / np.maximum(count_vals, 1)),
            np.nan,
        )

    mask = count_vals > min_count
    return W_centers, means, sigma_mean, mask


def mean_n90_from_hist(hists, plab_bin_idx):
    """Extract ⟨n₉₀⟩ ± σ_mean vs W for a given |p_lab| bin (primary obs.)."""
    return _mean_from_sum_count(
        hists, "sum_n90_vs_W", "count_n90_vs_W", plab_bin_idx)


def mean_mult_from_hist(hists, plab_bin_idx):
    """Extract ⟨N_charged⟩ ± σ_mean vs W for a given |p_lab| bin (secondary obs.)."""
    return _mean_from_sum_count(
        hists, "sum_N_vs_W", "count_vs_W", plab_bin_idx)


# ---------------------------------------------------------------------------
# Figure 1 – Primary FFS result: ⟨n₉₀⟩ vs W
# ---------------------------------------------------------------------------

def plot_ffs_main(hists, outdir):
    """⟨n₉₀⟩ vs W for each |p_lab| bin — primary FFS observable (arXiv:2308.10951)."""
    fig, ax = plt.subplots(figsize=(8, 6))

    plotted_any = False
    for k in range(len(P_LAB_BINS) - 1):
        W_c, means, errs, mask = mean_n90_from_hist(hists, k)
        if mask.sum() == 0:
            continue
        plotted_any = True
        label = (rf"$|p|_{{\rm lab}}\in"
                 rf"[{P_LAB_BINS[k]:.0f},{P_LAB_BINS[k+1]:.0f}]$ GeV")
        ax.errorbar(
            W_c[mask], means[mask], yerr=errs[mask],
            fmt=f"{MARKERS[k]}-",
            color=PALETTE[k],
            label=label,
            linewidth=2.2,
            markersize=8,
            capsize=4,
        )

    if not plotted_any:
        _plot_ffs_from_3d(hists, ax)

    ax.set_xlabel(r"$W = \sqrt{(P+q)^2}$  [GeV]", fontsize=FONT_LABEL)
    ax.set_ylabel(r"$\langle n_{90} \rangle$ per jet", fontsize=FONT_LABEL)
    ax.set_title(
        r"FFS Effect: $\langle n_{90}\rangle$ vs $W$ at Fixed Lab-Frame $|p|$",
        fontsize=FONT_TITLE, pad=12)
    ax.set_xlim(5, 55)
    ax.set_ylim(bottom=0)
    ax.xaxis.set_minor_locator(mticker.AutoMinorLocator())
    ax.yaxis.set_minor_locator(mticker.AutoMinorLocator())

    legend = ax.legend(
        title=r"Fixed lab-frame jet $|p|$",
        title_fontsize=FONT_LEGEND,
        loc="upper left",
        framealpha=0.93,
    )
    legend.get_frame().set_linewidth(0.8)

    eic_label(ax)
    _add_ffs_arrow(ax)

    fig.tight_layout()
    save_fig(fig, os.path.join(outdir, "ffs_main.pdf"))
    save_fig(fig, os.path.join(outdir, "ffs_main.png"))


def _add_ffs_arrow(ax):
    """Draw a subtle arrow hinting at the FFS trend."""
    ylim = ax.get_ylim()
    y_pos = ylim[0] + 0.08 * (ylim[1] - ylim[0])
    ax.annotate(
        r"FFS: $\langle n_{90}\rangle$ varies" "\n"
        r"with $W$ at fixed $|p|_{\rm lab}$",
        xy=(45, y_pos + 0.5 * (ylim[1] - ylim[0])),
        xytext=(35, y_pos + 0.15 * (ylim[1] - ylim[0])),
        arrowprops=dict(arrowstyle="->", color="#888888", lw=1.5),
        fontsize=FONT_ANNOT,
        color="#666666",
        ha="center",
    )


def _plot_ffs_from_3d(hists, ax):
    """
    Fall-back: extract ⟨n₉₀⟩ vs W from the 3D (W × p_lab × n₉₀) histogram
    when the sum/count histograms are empty.
    """
    h3 = hists["n90_3d"]
    W_edges    = h3.axes[0].edges
    PL_edges   = h3.axes[1].edges
    N90_centers = h3.axes[2].centers

    for k in range(len(PL_edges) - 1):
        values  = h3.values()[:, k, :]      # shape (nW, nN90)
        W_c     = 0.5 * (W_edges[:-1] + W_edges[1:])
        means, errs = [], []
        masks = []
        for iw in range(len(W_c)):
            counts = values[iw]
            tot = counts.sum()
            if tot > 10:
                mean = float(np.sum(counts * N90_centers) / tot)
                var  = float(np.sum(counts * (N90_centers - mean)**2) / tot)
                err  = float(np.sqrt(var / tot))
            else:
                mean = err = np.nan
            means.append(mean)
            errs.append(err)
            masks.append(tot > 10)

        means = np.array(means)
        errs  = np.array(errs)
        mask  = np.array(masks)
        if mask.sum() == 0:
            continue
        label = (rf"$|p|_{{\rm lab}}\in"
                 rf"[{PL_edges[k]:.0f},{PL_edges[k+1]:.0f}]$ GeV")
        ax.errorbar(
            W_c[mask], means[mask], yerr=errs[mask],
            fmt=f"{MARKERS[k]}-", color=PALETTE[k], label=label,
            linewidth=2.2, markersize=8, capsize=4,
        )


# ---------------------------------------------------------------------------
# Figure 2 – Ratio plot (FFS magnitude)
# ---------------------------------------------------------------------------

def plot_ffs_ratio(hists, outdir):
    """Ratio ⟨n₉₀⟩(W) / ⟨n₉₀⟩_ref vs W, where ref = lowest |p_lab| bin."""
    fig, (ax_main, ax_ratio) = plt.subplots(
        2, 1, figsize=(8, 7),
        gridspec_kw={"height_ratios": [2, 1], "hspace": 0.08},
        sharex=True,
    )

    ref_bin = 0    # lowest |p_lab| bin as reference

    all_means = {}
    all_errs  = {}
    all_masks = {}
    W_c_global = None

    for k in range(len(P_LAB_BINS) - 1):
        W_c, means, errs, mask = mean_n90_from_hist(hists, k)
        all_means[k] = means
        all_errs[k]  = errs
        all_masks[k] = mask
        W_c_global   = W_c

    if W_c_global is None:
        print("  [ratio plot] No data available; skipping.")
        plt.close(fig)
        return

    for k in range(len(P_LAB_BINS) - 1):
        means = all_means[k]
        errs  = all_errs[k]
        mask  = all_masks[k]
        label = (rf"$|p|_{{\rm lab}}\in"
                 rf"[{P_LAB_BINS[k]:.0f},{P_LAB_BINS[k+1]:.0f}]$ GeV")

        ax_main.errorbar(
            W_c_global[mask], means[mask], yerr=errs[mask],
            fmt=f"{MARKERS[k]}-", color=PALETTE[k], label=label,
            linewidth=2.2, markersize=7, capsize=3,
        )

        # Ratio to lowest |p_lab| bin
        if k != ref_bin:
            ref_means = all_means[ref_bin]
            with np.errstate(invalid="ignore", divide="ignore"):
                ratio = np.where(ref_means > 0, means / ref_means, np.nan)
                ratio_err = ratio * np.sqrt(
                    (errs / np.where(means > 0, means, 1))**2
                    + (all_errs[ref_bin] / np.where(ref_means > 0, ref_means, 1))**2
                )
            combined_mask = mask & all_masks[ref_bin]
            ax_ratio.errorbar(
                W_c_global[combined_mask], ratio[combined_mask],
                yerr=ratio_err[combined_mask],
                fmt=f"{MARKERS[k]}-", color=PALETTE[k],
                linewidth=2.0, markersize=6, capsize=3,
            )

    ax_ratio.axhline(1.0, color="#888888", linestyle="--", linewidth=1.2)
    ax_ratio.set_xlabel(r"$W$  [GeV]", fontsize=FONT_LABEL)
    ax_ratio.set_ylabel(r"Ratio to $|p|_{\rm ref}$", fontsize=FONT_LABEL - 1)
    ax_ratio.set_ylim(0.5, 2.5)
    ax_ratio.yaxis.set_minor_locator(mticker.AutoMinorLocator())

    ax_main.set_ylabel(r"$\langle n_{90}\rangle$ per jet",
                       fontsize=FONT_LABEL)
    ax_main.set_title(
        r"FFS Ratio: $\langle n_{90}\rangle$ Relative to Lowest $|p|_{\rm lab}$ Bin",
        fontsize=FONT_TITLE - 1, pad=10,
    )
    ax_main.legend(title=r"Fixed lab-frame $|p|$", title_fontsize=FONT_LEGEND,
                   loc="upper left")
    ax_main.set_ylim(bottom=0)
    ax_main.xaxis.set_minor_locator(mticker.AutoMinorLocator())
    ax_main.yaxis.set_minor_locator(mticker.AutoMinorLocator())
    eic_label(ax_main)

    fig.subplots_adjust(left=0.12, right=0.97, top=0.93, bottom=0.10)
    save_fig(fig, os.path.join(outdir, "ffs_ratio.pdf"))
    save_fig(fig, os.path.join(outdir, "ffs_ratio.png"))


# ---------------------------------------------------------------------------
# Figure 3 – DIS kinematic distributions
# ---------------------------------------------------------------------------

def plot_kinematics(hists, outdir):
    """Q² vs W plane + marginal distributions."""
    fig = plt.figure(figsize=(12, 9))
    gs  = GridSpec(2, 3, figure=fig, hspace=0.42, wspace=0.38)

    ax_plane = fig.add_subplot(gs[0, :2])   # Q² vs W (large)
    ax_Q2    = fig.add_subplot(gs[0, 2])
    ax_x     = fig.add_subplot(gs[1, 0])
    ax_y     = fig.add_subplot(gs[1, 1])
    ax_W     = fig.add_subplot(gs[1, 2])

    # ── Q² vs W 2D heatmap ────────────────────────────────────────────────
    h2 = hists["Q2_vs_W"]
    vals  = h2.values().T
    W_e   = h2.axes[0].edges
    Q2_e  = h2.axes[1].edges

    smooth = gaussian_filter(vals.astype(float), sigma=0.8)
    smooth = np.where(smooth < 0.5, np.nan, smooth)

    pcm = ax_plane.pcolormesh(
        W_e, Q2_e, smooth,
        norm=mcolors.LogNorm(vmin=1, vmax=smooth[np.isfinite(smooth)].max() + 1),
        cmap="YlOrRd", shading="auto",
    )
    cb = fig.colorbar(pcm, ax=ax_plane, pad=0.01)
    cb.set_label("Events (arb.)", fontsize=FONT_ANNOT)
    ax_plane.set_xlabel(r"$W$  [GeV]", fontsize=FONT_LABEL)
    ax_plane.set_ylabel(r"$Q^2$  [GeV$^2$]", fontsize=FONT_LABEL)
    ax_plane.set_yscale("log")
    ax_plane.set_title(r"DIS Kinematic Plane: $Q^2$ vs $W$",
                       fontsize=FONT_TITLE - 2)
    eic_label(ax_plane, loc="upper right")

    # ── 1D marginals ──────────────────────────────────────────────────────
    def _fill_1d(ax, h_name, xlabel, log_y=True, color=EIC_TEAL):
        h  = hists[h_name]
        ax.stairs(h.values(), h.axes[0].edges,
                  fill=True, alpha=0.25, color=color, linewidth=0)
        ax.stairs(h.values(), h.axes[0].edges,
                  color=color, linewidth=2.0)
        ax.set_xlabel(xlabel, fontsize=FONT_LABEL - 1)
        ax.set_ylabel("Events", fontsize=FONT_LABEL - 1)
        if log_y:
            ax.set_yscale("log")

    _fill_1d(ax_Q2, "Q2",  r"$Q^2$  [GeV$^2$]",  color=EIC_BLUE)
    _fill_1d(ax_x,  "x",   r"Bjorken $x$",        color=EIC_CORAL)
    _fill_1d(ax_y,  "y",   r"Inelasticity $y$",    color=EIC_GOLD,  log_y=False)
    _fill_1d(ax_W,  "W",   r"$W$  [GeV]",          color=EIC_LILAC, log_y=False)

    ax_Q2.set_xscale("log")
    ax_x.set_xscale("log")

    fig.suptitle("EIC DIS Event Kinematics  (Pythia8 NC-DIS, 10×100 GeV)",
                 fontsize=FONT_TITLE, y=1.01)
    save_fig(fig, os.path.join(outdir, "kinematics.pdf"))
    save_fig(fig, os.path.join(outdir, "kinematics.png"))


# ---------------------------------------------------------------------------
# Figure 4 – Jet landscape
# ---------------------------------------------------------------------------

def plot_jet_landscape(hists, outdir):
    """Jet η vs pT heat-map and jet multiplicity distribution."""
    fig, (ax_etapt, ax_njet) = plt.subplots(1, 2, figsize=(12, 5))

    # ── η–pT heat map ────────────────────────────────────────────────────
    h_ep = hists["jet_eta_pt"]
    vals  = h_ep.values()
    smooth = gaussian_filter(vals.astype(float) + 0.01, sigma=0.7)
    smooth = np.where(smooth < 0.5, np.nan, smooth)

    pcm = ax_etapt.pcolormesh(
        h_ep.axes[0].edges,
        h_ep.axes[1].edges,
        smooth.T,
        norm=mcolors.LogNorm(vmin=1,
                             vmax=smooth[np.isfinite(smooth)].max() + 1),
        cmap="Blues", shading="auto",
    )
    cb = fig.colorbar(pcm, ax=ax_etapt, pad=0.01)
    cb.set_label("Jets (arb.)", fontsize=FONT_ANNOT)

    # EIC detector acceptance overlay
    ax_etapt.axvline(-3.5, color=EIC_CORAL, linestyle="--",
                     linewidth=1.5, label="EIC acceptance  |η| < 3.5")
    ax_etapt.axvline( 3.5, color=EIC_CORAL, linestyle="--", linewidth=1.5)
    ax_etapt.legend(fontsize=FONT_ANNOT, loc="upper right")
    ax_etapt.set_xlabel(r"Jet $\eta_{\rm lab}$", fontsize=FONT_LABEL)
    ax_etapt.set_ylabel(r"Jet $p_T$  [GeV]", fontsize=FONT_LABEL)
    ax_etapt.set_title(r"Lab-frame Jet Landscape ($\eta$ vs $p_T$)",
                       fontsize=FONT_TITLE - 2)

    # ── N_jets distribution ──────────────────────────────────────────────
    h_nj = hists["n_jets"]
    centers = h_nj.axes[0].centers
    values  = h_nj.values()
    mask = values > 0
    ax_njet.bar(centers[mask], values[mask],
                width=0.85, color=EIC_TEAL, alpha=0.85, edgecolor="white",
                linewidth=0.6)
    ax_njet.set_xlabel(r"$N_{\rm jets}$ per event", fontsize=FONT_LABEL)
    ax_njet.set_ylabel("Events", fontsize=FONT_LABEL)
    ax_njet.set_title("Anti-$k_T$ Jet Multiplicity per Event",
                      fontsize=FONT_TITLE - 2)
    ax_njet.set_yscale("log")
    ax_njet.xaxis.set_major_locator(mticker.MaxNLocator(integer=True))
    eic_label(ax_njet, loc="upper right")

    fig.suptitle("EIC Jet Observables  (anti-$k_T$, $R=0.4$, $p_T > 2$ GeV)",
                 fontsize=FONT_TITLE, y=1.02)
    fig.tight_layout()
    save_fig(fig, os.path.join(outdir, "jet_landscape.pdf"))
    save_fig(fig, os.path.join(outdir, "jet_landscape.png"))


# ---------------------------------------------------------------------------
# Figure 5 – FFS heat-map: ⟨n₉₀⟩(|p_lab|, W)
# ---------------------------------------------------------------------------

def plot_ffs_heatmap(hists, outdir):
    """2D colour map of ⟨n₉₀⟩ as a function of (W, |p_lab|)."""
    h3 = hists["n90_3d"]

    W_edges  = h3.axes[0].edges
    PL_edges = h3.axes[1].edges
    N90_c    = h3.axes[2].centers

    vals = h3.values()    # shape (nW, nPL, nN90)

    nW  = len(W_edges) - 1
    nPL = len(PL_edges) - 1

    mean_grid  = np.full((nW, nPL), np.nan)
    count_grid = np.zeros((nW, nPL))

    for iw in range(nW):
        for ip in range(nPL):
            counts = vals[iw, ip, :]
            tot    = counts.sum()
            if tot > 10:
                mean_grid[iw, ip]  = float(np.sum(counts * N90_c) / tot)
                count_grid[iw, ip] = tot

    fig, ax = plt.subplots(figsize=(8, 5))

    # Smooth the grid (only filled cells)
    filled = np.isfinite(mean_grid)
    smooth_grid = np.where(filled, mean_grid, 0.0)
    # Light smoothing only if we have enough cells
    if filled.sum() >= 4:
        smooth_grid = gaussian_filter(smooth_grid, sigma=0.5)
        smooth_grid = np.where(filled, smooth_grid, np.nan)

    vmin = np.nanmin(smooth_grid) if np.any(np.isfinite(smooth_grid)) else 0
    vmax = np.nanmax(smooth_grid) if np.any(np.isfinite(smooth_grid)) else 1

    pcm = ax.pcolormesh(
        W_edges, PL_edges, smooth_grid.T,
        cmap="viridis",
        vmin=max(vmin - 0.5, 0),
        vmax=vmax + 0.5,
        shading="flat",
    )
    cb = fig.colorbar(pcm, ax=ax, pad=0.02)
    cb.set_label(r"$\langle n_{90}\rangle$ per jet", fontsize=FONT_LABEL)

    # Contour lines for aesthetics
    W_c  = 0.5 * (W_edges[:-1] + W_edges[1:])
    PL_c = 0.5 * (PL_edges[:-1] + PL_edges[1:])
    if np.any(np.isfinite(smooth_grid)):
        try:
            ax.contour(
                W_c, PL_c, smooth_grid.T,
                levels=5, colors="white", linewidths=0.8, alpha=0.6,
            )
        except Exception:
            pass

    ax.set_xlabel(r"$W$  [GeV]", fontsize=FONT_LABEL)
    ax.set_ylabel(r"Jet $|p|_{\rm lab}$  [GeV]", fontsize=FONT_LABEL)
    ax.set_title(
        r"FFS Heat Map: $\langle n_{90}\rangle$ vs ($W$, $|p|_{\rm lab}$)",
        fontsize=FONT_TITLE - 1, pad=10,
    )
    ax.xaxis.set_minor_locator(mticker.AutoMinorLocator())
    ax.yaxis.set_minor_locator(mticker.AutoMinorLocator())
    eic_label(ax, loc="upper left")
    arxiv_label(ax)

    fig.tight_layout()
    save_fig(fig, os.path.join(outdir, "ffs_heatmap.pdf"))
    save_fig(fig, os.path.join(outdir, "ffs_heatmap.png"))


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Make publication-quality FFS effect plots",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("input", type=str,
                   help="Input ROOT file from analyze_events.py")
    p.add_argument("--outdir", type=str, default="plots/",
                   help="Output directory for figures")
    p.add_argument("--format", type=str, default="pdf",
                   choices=["pdf", "png", "svg"],
                   help="Primary figure format (PNG is always also saved)")
    return p.parse_args()


def main():
    args = parse_args()
    apply_style()
    os.makedirs(args.outdir, exist_ok=True)

    print(f"Loading histograms from {args.input} …", flush=True)
    hists = {}
    with uproot.open(args.input) as f:
        for key in f.keys():
            name = key.rstrip(";1")
            try:
                hists[name] = f[key].to_hist()
            except Exception as exc:
                print(f"  Warning: could not load {key}: {exc}")

    required = {"n90_3d", "sum_n90_vs_W", "count_n90_vs_W",
                "mult_3d", "sum_N_vs_W", "count_vs_W",
                "Q2", "x", "y", "W", "Q2_vs_W", "jet_eta_pt", "n_jets"}
    missing = required - set(hists)
    if missing:
        sys.exit(f"ERROR: missing histograms in ROOT file: {missing}")

    print("Generating figures …", flush=True)
    plot_ffs_main(hists, args.outdir)
    plot_ffs_ratio(hists, args.outdir)
    plot_kinematics(hists, args.outdir)
    plot_jet_landscape(hists, args.outdir)
    plot_ffs_heatmap(hists, args.outdir)

    print(f"\nAll figures written to {args.outdir}/")


if __name__ == "__main__":
    main()
