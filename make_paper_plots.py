#!/usr/bin/env python3
"""
PRL figure stage for the EIC FFS study (ANALYSIS_DESIGN.md Sec. 3).

Reads results/results.json (from make_results.py) and produces:

  fig1_concept.pdf       boost factor + string-anchor mapping (reach)
  fig2_ffs_main.pdf      <n90_lab> vs W at fixed (|p|_lab, Q^2): the effect
  fig3_universality.pdf  lab-frame splay vs CM-frame collapse
  fig4_projection.pdf    H0-rejection significance vs integrated luminosity
  fig5_generators.pdf    string vs cluster hadronization comparison
  fig6_lhc.pdf           universal curve extrapolated to LHC scales:
                         the jet-individualism bias
  fig7_decomposition.pdf frame-dependence decomposition: matched-cell
                         contrasts prove structure = f(p_CM), not f(p_lab)

Typography: Tufte-style (high data-ink ratio, despined axes, frameless
legends, direct labels); text in EB Garamond, math in STIX serif.

Usage:  python make_paper_plots.py [--results results/results.json] [--outdir plots]
"""

import argparse
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── Tufte-style serif theme (EB Garamond) ───────────────────────────────────
# Register Garamond faces directly (robust against stale matplotlib caches).
import glob as _glob
import matplotlib.font_manager as _fm
for _f in _glob.glob("/usr/share/fonts/opentype/ebgaramond/*.otf"):
    try:
        _fm.fontManager.addfont(_f)
    except Exception:
        pass
_GARAMOND = next((n for n in ("EB Garamond", "EB Garamond 12")
                  if any(f.name == n for f in _fm.fontManager.ttflist)),
                 "serif")

# High data-ink ratio: no top/right spines, muted axes, frameless legends,
# direct labeling where possible, no titles (captions carry the prose).
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": [_GARAMOND, "Garamond", "DejaVu Serif"],
    # math in STIX: a serif face that pairs with Garamond and has full
    # Greek/symbol coverage (EB Garamond itself has no Greek glyphs)
    "mathtext.fontset": "stix",
    "figure.dpi": 150,
    "font.size": 13,
    "axes.labelsize": 14,
    "axes.titlesize": 13,
    "legend.fontsize": 11,
    "xtick.labelsize": 11.5,
    "ytick.labelsize": 11.5,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.7,
    "axes.edgecolor": "0.3",
    "xtick.color": "0.3",
    "ytick.color": "0.3",
    "xtick.labelcolor": "0.15",
    "ytick.labelcolor": "0.15",
    "axes.labelcolor": "0.1",
    "text.color": "0.1",
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.major.size": 3.2,
    "ytick.major.size": 3.2,
    "xtick.minor.size": 1.8,
    "ytick.minor.size": 1.8,
    "legend.frameon": False,
    "lines.solid_capstyle": "round",
    "axes.unicode_minus": False,     # EB Garamond lacks U+2212
})


def tufte(ax):
    """Lighten an axes beyond the rcParams: nudge spines outward."""
    for side in ("left", "bottom"):
        ax.spines[side].set_position(("outward", 4))

CONFIG_MARKERS = {"5x41": "o", "10x100": "s", "18x275": "^"}
CONFIG_LABELS = {"5x41": r"$5\times41$", "10x100": r"$10\times100$",
                 "18x275": r"$18\times275$"}
CONFIG_COLORS = {"5x41": "#4477AA", "10x100": "#228833", "18x275": "#AA3377"}

EIC_LABEL = "EIC NC-DIS, Pythia 8"
SEL_LABEL = (r"$25 < Q^2 < 100$ GeV$^2$, current-hemisphere jet"
             "\n" + r"ee-gen-$k_T$ $R=1.0$ in $\gamma^*p$ frame")


def _saveplot(fig, outdir, name):
    path = os.path.join(outdir, name)
    fig.savefig(path, bbox_inches="tight")
    fig.savefig(path.replace(".pdf", ".png"), bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {path}")


# ---------------------------------------------------------------------------
# Fig. 1 — concept & reach
# ---------------------------------------------------------------------------

def fig1(res, outdir):
    bmap = res["fig1"]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.2))

    for config, rows in bmap.items():
        W = [r["W_center"] for r in rows]
        g = [r["gamma_mean"] for r in rows]
        ax1.plot(W, g, marker=CONFIG_MARKERS[config], ms=3.5, lw=1.1,
                 color=CONFIG_COLORS[config])
        # direct label at the curve end (Tufte: no legend box)
        ax1.annotate(CONFIG_LABELS[config] + " GeV", (W[-1], g[-1]),
                     xytext=(6, 0), textcoords="offset points",
                     color=CONFIG_COLORS[config], fontsize=11, va="center")
    ax1.set_xlabel(r"$W$ [GeV]")
    ax1.set_ylabel(r"$\langle\gamma\rangle$ of $\gamma^*p$ frame in lab")
    ax1.set_yscale("log")
    ax1.set_yticks([2, 3, 5, 10])
    ax1.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax1.get_yaxis().set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax1.set_xlim(right=150)
    ax1.text(0.97, 0.92, "fragmentation frame is\nstrongly boosted in the lab",
             transform=ax1.transAxes, fontsize=11, style="italic", color="0.35",
             ha="right")
    tufte(ax1)

    for config, rows in bmap.items():
        W = np.array([r["W_center"] for r in rows])
        ax2.plot(W, [r["pcm_mean"] for r in rows], marker=CONFIG_MARKERS[config],
                 color=CONFIG_COLORS[config], ls="-", ms=3.5, lw=1.1)
        ax2.plot(W, [r["plab_mean"] for r in rows], marker=CONFIG_MARKERS[config],
                 color=CONFIG_COLORS[config], ls=":", alpha=0.45, mfc="none",
                 ms=3.5, lw=1.0)
    Wd = np.linspace(8, 135, 50)
    ax2.plot(Wd, Wd / 2, ls="--", lw=0.9, color="0.2")
    ax2.annotate(r"$W/2$ string anchor", (Wd[-1], Wd[-1] / 2), xytext=(6, 0),
                 textcoords="offset points", fontsize=11, va="center",
                 color="0.2")
    ax2.set_xlabel(r"$W$ [GeV]")
    ax2.set_ylabel(r"jet momentum [GeV]")
    ax2.set_xlim(right=150)
    ax2.text(0.04, 0.84, r"filled: $\langle p_{\rm CM}\rangle$"
             "\n" + r"open: $\langle |p|_{\rm lab}\rangle$",
             transform=ax2.transAxes, fontsize=11, color="0.35")
    ax2.text(0.40, 0.06,
             "CM momentum locked to $W/2$;\nlab momentum is a frame artifact",
             transform=ax2.transAxes, fontsize=11, style="italic", color="0.35")
    ax2.text(0.04, 0.95, EIC_LABEL, transform=ax2.transAxes, fontsize=10,
             color="0.45")
    tufte(ax2)
    _saveplot(fig, outdir, "fig1_concept.pdf")


# ---------------------------------------------------------------------------
# Fig. 2 — the FFS effect
# ---------------------------------------------------------------------------

def fig2(res, outdir):
    f2 = res["fig2"]
    anchors = res["h0_anchors"]
    # drop the sparse highest-|p| bin from the figure (kept in results.json)
    plab_keys = [f"{lo:g}-{hi:g}" for lo, hi in res["binnings"]["plab_bins"]][:4]
    cmap = plt.get_cmap("plasma")
    pcolors = {k: cmap(0.12 + 0.75 * i / max(len(plab_keys) - 1, 1))
               for i, k in enumerate(plab_keys)}

    # variation envelope per (plab bin, W bin): min/max across variations
    def envelope(pkey):
        per_w = {}
        for skey, entry in f2.items():
            config, variation, level = skey.split("|")
            if level != "truth" or variation == "baseline_mpi":
                continue
            for r in entry.get(pkey, []):
                per_w.setdefault(r["W_center"], []).append(r["mean"])
        Ws = sorted(per_w)
        return (np.array(Ws), np.array([min(per_w[w]) for w in Ws]),
                np.array([max(per_w[w]) for w in Ws]))

    fig, (ax, axr) = plt.subplots(
        2, 1, figsize=(8.5, 8.5), sharex=True,
        gridspec_kw={"height_ratios": [2.6, 1], "hspace": 0.05})

    shown_p = set()
    for pkey in plab_keys:
        if pkey not in anchors:
            continue
        anc = anchors[pkey]["value"]
        Ws, lo, hi = envelope(pkey)
        if len(Ws) >= 3:
            ax.fill_between(Ws, lo, hi, color=pcolors[pkey], alpha=0.12, lw=0)
            axr.fill_between(Ws, lo / anc, hi / anc, color=pcolors[pkey],
                             alpha=0.12, lw=0)
        for skey, entry in f2.items():
            config, variation, level = skey.split("|")
            if variation != "baseline" or level != "truth":
                continue
            rows = entry.get(pkey, [])
            if not rows:
                continue
            W = [r["W_center"] for r in rows]
            m = [r["mean"] for r in rows]
            e = [r["sem"] for r in rows]
            ax.errorbar(W, m, yerr=e, marker=CONFIG_MARKERS[config], ms=3.5,
                        color=pcolors[pkey], ls="-", lw=1.0, capsize=0,
                        elinewidth=0.8)
            axr.errorbar(W, np.array(m) / anc, yerr=np.array(e) / anc,
                         marker=CONFIG_MARKERS[config], ms=3.5,
                         color=pcolors[pkey], ls="-", lw=1.0, elinewidth=0.8)
            shown_p.add(pkey)
        # H0: frame-independent null
        ax.axhline(anc, color=pcolors[pkey], ls="--", lw=0.8, alpha=0.5)

    axr.axhline(1.0, color="0.2", ls="--", lw=0.8)
    ax.set_ylabel(r"$\langle n_{90}\rangle$  (lab frame)")
    axr.set_ylabel(r"ratio to $H_0$")
    axr.set_xlabel(r"$W$ [GeV]")
    ax.set_xscale("log")
    axr.set_xscale("log")
    axr.set_xticks([10, 20, 30, 50, 100])
    axr.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

    from matplotlib.lines import Line2D
    handles = [Line2D([], [], color=pcolors[k], marker="none", lw=2.5,
                      label=rf"$|p|_{{\rm lab}} \in [{k.replace('-', ',')})$ GeV")
               for k in plab_keys if k in shown_p]
    handles += [Line2D([], [], color="gray", marker=CONFIG_MARKERS[c], ls="",
                       label=CONFIG_LABELS[c] + " GeV") for c in CONFIG_MARKERS]
    handles += [Line2D([], [], color="gray", ls="--",
                       label=r"$H_0$: frame-indep. frag."),
                ]
    ax.legend(handles=handles, ncol=3, fontsize=10, handlelength=1.5,
              columnspacing=1.0, labelspacing=0.35, loc="lower left",
              bbox_to_anchor=(-0.02, 1.01, 1.04, 0.12), mode="expand",
              borderaxespad=0)
    lumi_note = ""
    li = res.get("fig4", {}).get("luminosity_inputs", {}).get("per_config", {})
    if li:
        mc_pb = min(v["mc_equiv_fb"] for v in li.values()) * 1e3
        ann = [v["annual_fb"] for v in li.values() if v.get("annual_fb")]
        lumi_note = (f"\nMC stats $\\approx$ first {mc_pb:.0f} pb$^{{-1}}$; "
                     f"1 EIC yr = {min(ann):.0f}--{max(ann):.0f} fb$^{{-1}}$")
    ax.text(0.98, 0.03, EIC_LABEL + "\n" + SEL_LABEL + lumi_note,
            transform=ax.transAxes, fontsize=9, ha="right", va="bottom",
            color="0.45")
    ax.set_ylim(top=8.4)
    tufte(ax)
    tufte(axr)
    _saveplot(fig, outdir, "fig2_ffs_main.pdf")


# ---------------------------------------------------------------------------
# Fig. 3 — universality: splay vs collapse
# ---------------------------------------------------------------------------

def fig3(res, outdir):
    f3 = res["fig3"]
    # color by W-slice center, common normalisation
    all_keys = sorted(set(f3["splay"]) | set(f3["collapse"]),
                      key=lambda k: float(k.split("|W")[1].split("-")[0]))

    def wcenter(key):
        lo, hi = key.split("|W")[1].split("-")
        return 0.5 * (float(lo) + float(hi))

    wvals = [wcenter(k) for k in all_keys]
    norm = matplotlib.colors.LogNorm(vmin=min(wvals) * 0.9,
                                     vmax=max(wvals) * 1.1)
    cmap = plt.get_cmap("viridis")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.5, 5.6), sharey=True)

    for ax, part, xlab, obs in (
            (ax1, "splay", r"$|p|_{\rm lab}$ [GeV]",
             r"$\langle n_{90}\rangle$ (lab frame)"),
            (ax2, "collapse", r"$p_{\rm CM}$ [GeV]",
             r"$\langle n_{90}\rangle$ (CM frame)")):
        for key in all_keys:
            rows = f3[part].get(key)
            if not rows:
                continue
            # same bulk filter as the universality metric (drop sparse tails)
            bulk = res["binnings"].get("bulk_fraction", 0.0)
            n_max = max(r["n"] for r in rows)
            rows = [r for r in rows if r["n"] >= bulk * n_max]
            config = key.split("|")[0]
            x = [r["x_mean"] for r in rows]
            m = [r["mean"] for r in rows]
            e = [r["sem"] for r in rows]
            ax.errorbar(x, m, yerr=e, marker=CONFIG_MARKERS[config], ms=4.5,
                        ls="-", lw=1.1, elinewidth=0.8,
                        color=cmap(norm(wcenter(key))))
        ax.set_xlabel(xlab)
        ax.set_xscale("log")
        ax.set_xticks([3, 5, 10, 20, 40])
        ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

    ax1.set_ylabel(r"$\langle n_{90}\rangle$")
    ax1.set_title("binned in the lab frame", fontsize=13, color="0.15",
                  loc="left")
    ax2.set_title(r"binned in the $\gamma^*p$ (color) frame", fontsize=13,
                  color="0.15", loc="left")
    ax1.text(0.04, 0.90, "universality broken:\n"
             rf"slice spread {100*f3['spread_lab']:.0f}% RMS",
             transform=ax1.transAxes, fontsize=12, style="italic", color="0.25")
    ax2.text(0.04, 0.90, "universality restored:\n"
             rf"slice spread {100*f3['spread_cm']:.1f}% RMS",
             transform=ax2.transAxes, fontsize=12, style="italic", color="0.25")

    from matplotlib.lines import Line2D
    handles = [Line2D([], [], color="0.4", marker=CONFIG_MARKERS[c], ls="",
                      ms=4.5, label=CONFIG_LABELS[c] + " GeV")
               for c in CONFIG_MARKERS]
    ax2.legend(handles=handles, loc="lower right", fontsize=10.5,
               labelspacing=0.35)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    cbar = fig.colorbar(sm, ax=[ax1, ax2], pad=0.012, aspect=40, shrink=0.95)
    cbar.set_label(r"$W$ slice center [GeV]")
    cbar.outline.set_visible(False)
    cbar.ax.tick_params(size=0)
    ax1.text(0.04, 0.05, EIC_LABEL, transform=ax1.transAxes, fontsize=10,
             color="0.45")
    tufte(ax1)
    tufte(ax2)
    _saveplot(fig, outdir, "fig3_universality.pdf")


# ---------------------------------------------------------------------------
# Fig. 4 — luminosity projection
# ---------------------------------------------------------------------------

def fig4(res, outdir):
    proj = res["fig4"]
    L = np.array(proj["lumi_grid"])
    fig, ax = plt.subplots(figsize=(8, 6))

    lumi_info = proj.get("luminosity_inputs", {}).get("per_config", {})
    for config, v in proj["configs"].items():
        ax.plot(L, v["significance_syst"], color=CONFIG_COLORS[config], lw=1.1,
                alpha=0.75, label=CONFIG_LABELS[config] + " GeV")
        ann = lumi_info.get(config, {}).get("annual_fb")
        if ann:
            for Lmark, mk, ms in ((ann / 365.25, "D", 5), (ann, "*", 12)):
                sig_m = np.interp(np.log10(Lmark), np.log10(L),
                                  v["significance_syst"])
                ax.plot(Lmark, sig_m, mk, color=CONFIG_COLORS[config],
                        ms=ms, mec="0.2", mew=0.5, zorder=5)
    syst = res["binnings"]["syst_floor"]
    cons = res["binnings"].get("syst_floor_cons", 0.05)
    ax.plot(L, proj["combined"]["significance_syst"], color="0.1", lw=2.0,
            label=f"combined ({100*syst:g}% syst floor)")
    ax.plot(L, proj["combined"]["significance_syst_cons"], color="0.1",
            ls="-.", lw=1.3, label=f"combined ({100*cons:g}% syst floor)")
    ax.plot(L, proj["combined"]["significance_stat"], color="0.1", ls="--",
            lw=1.1, label="combined (stat only)")

    ax.axhline(5, color="#9e2a2b", ls=":", lw=1.1, alpha=0.9)
    ax.text(L[-1] * 0.6, 5.5, r"$5\sigma$", color="#9e2a2b", fontsize=12,
            ha="right")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xticks([1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 100])
    ax.get_xaxis().set_major_formatter(
        matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    ax.set_yticks([10, 100, 1000])
    ax.get_yaxis().set_major_formatter(
        matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    ax.set_xlabel(r"integrated luminosity per configuration [fb$^{-1}$]")
    ax.set_ylabel(r"expected significance of $H_0$ rejection [$\sigma$]")
    from matplotlib.lines import Line2D
    handles, labels = ax.get_legend_handles_labels()
    handles += [Line2D([], [], marker="D", color="0.5", ls="", ms=5,
                       mec="0.2", label="1 day at design lumi."),
                Line2D([], [], marker="*", color="0.5", ls="", ms=11,
                       mec="0.2", label="1 year at design lumi.")]
    ax.legend(handles=handles, loc="upper left", fontsize=10,
              labelspacing=0.4,
              title="reco level: track jets, smeared,\nelectron-method $W$",
              title_fontsize=9.5)
    ax.text(0.97, 0.04, EIC_LABEL +
            "\ndesign lumi: EIC YR Tab. 10.1", transform=ax.transAxes,
            fontsize=9.5, ha="right", color="0.4")
    tufte(ax)
    _saveplot(fig, outdir, "fig4_projection.pdf")


# ---------------------------------------------------------------------------
# Fig. 5 — hadronization-model comparison (Lund string vs cluster)
# ---------------------------------------------------------------------------

def fig5(res, outdir):
    f2 = res["fig2"]
    have_herwig = any(k.split("|")[1] == "herwig" for k in f2)
    if not have_herwig:
        print("  (no Herwig samples; skipping fig5)")
        return

    pkeys = [f"{lo:g}-{hi:g}" for lo, hi in res["binnings"]["plab_bins"]][:4]
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), sharex=True)
    variants = [("baseline", "-", "o", "Pythia 8 (string, Monash)"),
                ("lund-soft", ":", "", "Pythia 8 string variations"),
                ("lund-hard", ":", "", None),
                ("herwig", "--", "s", "Herwig 7.3 (cluster)")]

    for ax, pkey in zip(axes.ravel(), pkeys):
        for variation, ls, mk, label in variants:
            shown = False
            for config in CONFIG_MARKERS:
                rows = f2.get(f"{config}|{variation}|truth", {}).get(pkey, [])
                if not rows:
                    continue
                color = ("#CC3311" if variation == "herwig" else
                         "#004488" if variation == "baseline" else "#88AAcc")
                ax.errorbar([r["W_center"] for r in rows],
                            [r["mean"] for r in rows],
                            yerr=[r["sem"] for r in rows],
                            ls=ls, marker=mk, ms=3, lw=1.1, elinewidth=0.8,
                            color=color,
                            label=(label if (label and not shown) else None))
                shown = True
        ax.set_xscale("log")
        ax.text(0.05, 0.90, rf"$|p|_{{\rm lab}} \in [{pkey.replace('-', ',')})$ GeV",
                transform=ax.transAxes, fontsize=12, color="0.15")
        tufte(ax)
    for ax in axes[1]:
        ax.set_xlabel(r"$W$ [GeV]")
        ax.set_xticks([10, 20, 50, 100])
        ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    for ax in axes[:, 0]:
        ax.set_ylabel(r"$\langle n_{90}\rangle$ (lab frame)")
    axes[0, 0].legend(fontsize=10.5, loc="lower right", labelspacing=0.4)
    axes[0, 1].text(0.95, 0.05, EIC_LABEL + "\n" + SEL_LABEL,
                    transform=axes[0, 1].transAxes, fontsize=9, ha="right",
                    color="0.4")
    _saveplot(fig, outdir, "fig5_generators.pdf")


# ---------------------------------------------------------------------------
# Fig. 6 — implication for the LHC: the jet-individualism bias
# ---------------------------------------------------------------------------

def fig6(res, outdir):
    """
    The EIC-measured universal color-frame curve <n90>(p_CM), with its
    logarithmic fit extrapolated to LHC scales.  Under jet individualism
    (substructure = f(lab pT)), a W-daughter jet and a QCD jet at the same
    lab pT are identical; the FFS prediction is that they sit at very
    different points of this curve (p_CM ~ m_V/2 vs p_CM ~ pT), and the
    vertical gap is the bias inherited by every multiplicity-sensitive
    tagger trained on one topology and applied to the other.
    """
    f3 = res["fig3"]
    bulk = res["binnings"].get("bulk_fraction", 0.0)

    # global weighted mean per p_CM bin over all (config, W-slice) profiles
    per_bin = {}
    for rows in f3["collapse"].values():
        n_max = max(r["n"] for r in rows)
        for r in rows:
            if r["n"] >= bulk * n_max:
                per_bin.setdefault(r["x_center"], []).append(r)
    pts = []
    for c, rs in sorted(per_bin.items()):
        w = np.array([1.0 / r["sem"]**2 for r in rs])
        m = np.array([r["mean"] for r in rs])
        x = np.array([r["x_mean"] for r in rs])
        pts.append((np.average(x, weights=w), np.average(m, weights=w),
                    1.0 / np.sqrt(w.sum())))
    x, y, ye = map(np.array, zip(*pts))

    # weighted fit  n90 = a + b ln(p_CM)
    A = np.vstack([np.ones_like(x), np.log(x)]).T
    Wm = np.diag(1.0 / ye**2)
    cov = np.linalg.inv(A.T @ Wm @ A)
    a, b = cov @ A.T @ Wm @ y

    fig, ax = plt.subplots(figsize=(8.6, 5.6))
    ax.errorbar(x, y, yerr=ye, ls="", marker="o", ms=4, color="0.2",
                elinewidth=0.8)
    xf = np.linspace(x.min(), x.max(), 50)
    ax.plot(xf, a + b * np.log(xf), color="#9e2a2b", lw=1.4)
    xe = np.geomspace(x.max(), 1500, 60)
    ax.plot(xe, a + b * np.log(xe), color="#9e2a2b", lw=1.2, ls="--",
            alpha=0.8)
    ax.text(7.5, a + b * np.log(7.5) + 0.45,
            "EIC measurement\n(this campaign)", fontsize=11, color="0.25",
            ha="center", style="italic")
    ax.text(0.97, 0.06,
            rf"fit: $\langle n_{{90}}\rangle = {a:.2f} + {b:.2f}\,\ln p_{{\rm CM}}$"
            "\n(dashed: extrapolation)", fontsize=11, color="#9e2a2b",
            ha="right", transform=ax.transAxes)

    # LHC reference points: color-frame momentum of the fragmenting system
    refs = [
        (45.0, (10, -16), "right",
         "$W/Z$-daughter jet\n($p_{\\rm CM}=m_V/2$, any lab $p_T$)"),
        (200.0, (-8, 6), "right", "QCD jet, $p_T=200$ GeV\n($p_{\\rm CM}\\approx p_T$)"),
        (1000.0, (-8, 2), "right", "QCD jet, $p_T=1$ TeV"),
    ]
    for pcm, (dx, dy), ha, label in refs:
        yv = a + b * np.log(pcm)
        ax.plot(pcm, yv, "s", ms=6, color="#1d3557", mec="0.15", mew=0.5,
                zorder=5)
        ax.annotate(label, (pcm, yv), xytext=(dx, dy),
                    textcoords="offset points", fontsize=10,
                    ha="left" if dx > 0 else "right",
                    va="bottom" if dy > 0 else "top", color="#1d3557")

    # the individualism bias: W-jet vs QCD jet at the same lab pT = 200 GeV
    y_w, y_q = a + b * np.log(45.0), a + b * np.log(200.0)
    ax.plot([45, 200], [y_w, y_w], color="0.6", lw=0.7, ls=":")
    ax.annotate("", xy=(200, y_q), xytext=(200, y_w),
                arrowprops=dict(arrowstyle="<->", color="0.25", lw=0.9))
    ax.text(178, 0.5 * (y_w + y_q),
            "jet-individualism bias\nat the same lab $p_T$:\n"
            rf"$\Delta\langle n_{{90}}\rangle = {y_q - y_w:.2f}"
            rf"\ \approx {100 * (y_q / y_w - 1):.0f}\%$",
            fontsize=11, color="0.2", va="center", ha="right")

    ax.set_xscale("log")
    ax.set_xlabel(r"color-frame jet momentum $p_{\rm CM}$ [GeV]")
    ax.set_ylabel(r"$\langle n_{90}\rangle$ (color frame)")
    ax.set_xticks([5, 10, 50, 100, 500, 1000])
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    tufte(ax)
    _saveplot(fig, outdir, "fig6_lhc.pdf")
    print(f"  fig6 fit: n90 = {a:.3f} + {b:.3f} ln(pcm); "
          f"W-jet vs QCD@200: delta = {y_q - y_w:.2f} ({100*(y_q/y_w-1):.0f}%); "
          f"QCD@1TeV vs W-jet: {100*((a + b*np.log(1000))/y_w - 1):.0f}%")



# ---------------------------------------------------------------------------
# Fig. 7 — frame-dependence decomposition: p_CM, not |p|_lab
# ---------------------------------------------------------------------------

def fig7(res, outdir):
    """
    Direct demonstration that jet structure depends on the color-frame
    momentum and not the lab-frame momentum.

    (a) <n90_lab> on the (p_CM, |p|_lab) plane at fixed Q^2: vertical
        banding -- iso-structure lines are lines of constant p_CM.
    (b) matched-cell contrast: within narrow (p_CM, Q^2) cells, jets split
        at the median |p|_lab (flat segments: nothing happens), mirrored
        by (|p|_lab, Q^2) cells split in p_CM (steep segments).  The
        inclusive |p|_lab trend (dotted) even has the opposite sign --
        Simpson's paradox: lab momentum has no intrinsic explanatory power.
    """
    f7 = res.get("fig7")
    if not f7:
        print("  (no fig7 data; skipping)")
        return

    pcm_e = np.array(f7["pcm_edges"])
    plab_e = np.array(f7["plab_edges"])
    cmap = plt.get_cmap("viridis")
    norm = matplotlib.colors.LogNorm(vmin=pcm_e[0], vmax=pcm_e[-1])

    fig = plt.figure(figsize=(14.2, 5.4))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.5, 1, 1], wspace=0.26)
    axm = fig.add_subplot(gs[0])
    axl = fig.add_subplot(gs[1])
    axc = fig.add_subplot(gs[2], sharey=axl)

    # ── (a) the landscape ──────────────────────────────────────────────────
    Z = np.full((len(plab_e) - 1, len(pcm_e) - 1), np.nan)
    for r in f7["map"]:
        Z[r["j"], r["i"]] = r["mean"]
    pc = axm.pcolormesh(pcm_e, plab_e, np.ma.masked_invalid(Z),
                        cmap="viridis", shading="flat", rasterized=True)
    axm.set_xscale("log")
    axm.set_yscale("log")
    axm.set_xlabel(r"$p_{\rm CM}$ [GeV]")
    axm.set_ylabel(r"$|p|_{\rm lab}$ [GeV]")
    for axis in (axm.get_xaxis(), axm.get_yaxis()):
        axis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
        axis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    axm.set_xticks([5, 10, 20, 40])
    axm.set_yticks([5, 10, 20, 40, 70])
    cb = fig.colorbar(pc, ax=axm, pad=0.02, aspect=30)
    cb.set_label(r"$\langle n_{90}\rangle$ (lab frame)")
    cb.outline.set_visible(False)
    cb.ax.tick_params(size=0)
    q2lo, q2hi = f7["q2_narrow"]
    axm.text(0.03, 0.96, rf"${q2lo:g} < Q^2 < {q2hi:g}$ GeV$^2$",
             transform=axm.transAxes, fontsize=10.5, color="0.3", va="top")
    axm.text(0.03, 0.04, "vertical banding:\nstructure set by $p_{\\rm CM}$ alone",
             transform=axm.transAxes, fontsize=11, color="white",
             style="italic", va="bottom")

    # ── (b) matched-cell contrasts ─────────────────────────────────────────
    con = f7["contrasts"]["baseline"]

    for seg in con["vary_plab"]["segments"]:
        c = cmap(norm(np.sqrt(seg["cell"][0] * seg["cell"][1])))
        axl.plot([seg["x_lo"], seg["x_hi"]], [seg["y_lo"], seg["y_hi"]],
                 "-", color=c, lw=1.6, marker="o", ms=3, mec="none",
                 alpha=0.95)
    rows = f7["inclusive"]["rows"]
    axl.plot([r["plab"] for r in rows], [r["mean"] for r in rows],
             ls=":", color="0.25", lw=1.4)
    axl.annotate("inclusive\n(p$_{\\rm CM}$ uncontrolled)",
                 (rows[-1]["plab"], rows[-1]["mean"]), xytext=(0, -26),
                 textcoords="offset points", fontsize=9.5, color="0.3",
                 ha="right")
    ag = con["vary_plab"]["aggregate"]
    axl.set_title("vary $|p|_{\\rm lab}$ at fixed $(p_{\\rm CM}, Q^2)$",
                  fontsize=12.5, color="0.15", loc="left")
    axl.text(0.04, 0.97, "nothing happens:\n"
             rf"slope $= {ag['mean']:+.3f} \pm {ag['stat_err']:.3f}$"
             rf"$\,\pm\,{ag['cell_rms']:.2f}$ (cells)",
             transform=axl.transAxes, fontsize=11, va="top", style="italic",
             color="0.25")
    axl.set_xscale("log")
    axl.set_xlabel(r"$|p|_{\rm lab}$ [GeV]")
    axl.set_ylabel(r"$\langle n_{90}\rangle$ (lab frame)")
    axl.set_xticks([5, 10, 20, 40])
    axl.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    axl.get_xaxis().set_minor_formatter(matplotlib.ticker.NullFormatter())

    for seg in con["vary_pcm"]["segments"]:
        c = cmap(norm(np.sqrt(seg["x_lo"] * seg["x_hi"])))
        axc.plot([seg["x_lo"], seg["x_hi"]], [seg["y_lo"], seg["y_hi"]],
                 "-", color=c, lw=1.6, marker="o", ms=3, mec="none",
                 alpha=0.95)
    ag = con["vary_pcm"]["aggregate"]
    axc.set_title("vary $p_{\\rm CM}$ at fixed $(|p|_{\\rm lab}, Q^2)$",
                  fontsize=12.5, color="0.15", loc="left")
    axc.text(0.04, 0.97, "structure shifts:\n"
             rf"slope $= {ag['mean']:+.2f} \pm {ag['stat_err']:.2f}$"
             rf"$\,\pm\,{ag['cell_rms']:.2f}$ (cells)",
             transform=axc.transAxes, fontsize=11, va="top", style="italic",
             color="0.25")
    axc.set_xscale("log")
    axc.set_xlabel(r"$p_{\rm CM}$ [GeV]")
    axc.set_xticks([5, 10, 20, 40])
    axc.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    axc.get_xaxis().set_minor_formatter(matplotlib.ticker.NullFormatter())
    plt.setp(axc.get_yticklabels(), visible=False)

    hw = f7["contrasts"].get("herwig")
    note = EIC_LABEL
    if hw and hw["vary_plab"]["aggregate"]:
        note += (f"\nHerwig 7: ${hw['vary_plab']['aggregate']['mean']:+.2f}$"
                 f" vs ${hw['vary_pcm']['aggregate']['mean']:+.2f}$")
    axc.text(0.96, 0.04, note, transform=axc.transAxes, fontsize=9.5,
             ha="right", color="0.45")

    tufte(axl)
    tufte(axc)
    _saveplot(fig, outdir, "fig7_decomposition.pdf")


# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="results/results.json")
    ap.add_argument("--outdir", default="plots")
    args = ap.parse_args()

    res = json.load(open(args.results))
    os.makedirs(args.outdir, exist_ok=True)
    fig1(res, args.outdir)
    fig2(res, args.outdir)
    fig3(res, args.outdir)
    fig4(res, args.outdir)
    fig5(res, args.outdir)
    fig6(res, args.outdir)
    fig7(res, args.outdir)
    print("Done.")


if __name__ == "__main__":
    main()
