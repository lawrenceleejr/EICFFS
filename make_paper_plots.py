#!/usr/bin/env python3
"""
PRL figure stage for the EIC FFS study (ANALYSIS_DESIGN.md Sec. 3).

Reads results/results.json (from make_results.py) and produces:

  fig1_concept.pdf       boost factor + string-anchor mapping (reach)
  fig2_ffs_main.pdf      <n90_lab> vs W at fixed (|p|_lab, Q^2): the effect
  fig3_universality.pdf  lab-frame splay vs CM-frame collapse
  fig4_projection.pdf    H0-rejection significance vs integrated luminosity

Usage:  python make_paper_plots.py [--results results/results.json] [--outdir plots]
"""

import argparse
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mplhep as hep

plt.style.use(hep.style.CMS)
plt.rcParams.update({
    "figure.dpi": 150,
    "font.size": 15,
    "legend.fontsize": 11,
    "axes.labelsize": 17,
})

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
        ax1.plot(W, [r["gamma_mean"] for r in rows],
                 marker=CONFIG_MARKERS[config], color=CONFIG_COLORS[config],
                 label=CONFIG_LABELS[config] + " GeV")
    ax1.set_xlabel(r"$W$ [GeV]")
    ax1.set_ylabel(r"$\langle\gamma\rangle$ of $\gamma^*p$ frame in lab")
    ax1.set_yscale("log")
    ax1.legend(title=r"$E_e\times E_p$")
    ax1.text(0.05, 0.06, "fragmentation frame is\nstrongly boosted in the lab",
             transform=ax1.transAxes, fontsize=12, style="italic")

    for config, rows in bmap.items():
        W = np.array([r["W_center"] for r in rows])
        ax2.plot(W, [r["pcm_mean"] for r in rows], marker=CONFIG_MARKERS[config],
                 color=CONFIG_COLORS[config], ls="-", ms=5)
        ax2.plot(W, [r["plab_mean"] for r in rows], marker=CONFIG_MARKERS[config],
                 color=CONFIG_COLORS[config], ls=":", alpha=0.55, mfc="none",
                 ms=5)
    Wd = np.linspace(8, 135, 50)
    ax2.plot(Wd, Wd / 2, "k--", lw=1.2)
    ax2.set_xlabel(r"$W$ [GeV]")
    ax2.set_ylabel(r"jet momentum [GeV]")
    from matplotlib.lines import Line2D
    ax2.legend(handles=[
        Line2D([], [], color="gray", ls="-", marker="o", ms=5,
               label=r"$\langle p_{\rm CM}\rangle$ (filled)"),
        Line2D([], [], color="gray", ls=":", marker="o", ms=5, mfc="none",
               label=r"$\langle |p|_{\rm lab}\rangle$ (open)"),
        Line2D([], [], color="k", ls="--", label=r"$W/2$ string anchor"),
    ], fontsize=11, loc="upper left")
    ax2.text(0.45, 0.06,
             "CM momentum locked to $W/2$;\nlab momentum is a frame artifact",
             transform=ax2.transAxes, fontsize=12, style="italic")

    fig.suptitle(EIC_LABEL, x=0.5, y=1.0, fontsize=14)
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
            ax.fill_between(Ws, lo, hi, color=pcolors[pkey], alpha=0.18, lw=0)
            axr.fill_between(Ws, lo / anc, hi / anc, color=pcolors[pkey],
                             alpha=0.18, lw=0)
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
            ax.errorbar(W, m, yerr=e, marker=CONFIG_MARKERS[config], ms=5,
                        color=pcolors[pkey], ls="-", lw=1.2, capsize=0)
            axr.errorbar(W, np.array(m) / anc, yerr=np.array(e) / anc,
                         marker=CONFIG_MARKERS[config], ms=5,
                         color=pcolors[pkey], ls="-", lw=1.2)
            shown_p.add(pkey)
        # H0: frame-independent null
        ax.axhline(anc, color=pcolors[pkey], ls="--", lw=1.0, alpha=0.8)

    axr.axhline(1.0, color="k", ls="--", lw=1.0)
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
    ax.legend(handles=handles, ncol=2, loc="upper left", fontsize=10.5)
    lumi_note = ""
    li = res.get("fig4", {}).get("luminosity_inputs", {}).get("per_config", {})
    if li:
        mc_pb = min(v["mc_equiv_fb"] for v in li.values()) * 1e3
        ann = [v["annual_fb"] for v in li.values() if v.get("annual_fb")]
        lumi_note = (f"\nMC stats $\\approx$ first {mc_pb:.0f} pb$^{{-1}}$; "
                     f"1 EIC yr = {min(ann):.0f}--{max(ann):.0f} fb$^{{-1}}$")
    ax.text(0.97, 0.05, EIC_LABEL + "\n" + SEL_LABEL + lumi_note,
            transform=ax.transAxes, fontsize=10, ha="right")
    ax.set_title("Frame-dependent fragmentation shift at the EIC", fontsize=14)
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
            ax.errorbar(x, m, yerr=e, marker=CONFIG_MARKERS[config], ms=6,
                        ls="-", lw=1.3, color=cmap(norm(wcenter(key))))
        ax.set_xlabel(xlab)
        ax.set_xscale("log")
        ax.set_xticks([3, 5, 10, 20, 40])
        ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

    ax1.set_ylabel(r"$\langle n_{90}\rangle$")
    ax1.set_title("binned in the lab frame", fontsize=14)
    ax2.set_title(r"binned in the $\gamma^*p$ (color) frame", fontsize=14)
    ax1.text(0.04, 0.90, "universality broken:\n"
             rf"slice spread {100*f3['spread_lab']:.0f}% RMS",
             transform=ax1.transAxes, fontsize=13)
    ax2.text(0.04, 0.90, "universality restored:\n"
             rf"slice spread {100*f3['spread_cm']:.1f}% RMS",
             transform=ax2.transAxes, fontsize=13)

    from matplotlib.lines import Line2D
    handles = [Line2D([], [], color="gray", marker=CONFIG_MARKERS[c], ls="",
                      label=CONFIG_LABELS[c] + " GeV") for c in CONFIG_MARKERS]
    ax2.legend(handles=handles, loc="lower right", fontsize=11)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    cbar = fig.colorbar(sm, ax=[ax1, ax2], pad=0.012, aspect=28)
    cbar.set_label(r"$W$ slice center [GeV]")
    ax1.text(0.04, 0.05, EIC_LABEL, transform=ax1.transAxes, fontsize=11)
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
        ax.plot(L, v["significance_syst"], color=CONFIG_COLORS[config], lw=1.4,
                alpha=0.8, label=CONFIG_LABELS[config] + " GeV")
        ann = lumi_info.get(config, {}).get("annual_fb")
        if ann:
            for Lmark, mk, ms in ((ann / 365.25, "D", 7), (ann, "*", 16)):
                sig_m = np.interp(np.log10(Lmark), np.log10(L),
                                  v["significance_syst"])
                ax.plot(Lmark, sig_m, mk, color=CONFIG_COLORS[config],
                        ms=ms, mec="k", mew=0.6, zorder=5)
    syst = res["binnings"]["syst_floor"]
    cons = res["binnings"].get("syst_floor_cons", 0.05)
    ax.plot(L, proj["combined"]["significance_syst"], "k-", lw=2.6,
            label=f"combined ({100*syst:g}% syst floor)")
    ax.plot(L, proj["combined"]["significance_syst_cons"], "k-.", lw=1.8,
            label=f"combined ({100*cons:g}% syst floor)")
    ax.plot(L, proj["combined"]["significance_stat"], "k--", lw=1.6,
            label="combined (stat only)")

    ax.axhline(5, color="crimson", ls=":", lw=1.6)
    ax.text(L[-1] * 0.6, 5.5, r"$5\sigma$", color="crimson", fontsize=13,
            ha="right")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"integrated luminosity per configuration [fb$^{-1}$]")
    ax.set_ylabel(r"expected significance of $H_0$ rejection [$\sigma$]")
    from matplotlib.lines import Line2D
    handles, labels = ax.get_legend_handles_labels()
    handles += [Line2D([], [], marker="D", color="gray", ls="", ms=7,
                       mec="k", label="1 day at design lumi."),
                Line2D([], [], marker="*", color="gray", ls="", ms=14,
                       mec="k", label="1 year at design lumi.")]
    ax.legend(handles=handles, loc="upper left", fontsize=11,
              title="reco level: track jets, smeared,\nelectron-method $W$",
              title_fontsize=10)
    ax.set_title("EIC discovery reach for frame-dependent fragmentation",
                 fontsize=14)
    syst = res["binnings"]["syst_floor"]
    ax.text(0.97, 0.04, EIC_LABEL +
            f"\nsyst floor {100*syst:.1f}% per bin"
            "\ndesign lumi: EIC YR Tab. 10.1", transform=ax.transAxes,
            fontsize=10, ha="right")
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
                            ls=ls, marker=mk, ms=4, lw=1.4, color=color,
                            label=(label if (label and not shown) else None))
                shown = True
        ax.set_xscale("log")
        ax.text(0.05, 0.88, rf"$|p|_{{\rm lab}} \in [{pkey.replace('-', ',')})$ GeV",
                transform=ax.transAxes, fontsize=12)
    for ax in axes[1]:
        ax.set_xlabel(r"$W$ [GeV]")
        ax.set_xticks([10, 20, 50, 100])
        ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    for ax in axes[:, 0]:
        ax.set_ylabel(r"$\langle n_{90}\rangle$ (lab frame)")
    axes[0, 0].legend(fontsize=11, loc="lower right")
    axes[0, 1].text(0.95, 0.05, EIC_LABEL + "\n" + SEL_LABEL,
                    transform=axes[0, 1].transAxes, fontsize=9, ha="right")
    fig.suptitle("FFS trend: string vs cluster hadronization", fontsize=15)
    _saveplot(fig, outdir, "fig5_generators.pdf")


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
    print("Done.")


if __name__ == "__main__":
    main()
