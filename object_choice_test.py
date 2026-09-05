#!/usr/bin/env python3
"""
Which current-region definition to measure, and what breaks its frame independence.

A cone jet keeps a boost-dependent share of the shower (see slope_vs_radius.pdf).
Taking the whole current region instead removes that, but the current region has
its own choices: which frame defines it, and what particle-level thresholds are
applied.  This runs the beam-energy frame-independence test over those choices.

For each definition, n90 is computed from colour-frame momenta and the exponent
d ln<n90> / d ln|p|_lab is fitted across the 5x41, 10x100 and 18x275 GeV
configurations in cells of fixed (W, Q).  Zero means frame independent.

Usage
-----
    python object_choice_test.py --outdir figures/ \
        --beams 5x41=data/e5p41_events_*.parquet 10x100=data/events_*.parquet \
                18x275=data/e18p275_events_*.parquet
"""

import argparse
import glob
import os
import sys

import numpy as np
import awkward as ak
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
plt.style.use(os.path.join(HERE, "utils", "tufte.mplstyle"))

from utils.dis_kinematics import hcm_boost_matrix, breit_boost_matrix, apply_boost
from analyze_events import n_x_segments

INK, MUTED, FAINT = "#1f1f1f", "#8a8a8a", "#d9d9d9"
GOOD, BAD = "#2f6b4f", "#c44e52"
W_CELLS = [(10, 15), (15, 22), (22, 28)]
Q_CELLS = [(2.2, 3.3), (3.3, 5.0), (5.0, 7.5)]
MIN_ENTRIES = 200

# (label, selector, colour key).  Selectors receive the per-beam arrays.
DEFINITIONS = [
    ("Breit current hemisphere", lambda s: s["cur_breit"], "base"),
    ("  with lab $|\\vec p| > 0.5$ GeV", lambda s: s["cur_breit"] & (s["p"] > 0.5), "bad"),
    ("  with lab $p_T > 0.15$ GeV", lambda s: s["cur_breit"] & (s["pt"] > 0.15), "base"),
    ("$\\gamma^*p$ current region, $y^* > 0$", lambda s: s["ystar"] > 0, "good"),
    ("  with lab $|\\vec p| > 0.5$ GeV", lambda s: (s["ystar"] > 0) & (s["p"] > 0.5), "bad"),
    ("  with lab $p_T > 0.15$ GeV", lambda s: (s["ystar"] > 0) & (s["pt"] > 0.15), "good"),
]


def load_beam(pattern):
    files = sorted(glob.glob(pattern)) or [pattern]
    ev = ak.concatenate([ak.from_parquet(f) for f in files])
    W = ak.to_numpy(ev.W).astype(float)
    Q = np.sqrt(ak.to_numpy(ev.Q2).astype(float))
    P = np.stack([ak.to_numpy(ev.proton[c]) for c in ("px", "py", "pz", "e")], 1).astype(float)
    q = np.stack([ak.to_numpy(ev.q[c]) for c in ("px", "py", "pz", "e")], 1).astype(float)
    L, LB = hcm_boost_matrix(P, q), breit_boost_matrix(P, q)
    qB = apply_boost(LB, q)
    qhb = qB[:, :3] / np.linalg.norm(qB[:, :3], axis=1)[:, None]
    qc = apply_boost(L, q)
    qhc = qc[:, :3] / np.linalg.norm(qc[:, :3], axis=1)[:, None]
    n = ak.to_numpy(ak.num(ev.particles.px))
    evi = np.repeat(np.arange(len(W)), n)
    C = np.stack([ak.to_numpy(ak.flatten(ev.particles[c])) for c in ("px", "py", "pz", "e")],
                 1).astype(float)
    CB = np.einsum("nij,nj->ni", LB[evi], C)
    Ccm = np.einsum("nij,nj->ni", L[evi], C)
    pl = np.einsum("ni,ni->n", Ccm[:, :3], qhc[evi])
    ystar = 0.5 * np.log(np.clip((Ccm[:, 3] + pl) / np.clip(Ccm[:, 3] - pl, 1e-9, None), 1e-12, None))
    return dict(W=W, Q=Q, evi=evi, C=C, Ccm=Ccm, nev=len(W),
                cur_breit=np.einsum("ni,ni->n", CB[:, :3], qhb[evi]) > 0,
                ystar=ystar, p=np.linalg.norm(C[:, :3], axis=1),
                pt=np.hypot(C[:, 0], C[:, 1]))


def evaluate(beams, selector):
    """Per-cell exponents across beams, plus the mean constituent multiplicity."""
    per = {}
    for lab, s in beams.items():
        sel = selector(s)
        H = np.zeros((s["nev"], 4))
        np.add.at(H, s["evi"][sel], s["C"][sel])
        per[lab] = dict(
            plab=np.linalg.norm(H[:, :3], axis=1),
            n90=n_x_segments(np.linalg.norm(s["Ccm"][sel][:, :3], axis=1), s["evi"][sel], s["nev"]),
            nc=np.bincount(s["evi"][sel], minlength=s["nev"]))
    slopes = []
    for wlo, whi in W_CELLS:
        for qlo, qhi in Q_CELLS:
            xs, ys = [], []
            for lab, s in beams.items():
                d = per[lab]
                m = ((s["W"] >= wlo) & (s["W"] < whi) & (s["Q"] >= qlo) & (s["Q"] < qhi)
                     & (d["plab"] > 1) & np.isfinite(d["n90"]))
                if m.sum() < MIN_ENTRIES:
                    xs = []
                    break
                xs.append(np.median(d["plab"][m]))
                ys.append(d["n90"][m].mean())
            if xs and min(ys) > 0:
                slopes.append(np.polyfit(np.log(xs), np.log(ys), 1)[0])
    nc = np.mean([per[l]["nc"][per[l]["plab"] > 1].mean() for l in per])
    empty = np.mean([np.mean(per[l]["nc"][per[l]["plab"] > 1] <= 1) for l in per])
    return np.array(slopes), nc, empty


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--beams", nargs="+", required=True, metavar="LABEL=GLOB")
    p.add_argument("--outdir", default="figures/")
    args = p.parse_args()

    beams = {}
    for item in args.beams:
        lab, pat = item.split("=", 1)
        beams[lab] = load_beam(pat)
        print(f"{lab}: {beams[lab]['nev']:,} events")

    rows = []
    print(f"\n{'definition':<40}{'median':>9}{'|max|':>8}{'<N>':>7}{'N<=1':>7}{'cells':>7}")
    for name, sel, kind in DEFINITIONS:
        sl, nc, empty = evaluate(beams, sel)
        rows.append((name, np.median(sl), np.max(np.abs(sl)), nc, empty, kind, sl))
        plain = name.replace("$", "").replace("\\gamma^*p", "gamma*p").replace("\\vec p", "|p|")
        print(f"{plain:<40}{np.median(sl):>+9.3f}{np.max(np.abs(sl)):>8.3f}"
              f"{nc:>7.2f}{100*empty:>6.0f}%{len(sl):>7}")

    fig, ax = plt.subplots(figsize=(5.2, 3.4))
    ax.axvline(0.0, color=FAINT, lw=0.9, zorder=0)
    ys = np.arange(len(rows))[::-1]
    cols = {"base": INK, "good": GOOD, "bad": BAD}
    for y, (name, med, mx, nc, empty, kind, sl) in zip(ys, rows):
        c = cols[kind]
        ax.plot([sl.min(), sl.max()], [y, y], color=c, lw=1.0, alpha=0.35, solid_capstyle="butt")
        ax.plot([med], [y], marker="o", ms=5, color=c, mec="white", mew=0.6, zorder=3)
        ax.annotate(f"{med:+.3f}", (med, y), xytext=(0, 8), textcoords="offset points",
                    ha="center", fontsize=7.5, color=c)
        ax.annotate(f"{nc:.1f} particles", (0.30, y), xytext=(0, 0), textcoords="offset points",
                    ha="left", va="center", fontsize=7, color=MUTED)
    ax.set_yticks(ys)
    ax.set_yticklabels([r[0] for r in rows], fontsize=8)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_ylim(-0.6, len(rows) - 0.4)
    ax.set_xlim(-0.10, 0.42)
    ax.set_xticks([-0.1, 0.0, 0.1, 0.2])
    ax.set_xlabel(r"lab-frame dependence,  $\mathrm{d}\ln\langle n_{90}\rangle\,/\,\mathrm{d}\ln|\vec p|_{\rm lab}$")
    ax.spines["bottom"].set_bounds(-0.1, 0.2)
    t = ax.text(0.0, -0.30,
                "Taking the whole current region removes the cone entirely, and a detector's "
                "angular acceptance does not spoil it: 99.7% of the current hemisphere's momentum "
                "lies inside $|\\eta| < 3.5$.  What does spoil it is a threshold on total momentum, "
                "which is not invariant under a boost along the axis.  A transverse-momentum "
                "threshold is, and costs nothing.  The $\\gamma^*p$ region also holds four times "
                "the particles of the Breit hemisphere.",
                transform=ax.transAxes, fontsize=7.5, color=MUTED, va="top", ha="left", wrap=True)
    t._is_caption = True
    os.makedirs(args.outdir, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(args.outdir, f"object_choice.{ext}"), bbox_inches="tight",
                    dpi=220 if ext == "png" else None)
    print(f"  -> {args.outdir}/object_choice.pdf")


if __name__ == "__main__":
    main()
