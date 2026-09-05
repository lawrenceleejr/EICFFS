#!/usr/bin/env python3
"""
Infrared and collinear safety of n90 against iterated soft-drop multiplicity.

Takes the Breit current hemisphere of real events, applies two deformations
that any infrared- and collinear-safe observable must be insensitive to in the
limit, and measures how each observable responds:

  collinear   every constituent is split into two equal halves separated by an
              angle delta;  the observable must return to its unsplit value as
              delta -> 0
  soft        three particles each carrying a fraction eps of the object's
              momentum are added at random angles;  the observable must return
              to its unmodified value as eps -> 0

Usage
-----
    python irc_safety_test.py data/events_1.parquet --outdir figures/
"""

import argparse
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
from analyze_events import compute_n_x, isd_multiplicity

INK, MUTED, FAINT, ACCENT = "#1f1f1f", "#8a8a8a", "#d9d9d9", "#c44e52"
DELTAS = [0.3, 0.1, 0.03, 0.01, 0.003, 0.001]
EPSILONS = [0.05, 0.01, 0.003, 1e-3, 1e-4]


def hemispheres(path, n_events, min_const):
    """Constituent four-vectors of the Breit current hemisphere, per event."""
    ev = ak.from_parquet(path)[:n_events]
    W = ak.to_numpy(ev.W).astype(float)
    P = np.stack([ak.to_numpy(ev.proton[c]) for c in ("px", "py", "pz", "e")], 1).astype(float)
    q = np.stack([ak.to_numpy(ev.q[c]) for c in ("px", "py", "pz", "e")], 1).astype(float)
    LB = breit_boost_matrix(P, q)
    qB = apply_boost(LB, q)
    qhat = qB[:, :3] / np.linalg.norm(qB[:, :3], axis=1)[:, None]
    n = ak.to_numpy(ak.num(ev.particles.px))
    evi = np.repeat(np.arange(len(W)), n)
    C = np.stack([ak.to_numpy(ak.flatten(ev.particles[c])) for c in ("px", "py", "pz", "e")], 1).astype(float)
    CB = np.einsum("nij,nj->ni", LB[evi], C)
    cur = np.einsum("ni,ni->n", CB[:, :3], qhat[evi]) > 0
    e, v = evi[cur], C[cur]
    order = np.argsort(e, kind="stable")
    e, v = e[order], v[order]
    cnt = np.bincount(e, minlength=len(W))
    starts = np.concatenate([[0], np.cumsum(cnt)[:-1]])
    return [v[starts[i]:starts[i] + cnt[i]] for i in np.where(cnt >= min_const)[0]]


def split_collinear(p4, delta, rng):
    out = []
    for v in p4:
        p = v[:3]
        nrm = np.linalg.norm(p)
        if nrm <= 0:
            out.append(v)
            continue
        u = p / nrm
        a = np.cross(u, rng.normal(size=3))
        a /= max(np.linalg.norm(a), 1e-12)
        for s in (+1, -1):
            d = u * np.cos(delta / 2) + s * a * np.sin(delta / 2)
            out.append(np.concatenate([0.5 * nrm * d, [0.5 * v[3]]]))
    return np.array(out)


def add_soft(p4, eps, rng, m=3):
    tot = np.linalg.norm(p4[:, :3], axis=1).sum()
    extra = []
    for _ in range(m):
        d = rng.normal(size=3)
        d /= np.linalg.norm(d)
        e = eps * tot
        extra.append(np.concatenate([e * d, [e]]))
    return np.vstack([p4, np.array(extra)])


def measure(objs, transform):
    a = [compute_n_x(np.linalg.norm(transform(o)[:, :3], axis=1)) for o in objs]
    b = [isd_multiplicity(transform(o)) for o in objs]
    return float(np.nanmean(a)), float(np.mean(b))


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("input", help="event Parquet file")
    p.add_argument("--outdir", default="figures/")
    p.add_argument("--n-objects", type=int, default=1500)
    p.add_argument("--n-events", type=int, default=60000)
    p.add_argument("--min-const", type=int, default=3)
    args = p.parse_args()

    rng = np.random.default_rng(7)
    objs = hemispheres(args.input, args.n_events, args.min_const)[:args.n_objects]
    print(f"{len(objs)} current hemispheres with >= {args.min_const} constituents")
    base = measure(objs, lambda v: v)
    print(f"unmodified:  <n90> = {base[0]:.3f}   <n_SD> = {base[1]:.3f}")

    coll = [measure(objs, lambda v, d=d: split_collinear(v, d, rng)) for d in DELTAS]
    soft = [measure(objs, lambda v, e=e: add_soft(v, e, rng)) for e in EPSILONS]
    for d, m in zip(DELTAS, coll):
        print(f"  collinear delta={d:<7g}  n90 {m[0]:.3f} ({m[0]-base[0]:+.3f})   "
              f"n_SD {m[1]:.3f} ({m[1]-base[1]:+.3f})")
    for e, m in zip(EPSILONS, soft):
        print(f"  soft      eps={e:<9g}  n90 {m[0]:.3f} ({m[0]-base[0]:+.3f})   "
              f"n_SD {m[1]:.3f} ({m[1]-base[1]:+.3f})")

    fig, ax = plt.subplots(figsize=(4.8, 3.3))
    ax.axhline(0.0, color=FAINT, lw=0.9, zorder=0)
    d = np.array(DELTAS)
    e = np.array(EPSILONS)
    rel = lambda vals, i, b: 100 * (np.array([v[i] for v in vals]) - b) / b
    ax.plot(d, rel(coll, 0, base[0]), color=ACCENT, lw=1.3, marker="o", ms=4, mec="white", mew=0.5)
    ax.plot(d, rel(coll, 1, base[1]), color=INK, lw=1.3, marker="o", ms=4, mec="white", mew=0.5)
    ax.plot(e, rel(soft, 0, base[0]), color=ACCENT, lw=1.1, ls=(0, (3, 2)), marker="s", ms=3.5,
            mec="white", mew=0.5)
    ax.plot(e, rel(soft, 1, base[1]), color=INK, lw=1.1, ls=(0, (3, 2)), marker="s", ms=3.5,
            mec="white", mew=0.5)
    ax.annotate(r"$n_{90}$, collinear split", (d[-1], rel(coll, 0, base[0])[-1]),
                xytext=(7, 0), textcoords="offset points", fontsize=8, color=ACCENT, va="center")
    ax.annotate(r"$n_{\rm SD}$, collinear split: zero below $\theta_{\rm cut}$",
                (d[0], rel(coll, 1, base[1])[0]), xytext=(8, 0), textcoords="offset points",
                fontsize=8, color=INK, va="center")
    ax.annotate("soft additions (dashed):\nboth vanish in the limit", (e[-1], 6),
                xytext=(7, 0), textcoords="offset points", fontsize=7.5, color=MUTED, va="center")
    from analyze_events import SD_THETA_CUT
    ax.axvline(SD_THETA_CUT, color=FAINT, lw=0.9, ls=(0, (3, 3)), zorder=0)
    ax.annotate(r"$\theta_{\rm cut}$", (SD_THETA_CUT, 104), ha="center", va="bottom",
                fontsize=7.5, color=MUTED)
    ax.set_xscale("log")
    ax.invert_xaxis()
    ax.set_xlabel(r"splitting angle $\delta$ [rad],  or soft fraction $\epsilon$   $\longrightarrow$ limit")
    ax.set_ylabel("change in the observable  [%]")
    ax.set_ylim(-8, 118)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.text(0.0, -0.24,
            "Every constituent split into two equal halves at angle $\\delta$ (solid), or three "
            "particles of momentum fraction $\\epsilon$ added (dashed).  A collinear-safe observable "
            "returns to zero as $\\delta\\to0$.  $n_{90}$ does not: a democratic split doubles it, "
            "however small the angle.  Both observables are infrared safe.",
            transform=ax.transAxes, fontsize=7.5, color=MUTED, va="top", ha="left", wrap=True)
    os.makedirs(args.outdir, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(args.outdir, f"irc_safety.{ext}"), bbox_inches="tight",
                    dpi=220 if ext == "png" else None)
    print(f"  -> {args.outdir}/irc_safety.pdf")


if __name__ == "__main__":
    main()
