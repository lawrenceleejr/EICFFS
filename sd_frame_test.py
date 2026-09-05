#!/usr/bin/env python3
"""
Which frame, and which variables, should soft-drop multiplicity be computed in?

Soft drop is normally written for pp with transverse-momentum fractions and an
angular distance in rapidity-azimuth (Larkoski, Marzani, Soyez, Thaler,
arXiv:1402.2657), which makes it invariant under boosts along the *beam*.  The
e+e- form uses energy fractions and opening angles and is used in the CM frame.
Iterated soft drop (arXiv:1704.06266) inherits the pp convention.

Neither is automatically safe in DIS, because the boost that relates the
laboratory to the colour rest frame runs along P + q, not along the beam.  This
compares three prescriptions on the beam-energy frame-independence test:

  e+e- variables in the laboratory        the naive choice
  pp variables about the P + q axis       standard prescription, right axis
  e+e- variables in the object rest frame boost the object to rest first

Usage
-----
    python sd_frame_test.py --outdir figures/ \
        --beams 5x41=data/e5p41_events_[12].parquet 10x100=data/events_[12].parquet \
                18x275=data/e18p275_events_[12].parquet
"""

import argparse
import glob
import os
import sys

import numpy as np
import awkward as ak
import fastjet as fj
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
plt.style.use(os.path.join(HERE, "utils", "tufte.mplstyle"))

from utils.dis_kinematics import (hcm_boost_matrix, breit_boost_matrix, apply_boost,
                                  rest_frame_boost_matrix)
from analyze_events import isd_multiplicity, SD_ZCUT, SD_BETA, SD_THETA_CUT, SD_R0

INK, MUTED, FAINT, GOOD, BAD = "#1f1f1f", "#8a8a8a", "#d9d9d9", "#2f6b4f", "#c44e52"
W_CELLS = [(10, 15), (15, 22), (22, 28)]
Q_CELLS = [(2.2, 3.3), (3.3, 5.0), (5.0, 7.5)]


def rotate_to_axis(V, axis):
    """Rotate four-vectors so that ``axis`` becomes +z."""
    z = axis / np.linalg.norm(axis)
    a = np.array([0.0, 0.0, 1.0]) if abs(z[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    x = np.cross(a, z)
    x /= np.linalg.norm(x)
    R = np.stack([x, np.cross(z, x), z])
    out = V.copy()
    out[:, :3] = V[:, :3] @ R.T
    return out


def isd_pp(p4, axis, zcut=SD_ZCUT, beta=SD_BETA, dr_cut=SD_THETA_CUT, R0=SD_R0):
    """
    Iterated soft drop in pp variables about ``axis``: Cambridge/Aachen in
    rapidity-azimuth, z from transverse momenta, distance ΔR.  Invariant under
    boosts along ``axis``.
    """
    n = len(p4)
    if n < 2:
        return 0
    V = rotate_to_axis(p4, axis)
    pj = [fj.PseudoJet(*map(float, v)) for v in V]
    cs = fj.ClusterSequence(pj, fj.JetDefinition(fj.cambridge_algorithm, 100.0))
    node = cs.exclusive_jets(1)[0]
    if len(node.constituents()) != n:
        raise RuntimeError("reclustering lost constituents")
    p1, p2 = fj.PseudoJet(), fj.PseudoJet()
    count = 0
    while node.has_parents(p1, p2):
        a, b = fj.PseudoJet(p1), fj.PseudoJet(p2)
        if a.pt() <= 0 or b.pt() <= 0:
            break
        dr = a.delta_R(b)
        z = min(a.pt(), b.pt()) / (a.pt() + b.pt())
        if dr > dr_cut and z > zcut * (dr / R0) ** beta:
            count += 1
        node = a if a.pt() >= b.pt() else b
    return count


def load_beam(pattern):
    files = sorted(glob.glob(pattern)) or [pattern]
    ev = ak.concatenate([ak.from_parquet(f) for f in files])
    W = ak.to_numpy(ev.W).astype(float)
    Q = np.sqrt(ak.to_numpy(ev.Q2).astype(float))
    P = np.stack([ak.to_numpy(ev.proton[c]) for c in ("px", "py", "pz", "e")], 1).astype(float)
    q = np.stack([ak.to_numpy(ev.q[c]) for c in ("px", "py", "pz", "e")], 1).astype(float)
    L, LB = hcm_boost_matrix(P, q), breit_boost_matrix(P, q)
    qB = apply_boost(LB, q)
    qh = qB[:, :3] / np.linalg.norm(qB[:, :3], axis=1)[:, None]
    n = ak.to_numpy(ak.num(ev.particles.px))
    evi = np.repeat(np.arange(len(W)), n)
    C = np.stack([ak.to_numpy(ak.flatten(ev.particles[c])) for c in ("px", "py", "pz", "e")],
                 1).astype(float)
    CB = np.einsum("nij,nj->ni", LB[evi], C)
    cur = np.einsum("ni,ni->n", CB[:, :3], qh[evi]) > 0
    Ccm = np.einsum("nij,nj->ni", L[evi], C)
    boost_axis = (P + q)[:, :3]
    beam_axis = np.tile(np.array([0.0, 0.0, 1.0]), (len(W), 1))

    e = evi[cur]
    o = np.argsort(e, kind="stable")
    e = e[o]
    Vlab, Vcm = C[cur][o], Ccm[cur][o]
    cnt = np.bincount(e, minlength=len(W))
    st = np.concatenate([[0], np.cumsum(cnt)[:-1]])
    H = np.zeros((len(W), 4))
    np.add.at(H, evi[cur], C[cur])
    res = {k: np.zeros(len(W)) for k in ("lab_ee", "beam_pp", "axis_pp", "rest_ee")}
    for i in np.where(cnt >= 2)[0]:
        vl, vc = Vlab[st[i]:st[i] + cnt[i]], Vcm[st[i]:st[i] + cnt[i]]
        res["lab_ee"][i] = isd_multiplicity(vl)
        res["beam_pp"][i] = isd_pp(vl, beam_axis[i])
        res["axis_pp"][i] = isd_pp(vl, boost_axis[i])
        Hi = vc.sum(axis=0)
        if Hi[3]**2 - (Hi[:3]**2).sum() > 1e-9:
            Lr = rest_frame_boost_matrix(Hi[None, :])[0]
            res["rest_ee"][i] = isd_multiplicity((Lr @ vc.T).T)
    return dict(W=W, Q=Q, plab=np.linalg.norm(H[:, :3], axis=1), nc=cnt, **res)


def slopes(beams, key):
    out = []
    for wlo, whi in W_CELLS:
        for qlo, qhi in Q_CELLS:
            xs, ys = [], []
            for d in beams.values():
                m = ((d["W"] >= wlo) & (d["W"] < whi) & (d["Q"] >= qlo) & (d["Q"] < qhi)
                     & (d["plab"] > 1) & (d["nc"] >= 2))
                if m.sum() < 200:
                    xs = []
                    break
                xs.append(np.median(d["plab"][m]))
                ys.append(d[key][m].mean())
            if xs and min(ys) > 0:
                out.append(np.polyfit(np.log(xs), np.log(ys), 1)[0])
    return np.array(out)


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
        print(f"{lab}: {len(beams[lab]['W']):,} events")

    variants = [
        (r"$e^+e^-$ variables ($E$, $\theta$)" "\n" "in the laboratory", "lab_ee", BAD),
        (r"$pp$ variables ($p_T$, $\Delta R$)" "\n" r"about the beam axis, in the lab", "beam_pp", INK),
        (r"$pp$ variables ($p_T$, $\Delta R$)" "\n" r"about the $P+q$ axis, in the lab", "axis_pp", GOOD),
        (r"$e^+e^-$ variables ($E$, $\theta$)" "\n" "in the object's rest frame", "rest_ee", GOOD),
    ]
    rows = []
    print(f"\n{'variant':<46}{'median':>9}{'|max|':>8}{'mean':>7}")
    for name, key, col in variants:
        sl = slopes(beams, key)
        mean = np.mean([d[key][d["nc"] >= 2].mean() for d in beams.values()])
        rows.append((name, sl, mean, col))
        print(f"{name.replace(chr(10), ' '):<46}{np.median(sl):>+9.3f}"
              f"{np.max(np.abs(sl)):>8.3f}{mean:>7.2f}")

    fig, ax = plt.subplots(figsize=(5.0, 3.4))
    ax.axvline(0.0, color=FAINT, lw=0.9, zorder=0)
    ys = np.arange(len(rows))[::-1]
    for y, (name, sl, mean, col) in zip(ys, rows):
        ax.plot([sl.min(), sl.max()], [y, y], color=col, lw=1.0, alpha=0.35, solid_capstyle="butt")
        ax.plot([np.median(sl)], [y], marker="o", ms=5, color=col, mec="white", mew=0.6, zorder=3)
        ax.annotate(f"{np.median(sl):+.3f}", (np.median(sl), y), xytext=(0, 8),
                    textcoords="offset points", ha="center", fontsize=7.5, color=col)
    ax.set_yticks(ys)
    ax.set_yticklabels([r[0] for r in rows], fontsize=7.5)
    ax.tick_params(axis="y", length=0)
    for side in ("left", "top", "right"):
        ax.spines[side].set_visible(False)
    ax.set_ylim(-0.6, len(rows) - 0.4)
    ax.set_xlim(-0.62, 0.12)
    ax.set_xticks([-0.6, -0.4, -0.2, 0.0])
    ax.spines["bottom"].set_bounds(-0.6, 0.0)
    ax.set_xlabel(r"lab-frame dependence,  $\mathrm{d}\ln\langle n_{\rm SD}\rangle\,/\,\mathrm{d}\ln|\vec p|_{\rm lab}$")
    t = ax.text(0.0, -0.34,
                "Soft drop is normally written with transverse-momentum fractions and a "
                "rapidity-azimuth distance, which is invariant under boosts along the beam.  In DIS "
                "the boost that matters runs along $P+q$: measure the same standard variables about "
                "that axis and the observable is frame independent without boosting anything.  "
                "Going to the object's rest frame works equally well but is not the standard "
                "convention and destroys the angle for two-body objects.",
                transform=ax.transAxes, fontsize=7.5, color=MUTED, va="top", ha="left", wrap=True)
    t._is_caption = True
    os.makedirs(args.outdir, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(args.outdir, f"sd_frame_choice.{ext}"), bbox_inches="tight",
                    dpi=220 if ext == "png" else None)
    print(f"  -> {args.outdir}/sd_frame_choice.pdf")


if __name__ == "__main__":
    main()
