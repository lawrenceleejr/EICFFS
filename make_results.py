#!/usr/bin/env python3
"""
FFS statistics stage: per-jet tables -> results.json

Implements ANALYSIS_DESIGN.md Secs. 6 and 8:
  * Fig.2 profiles: <n90_lab> vs W in fixed (|p|_lab, Q^2) bins, per
    config / variation / level, with the frame-independent null H0
    anchored at the lowest populated W bin;
  * Fig.3 data: lab-frame splay vs CM-frame collapse, with the
    universality chi^2 metric in both binnings;
  * Q^2-binned cross-check of the W trend (DGLAP control);
  * Fig.4: projected significance of H0 rejection vs integrated
    luminosity at reco level (smeared, track jets, electron-method W),
    per config and combined, with and without a systematics floor.

Usage:  python make_results.py [--datadir data] [--output results/results.json]
"""

import argparse
import glob
import json
import os

import numpy as np
import awkward as ak

# ---------------------------------------------------------------------------
# Binnings (ANALYSIS_DESIGN.md Sec. 6)
# ---------------------------------------------------------------------------

Q2_BINS = [(25., 50.), (50., 100.), (100., 250.), (250., 1000.)]
Q2_MAIN = (25., 100.)

PLAB_BINS = [(4., 6.), (6., 9.), (9., 13.), (13., 20.), (20., 30.)]

W_EDGES = np.array([8, 12, 16, 20, 25, 30, 36, 43, 51, 60, 75, 95, 135.])

# W slices for the universality (splay/collapse) figure
W_SLICES = {
    "5x41":   [(10., 18.), (18., 26.)],
    "10x100": [(18., 30.), (30., 42.), (42., 56.)],
    "18x275": [(40., 60.), (60., 85.), (85., 130.)],
}

PCM_EDGES = np.array([2, 4, 6, 8, 10, 12, 15, 18, 22, 27, 33, 40, 50, 65.])
PLAB_FINE = np.array([3, 4.5, 6, 8, 10, 13, 17, 22, 28, 36, 46, 60.])

MIN_JETS_PER_BIN = 100
BULK_FRACTION = 0.10        # universality metric: keep bins with
                            # n >= BULK_FRACTION * (slice's most populated bin),
                            # excluding sparse radiative-tail populations
SYST_FLOOR = 0.015          # fractional uncertainty floor on <n90> per bin
SYST_FLOOR_CONS = 0.05      # very conservative floor (robustness curve)
LUMI_GRID = np.logspace(-4, 2, 49)      # fb^-1

# ── EIC design luminosities ────────────────────────────────────────────────
# Peak luminosity per e-p beam configuration, EIC Yellow Report
# (arXiv:2103.05419) Table 10.1 / EIC CDR, high-divergence configuration.
PEAK_LUMI_1E33 = {"5x41": 0.44, "10x100": 4.48, "18x275": 1.54}  # 1e33/cm^2/s
ANNUAL_SECONDS = 1.0e7      # CDR operations year (~60% duty factor)
# 1e33 cm^-2 s^-1 x 1e7 s = 1e40 cm^-2 = 10 fb^-1
ANNUAL_FB = {c: v * 10.0 for c, v in PEAK_LUMI_1E33.items()}


# ---------------------------------------------------------------------------
# Sample loading
# ---------------------------------------------------------------------------

def load_samples(datadir):
    samples = {}
    for path in sorted(glob.glob(os.path.join(datadir, "jets_*.parquet"))):
        meta = json.load(open(path.replace(".parquet", ".json")))
        key = (meta["config"], meta["variation"] + ("_mpi" if meta.get("mpi") else ""),
               meta["level"])
        arr = ak.from_parquet(path)
        samples[key] = ({f: np.asarray(arr[f], dtype=float) for f in arr.fields},
                        meta)
        print(f"  loaded {key}: {len(arr):,} jets")
    return samples


def profile(values, mask):
    """Return (mean, sem, rms, n) of values[mask]."""
    v = values[mask]
    n = len(v)
    if n == 0:
        return np.nan, np.nan, np.nan, 0
    rms = float(np.std(v))
    return float(np.mean(v)), rms / np.sqrt(n), rms, int(n)


# ---------------------------------------------------------------------------
# Fig. 2: <n90_lab> vs W at fixed (|p|_lab, Q^2)
# ---------------------------------------------------------------------------

def fig2_profiles(samples, q2bin=Q2_MAIN, obs="n90lab"):
    out = {}
    for (config, variation, level), (d, meta) in samples.items():
        q2m = (d["Q2"] >= q2bin[0]) & (d["Q2"] < q2bin[1])
        entry = {}
        for plo, phi in PLAB_BINS:
            pm = q2m & (d["plab"] >= plo) & (d["plab"] < phi)
            rows = []
            for wlo, whi in zip(W_EDGES[:-1], W_EDGES[1:]):
                m = pm & (d["W"] >= wlo) & (d["W"] < whi)
                mean, sem, rms, n = profile(d[obs], m)
                if n >= MIN_JETS_PER_BIN:
                    rows.append({"W_lo": wlo, "W_hi": whi,
                                 "W_center": 0.5 * (wlo + whi),
                                 "mean": mean, "sem": sem, "rms": rms, "n": n})
            entry[f"{plo:g}-{phi:g}"] = rows
        out[f"{config}|{variation}|{level}"] = entry
    return out


def h0_anchor(fig2, level="truth", variation="baseline"):
    """
    Frame-independent null: for each |p|_lab bin, the <n90_lab> of the
    lowest populated W bin across all configs at the given level/variation.
    """
    anchors = {}
    for pkey in [f"{lo:g}-{hi:g}" for lo, hi in PLAB_BINS]:
        best = None
        for skey, entry in fig2.items():
            config, var, lev = skey.split("|")
            if var != variation or lev != level:
                continue
            for row in entry.get(pkey, []):
                if best is None or row["W_lo"] < best["W_lo"]:
                    best = row
        if best is not None:
            anchors[pkey] = {"value": best["mean"], "sem": best["sem"],
                             "W_lo": best["W_lo"], "W_hi": best["W_hi"]}
    return anchors


def effect_sizes(fig2, anchors, level="truth", variation="baseline"):
    """Relative shift of <n90_lab> at the highest vs lowest populated W."""
    out = {}
    for pkey, anc in anchors.items():
        hi = None
        for skey, entry in fig2.items():
            config, var, lev = skey.split("|")
            if var != variation or lev != level:
                continue
            for row in entry.get(pkey, []):
                if hi is None or row["W_lo"] > hi["W_lo"]:
                    hi = row
        if hi is None or hi["W_lo"] <= anc["W_lo"]:
            continue
        shift = hi["mean"] / anc["value"] - 1.0
        err = shift_err(hi, anc)
        out[pkey] = {"shift": shift, "err": err,
                     "W_low": [anc["W_lo"], anc["W_hi"]],
                     "W_high": [hi["W_lo"], hi["W_hi"]],
                     "n90_low": anc["value"], "n90_high": hi["mean"]}
    return out


def shift_err(hi, anc):
    r = hi["mean"] / anc["value"]
    return float(r * np.sqrt((hi["sem"] / hi["mean"])**2 +
                             (anc["sem"] / anc["value"])**2))


# ---------------------------------------------------------------------------
# Q^2-binned cross-check (DGLAP control)
# ---------------------------------------------------------------------------

def q2_crosscheck(samples):
    out = {}
    for q2bin in Q2_BINS:
        f2 = fig2_profiles(samples, q2bin=q2bin)
        anc = h0_anchor(f2)
        out[f"{q2bin[0]:g}-{q2bin[1]:g}"] = effect_sizes(f2, anc)
    return out


# ---------------------------------------------------------------------------
# Fig. 3: splay (lab) vs collapse (CM)
# ---------------------------------------------------------------------------

def slice_profiles(samples, var_x, var_y, edges, q2bin=Q2_MAIN,
                   variation="baseline"):
    """Profiles of var_y vs var_x per (config, W slice), truth level."""
    out = {}
    for (config, var, level), (d, meta) in samples.items():
        if var != variation or level != "truth":
            continue
        q2m = (d["Q2"] >= q2bin[0]) & (d["Q2"] < q2bin[1])
        for wlo, whi in W_SLICES[config]:
            wm = q2m & (d["W"] >= wlo) & (d["W"] < whi)
            rows = []
            for xlo, xhi in zip(edges[:-1], edges[1:]):
                m = wm & (d[var_x] >= xlo) & (d[var_x] < xhi)
                mean, sem, rms, n = profile(d[var_y], m)
                if n >= MIN_JETS_PER_BIN:
                    rows.append({"x_lo": xlo, "x_hi": xhi,
                                 "x_center": 0.5 * (xlo + xhi),
                                 "x_mean": float(np.mean(d[var_x][m])),
                                 "mean": mean, "sem": sem, "n": n})
            if rows:
                out[f"{config}|W{wlo:g}-{whi:g}"] = rows
    return out


def collapse_chi2(slices, edges):
    """
    Universality metrics over bins where >= 2 (config, W-slice) profiles
    contribute:
      * chi^2/ndf against the bin-by-bin global weighted mean
        (statistics-sensitive: inflates with MC sample size);
      * RMS fractional spread of slice means around the global mean
        (statistics-independent physics metric quoted in the paper).
    """
    centers = 0.5 * (np.asarray(edges[:-1]) + np.asarray(edges[1:]))
    per_bin = {c: [] for c in centers}
    for rows in slices.values():
        n_max = max(r["n"] for r in rows)
        for r in rows:
            if r["n"] >= BULK_FRACTION * n_max:
                per_bin[r["x_center"]].append((r["mean"], r["sem"]))
    chi2 = 0.0
    ndf = 0
    resid = []
    for c, pts in per_bin.items():
        if len(pts) < 2:
            continue
        means = np.array([p[0] for p in pts])
        sems = np.array([p[1] for p in pts])
        w = 1.0 / sems**2
        gmean = np.sum(w * means) / np.sum(w)
        chi2 += float(np.sum((means - gmean)**2 * w))
        ndf += len(pts) - 1
        resid.extend((means / gmean - 1.0).tolist())
    rms_spread = float(np.sqrt(np.mean(np.square(resid)))) if resid else np.nan
    return chi2, ndf, rms_spread


# ---------------------------------------------------------------------------
# Fig. 4: significance projection vs luminosity (reco level)
# ---------------------------------------------------------------------------

def luminosity_projection(samples, q2bin=Q2_MAIN):
    """
    Expected significance of rejecting the frame-independent null H0
    with EIC data, vs integrated luminosity per beam configuration.

    Uses reco-level (smeared, track-jet) baseline samples.  Bin counts
    are scaled from the generator cross section; H0 is anchored at the
    lowest populated W bin per |p|_lab bin (combined configs); the H0
    anchor uncertainty and an optional per-bin systematics floor are
    included.
    """
    f2 = fig2_profiles(samples, q2bin=q2bin)
    anchors = h0_anchor(f2, level="reco")

    # jets per fb^-1 per bin
    rates = {}
    sample_info = {}
    for skey, entry in f2.items():
        config, variation, level = skey.split("|")
        if variation != "baseline" or level != "reco":
            continue
        d, meta = samples[(config, "baseline", "reco")]
        sigma_fb = meta["sigma_gen_mb"] * 1e12
        per_event = sigma_fb / meta["n_tried"]
        rates[config] = {
            pkey: [dict(r, rate=r["n"] * per_event) for r in rows]
            for pkey, rows in entry.items()
        }
        sample_info[config] = {
            "sigma_fb": sigma_fb,
            "events_per_fb": sigma_fb,                       # DIS events
            "jets_per_fb": meta["n_jets"] * per_event,       # selected jets
            "annual_fb": ANNUAL_FB.get(config),
            "peak_lumi_1e33": PEAK_LUMI_1E33.get(config),
            "mc_equiv_fb": meta["n_tried"] / sigma_fb,       # MC sample size
        }

    results = {"lumi_grid": LUMI_GRID.tolist(), "configs": {}}
    chi2_grid = {}      # config -> (nL,) chi2 arrays: stat, +syst, +syst_cons
    for config, entry in rates.items():
        c_stat = np.zeros(len(LUMI_GRID))
        c_syst = np.zeros(len(LUMI_GRID))
        c_cons = np.zeros(len(LUMI_GRID))
        for pkey, rows in entry.items():
            if pkey not in anchors:
                continue
            anc = anchors[pkey]
            for r in rows:
                # exclude the anchor bin itself
                if r["W_lo"] == anc["W_lo"]:
                    continue
                delta = r["mean"] - anc["value"]
                for iL, L in enumerate(LUMI_GRID):
                    n_exp = r["rate"] * L
                    if n_exp < 25:
                        continue
                    stat2 = r["rms"]**2 / n_exp + anc["sem"]**2
                    c_stat[iL] += delta**2 / stat2
                    for floor, acc in ((SYST_FLOOR, c_syst),
                                       (SYST_FLOOR_CONS, c_cons)):
                        syst2 = stat2 + (floor * r["mean"])**2 \
                            + (floor * anc["value"])**2
                        acc[iL] += delta**2 / syst2
        chi2_grid[config] = (c_stat, c_syst, c_cons)
        results["configs"][config] = {
            "significance_stat": np.sqrt(c_stat).tolist(),
            "significance_syst": np.sqrt(c_syst).tolist(),
            "significance_syst_cons": np.sqrt(c_cons).tolist(),
        }

    if chi2_grid:
        comb = [np.sum([v[i] for v in chi2_grid.values()], axis=0)
                for i in range(3)]
    else:
        comb = [np.zeros(len(LUMI_GRID))] * 3
    results["combined"] = {
        "significance_stat": np.sqrt(comb[0]).tolist(),
        "significance_syst": np.sqrt(comb[1]).tolist(),
        "significance_syst_cons": np.sqrt(comb[2]).tolist(),
    }

    results["luminosity_inputs"] = {
        "source": "EIC Yellow Report (arXiv:2103.05419) Table 10.1 / CDR, "
                  "high-divergence configuration; "
                  f"operations year = {ANNUAL_SECONDS:.0e} s",
        "per_config": sample_info,
    }

    def l_5sigma(sig):
        sig = np.asarray(sig)
        above = np.where(sig >= 5.0)[0]
        return float(LUMI_GRID[above[0]]) if len(above) else None

    results["L_5sigma_fb"] = {
        "combined_stat": l_5sigma(results["combined"]["significance_stat"]),
        "combined_syst": l_5sigma(results["combined"]["significance_syst"]),
        "combined_syst_cons":
            l_5sigma(results["combined"]["significance_syst_cons"]),
        **{f"{c}_syst": l_5sigma(v["significance_syst"])
           for c, v in results["configs"].items()},
    }
    # running time to 5 sigma at design luminosity, per configuration
    results["time_to_5sigma_hours"] = {}
    for c, v in results["configs"].items():
        L5 = l_5sigma(v["significance_syst"])
        if L5 is not None and c in ANNUAL_FB:
            results["time_to_5sigma_hours"][c] = \
                L5 / ANNUAL_FB[c] * ANNUAL_SECONDS / 3600.0
    return results


# ---------------------------------------------------------------------------
# Fig. 1 support: boost factor and p_CM mapping
# ---------------------------------------------------------------------------

def boost_map(samples):
    out = {}
    for (config, variation, level), (d, meta) in samples.items():
        if variation != "baseline" or level != "truth":
            continue
        rows = []
        for wlo, whi in zip(W_EDGES[:-1], W_EDGES[1:]):
            m = (d["W"] >= wlo) & (d["W"] < whi)
            if m.sum() < MIN_JETS_PER_BIN:
                continue
            rows.append({"W_center": 0.5 * (wlo + whi),
                         "gamma_mean": float(np.mean(d["gboost"][m])),
                         "pcm_mean": float(np.mean(d["pcm"][m])),
                         "plab_mean": float(np.mean(d["plab"][m])),
                         "n": int(m.sum())})
        out[config] = rows
    return out


# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--datadir", default="data")
    ap.add_argument("--output", default="results/results.json")
    args = ap.parse_args()

    print("Loading jet tables …")
    samples = load_samples(args.datadir)

    print("Fig.2 profiles …")
    f2_truth = fig2_profiles(samples)
    anchors = h0_anchor(f2_truth)
    eff = effect_sizes(f2_truth, anchors)
    eff_reco = effect_sizes(f2_truth, h0_anchor(f2_truth, "reco"), "reco")
    eff_herwig = effect_sizes(
        f2_truth, h0_anchor(f2_truth, variation="herwig"),
        variation="herwig")

    print("Q^2 cross-check …")
    q2x = q2_crosscheck(samples)

    print("Fig.3 splay/collapse …")
    splay = slice_profiles(samples, "plab", "n90lab", PLAB_FINE)
    collapse = slice_profiles(samples, "pcm", "n90cm", PCM_EDGES)
    chi2_lab, ndf_lab, spread_lab = collapse_chi2(splay, PLAB_FINE)
    chi2_cm, ndf_cm, spread_cm = collapse_chi2(collapse, PCM_EDGES)

    # same metric for the cluster-hadronization (Herwig) samples, if present
    herwig_univ = None
    splay_h = slice_profiles(samples, "plab", "n90lab", PLAB_FINE,
                             variation="herwig")
    if splay_h:
        collapse_h = slice_profiles(samples, "pcm", "n90cm", PCM_EDGES,
                                    variation="herwig")
        _, _, sp_lab_h = collapse_chi2(splay_h, PLAB_FINE)
        _, _, sp_cm_h = collapse_chi2(collapse_h, PCM_EDGES)
        herwig_univ = {"spread_lab": sp_lab_h, "spread_cm": sp_cm_h,
                       "restoration_factor": sp_lab_h / max(sp_cm_h, 1e-12)}

    print("Fig.4 luminosity projection …")
    proj = luminosity_projection(samples)

    print("Fig.1 boost map …")
    bmap = boost_map(samples)

    results = {
        "binnings": {
            "Q2_main": Q2_MAIN, "Q2_bins": Q2_BINS,
            "plab_bins": PLAB_BINS, "W_edges": W_EDGES.tolist(),
            "W_slices": W_SLICES, "pcm_edges": PCM_EDGES.tolist(),
            "plab_fine": PLAB_FINE.tolist(),
            "syst_floor": SYST_FLOOR,
            "syst_floor_cons": SYST_FLOOR_CONS,
            "bulk_fraction": BULK_FRACTION,
        },
        "samples": {f"{k[0]}|{k[1]}|{k[2]}": m["n_jets"]
                    for k, (_, m) in samples.items()},
        "fig2": f2_truth,
        "h0_anchors": anchors,
        "effect_sizes_truth": eff,
        "effect_sizes_reco": eff_reco,
        "effect_sizes_herwig": eff_herwig,
        "q2_crosscheck": q2x,
        "fig3": {
            "splay": splay, "collapse": collapse,
            "chi2_lab": chi2_lab, "ndf_lab": ndf_lab,
            "chi2ndf_lab": chi2_lab / max(ndf_lab, 1),
            "chi2_cm": chi2_cm, "ndf_cm": ndf_cm,
            "chi2ndf_cm": chi2_cm / max(ndf_cm, 1),
            "spread_lab": spread_lab, "spread_cm": spread_cm,
            "restoration_factor": spread_lab / max(spread_cm, 1e-12),
            "herwig": herwig_univ,
        },
        "fig4": proj,
        "fig1": bmap,
    }

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    with open(args.output, "w") as f:
        json.dump(results, f, indent=1)
    print(f"\nWrote {args.output}")

    # ── Headline numbers ───────────────────────────────────────────────────
    print("\n=== HEADLINE RESULTS ===")
    for pkey, e in eff.items():
        print(f"  FFS shift, plab {pkey} GeV: "
              f"{100*e['shift']:+.1f} +- {100*e['err']:.1f} %   "
              f"(W {e['W_low'][0]:g}-{e['W_low'][1]:g} -> "
              f"{e['W_high'][0]:g}-{e['W_high'][1]:g} GeV)")
    print(f"  Universality RMS spread: lab = {100*results['fig3']['spread_lab']:.1f}%, "
          f"CM = {100*results['fig3']['spread_cm']:.1f}%  "
          f"(restoration x{results['fig3']['restoration_factor']:.1f})")
    print(f"  [chi2/ndf at MC stats: lab = {results['fig3']['chi2ndf_lab']:.1f}, "
          f"CM = {results['fig3']['chi2ndf_cm']:.1f}]")
    if herwig_univ:
        print(f"  Herwig universality RMS spread: "
              f"lab = {100*herwig_univ['spread_lab']:.1f}%, "
              f"CM = {100*herwig_univ['spread_cm']:.1f}%")
    for pkey, e in eff_herwig.items():
        print(f"  Herwig FFS shift, plab {pkey} GeV: "
              f"{100*e['shift']:+.1f} +- {100*e['err']:.1f} %")
    for k, v in proj["L_5sigma_fb"].items():
        print(f"  L(5 sigma) {k}: {v if v is not None else '> grid'} fb^-1")
    for c, h in proj.get("time_to_5sigma_hours", {}).items():
        ann = ANNUAL_FB.get(c)
        print(f"  time to 5 sigma, {c} at design lumi ({ann:g} fb^-1/yr): "
              f"{h:.2g} h")
    for c, si in proj.get("luminosity_inputs", {}).get("per_config", {}).items():
        print(f"  {c}: sigma = {si['sigma_fb']/1e6:.1f} nb, "
              f"{si['jets_per_fb']:.3g} jets/fb^-1, "
              f"MC sample = {si['mc_equiv_fb']*1000:.1f} pb^-1 equivalent")


if __name__ == "__main__":
    main()
