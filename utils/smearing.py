"""
Parametric detector smearing for the EIC FFS study (ANALYSIS_DESIGN.md Sec. 7).

Implements an ePIC-like fast simulation:
  * charged-particle tracking: momentum resolution and efficiency,
    eta-dependent (barrel / forward / far-forward);
  * scattered-electron EM calorimetry: energy resolution (angles assumed
    measured precisely by tracking);
  * acceptance: |eta| < 3.5, pT > 0.2 GeV for tracks.

The smeared ("reco-level") analysis uses charged tracks only (track-jets),
which is how a first-data EIC measurement of n90 would actually be done.
"""

import numpy as np

TRACK_PT_MIN = 0.2      # GeV
TRACK_ETA_MAX = 3.5
TRACK_EFF = 0.95

# sigma(p)/p = a*p (+) b, per eta region  [ePIC Yellow-Report-like]
_TRACK_RES = (
    (1.0, 0.0005, 0.005),   # |eta| < 1.0
    (2.5, 0.0005, 0.010),   # 1.0 < |eta| < 2.5
    (3.5, 0.0010, 0.020),   # 2.5 < |eta| < 3.5
)

# Electron EMCal: sigma(E)/E = a/sqrt(E) (+) b
ECAL_A = 0.02
ECAL_B = 0.01


def smear_tracks(px, py, pz, e, charge, rng):
    """
    Apply tracking efficiency, acceptance, and momentum smearing to
    charged particles.  Neutral particles are dropped (track-jet analysis).

    Returns (px, py, pz, e, charge) arrays of the surviving smeared tracks.
    """
    px = np.asarray(px, dtype=float)
    py = np.asarray(py, dtype=float)
    pz = np.asarray(pz, dtype=float)
    e = np.asarray(e, dtype=float)
    charge = np.asarray(charge, dtype=int)

    sel = charge != 0
    px, py, pz, e = px[sel], py[sel], pz[sel], e[sel]
    charge = charge[sel]
    if len(px) == 0:
        return px, py, pz, e, charge

    p = np.sqrt(px**2 + py**2 + pz**2)
    pt = np.sqrt(px**2 + py**2)
    with np.errstate(divide="ignore", invalid="ignore"):
        eta = np.arctanh(np.clip(pz / np.maximum(p, 1e-12), -1 + 1e-9, 1 - 1e-9))

    # Acceptance + efficiency
    keep = (pt > TRACK_PT_MIN) & (np.abs(eta) < TRACK_ETA_MAX)
    keep &= rng.random(len(px)) < TRACK_EFF
    px, py, pz, e, charge, p = px[keep], py[keep], pz[keep], e[keep], \
        charge[keep], p[keep]
    eta = eta[keep]
    if len(px) == 0:
        return px, py, pz, e, charge

    # Momentum resolution
    sigma_rel = np.empty(len(px))
    abs_eta = np.abs(eta)
    lo = 0.0
    for hi, a, b in _TRACK_RES:
        m = (abs_eta >= lo) & (abs_eta < hi)
        sigma_rel[m] = np.sqrt((a * p[m])**2 + b**2)
        lo = hi

    scale = 1.0 + sigma_rel * rng.standard_normal(len(px))
    scale = np.maximum(scale, 0.05)
    m2 = np.maximum(e**2 - p**2, 0.0)        # preserve particle mass
    px, py, pz = px * scale, py * scale, pz * scale
    e = np.sqrt((p * scale)**2 + m2)
    return px, py, pz, e, charge


def smear_electron(k_out, rng):
    """
    Smear the scattered-electron energy (EMCal resolution); direction
    is kept (tracking measures angles much better than the energy).
    """
    k = np.asarray(k_out, dtype=float)
    E = k[3]
    if E <= 0:
        return k
    sigma_rel = np.sqrt(ECAL_A**2 / E + ECAL_B**2)
    scale = max(1.0 + sigma_rel * rng.standard_normal(), 0.05)
    return k * scale          # massless scaling of the whole 4-vector
