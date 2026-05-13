"""
DIS kinematics utilities for the EIC FFS (Frame-dependent Fragmentation Shift) study.

Implements 4-vector operations, DIS invariant reconstruction, and frame boosts
(lab → Breit, lab → photon-proton CM) used to study the FFS effect.

Reference: arXiv:2308.10951
"""

import numpy as np
from typing import Union


# Proton mass (GeV)
M_PROTON = 0.938272


def four_dot(a: np.ndarray, b: np.ndarray) -> Union[float, np.floating]:
    """
    Minkowski inner product using (+,-,-,-) metric.

    Parameters
    ----------
    a, b : array_like, shape (4,)
        Four-vectors in (px, py, pz, E) convention.

    Returns
    -------
    float
        a·b = E_a·E_b - px_a·px_b - py_a·py_b - pz_a·pz_b
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    return a[3] * b[3] - a[0] * b[0] - a[1] * b[1] - a[2] * b[2]


def invariant_mass2(p: np.ndarray) -> Union[float, np.floating]:
    """Return m² = p·p for a 4-vector (px, py, pz, E)."""
    return four_dot(p, p)


def lorentz_boost_to_rest(p_boost: np.ndarray, p_target: np.ndarray) -> np.ndarray:
    """
    Boost a 4-vector p_target into the rest frame of p_boost.

    Parameters
    ----------
    p_boost : array_like, shape (4,)
        The 4-vector whose rest frame we want (px, py, pz, E).
    p_target : array_like, shape (4,)
        The 4-vector to be boosted.

    Returns
    -------
    np.ndarray, shape (4,)
        p_target in the rest frame of p_boost.
    """
    p_boost = np.asarray(p_boost, dtype=float)
    p_target = np.asarray(p_target, dtype=float)

    E_b = p_boost[3]
    pvec_b = p_boost[:3]
    mass_b = np.sqrt(max(invariant_mass2(p_boost), 1e-12))

    beta = pvec_b / E_b
    beta2 = np.dot(beta, beta)
    gamma = E_b / mass_b

    # Lorentz boost formula
    beta_dot_p = np.dot(beta, p_target[:3])
    factor = gamma * (gamma * beta_dot_p / (gamma + 1) - p_target[3])

    out = np.empty(4)
    out[:3] = p_target[:3] + factor * beta
    out[3] = gamma * (p_target[3] - beta_dot_p)
    return out


def boost_to_breit_frame(p_proton: np.ndarray, p_photon: np.ndarray,
                          p_target: np.ndarray) -> np.ndarray:
    """
    Boost a 4-vector to the Breit (brick-wall) frame.

    In the Breit frame: q = (0, 0, -Q, 0) (photon has no energy,
    carries only 3-momentum −Q in the z-direction).  The struck quark
    comes in with pz = +Q/2 and is reflected to pz = -Q/2.

    Parameters
    ----------
    p_proton : array_like, shape (4,)
        Incoming proton 4-vector in the lab frame.
    p_photon : array_like, shape (4,)
        Virtual photon 4-vector q = k - k' in the lab frame.
    p_target : array_like, shape (4,)
        Particle 4-vector to be boosted.

    Returns
    -------
    np.ndarray, shape (4,)
        4-vector in the Breit frame.
    """
    p_proton = np.asarray(p_proton, dtype=float)
    p_photon = np.asarray(p_photon, dtype=float)
    p_target = np.asarray(p_target, dtype=float)

    # The Breit-frame boost is defined by the condition that in that frame
    # 2x*P + q = 0, i.e. the proton 4-momentum (scaled by 2x) + photon = 0.
    # The system to boost to rest is: q + 2x*P
    Q2 = -four_dot(p_photon, p_photon)
    Pdotq = four_dot(p_proton, p_photon)
    if Pdotq == 0:
        return p_target.copy()
    x = Q2 / (2.0 * Pdotq)

    # 4-vector whose rest frame is the Breit frame
    p_breit_system = p_photon + 2.0 * x * p_proton
    return lorentz_boost_to_rest(p_breit_system, p_target)


class DISKinematics:
    """
    Reconstructed DIS kinematics from a single event.

    Attributes
    ----------
    k_in   : 4-vector of incoming lepton (px,py,pz,E)
    k_out  : 4-vector of scattered lepton
    P_in   : 4-vector of incoming proton
    q      : 4-vector of virtual photon  (q = k_in - k_out)
    Q2     : photon virtuality  -q²  (GeV²)
    W      : invariant mass of photon-proton system  (GeV)
    x      : Bjorken-x
    y      : inelasticity
    nu     : lepton energy loss (lab frame)  (GeV)
    valid  : True if all invariants are physical
    """

    def __init__(self, k_in, k_out, P_in):
        self.k_in = np.asarray(k_in, dtype=float)
        self.k_out = np.asarray(k_out, dtype=float)
        self.P_in = np.asarray(P_in, dtype=float)
        self.valid = False
        self._compute()

    def _compute(self):
        """Compute all DIS invariants."""
        self.q = self.k_in - self.k_out

        self.Q2 = -four_dot(self.q, self.q)
        Pdotq = four_dot(self.P_in, self.q)
        Pdotk = four_dot(self.P_in, self.k_in)
        Mp2 = four_dot(self.P_in, self.P_in)

        if self.Q2 <= 0 or Pdotq <= 0 or Pdotk <= 0:
            return

        self.y = Pdotq / Pdotk
        self.x = self.Q2 / (2.0 * Pdotq)
        W2 = Mp2 + 2.0 * Pdotq - self.Q2
        if W2 <= 0:
            return
        self.W = float(np.sqrt(W2))

        # Lepton energy loss in the proton rest frame
        Mp = float(np.sqrt(max(Mp2, 0.0)))
        self.nu = Pdotq / Mp if Mp > 0 else 0.0

        # Validity checks
        if not (0 < self.x <= 1) or not (0 < self.y < 1):
            return
        if self.Q2 < 0.5 or self.W < 1.0:
            return

        self.valid = True

    def boost_to_breit(self, p_target):
        """Boost p_target to the Breit frame."""
        return boost_to_breit_frame(self.P_in, self.q, p_target)

    def boost_to_gamma_p_cm(self, p_target):
        """Boost p_target to the photon-proton CM frame."""
        p_system = self.P_in + self.q
        return lorentz_boost_to_rest(p_system, p_target)

    def __repr__(self):
        if self.valid:
            return (f"DISKinematics(Q2={self.Q2:.2f} GeV², "
                    f"W={self.W:.2f} GeV, x={self.x:.4f}, y={self.y:.3f})")
        return "DISKinematics(invalid)"
