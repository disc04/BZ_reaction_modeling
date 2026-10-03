"""Oregonator with slow reactant depletion (closed batch reactor).

The classic Oregonator holds bromate (A) and malonic acid (B) constant, so it
cannot age. Here both are dynamic:

    dA/dt = s * (-r1 - r3 + r4)      bromate used in O1, O3; returned in O4
    dB/dt = s * (-r5)                malonic acid used in O5

k_c follows [BrO3-] (k_c = k_c0 * A / A0), the hypothesis that best fitted
the acid series in Test 1.

s (< 1) is a depletion scale factor. It is needed because the Oregonator does
not conserve catalyst: its ferriin/Ce4+ excursions (~10 mM) exceed the 2 mM of
ferroin actually present, so the raw consumption steps run far too fast.
s is calibrated against the observed oscillation lifetime.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks

from .recipe import params_from_recipe


def simulate_depleting(
    kc0: float,
    s: float,
    h2so4: float = 0.5,
    bromate: float = 0.47,
    malonic: float = 0.18,
    f: float = 1.0,
    t_end: float = 6000.0,
    dt: float = 0.05,
):
    """Integrate (x, y, z, A, B). Returns t (s) and the (5, n) state array."""
    p = params_from_recipe(h2so4, bromate, malonic, kc0, f)
    A0 = bromate

    def rhs(t, u):
        x, y, z, A, B = u
        kc = kc0 * A / A0
        r1 = p.k1 * A * y
        r2 = p.k2 * x * y
        r3 = p.k3 * A * x
        r4 = p.k4 * x * x
        r5 = kc * B * z
        return [r1 - r2 + r3 - 2 * r4,
                -r1 - r2 + 0.5 * f * r5,
                2 * r3 - r5,
                s * (-r1 - r3 + r4),
                -s * r5]

    t = np.arange(0.0, t_end + dt / 2, dt)
    sol = solve_ivp(rhs, (0.0, t_end), [1e-10, 1e-6, 1e-6, bromate, malonic],
                    method="LSODA", t_eval=t, rtol=1e-8, atol=1e-14)
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol.t, sol.y


def oscillation_peaks(t: np.ndarray, u: np.ndarray):
    """Peak times (s) and Ce4+ peak heights (M) of each oscillation."""
    pk, _ = find_peaks(np.log10(np.clip(u[0], 1e-30, None)), prominence=1.0)
    # catalyst peak: maximum of z within each cycle (z lags the HBrO2 spike)
    bounds = list(pk) + [len(t) - 1]
    zmax = np.array([u[2, a:b].max() for a, b in zip(bounds[:-1], bounds[1:])])
    return t[pk], zmax
