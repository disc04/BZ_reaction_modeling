"""Classic three-variable Oregonator (Field & Noyes, 1974).

Mechanism (A = BrO3-, B = organic substrate, e.g. malonic acid):

    O1:  A + Y  -> X + P        rate k1*A*Y
    O2:  X + Y  -> 2P           rate k2*X*Y
    O3:  A + X  -> 2X + 2Z      rate k3*A*X   (autocatalysis)
    O4:  2X     -> A + P        rate k4*X^2
    O5:  B + Z  -> (f/2) Y      rate kc*B*Z   (catalyst reset, Br- regeneration)

Dynamic variables (molar concentrations):
    x = [HBrO2]   activator
    y = [Br-]     inhibitor
    z = [Ce(IV)]  oxidised catalyst (or ferriin with ferroin)

A and B are treated as constant (pool chemical approximation).

Conventions follow Field & Noyes (1974). O3 is the net of
    BrO3- + HBrO2 + H+       -> 2 BrO2* + H2O
    2 BrO2* + 2 Ce3+ + 2 H+  -> 2 HBrO2 + 2 Ce4+
so bromate (A) is the reactant and two Ce4+ are formed per event.
Vilcu & Bala (2004) write B + X -> 2X + Z and Z -> fY; with the f/2 in O5
both forms give identical x, y dynamics and differ only in the scale of z.
Vilcu & Bala's eq. 5.2 also misprint dz/dt = k3*b*z - k5*z (should be in x).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp


@dataclass(frozen=True)
class OregonatorParams:
    """Rate constants and pool concentrations (units: M, s).

    Defaults are the classic Field-Noyes values with [H+] = 0.8 M folded
    into the rate constants. kc is a lumped, empirical constant for the
    complex catalyst-reduction chemistry; together with B it sets the
    oscillation period (period ~ 8 / (kc*B)).
    """

    k1: float = 1.28      # M^-1 s^-1
    k2: float = 2.4e6     # M^-1 s^-1
    k3: float = 33.6      # M^-1 s^-1
    k4: float = 2.4e3     # M^-1 s^-1
    kc: float = 1.0       # M^-1 s^-1
    A: float = 0.06       # [BrO3-], M
    B: float = 0.02       # [organic substrate], M
    f: float = 1.0        # Br- produced per 2 Ce4+ reduced; oscillates for ~0.5 < f < 2.4


def oregonator_rhs(t: float, c: np.ndarray, p: OregonatorParams) -> list[float]:
    """Right-hand side of the Oregonator rate equations."""
    x, y, z = c
    r1 = p.k1 * p.A * y
    r2 = p.k2 * x * y
    r3 = p.k3 * p.A * x
    r4 = p.k4 * x * x
    r5 = p.kc * p.B * z
    dx = r1 - r2 + r3 - 2.0 * r4
    dy = -r1 - r2 + 0.5 * p.f * r5
    dz = 2.0 * r3 - r5
    return [dx, dy, dz]


def simulate(
    params: OregonatorParams | None = None,
    t_end: float = 1500.0,
    n_points: int = 15001,
    c0: tuple[float, float, float] = (1e-10, 1e-6, 1e-6),
) -> tuple[np.ndarray, np.ndarray]:
    """Integrate the Oregonator.

    The system is stiff (rate constants span ~8 orders of magnitude),
    so LSODA with tight tolerances is used.

    Returns
    -------
    t : (n_points,) array of times in s
    c : (3, n_points) array of concentrations [x, y, z] in M
    """
    p = params or OregonatorParams()
    t_eval = np.linspace(0.0, t_end, n_points)
    sol = solve_ivp(
        oregonator_rhs,
        (0.0, t_end),
        c0,
        args=(p,),
        method="LSODA",
        t_eval=t_eval,
        rtol=1e-8,
        atol=1e-14,
    )
    if not sol.success:
        raise RuntimeError(f"Integration failed: {sol.message}")
    return sol.t, sol.y


def estimate_period(t: np.ndarray, x: np.ndarray, discard: float = 0.25) -> float:
    """Mean oscillation period from peaks of log10(x).

    The first `discard` fraction of the trace is skipped to remove the
    initial transient. Returns NaN if fewer than two peaks are found
    (i.e. the system has settled to a steady state).
    """
    from scipy.signal import find_peaks

    start = int(len(t) * discard)
    peaks, _ = find_peaks(np.log10(np.clip(x[start:], 1e-30, None)), prominence=1.0)
    if len(peaks) < 2:
        return float("nan")
    return float(np.mean(np.diff(t[start:][peaks])))


def jacobian(c: np.ndarray, p: OregonatorParams) -> np.ndarray:
    """Analytical Jacobian d(rhs)/d(x, y, z) at state c."""
    x, y, _ = c
    return np.array([
        [-p.k2 * y + p.k3 * p.A - 4.0 * p.k4 * x, p.k1 * p.A - p.k2 * x, 0.0],
        [-p.k2 * y, -p.k1 * p.A - p.k2 * x, 0.5 * p.f * p.kc * p.B],
        [2.0 * p.k3 * p.A, 0.0, -p.kc * p.B],
    ])


def steady_state(p: OregonatorParams | None = None) -> np.ndarray:
    """Positive steady state (x*, y*, z*).

    Setting dz/dt = dy/dt = 0 gives z and y as functions of x:
        z = 2*k3*A*x / (kc*B)
        y = f*k3*A*x / (k1*A + k2*x)
    Substituting into dx/dt = 0 leaves one equation in x, solved by bisection.
    """
    from scipy.optimize import brentq

    p = p or OregonatorParams()

    def y_of(x):
        return p.f * p.k3 * p.A * x / (p.k1 * p.A + p.k2 * x)

    def g(x):  # dx/dt divided by x
        return (p.k1 * p.A - p.k2 * x) * y_of(x) / x + p.k3 * p.A - 2.0 * p.k4 * x

    x = brentq(g, 1e-15, p.k3 * p.A / (2.0 * p.k4), xtol=1e-20, rtol=1e-12)
    return np.array([x, y_of(x), 2.0 * p.k3 * p.A * x / (p.kc * p.B)])
