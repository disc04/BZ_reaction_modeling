"""Reduced two-variable Oregonator (Tyson & Fife, 1980).

Starting from the Field-Noyes Oregonator, concentrations and time are
scaled (Tyson, 1985):

    x = [HBrO2] / X0,   X0 = k3*A / (2*k4)
    y = [Br-]   / Y0,   Y0 = k3*A / k2
    z = [M_ox]  / Z0,   Z0 = (k3*A)^2 / (k4*kc*B)
    tau = kc*B * t

giving the dimensionless three-variable model

    eps  dx/dtau = q*y - x*y + x*(1 - x)
    eps' dy/dtau = -q*y - x*y + f*z
         dz/dtau = x - z

with eps = kc*B/(k3*A), eps' = 2*k4*kc*B/(k2*k3*A), q = 2*k1*k4/(k2*k3).

Because eps' << eps (Br- is by far the fastest variable), y follows its
quasi-steady state y = f*z / (q + x). Substituting gives two variables:

    eps dx/dtau = x*(1 - x) - f*z*(x - q)/(x + q)
        dz/dtau = x - z

The x-nullcline is N-shaped; the z-nullcline is the line z = x.
Oscillations occur when the steady state sits on the middle (repelling)
branch of the N, between the two Hopf points.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from .oregonator import OregonatorParams


@dataclass(frozen=True)
class ReducedParams:
    """Dimensionless parameters of the two-variable Oregonator."""

    eps: float = 9.92e-3
    q: float = 7.62e-5
    f: float = 1.0
    time_scale: float = 0.02   # kc*B in 1/s, converts tau to seconds


def from_oregonator(p: OregonatorParams | None = None) -> ReducedParams:
    """Dimensionless parameters from Field-Noyes rate constants."""
    p = p or OregonatorParams()
    return ReducedParams(
        eps=p.kc * p.B / (p.k3 * p.A),
        q=2.0 * p.k1 * p.k4 / (p.k2 * p.k3),
        f=p.f,
        time_scale=p.kc * p.B,
    )


def reduced_rhs(tau: float, s: np.ndarray, r: ReducedParams) -> list[float]:
    """Right-hand side of the two-variable Oregonator."""
    x, z = s
    dx = (x * (1.0 - x) - r.f * z * (x - r.q) / (x + r.q)) / r.eps
    dz = x - z
    return [dx, dz]


def simulate_reduced(
    r: ReducedParams | None = None,
    tau_end: float = 30.0,
    n_points: int = 30001,
    s0: tuple[float, float] = (0.01, 0.01),
) -> tuple[np.ndarray, np.ndarray]:
    """Integrate the reduced model. Returns (tau, [x, z])."""
    r = r or ReducedParams()
    tau = np.linspace(0.0, tau_end, n_points)
    sol = solve_ivp(reduced_rhs, (0.0, tau_end), s0, args=(r,), method="LSODA",
                    t_eval=tau, rtol=1e-9, atol=1e-14)
    if not sol.success:
        raise RuntimeError(f"Integration failed: {sol.message}")
    return sol.t, sol.y


def x_nullcline(x: np.ndarray, r: ReducedParams) -> np.ndarray:
    """z on the x-nullcline (dx/dtau = 0), defined for x > q."""
    return x * (1.0 - x) * (x + r.q) / (r.f * (x - r.q))


def steady_state_reduced(r: ReducedParams) -> np.ndarray:
    """Positive steady state (x*, z*): intersection of z = x with the x-nullcline.

    With z = x the condition reduces to the quadratic
        x^2 + (f - 1 + q) x - q (1 + f) = 0.
    """
    b = r.f - 1.0 + r.q
    x = 0.5 * (-b + np.sqrt(b * b + 4.0 * r.q * (1.0 + r.f)))
    return np.array([x, x])


def jacobian_reduced(s: np.ndarray, r: ReducedParams) -> np.ndarray:
    """Analytical Jacobian of the reduced model."""
    x, z = s
    g = (x - r.q) / (x + r.q)
    dg = 2.0 * r.q / (x + r.q) ** 2
    return np.array([
        [(1.0 - 2.0 * x - r.f * z * dg) / r.eps, -r.f * g / r.eps],
        [1.0, -1.0],
    ])


def hopf_points(r: ReducedParams, f_range: tuple[float, float] = (0.3, 3.0)) -> list[float]:
    """Values of f where the steady state changes stability (trace J = 0)."""
    from dataclasses import replace

    def trace(f):
        rf = replace(r, f=f)
        return np.trace(jacobian_reduced(steady_state_reduced(rf), rf))

    grid = np.linspace(*f_range, 2001)
    tr = np.array([trace(f) for f in grid])
    roots = []
    for i in np.where(np.sign(tr[:-1]) != np.sign(tr[1:]))[0]:
        roots.append(brentq(trace, grid[i], grid[i + 1], xtol=1e-12))
    return roots
