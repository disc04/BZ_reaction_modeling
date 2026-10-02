"""Reaction-diffusion Oregonator in 2D: target and spiral waves.

The reduced two-variable Oregonator (see reduced.py) with diffusion:

    dx/dt = [x(1 - x) - f z (x - q)/(x + q)] / eps + Dx * lap(x)
    dz/dt = x - z                                   + Dz * lap(z)

x = scaled [HBrO2] (activator), z = scaled oxidised catalyst (ferriin / Ce4+).
Time and space are dimensionless. The medium is *excitable* when f lies
above the upper Hopf point: the rest state is stable, but a supra-threshold
kick triggers a full excursion that spreads to neighbours by diffusion of x.

f may be a 2D array, which allows heterogeneous media (e.g. a pacemaker
region with oscillatory kinetics embedded in an excitable medium).

Numerics: isotropic 9-point Laplacian, no-flux (Neumann) boundaries, explicit Euler.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .reduced import ReducedParams, steady_state_reduced


@dataclass(frozen=True)
class RDParams:
    eps: float = 0.05
    q: float = 0.002
    f: float = 2.5          # excitable (upper Hopf at f ~ 2.31 for these eps, q)
    Dx: float = 1.0
    Dz: float = 0.6
    n: int = 400            # grid points per side
    dx: float = 0.25        # grid spacing
    dt: float = 0.01        # time step (explicit; must satisfy dt < dx^2 / (4 Dx))

    def __post_init__(self):
        if self.dt > self.dx ** 2 / (4.0 * max(self.Dx, self.Dz)):
            raise ValueError("dt too large for explicit diffusion stability")


def laplacian(u: np.ndarray, dx: float) -> np.ndarray:
    """Isotropic 9-point Laplacian with no-flux boundaries.

    Weights 4 (edge neighbours), 1 (corner neighbours), -20 (centre), / 6 dx^2.
    Compared with the 5-point stencil, this removes the grid-aligned
    (square-looking) distortion of circular and spiral waves.
    """
    up = np.pad(u, 1, mode="edge")
    edges = up[:-2, 1:-1] + up[2:, 1:-1] + up[1:-1, :-2] + up[1:-1, 2:]
    corners = up[:-2, :-2] + up[:-2, 2:] + up[2:, :-2] + up[2:, 2:]
    return (4.0 * edges + corners - 20.0 * u) / (6.0 * dx**2)


def rest_state(p: RDParams) -> tuple[float, float]:
    """Homogeneous steady state of the excitable medium."""
    x, z = steady_state_reduced(ReducedParams(eps=p.eps, q=p.q, f=p.f))
    return float(x), float(z)


def step(x: np.ndarray, z: np.ndarray, p: RDParams, f: np.ndarray | float) -> None:
    """Advance (x, z) in place by one explicit Euler step."""
    react_x = (x * (1.0 - x) - f * z * (x - p.q) / (x + p.q)) / p.eps
    dx_ = react_x + p.Dx * laplacian(x, p.dx)
    dz_ = (x - z) + p.Dz * laplacian(z, p.dx)
    x += p.dt * dx_
    z += p.dt * dz_
    np.clip(x, p.q * 0.5, 1.0, out=x)


def run(
    x: np.ndarray,
    z: np.ndarray,
    p: RDParams,
    t_end: float,
    frame_every: float,
    f: np.ndarray | float | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Integrate to t_end, returning (times, z frames)."""
    f = p.f if f is None else f
    n_steps = int(round(t_end / p.dt))
    every = max(1, int(round(frame_every / p.dt)))
    times, frames = [], []
    for i in range(n_steps + 1):
        if i % every == 0:
            times.append(i * p.dt)
            frames.append(z.astype(np.float32))
        if i < n_steps:
            step(x, z, p, f)
    return np.array(times), np.array(frames)


# --- initial conditions -------------------------------------------------------

def _coords(p: RDParams):
    c = (np.arange(p.n) - p.n / 2) * p.dx
    return np.meshgrid(c, c, indexing="xy")


def target_setup(
    p: RDParams,
    pacemakers: tuple[tuple[float, float], ...] = ((-0.2, -0.15), (0.22, 0.18)),
    radius: float = 2.0,
    f_pacemaker: float = 1.0,
):
    """Excitable medium at rest with small oscillatory pacemaker regions.

    Inside each pacemaker disc f is lowered into the oscillatory window, so
    the region oscillates on its own and sends out concentric (target) waves,
    like a dust particle or bubble in the dish. Pacemaker centres are given
    as fractions of the domain size, measured from its centre.
    """
    L = p.n * p.dx
    X, Y = _coords(p)
    xr, zr = rest_state(p)
    x = np.full((p.n, p.n), xr)
    z = np.full((p.n, p.n), zr)
    f = np.full((p.n, p.n), p.f)
    for fx, fy in pacemakers:
        cx, cy = fx * L, fy * L
        f[(X - cx) ** 2 + (Y - cy) ** 2 < radius**2] = f_pacemaker
    return x, z, f


def spiral_setup(p: RDParams, width: float = 1.0):
    """Broken planar wave: its free end curls into a spiral.

    A vertical excited stripe (high x) with a refractory band (high z) on its
    left travels to the right. It is cut off in the upper half of the domain,
    leaving a free end in the centre.
    """
    X, Y = _coords(p)
    xr, zr = rest_state(p)
    x = np.full((p.n, p.n), xr)
    z = np.full((p.n, p.n), zr)
    front = (X > -width) & (X < width) & (Y < 0)
    back = (X > -6 * width) & (X <= -width) & (Y < 0)
    x[front] = 0.8
    z[back] = 0.4
    return x, z
