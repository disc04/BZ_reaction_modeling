"""Plotting utilities for BZ model output."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

# Fixed colour per species (validated categorical palette, CVD-safe).
SPECIES = (
    ("x", r"HBrO$_2$ (X)", "#eb6834"),
    ("y", r"Br$^-$ (Y)", "#2a78d6"),
    ("z", r"Ce$^{4+}$ (Z)", "#1baf7a"),
)

TEXT = "#3d3d3a"
GRID = "#e4e3dc"


def style_axes(ax: plt.Axes) -> None:
    """Recessive axes: light grid, no top/right spines."""
    ax.grid(color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color("#b5b4ac")
    ax.tick_params(colors=TEXT, labelsize=9)
    ax.xaxis.label.set_color(TEXT)
    ax.yaxis.label.set_color(TEXT)


def plot_timeseries(
    t: np.ndarray,
    c: np.ndarray,
    ax: plt.Axes | None = None,
    title: str | None = None,
    legend: bool = True,
    direct_labels: bool = True,
    ylim: tuple[float, float] = (1e-9, 1e-1),
) -> plt.Axes:
    """Plot all three Oregonator species on one log-scaled axis.

    Parameters
    ----------
    t : times (s)
    c : (3, n) concentrations [x, y, z] (M)
    ax : existing axes to draw on; a new figure is created if None
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(8, 4))

    for data, (_, label, color) in zip(c, SPECIES):
        ax.semilogy(t, data, color=color, lw=1.6, label=label)
        if direct_labels:
            ax.annotate(
                label, xy=(t[-1], data[-1]), xytext=(4, 0),
                textcoords="offset points", va="center", fontsize=9, color=TEXT,
                annotation_clip=False,
            )

    ax.set_ylim(*ylim)
    ax.set_xlim(t[0], t[-1])
    ax.set_xlabel("time (s)")
    ax.set_ylabel("concentration (M)")
    if title:
        ax.set_title(title, loc="left", fontsize=11, color=TEXT)
    if legend:
        ax.legend(loc="upper left", ncol=3, frameon=False, fontsize=9,
                  bbox_to_anchor=(0, 1.0))
    style_axes(ax)
    return ax


CYCLE = "#2a78d6"
TRANSIENT = "#a3a29a"
AXIS_LABELS = {
    "x": r"log$_{10}$ [HBrO$_2$]",
    "y": r"log$_{10}$ [Br$^-$]",
    "z": r"log$_{10}$ [Ce$^{4+}$]",
}


def plot_phase_3d(
    c: np.ndarray,
    ax=None,
    color: str = CYCLE,
    lw: float = 1.8,
    alpha: float = 1.0,
    label: str | None = None,
):
    """Trajectory in (log x, log y, log z) phase space on a 3D axes.

    Log coordinates are used because concentrations span ~6 decades;
    matplotlib's 3D axes do not support log scaling, so log10 values
    are plotted directly.
    """
    if ax is None:
        fig = plt.figure(figsize=(6, 5))
        ax = fig.add_subplot(projection="3d")
    lx, ly, lz = np.log10(np.clip(c, 1e-30, None))
    ax.plot(lx, ly, lz, color=color, lw=lw, alpha=alpha, label=label)
    ax.set_xlabel(AXIS_LABELS["x"], fontsize=9, color=TEXT, labelpad=6)
    ax.set_ylabel(AXIS_LABELS["y"], fontsize=9, color=TEXT, labelpad=6)
    ax.set_zlabel(AXIS_LABELS["z"], fontsize=9, color=TEXT, labelpad=6)
    ax.tick_params(colors=TEXT, labelsize=8)
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.set_pane_color((1, 1, 1, 0))
        axis._axinfo["grid"]["color"] = GRID
    return ax


def plot_phase_2d(
    c: np.ndarray,
    ax: plt.Axes | None = None,
    xvar: str = "z",
    yvar: str = "y",
    color: str = CYCLE,
    lw: float = 1.8,
    alpha: float = 1.0,
    label: str | None = None,
) -> plt.Axes:
    """2D projection of a trajectory on log-log axes (e.g. Br- vs Ce4+)."""
    if ax is None:
        _, ax = plt.subplots(figsize=(5, 5))
    idx = {"x": 0, "y": 1, "z": 2}
    names = {k: v.replace(r"log$_{10}$ ", "") for k, v in AXIS_LABELS.items()}
    ax.loglog(c[idx[xvar]], c[idx[yvar]], color=color, lw=lw, alpha=alpha, label=label)
    ax.set_xlabel(f"{names[xvar]} (M)")
    ax.set_ylabel(f"{names[yvar]} (M)")
    style_axes(ax)
    return ax


def plot_phase_plane(
    r,
    ax: plt.Axes | None = None,
    trajectory: np.ndarray | None = None,
    x_range: tuple[float, float] = (1e-5, 1.0),
    z_range: tuple[float, float] = (1e-5, 10.0),
) -> plt.Axes:
    """Nullclines of the reduced Oregonator on log-log axes.

    Parameters
    ----------
    r : ReducedParams
    trajectory : optional (2, n) array [x, z] drawn on top (e.g. the limit cycle)
    """
    from .reduced import jacobian_reduced, steady_state_reduced, x_nullcline

    if ax is None:
        _, ax = plt.subplots(figsize=(5, 5))

    x = np.logspace(np.log10(r.q * 1.0005), 0, 4000)
    zn = x_nullcline(x, r)
    ax.loglog(x, zn, color=SPECIES[0][2], lw=1.8, label=r"$\dot{x}=0$ (N-shaped)")
    xl = np.logspace(np.log10(x_range[0]), np.log10(x_range[1]), 50)
    ax.loglog(xl, xl, color=SPECIES[2][2], lw=1.8, label=r"$\dot{z}=0$ ($z=x$)")

    if trajectory is not None:
        ax.loglog(trajectory[0], trajectory[1], color=CYCLE, lw=2.0, label="trajectory")

    ss = steady_state_reduced(r)
    stable = np.all(np.linalg.eigvals(jacobian_reduced(ss, r)).real < 0)
    ax.plot(ss[0], ss[1], "o", ms=8, mew=1.8, mec=TEXT,
            mfc=TEXT if stable else "white",
            label="stable steady state" if stable else "unstable steady state")

    ax.set_xlim(*x_range)
    ax.set_ylim(*z_range)
    ax.set_xlabel(r"x  (scaled [HBrO$_2$])")
    ax.set_ylabel(r"z  (scaled [Ce$^{4+}$])")
    style_axes(ax)
    return ax


# --- spatial patterns ---------------------------------------------------------

def ferroin_cmap():
    """Colour map mimicking the ferroin indicator: red (reduced) -> blue (oxidised).

    This is deliberately a physical-colour mimic rather than a perceptual
    data ramp: it reproduces what is seen in the Petri dish.
    """
    from matplotlib.colors import LinearSegmentedColormap

    return LinearSegmentedColormap.from_list(
        "ferroin", ["#b83a1e", "#e0623a", "#f2e6dc", "#5b8fd9", "#1f3f9e"]
    )


def _norm(vmax: float = 0.3):
    from matplotlib.colors import PowerNorm

    return PowerNorm(gamma=0.5, vmin=0.0, vmax=vmax)


def plot_pattern(z: np.ndarray, ax: plt.Axes | None = None, vmax: float = 0.3,
                 title: str | None = None):
    """Show one 2D catalyst field in ferroin colours."""
    if ax is None:
        _, ax = plt.subplots(figsize=(4, 4))
    im = ax.imshow(z, origin="lower", cmap=ferroin_cmap(), norm=_norm(vmax),
                   interpolation="bilinear")
    ax.set_xticks([])
    ax.set_yticks([])
    for side in ax.spines.values():
        side.set_visible(False)
    if title:
        ax.set_title(title, loc="left", fontsize=10, color=TEXT)
    return im


def save_gif(frames: np.ndarray, path, vmax: float = 0.3, frame_ms: int = 200,
             size: int = 300) -> None:
    """Write z frames as an animated GIF in ferroin colours.

    frame_ms is the display time of each frame in milliseconds.
    """
    from PIL import Image

    cmap, norm = ferroin_cmap(), _norm(vmax)
    images = []
    for fr in frames:
        rgb = (cmap(norm(np.flipud(fr)))[:, :, :3] * 255).astype(np.uint8)
        img = Image.fromarray(rgb).resize((size, size), Image.BILINEAR)
        images.append(img.convert("P", palette=Image.ADAPTIVE, colors=32))
    images[0].save(path, save_all=True, append_images=images[1:],
                   duration=frame_ms, loop=0, optimize=True)
