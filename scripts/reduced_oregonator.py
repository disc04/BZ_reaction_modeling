"""Reduced two-variable Oregonator: nullclines and Hopf bifurcations.

Produces:
  figures/reduced_nullclines.png    phase plane for f below / inside / above the window
  figures/reduced_bifurcation.png   bifurcation diagram in f with zooms on both edges

Usage (from the repository root):
    python scripts/reduced_oregonator.py
"""

from dataclasses import replace
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from bz_models import OregonatorParams, estimate_period, simulate
from bz_models.plotting import CYCLE, TEXT, plot_phase_plane, style_axes
from bz_models.reduced import (
    from_oregonator, hopf_points, jacobian_reduced, simulate_reduced,
    steady_state_reduced,
)

ROOT = Path(__file__).resolve().parents[1]  # repository root, for figure paths

FIG_DIR = ROOT / "figures"
R = from_oregonator(OregonatorParams())


def is_stable(r) -> bool:
    ss = steady_state_reduced(r)
    return bool(np.all(np.linalg.eigvals(jacobian_reduced(ss, r)).real < 0))


def attractor_range(r, kick: float, tau_end: float = 80.0) -> tuple[float, float]:
    """Min and max of x on the long-time attractor, starting from a kicked steady state."""
    ss = steady_state_reduced(r)
    tau, s = simulate_reduced(r, tau_end=tau_end, n_points=40001, s0=tuple(ss * (1 + kick)))
    x = s[0][tau > tau_end - 15]
    return x.min(), x.max()


def compare_with_full_model() -> None:
    t, c = simulate(OregonatorParams(), t_end=3000.0, n_points=60001)
    tau, s = simulate_reduced(R, tau_end=60.0, n_points=60001)
    p_full = estimate_period(t, c[0])
    p_red = estimate_period(tau, s[0]) / R.time_scale
    print(f"Scaled parameters: eps = {R.eps:.3g}, q = {R.q:.3g}")
    print(f"Period: full model {p_full:.0f} s, reduced model {p_red:.0f} s")


def nullcline_figure() -> None:
    cases = [(0.4, "f = 0.4: oxidised steady state"),
             (1.0, "f = 1.0: limit cycle"),
             (2.6, "f = 2.6: reduced steady state")]
    fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=True)
    for i, (ax, (f, title)) in enumerate(zip(axes, cases)):
        r = replace(R, f=f)
        tau, s = simulate_reduced(r, tau_end=40.0, n_points=40001, s0=(0.3, 0.3))
        traj = s[:, tau > 25.0] if f == 1.0 else s
        plot_phase_plane(r, ax=ax, trajectory=traj)
        ax.set_title(title, loc="left", fontsize=11, color=TEXT)
        if i > 0:
            ax.set_ylabel("")
        if i == 1:
            ax.legend(loc="lower right", frameon=False, fontsize=9)
    fig.tight_layout()
    out = FIG_DIR / "reduced_nullclines.png"
    fig.savefig(out, dpi=150)
    print(f"Saved {out.relative_to(ROOT)}")


def bifurcation_figure() -> None:
    f_lo, f_hi = hopf_points(R)
    print(f"Hopf points: f = {f_lo:.4f} and f = {f_hi:.4f}")

    def branches(f_grid):
        ss, stab, cmin, cmax = [], [], [], []
        for f in f_grid:
            r = replace(R, f=f)
            ss.append(steady_state_reduced(r)[0])
            stab.append(is_stable(r))
            lo, hi = attractor_range(r, kick=0.1)
            osc = hi / lo > 1.01
            cmin.append(lo if osc else np.nan)
            cmax.append(hi if osc else np.nan)
        return map(np.array, (ss, stab, cmin, cmax))

    def draw(ax, f_grid, ss, stab, cmin, cmax, label=True):
        ax.semilogy(np.where(stab, f_grid, np.nan), np.where(stab, ss, np.nan),
                    color=TEXT, lw=2, label="stable steady state" if label else None)
        ax.semilogy(np.where(~stab, f_grid, np.nan), np.where(~stab, ss, np.nan),
                    color=TEXT, lw=1.6, ls="--", label="unstable steady state" if label else None)
        ax.fill_between(f_grid, cmin, cmax, color=CYCLE, alpha=0.15, lw=0)
        ax.semilogy(f_grid, cmax, color=CYCLE, lw=2, label="limit cycle (max / min of x)" if label else None)
        ax.semilogy(f_grid, cmin, color=CYCLE, lw=2)
        for fh in (f_lo, f_hi):
            if f_grid[0] <= fh <= f_grid[-1]:
                xs = steady_state_reduced(replace(R, f=fh))[0]
                ax.plot(fh, xs, "D", ms=7, color=TEXT, zorder=5)
        ax.set_xlabel("stoichiometric factor f")
        ax.set_ylabel("x (scaled [HBrO$_2$])")
        style_axes(ax)

    fig = plt.figure(figsize=(13, 9))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.25, 1])

    # Full range
    ax = fig.add_subplot(gs[0, :])
    f_grid = np.unique(np.concatenate([np.linspace(0.3, 2.7, 121),
                                       np.linspace(0.50, 0.51, 11),
                                       np.linspace(2.38, 2.40, 11)]))
    draw(ax, f_grid, *branches(f_grid))
    ax.set_title("Bifurcation diagram of the reduced Oregonator (diamonds: Hopf points)",
                 loc="left", fontsize=11, color=TEXT)
    ax.legend(loc="center right", bbox_to_anchor=(0.80, 0.6), frameon=False, fontsize=9)
    ax.annotate(f"subcritical Hopf\nf = {f_lo:.3f}", xy=(f_lo, 0.5), xytext=(18, -110),
                textcoords="offset points", fontsize=9, color=TEXT,
                arrowprops=dict(arrowstyle="-", color=TEXT, lw=0.8, shrinkA=2, shrinkB=6))
    ax.annotate(f"supercritical Hopf\nf = {f_hi:.3f}", xy=(f_hi, 2e-4), xytext=(-12, 30),
                textcoords="offset points", fontsize=9, color=TEXT, ha="right")

    # Zoom: lower edge, bistability
    ax1 = fig.add_subplot(gs[1, 0])
    g1 = np.linspace(0.498, 0.512, 57)
    ss, stab, cmin, cmax = branches(g1)
    draw(ax1, g1, ss, stab, cmin, cmax, label=False)
    both = stab & np.isfinite(cmax)
    if both.any():
        ax1.axvspan(g1[both].min(), g1[both].max(), color="#eda100", alpha=0.12, lw=0)
        ax1.annotate("bistable:\nsteady state and\ncycle coexist",
                     xy=(g1[both].mean(), 3e-3), ha="center", fontsize=9, color=TEXT)
    ax1.set_title("Lower edge: oscillations appear abruptly", loc="left", fontsize=11, color=TEXT)

    # Zoom: upper edge, small cycle then canard explosion
    ax2 = fig.add_subplot(gs[1, 1])
    g2 = np.linspace(2.375, 2.405, 61)
    draw(ax2, g2, *branches(g2), label=False)
    ax2.annotate("canard explosion:\namplitude jumps\n~4 decades", xy=(2.3875, 3e-3),
                 ha="left", fontsize=9, color=TEXT, xytext=(8, 0), textcoords="offset points")
    ax2.set_title("Upper edge: small cycle, then explosion", loc="left", fontsize=11, color=TEXT)

    fig.tight_layout()
    out = FIG_DIR / "reduced_bifurcation.png"
    fig.savefig(out, dpi=150)
    print(f"Saved {out.relative_to(ROOT)}")


if __name__ == "__main__":
    FIG_DIR.mkdir(exist_ok=True)
    compare_with_full_model()
    nullcline_figure()
    bifurcation_figure()
