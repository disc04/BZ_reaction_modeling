"""Effect of the stoichiometric factor f on Oregonator dynamics.

Produces two figures:
  figures/oregonator_f_timeseries.png  time series for selected f values
  figures/oregonator_f_period.png      oscillation period vs f

Usage (from the repository root):
    python scripts/explore_f.py
"""

from dataclasses import replace
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from bz_models import OregonatorParams, estimate_period, plot_timeseries, simulate
from bz_models.plotting import TEXT, style_axes

ROOT = Path(__file__).resolve().parents[1]  # repository root, for figure paths

BASE = OregonatorParams()
FIG_DIR = ROOT / "figures"


def timeseries_panels(f_values=(0.4, 0.6, 1.0, 2.2, 2.6), t_end=1500.0) -> None:
    fig, axes = plt.subplots(len(f_values), 1, figsize=(9, 2.3 * len(f_values)),
                             sharex=True)
    for i, (ax, f) in enumerate(zip(axes, f_values)):
        t, c = simulate(replace(BASE, f=f), t_end=t_end)
        period = estimate_period(t, c[0])
        state = f"period ≈ {period:.0f} s" if np.isfinite(period) else "steady state"
        plot_timeseries(t, c, ax=ax, title=f"f = {f}: {state}",
                        legend=(i == 0), direct_labels=False)
        if i < len(f_values) - 1:
            ax.set_xlabel("")
    fig.tight_layout()
    out = FIG_DIR / "oregonator_f_timeseries.png"
    fig.savefig(out, dpi=150)
    print(f"Saved {out.relative_to(ROOT)}")


def period_scan(f_min=0.3, f_max=2.7, step=0.025, t_end=4000.0) -> None:
    f_grid = np.round(np.arange(f_min, f_max + 1e-9, step), 4)
    periods = []
    for f in f_grid:
        t, c = simulate(replace(BASE, f=f), t_end=t_end, n_points=40001)
        periods.append(estimate_period(t, c[0]))
    periods = np.array(periods)

    osc = np.isfinite(periods)
    lo, hi = f_grid[osc].min(), f_grid[osc].max()

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.axvspan(lo, hi, color="#2a78d6", alpha=0.08, lw=0)
    ax.plot(f_grid[osc], periods[osc], color="#2a78d6", lw=2, marker="o", ms=4)
    ax.set_xlim(f_min, f_max)
    ax.set_xlabel("stoichiometric factor f")
    ax.set_ylabel("period (s)")
    ax.set_title("Oscillation period vs f (shaded: oscillatory window)",
                 loc="left", fontsize=11, color=TEXT)
    ymin = np.nanmin(periods)
    for xpos, text, ha in ((lo, "← oxidised\nsteady state", "right"),
                           (hi, "→ reduced\nsteady state", "left")):
        ax.annotate(text, xy=(xpos, ymin), xytext=(-6 if ha == "right" else 6, 0),
                    textcoords="offset points", ha=ha, va="bottom",
                    fontsize=9, color=TEXT)
    style_axes(ax)
    fig.tight_layout()
    out = FIG_DIR / "oregonator_f_period.png"
    fig.savefig(out, dpi=150)
    print(f"Oscillations for {lo:.3f} <= f <= {hi:.3f}. Saved {out.relative_to(ROOT)}")


if __name__ == "__main__":
    FIG_DIR.mkdir(exist_ok=True)
    timeseries_panels()
    period_scan()
