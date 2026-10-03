"""Part 2, test 1: oscillation frequency vs [H2SO4] and [NaBrO3].

Compares the Oregonator with the droplet data of Chang et al. (2018).
One free parameter (k_c) is calibrated on the base recipe; everything else is
a prediction. Two hypotheses for k_c are tested:

  A  k_c constant                 (catalyst reduction independent of recipe)
  B  k_c proportional to [H+][BrO3-]  (catalyst reduction speeds up with acid
                                   and bromate, keeping eps = k_c B / k3 A fixed)

Produces figures/validation_acid_bromate.png.

Usage (from the repository root; needs Chang_BZ_data/, see README):
    python scripts/validate_acid_bromate.py
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from bz_models.data import BASE_RECIPE, load_series
from bz_models.oregonator import oscillation_frequency
from bz_models.plotting import TEXT, TRANSIENT, style_axes
from bz_models.recipe import hplus, params_from_recipe
from bz_models.validation import base_frequency, calibrate_kc

ROOT = Path(__file__).resolve().parents[1]  # repository root, for data and figure paths

COLOR_A = "#eb6834"   # model A
COLOR_B = "#2a78d6"   # model B
DATA = "#3d3d3a"

H0 = float(hplus(BASE_RECIPE["h2so4"]))
A0 = BASE_RECIPE["bromate"]


def model_frequency(kc_base: float, hypothesis: str, h2so4: float, bromate: float) -> float:
    kc = kc_base
    if hypothesis == "B":
        kc = kc_base * float(hplus(h2so4)) * bromate / (H0 * A0)
    p = params_from_recipe(h2so4, bromate, BASE_RECIPE["malonic"], kc)
    return oscillation_frequency(p)


def loglog_slope(x, y) -> float:
    ok = y > 0
    return float(np.polyfit(np.log(x[ok]), np.log(y[ok]), 1)[0])


def main() -> None:
    series = {c: load_series(ROOT, c) for c in ("h2so4", "bromate", "malonic", "ferroin")}

    target = base_frequency(ROOT)
    kc_base = calibrate_kc(target)
    print(f"Base-recipe frequency (data): {target:.4f} Hz  ->  calibrated k_c = {kc_base:.1f} M^-1 s^-1")

    panels = [
        ("h2so4", "[H$_2$SO$_4$] (M)", lambda c: (c, A0)),
        ("bromate", "[NaBrO$_3$] (M)", lambda c: (BASE_RECIPE["h2so4"], c)),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    for ax, (comp, xlabel, recipe) in zip(axes, panels):
        s = series[comp]
        grid = np.geomspace(s.conc.min() * 0.95, s.conc.max() * 1.05, 25)
        fa = np.array([model_frequency(kc_base, "A", *recipe(c)) for c in grid])
        fb = np.array([model_frequency(kc_base, "B", *recipe(c)) for c in grid])

        # data
        x, y = s.pooled()
        ax.plot(x, y, "o", ms=4, color=TRANSIENT, alpha=0.7, label="droplets (repeats)")
        ax.plot(s.conc, s.mean, "o", ms=8, color=DATA, label="data mean")

        # model A: show the no-oscillation range explicitly
        osc = fa > 0
        ax.plot(grid[osc], fa[osc], color=COLOR_A, lw=2.2, label="model A: $k_c$ constant")
        if (~osc).any():
            ax.axvspan(grid[~osc].min(), grid[~osc].max(), color=COLOR_A, alpha=0.08, lw=0)
            ax.text(np.sqrt(grid[~osc].min() * grid[~osc].max()), 4.6e-3,
                    "model A:\nno oscillations", ha="center", va="bottom", fontsize=9, color=TEXT)
        ax.plot(grid, fb, color=COLOR_B, lw=2.2,
                label=r"model B: $k_c \propto [\mathrm{H^+}][\mathrm{BrO_3^-}]$")

        sd, sb = loglog_slope(x, y), loglog_slope(grid, fb)
        ax.text(0.97, 0.05, f"log-log slope\ndata {sd:.2f}\nmodel B {sb:.2f}",
                transform=ax.transAxes, ha="right", va="bottom", fontsize=9, color=TEXT)
        print(f"{comp:8s} slope: data {sd:.2f}, model B {sb:.2f}; "
              f"model A oscillates only above {grid[osc].min():.2f} M")

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_ylim(4e-3, 0.3)
        ax.set_xlabel(xlabel)
        ax.set_title(f"Frequency vs {xlabel.split(' ')[0]}", loc="left", fontsize=11, color=TEXT)
        style_axes(ax)
    axes[0].set_ylabel("oscillation frequency (Hz)")
    axes[0].legend(loc="upper left", frameon=False, fontsize=9)
    fig.suptitle("Test 1: acid and bromate (droplet data: Chang et al. 2018)",
                 x=0.01, ha="left", fontsize=12, color=TEXT)
    fig.tight_layout()
    out = ROOT / "figures" / "validation_acid_bromate.png"
    fig.savefig(out, dpi=150)
    print(f"Saved {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
