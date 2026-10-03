"""Part 2, test 2: oscillation frequency vs [malonic acid], and the role of f.

In the Oregonator, malonic acid (B) enters only through the catalyst-reset
rate k_c*B (step O5). The test asks:
  (a) What does the calibrated model predict as [MA] varies?
        A  k_c constant         -> k_c*B grows with [MA]
        B  k_c*B constant       -> MA has no effect
  (b) What k_c*B would each measured frequency require (f = 1)?
  (c) Could the bromide yield f explain the data instead? A frequency map over
      (k_c*B, f) shows how much each parameter moves the frequency.

Produces figures/validation_malonic_acid.png.

Usage (from the repository root; needs Chang_BZ_data/, see README):
    python scripts/validate_malonic_acid.py        (~3-5 min, mostly the map)
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from scipy.optimize import brentq

from bz_models.data import BASE_RECIPE, load_series
from bz_models.oregonator import oscillation_frequency
from bz_models.plotting import TEXT, TRANSIENT, style_axes
from bz_models.recipe import params_from_recipe
from bz_models.validation import base_frequency, calibrate_kc

ROOT = Path(__file__).resolve().parents[1]  # repository root, for data and figure paths

COLOR_A = "#eb6834"
COLOR_B = "#2a78d6"
DATA = "#3d3d3a"
SEQ = LinearSegmentedColormap.from_list("blues", ["#cde2fb", "#6da7ec", "#2a78d6", "#184f95", "#0d366b"])

H, A, B0 = BASE_RECIPE["h2so4"], BASE_RECIPE["bromate"], BASE_RECIPE["malonic"]


def freq_kcB(kcB: float, f: float = 1.0, **kw) -> float:
    """Model frequency at base acid/bromate for a given reset rate k_c*B (1/s)."""
    return oscillation_frequency(params_from_recipe(H, A, B0, kcB / B0, f), **kw)


def required_kcB(target: float, hi: float) -> float:
    """k_c*B (f = 1) that reproduces a measured frequency; frequency rises with k_c*B."""
    return brentq(lambda k: freq_kcB(k) - target, 0.02, hi, xtol=1e-3)


def main() -> None:
    s = load_series(ROOT, "malonic")
    target = base_frequency(ROOT)
    kc = calibrate_kc(target)
    kcB_base = kc * B0
    print(f"Calibrated k_c = {kc:.1f} M^-1 s^-1 (k_c*B = {kcB_base:.2f} s^-1 at base recipe)")

    fig, axes = plt.subplots(1, 3, figsize=(17, 5.2), gridspec_kw={"width_ratios": [1, 1, 1.15]})

    # --- (a) predictions vs data ---------------------------------------------
    ax = axes[0]
    grid = np.geomspace(0.15, 0.6, 30)
    fa = np.array([oscillation_frequency(params_from_recipe(H, A, b, kc)) for b in grid])
    x, y = s.pooled()
    ax.plot(x, y, "o", ms=4, color=TRANSIENT, alpha=0.7, label="droplets (repeats)")
    ax.plot(s.conc, s.mean, "o", ms=8, color=DATA, label="data mean")
    osc = fa > 0
    ax.plot(grid[osc], fa[osc], color=COLOR_A, lw=2.2, label="model A: $k_c$ constant")
    ax.axvspan(grid[~osc].min(), grid[~osc].max(), color=COLOR_A, alpha=0.08, lw=0)
    ax.text(np.sqrt(grid[~osc].min() * grid[~osc].max()), 0.21, "model A:\nno oscillations",
            ha="center", va="top", fontsize=9, color=TEXT)
    ax.axhline(target, color=COLOR_B, lw=2.2, ls="--", label="model B: $k_c B$ constant")
    print(f"Model A stops oscillating above [MA] = {grid[osc].max():.2f} M")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim(0.02, 0.25)
    ax.set_xlabel("[malonic acid] (M)")
    ax.set_ylabel("oscillation frequency (Hz)")
    ax.set_title("(a) Predictions vs data", loc="left", fontsize=11, color=TEXT)
    ax.legend(loc="lower left", frameon=False, fontsize=9)
    style_axes(ax)

    # --- (b) k_c*B the data would require -------------------------------------
    ax = axes[1]
    hi = 9.5  # just below the upper Hopf point in k_c*B
    req = np.array([required_kcB(m, hi) for m in s.mean])
    lo_err = req - np.array([required_kcB(r.min(), hi) for r in s.repeats])
    hi_err = np.array([required_kcB(min(r.max(), 0.17), hi) for r in s.repeats]) - req
    ax.errorbar(s.conc, req, yerr=[lo_err, hi_err], fmt="o", ms=8, color=DATA, ecolor=TRANSIENT,
                elinewidth=2, capsize=0, label="required by data (range: repeats)")
    ax.plot(grid, kc * grid, color=COLOR_A, lw=2.2, label="model A: $k_c B$ with $k_c$ constant")
    slope = np.polyfit(np.log(s.conc), np.log(req), 1)[0]
    print("Required k_c*B (1/s):", np.round(req, 2), f"log-log slope {slope:.2f}")
    ax.text(0.03, 0.95, f"required: log-log slope {slope:.1f}\n(~{req[0] / req[-1]:.0f}x drop)",
            transform=ax.transAxes, ha="left", va="top", fontsize=9, color=TEXT)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("[malonic acid] (M)")
    ax.set_ylabel(r"catalyst-reset rate $k_c B$ (s$^{-1}$)")
    ax.set_title("(b) Reset rate needed to match the data", loc="left", fontsize=11, color=TEXT)
    ax.legend(loc="lower left", frameon=False, fontsize=9)
    style_axes(ax)

    # --- (c) frequency map over (k_c*B, f) --------------------------------------
    ax = axes[2]
    kcb_grid = np.geomspace(0.2, 14, 26)
    f_grid = np.linspace(0.4, 2.6, 23)
    F = np.array([[freq_kcB(k, f, t_end=400.0, max_t_end=1600.0) for k in kcb_grid] for f in f_grid])
    Fm = np.ma.masked_less_equal(F, 0)
    mesh = ax.pcolormesh(kcb_grid, f_grid, Fm, cmap=SEQ, shading="nearest")
    ax.set_facecolor("#f4f3ee")
    cs = ax.contour(kcb_grid, f_grid, Fm, levels=[0.04, 0.06, 0.08], colors="white", linewidths=1.2)
    ax.clabel(cs, fmt=lambda v: f"{v:.2f} Hz", fontsize=8)
    ax.plot(req, np.ones_like(req), "o", ms=7, mfc="white", mec=DATA, mew=1.5,
            label="required points (f = 1)")
    ax.plot(kcB_base, 1.0, "D", ms=8, color=DATA, label="calibrated base recipe")
    ax.text(0.97, 0.95, "grey: steady state\n(no oscillations)", transform=ax.transAxes,
            ha="right", va="top", fontsize=9, color=TEXT)
    ax.set_xscale("log")
    ax.set_xlabel(r"catalyst-reset rate $k_c B$ (s$^{-1}$)")
    ax.set_ylabel("bromide yield f")
    ax.set_title("(c) Model frequency over ($k_c B$, f)", loc="left", fontsize=11, color=TEXT)
    ax.legend(loc="lower left", frameon=False, fontsize=9)
    cb = fig.colorbar(mesh, ax=ax, fraction=0.05, pad=0.02)
    cb.set_label("frequency (Hz)", color=TEXT)
    cb.outline.set_visible(False)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.tick_params(colors=TEXT, labelsize=9)

    fig.suptitle("Test 2: malonic acid (droplet data: Chang et al. 2018)",
                 x=0.01, ha="left", fontsize=12, color=TEXT)
    fig.tight_layout()
    out = ROOT / "figures" / "validation_malonic_acid.png"
    fig.savefig(out, dpi=150)
    print(f"Saved {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
