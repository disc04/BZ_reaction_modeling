"""Part 2, test 3: full time series of one droplet vs the depleting Oregonator.

The droplet (base recipe, Chang et al. 2018, Figure 3) oscillates for ~70 min
in three phases: initial (small, fast), main (large, slowly lengthening period)
and late (shorter period, collapsing amplitude) before it is exhausted.

The model is the Oregonator with slow consumption of bromate and malonic acid
(bz_models.depletion). Two parameters are calibrated:
    k_c0  so that the model period at 10-20 min equals the measured one
    s     (depletion scale) so that oscillations end when the droplet's do
Everything else (shape of the period drift, amplitude, phases, waveform) is a
test.

Produces figures/validation_timeseries.png.

Usage (from the repository root; needs Chang_BZ_data/, see README):
    python scripts/validate_timeseries.py        (~3 min)
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq

from bz_models.data import load_peaks, load_timeseries
from bz_models.depletion import oscillation_peaks, simulate_depleting
from bz_models.plotting import TEXT, TRANSIENT, style_axes

ROOT = Path(__file__).resolve().parents[1]  # repository root, for data and figure paths

DATA = "#3d3d3a"
MODEL = "#2a78d6"
PHASE = "#eda100"
REF_WINDOW = (600.0, 1200.0)   # 10-20 min: start of the main phase


def median_period(tp: np.ndarray, window) -> float:
    per = np.diff(tp)
    m = (tp[:-1] >= window[0]) & (tp[:-1] < window[1])
    return float(np.median(per[m])) if m.sum() >= 3 else np.nan


def model_stats(kc0: float, s: float, t_end: float = 6000.0):
    t, u = simulate_depleting(kc0, s, t_end=t_end, dt=0.1)
    tp, zmax = oscillation_peaks(t, u)
    return tp, zmax


def calibrate(ref_period: float, lifetime: float) -> tuple[float, float]:
    kc0, s = 11.0, 0.04
    for _ in range(2):   # alternate: k_c0 sets the period, s sets the lifetime
        kc0 = brentq(lambda k: median_period(model_stats(k, s)[0], REF_WINDOW) - ref_period,
                     5.0, 25.0, xtol=0.05)
        s = brentq(lambda v: model_stats(kc0, v)[0][-1] - lifetime, 0.02, 0.08, xtol=1e-4)
    return kc0, s


def phases(pk) -> tuple[float, float]:
    """Data-driven phase boundaries: amplitude above half the main-phase median."""
    big = pk.height > 0.5 * np.median(pk.height)
    return float(pk.t[big][0]), float(pk.t[big][-1])


def rolling_median(x, y, width):
    out = np.full_like(y, np.nan)
    for i, xi in enumerate(x):
        m = np.abs(x - xi) <= width / 2
        out[i] = np.median(y[m])
    return out


def main() -> None:
    t_d, I_d = load_timeseries(ROOT)
    pk = load_peaks(ROOT)
    t_main, t_late = phases(pk)
    ref_period = median_period(pk.t, REF_WINDOW)
    lifetime = float(pk.t[-1])
    print(f"Data: main phase {t_main / 60:.1f}-{t_late / 60:.1f} min, lifetime {lifetime / 60:.1f} min, "
          f"period at 10-20 min {ref_period:.1f} s")

    kc0, s = calibrate(ref_period, lifetime)
    print(f"Calibrated: k_c0 = {kc0:.1f} M^-1 s^-1, depletion scale s = {s:.4f}")
    t_m, u_m = simulate_depleting(kc0, s, t_end=lifetime + 600.0, dt=0.05)
    tp_m, zmax_m = oscillation_peaks(t_m, u_m)
    # After malonic acid runs out nothing reduces the catalyst, and the
    # non-conserving model lets ferriin grow without bound. Keep only real cycles.
    zmax_m = zmax_m[:-1]   # aligned with tp_m[:-1]
    t_stop = tp_m[-1] + 20.0
    print(f"Model: lifetime {tp_m[-1] / 60:.1f} min; at the end A = {u_m[3, -1]:.3f} M, B = {u_m[4, -1]:.4f} M")

    fig = plt.figure(figsize=(15, 10.5))
    gs = fig.add_gridspec(3, 3, height_ratios=[1, 1, 1.35], hspace=0.45, wspace=0.28, top=0.93)
    mins = lambda v: np.asarray(v) / 60.0  # noqa: E731
    xlim = (0, (lifetime + 300) / 60)

    # --- full traces ---------------------------------------------------------
    ax = fig.add_subplot(gs[0, :])
    ax.plot(mins(t_d), I_d, color=DATA, lw=0.5)
    for a, b, name in ((0, t_main, "initial"), (t_main, t_late, "main"), (t_late, lifetime, "late")):
        if name != "main":
            ax.axvspan(a / 60, b / 60, color=PHASE, alpha=0.12, lw=0)
        ax.text((a + b) / 120, 0.97, name, transform=ax.get_xaxis_transform(),
                ha="center", va="top", fontsize=9, color=TEXT)
    ax.set_xlim(*xlim)
    ax.set_ylabel("intensity (a.u.)")
    ax.set_title("Droplet: blue-channel intensity (ferriin), base recipe", loc="left",
                 fontsize=11, color=TEXT)
    style_axes(ax)

    ax = fig.add_subplot(gs[1, :], sharex=ax)
    keep = t_m <= t_stop
    ax.plot(mins(t_m[keep]), u_m[2, keep] * 1e3, color=MODEL, lw=0.5)
    ax.axvline(t_stop / 60, color=TEXT, lw=1, ls=":")
    ax.text(t_stop / 60, 0.97, " malonic acid\n exhausted", transform=ax.get_xaxis_transform(),
            ha="left", va="top", fontsize=9, color=TEXT)
    ax.set_ylim(0, 1.15 * zmax_m.max() * 1e3)
    ax.set_xlim(*xlim)
    ax.set_ylabel("[ferriin] (mM)")
    ax.set_xlabel("time after mixing (min)")
    ax.set_title(f"Model: Oregonator with reactant depletion (k$_{{c0}}$ = {kc0:.1f}, s = {s:.3f})",
                 loc="left", fontsize=11, color=TEXT)
    style_axes(ax)

    # --- period vs time --------------------------------------------------------
    ax = fig.add_subplot(gs[2, 0])
    per_d = np.diff(pk.t)
    ax.plot(mins(pk.t[:-1]), per_d, "o", ms=3, color=TRANSIENT, alpha=0.6, label="droplet (each cycle)")
    ax.plot(mins(pk.t[:-1]), rolling_median(pk.t[:-1], per_d, 300.0), color=DATA, lw=2,
            label="droplet (5-min median)")
    ax.plot(mins(tp_m[:-1]), np.diff(tp_m), color=MODEL, lw=2, label="model")
    ax.axvspan(*mins(REF_WINDOW), color=MODEL, alpha=0.08, lw=0)
    ax.set_ylim(0, 45)
    ax.set_xlim(*xlim)
    ax.set_xlabel("time after mixing (min)")
    ax.set_ylabel("period (s)")
    ax.set_title("Period", loc="left", fontsize=11, color=TEXT)
    ax.legend(loc="upper left", frameon=False, fontsize=8.5)
    style_axes(ax)

    # --- amplitude vs time -----------------------------------------------------
    ax = fig.add_subplot(gs[2, 1])
    main_d = (pk.t > t_main) & (pk.t < t_late)
    amp_d = pk.height / np.median(pk.height[main_d])
    main_m = (tp_m[:-1] > REF_WINDOW[0]) & (tp_m[:-1] < t_late)
    amp_m = zmax_m / np.median(zmax_m[main_m])
    ax.plot(mins(pk.t), amp_d, "o", ms=3, color=TRANSIENT, alpha=0.6)
    ax.plot(mins(pk.t), rolling_median(pk.t, amp_d, 300.0), color=DATA, lw=2, label="droplet")
    ax.plot(mins(tp_m[:-1]), amp_m, color=MODEL, lw=2, label="model")
    ax.set_xlim(*xlim)
    ax.set_xlabel("time after mixing (min)")
    ax.set_ylabel("amplitude / main-phase median")
    ax.set_title("Amplitude", loc="left", fontsize=11, color=TEXT)
    ax.legend(loc="upper left", frameon=False, fontsize=8.5)
    style_axes(ax)

    # --- waveform zoom -----------------------------------------------------------
    ax = fig.add_subplot(gs[2, 2])
    w0, w1 = 30 * 60.0, 30 * 60.0 + 75.0
    md = (t_d >= w0) & (t_d <= w1)
    yd = I_d[md]
    yd = (yd - np.nanmin(yd)) / (np.nanmax(yd) - np.nanmin(yd))
    # align the first model peak in the window with the first data peak
    pd0 = pk.t[pk.t >= w0][0]
    pm0 = tp_m[tp_m >= w0][0]
    shift = pm0 - pd0
    mm = (t_m >= w0 + shift) & (t_m <= w1 + shift)
    zm = u_m[2, mm]
    zm = (zm - zm.min()) / (zm.max() - zm.min())
    ax.plot(t_d[md] - w0, yd, color=DATA, lw=2, marker="o", ms=4, label="droplet (0.4 frames/s)")
    ax.plot(t_m[mm] - shift - w0, zm, color=MODEL, lw=1.8, label="model")
    ax.set_xlabel("time from 30 min (s)")
    ax.set_ylabel("normalised signal")
    ax.set_ylim(-0.05, 1.35)
    ax.set_title("Waveform (main phase)", loc="left", fontsize=11, color=TEXT)
    ax.legend(loc="upper right", frameon=False, fontsize=8.5, ncol=2)
    style_axes(ax)

    fig.suptitle("Test 3: time series of one droplet (Chang et al. 2018, Fig. 3)",
                 x=0.01, ha="left", fontsize=12, color=TEXT)
    out = ROOT / "figures" / "validation_timeseries.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"Saved {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
