"""Simulate the Oregonator and plot BZ oscillations.

Usage (from the repository root):
    python scripts/run_oregonator.py
"""

from pathlib import Path

import matplotlib.pyplot as plt

from bz_models import OregonatorParams, estimate_period, plot_timeseries, simulate

ROOT = Path(__file__).resolve().parents[1]  # repository root, for figure paths


def main() -> None:
    params = OregonatorParams()
    t, c = simulate(params, t_end=1500.0)
    period = estimate_period(t, c[0])

    fig, ax = plt.subplots(figsize=(9, 4.5))
    plot_timeseries(t, c, ax=ax,
                    title=f"Oregonator, f = {params.f}  (period ≈ {period:.0f} s)")
    fig.tight_layout()

    out = ROOT / "figures" / "oregonator_timeseries.png"
    out.parent.mkdir(exist_ok=True)
    fig.savefig(out, dpi=150)
    print(f"Period ≈ {period:.1f} s. Saved {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
