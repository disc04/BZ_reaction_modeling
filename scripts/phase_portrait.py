"""Phase portrait of the Oregonator: the limit cycle in (x, y, z) space.

Shows that trajectories from different initial conditions, including a
small perturbation of the (unstable) steady state, all converge to the
same closed orbit: a stable limit cycle.

Usage (from the repository root):
    python scripts/phase_portrait.py
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from bz_models import OregonatorParams, estimate_period, simulate
from bz_models.oregonator import jacobian, steady_state
from bz_models.plotting import (
    CYCLE, TEXT, TRANSIENT, plot_phase_2d, plot_phase_3d,
)

ROOT = Path(__file__).resolve().parents[1]  # repository root, for figure paths


def main() -> None:
    p = OregonatorParams()
    ss = steady_state(p)
    eig = np.linalg.eigvals(jacobian(ss, p))
    print("Steady state [x, y, z] (M):", ss)
    print("Jacobian eigenvalues (1/s):", eig)

    # Limit cycle: long run, keep only the settled part.
    t, c = simulate(p, t_end=3000.0, n_points=60001)
    period = estimate_period(t, c[0])
    cyc = c[:, t >= 3000.0 - 1.05 * period]

    # Transients from several starting points.
    starts = [
        (1e-10, 1e-6, 1e-6),
        (1e-5, 1e-3, 1e-2),
        (1e-9, 1e-4, 1e-5),
        tuple(ss * 1.001),  # tiny kick away from the steady state
    ]
    transients = [simulate(p, t_end=600.0, n_points=12001, c0=s)[1] for s in starts]

    fig = plt.figure(figsize=(13, 5.6))

    # --- 3D portrait ---------------------------------------------------------
    ax3 = fig.add_subplot(1, 2, 1, projection="3d")
    for i, tr in enumerate(transients):
        plot_phase_3d(tr, ax=ax3, color=TRANSIENT, lw=0.9, alpha=0.8,
                      label="transients from 4 initial conditions" if i == 0 else None)
    plot_phase_3d(cyc, ax=ax3, color=CYCLE, lw=2.4, label="limit cycle")
    ax3.scatter(*np.log10(ss), s=50, facecolor="white", edgecolor=TEXT,
                linewidth=1.5, depthshade=False, label="unstable steady state")
    for s in starts[:3]:
        ax3.scatter(*np.log10(s), s=18, color=TRANSIENT, depthshade=False)
    ax3.view_init(elev=22, azim=-58)
    ax3.legend(loc="upper left", frameon=False, fontsize=9)
    ax3.set_title("Phase space (log concentrations)", loc="left", fontsize=11, color=TEXT)

    # --- 2D projection: Br- vs oxidised catalyst -----------------------------
    ax2 = fig.add_subplot(1, 2, 2)
    plot_phase_2d(cyc, ax=ax2, xvar="z", yvar="y", lw=2.2)
    ax2.plot(ss[2], ss[1], "o", ms=7, mfc="white", mec=TEXT, mew=1.5)
    ax2.annotate("unstable\nsteady state", xy=(ss[2], ss[1]), xytext=(10, -4),
                 textcoords="offset points", fontsize=9, color=TEXT, va="top")

    # Label the FKN processes on the cycle; segments are located from the data.
    y, z = cyc[1], cyc[2]
    gy = np.sqrt(y.min() * y.max())  # geometric middle of the Br- range
    gz = np.sqrt(z.min() * z.max())  # geometric middle of the M_ox range
    low_z, low_y, high_z = z < 3 * z.min(), y < 3 * y.min(), z > 0.5 * z.max()
    upper = y > gy

    def closest(mask, target_z=None, target_y=None):
        cand = np.where(mask)[0]
        d = 0.0
        if target_z is not None:
            d = d + np.abs(np.log10(z[cand] / target_z))
        if target_y is not None:
            d = d + np.abs(np.log10(y[cand] / target_y))
        return int(cand[np.argmin(d)])

    labels = [
        (closest(low_z, target_y=1e-7), "A: Br⁻ consumed\nbelow threshold", (10, 0), "left"),
        (closest(low_y, target_z=gz), "B: autocatalysis,\nM$_\\mathrm{red}$ → M$_\\mathrm{ox}$", (0, 22), "center"),
        (closest(high_z, target_y=gy), "C: catalyst reduction\nreleases Br⁻", (-10, 0), "right"),
        (closest(upper, target_z=gz), "C: slow decay of\nM$_\\mathrm{ox}$ and Br⁻", (8, -16), "left"),
    ]
    for i, text, off, ha in labels:
        ax2.plot(z[i], y[i], "o", ms=5, color=CYCLE)
        ax2.annotate(text, xy=(z[i], y[i]), xytext=off, textcoords="offset points",
                     fontsize=9, color=TEXT, ha=ha, va="center")

    # Direction arrows along the cycle (shrink=0 so short arrows keep their heading).
    n = cyc.shape[1]
    for k in (int(0.25 * n), int(0.6 * n)):
        ax2.annotate("", xy=(z[k + 40], y[k + 40]), xytext=(z[k], y[k]),
                     arrowprops=dict(arrowstyle="-|>", color=CYCLE, lw=1.5,
                                     mutation_scale=16, shrinkA=0, shrinkB=0))
    ax2.set_title(f"Projection: Br⁻ vs M$_\\mathrm{{ox}}$ (one cycle, period ≈ {period:.0f} s)",
                  loc="left", fontsize=11, color=TEXT)

    fig.tight_layout()
    out = ROOT / "figures" / "oregonator_phase_portrait.png"
    fig.savefig(out, dpi=150)
    print(f"Saved {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
