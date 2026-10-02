"""Target and spiral waves in an excitable BZ medium (2D reaction-diffusion).

Produces:
  figures/bz_targets.png, figures/bz_targets.gif   pacemaker-driven target waves
  figures/bz_spiral.png,  figures/bz_spiral.gif    spiral from a broken wave front

Runtime: roughly 2-3 minutes on a laptop (pure NumPy).

Usage (from the repository root):
    python scripts/spatial_patterns.py
"""

from pathlib import Path

import matplotlib.pyplot as plt

from bz_models.plotting import TEXT, plot_pattern, save_gif
from bz_models.spatial import RDParams, run, spiral_setup, target_setup

ROOT = Path(__file__).resolve().parents[1]  # repository root, for figure paths

FIG_DIR = ROOT / "figures"
SNAPSHOT_TIMES = (6.0, 18.0, 36.0, 72.0)
FRAME_DT = 0.25          # simulated time between GIF frames
PLAYBACK = 1.2           # simulated time units shown per second of animation
FRAME_MS = int(round(1000 * FRAME_DT / PLAYBACK))


def snapshot_figure(times, frames, title, out):
    fig, axes = plt.subplots(1, len(SNAPSHOT_TIMES), figsize=(14, 3.9), layout="constrained")
    for ax, ts in zip(axes, SNAPSHOT_TIMES):
        i = int(abs(times - ts).argmin())
        im = plot_pattern(frames[i], ax=ax, title=f"t = {times[i]:.0f}")
    cbar = fig.colorbar(im, ax=axes, fraction=0.012, pad=0.01, shrink=0.9)
    cbar.set_ticks([0.0, 0.3])
    cbar.set_ticklabels(["reduced\n(ferroin, red)", "oxidised\n(ferriin, blue)"])
    cbar.ax.tick_params(labelsize=8, colors=TEXT)
    cbar.outline.set_visible(False)
    fig.suptitle(title, x=0.01, ha="left", fontsize=12, color=TEXT)
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"Saved {out.relative_to(ROOT)}")


def main() -> None:
    FIG_DIR.mkdir(exist_ok=True)
    p = RDParams()
    t_end = max(SNAPSHOT_TIMES)

    # Targets: excitable medium with two oscillatory pacemaker spots.
    x, z, f = target_setup(p)
    times, frames = run(x, z, p, t_end=t_end, frame_every=FRAME_DT, f=f)
    snapshot_figure(times, frames, "Target waves from two pacemakers: rings annihilate where they meet",
                    FIG_DIR / "bz_targets.png")
    save_gif(frames, FIG_DIR / "bz_targets.gif", frame_ms=FRAME_MS)
    print("Saved figures/bz_targets.gif")

    # Spiral: a broken wave front curls around its free end.
    x, z = spiral_setup(p)
    times, frames = run(x, z, p, t_end=t_end, frame_every=FRAME_DT)
    snapshot_figure(times, frames, "Spiral wave from a broken wave front",
                    FIG_DIR / "bz_spiral.png")
    save_gif(frames, FIG_DIR / "bz_spiral.gif", frame_ms=FRAME_MS)
    print("Saved figures/bz_spiral.gif")


if __name__ == "__main__":
    main()
