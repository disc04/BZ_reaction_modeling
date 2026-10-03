"""Loaders for the Chang, de Planque & Zauner (2018) BZ droplet dataset.

Data: doi:10.5258/SOTON/D0363 (CC-BY). Download and unpack into
`Chang_BZ_data/` at the repository root (the folder is git-ignored).
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np

DATA_DOI = "https://doi.org/10.5258/SOTON/D0363"

# Base recipe of the isolated-droplet experiments (Chang et al. 2018, Methods)
BASE_RECIPE = {"h2so4": 0.5, "bromate": 0.47, "malonic": 0.18, "ferroin": 0.002}

# Figure 4 (isolated droplets): folder and file prefix per varied component
FIGURE4_SERIES = {
    "h2so4": ("1_H2SO4", "H2SO4"),
    "bromate": ("2_NaBrO3", "NaBrO3"),
    "malonic": ("3_Malonic_acid", "MA"),
    "ferroin": ("4_Ferroin", "Ferroin"),
}


@dataclass
class Series:
    """One concentration series: repeats per concentration plus their mean."""
    component: str
    quantity: str
    conc: np.ndarray            # (n,) concentration of the varied component, M
    repeats: list[np.ndarray]   # n arrays of repeat values (missing repeats dropped)
    mean: np.ndarray            # (n,) mean as reported in the file

    def pooled(self) -> tuple[np.ndarray, np.ndarray]:
        """All individual repeats as flat (concentration, value) arrays."""
        x = np.concatenate([[c] * len(r) for c, r in zip(self.conc, self.repeats)])
        return x, np.concatenate(self.repeats)


def data_dir(root: Path) -> Path:
    d = Path(root) / "Chang_BZ_data"
    if not d.exists():
        raise FileNotFoundError(
            f"{d} not found. Download the dataset from {DATA_DOI} and unpack it there."
        )
    return d


def _num(v: str) -> float | None:
    v = v.strip()
    if v in ("", "--"):
        return None
    return float(v)


def load_series(root: Path, component: str, quantity: str = "Frequency") -> Series:
    """Load a Figure 4 series, e.g. load_series(ROOT, "h2so4", "Frequency")."""
    folder, prefix = FIGURE4_SERIES[component]
    path = data_dir(root) / "Figure_4" / folder / f"{prefix}_{quantity}.csv"
    rows = list(csv.reader(open(path, newline="")))[2:]   # 2 header rows
    conc, repeats, mean = [], [], []
    for r in rows:
        if not r or not r[0].strip():
            continue
        conc.append(float(r[0]))
        repeats.append(np.array([v for v in map(_num, r[2:-2]) if v is not None]))
        mean.append(_num(r[-1]))
    return Series(component, quantity, np.array(conc), repeats, np.array(mean, dtype=float))


# --- Figure 3: full time series of one droplet (base recipe) -----------------

PX_PER_SECOND = 6.0   # space-time plot: 15 px = 2.5 s (one video frame at 0.4 fps)


@dataclass
class Peaks:
    """Curated peaks of the Figure 3 droplet (Origin peak analysis)."""
    t: np.ndarray            # peak centre, s after mixing
    height: np.ndarray       # peak height above baseline, a.u.
    fwhm: np.ndarray         # full width at half maximum, s
    left: np.ndarray         # left (rise) half width, s
    right: np.ndarray        # right (decay) half width, s


def load_timeseries(root: Path) -> tuple[np.ndarray, np.ndarray]:
    """Flattened blue-channel intensity of the Figure 3 droplet: (t in s, intensity a.u.).

    Only times after mixing are returned; filtered (missing) samples are NaN.
    """
    path = data_dir(root) / "Figure_3" / "4_Hidden_Graph.csv"
    rows = list(csv.reader(open(path, newline="")))[3:]   # 3 header rows
    px, val = [], []
    for r in rows:
        x = _num(r[0]) if r else None
        if x is None or x < 0:
            continue
        v = _num(r[8]) if len(r) > 8 else None   # column I: flattened, height-calibrated
        px.append(x)
        val.append(np.nan if v is None else v)
    return np.array(px) / PX_PER_SECOND, np.array(val)


def load_peaks(root: Path) -> Peaks:
    path = data_dir(root) / "Figure_3" / "6_Peak_Parameters.csv"
    rows = [r for r in list(csv.reader(open(path, newline="")))[2:] if r and r[0].strip()]
    col = {name: i for i, name in enumerate(
        ["Index", "Area", "AreaIntgP", "CurveArea", "Row", "Begin", "End", "FWHM",
         "Left", "Right", "Center", "Height", "Centroid"])}
    get = lambda name: np.array([float(r[col[name]]) for r in rows])  # noqa: E731
    s = PX_PER_SECOND
    return Peaks(t=get("Center") / s, height=get("Height"), fwhm=get("FWHM") / s,
                 left=get("Left") / s, right=get("Right") / s)
