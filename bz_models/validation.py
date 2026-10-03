"""Shared helpers for Part 2: calibration of the Oregonator on droplet data."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from .data import BASE_RECIPE, load_series
from .oregonator import oscillation_frequency
from .recipe import params_from_recipe


def base_frequency(root: Path) -> float:
    """Measured frequency at the base recipe, pooled over the four Figure 4 series."""
    vals = []
    for comp in ("h2so4", "bromate", "malonic", "ferroin"):
        s = load_series(root, comp)
        vals.append(s.mean[int(np.argmin(abs(s.conc - BASE_RECIPE[comp])))])
    return float(np.mean(vals))


def base_frequency_model(kc: float, f: float = 1.0) -> float:
    r = BASE_RECIPE
    return oscillation_frequency(params_from_recipe(r["h2so4"], r["bromate"], r["malonic"], kc, f))


def calibrate_kc(target: float, bracket: tuple[float, float] = (20.0, 45.0)) -> float:
    """k_c (M^-1 s^-1) that reproduces `target` Hz at the base recipe (f = 1)."""
    return brentq(lambda kc: base_frequency_model(kc) - target, *bracket, xtol=0.01)
