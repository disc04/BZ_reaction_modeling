"""Oregonator parameters from a BZ recipe (concentrations of H2SO4, NaBrO3, MA).

Field-Noyes rate constants depend on acidity through the elementary steps:

    O1  BrO3- + Br-    + 2H+  ->  k1 = 2      * h^2   (M^-1 s^-1 after folding in h)
    O2  HBrO2 + Br-    +  H+  ->  k2 = 3.0e6  * h
    O3  BrO3- + HBrO2  +  H+  ->  k3 = 42     * h
    O4  2 HBrO2               ->  k4 = 3000           (no acid dependence)

(Field & Foersterling, 1986; at h = 0.8 M these give the classic k1, k2, k3.)
k_c (catalyst reduction by the organic substrate) has no first-principles value
and is calibrated against data.
"""

from __future__ import annotations

import numpy as np

from .oregonator import OregonatorParams

KA2_HSO4 = 0.0102  # second dissociation constant of H2SO4 at 25 C (M)


def hplus(h2so4: float | np.ndarray, ka2: float = KA2_HSO4):
    """[H+] from analytical H2SO4: first proton complete, second via Ka2."""
    c = np.asarray(h2so4, dtype=float)
    b = c + ka2
    x = 0.5 * (-b + np.sqrt(b * b + 4.0 * ka2 * c))
    return c + x


def params_from_recipe(
    h2so4: float,
    bromate: float,
    malonic: float,
    kc: float,
    f: float = 1.0,
) -> OregonatorParams:
    """OregonatorParams for given molar concentrations of H2SO4, NaBrO3 and malonic acid."""
    h = float(hplus(h2so4))
    return OregonatorParams(
        k1=2.0 * h * h, k2=3.0e6 * h, k3=42.0 * h, k4=3000.0,
        kc=kc, A=bromate, B=malonic, f=f,
    )
