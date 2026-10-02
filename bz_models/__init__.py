"""Kinetic models of the Belousov-Zhabotinsky (BZ) reaction."""

from .oregonator import (
    OregonatorParams,
    estimate_period,
    jacobian,
    oregonator_rhs,
    simulate,
    steady_state,
)
from .plotting import plot_phase_2d, plot_phase_3d, plot_phase_plane, plot_timeseries
from .reduced import (
    ReducedParams,
    from_oregonator,
    hopf_points,
    simulate_reduced,
    steady_state_reduced,
)

__all__ = [
    "OregonatorParams",
    "ReducedParams",
    "from_oregonator",
    "hopf_points",
    "simulate_reduced",
    "steady_state_reduced",
    "plot_phase_plane",
    "estimate_period",
    "jacobian",
    "oregonator_rhs",
    "plot_phase_2d",
    "plot_phase_3d",
    "plot_timeseries",
    "simulate",
    "steady_state",
]
