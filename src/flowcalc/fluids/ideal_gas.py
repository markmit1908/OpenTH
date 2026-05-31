"""Ideal / real-gas equation of state: p = s * rho * R * T  (Greyvenstein eq. 12).

This is the variant fully described in the paper; the worked examples use helium.
The compressibility factor ``s`` defaults to 1 (ideal gas) but is exposed so a real-gas
correlation can be substituted later.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from .base import FluidModel


@dataclass
class IdealGas(FluidModel):
    """Calorically-perfect gas.

    Parameters
    ----------
    name : str
        Fluid label (e.g. ``"helium"``).
    R : float
        Specific gas constant [J/(kg.K)].
    gamma : float
        Ratio of specific heats cp/cv.
    s : float
        Compressibility factor (1.0 for an ideal gas).
    """

    name: str
    R: float
    gamma: float
    s: float = 1.0

    @property
    def cp(self) -> float:
        """Constant-pressure specific heat from R and gamma."""
        return self.gamma * self.R / (self.gamma - 1.0)

    @property
    def cv(self) -> float:
        """Constant-volume specific heat."""
        return self.R / (self.gamma - 1.0)

    def density(self, p: float, T: float) -> float:
        # eq. (12): p = s * rho * R * T  ->  rho = p / (s R T)
        return p / (self.s * self.R * T)

    def enthalpy(self, T: float) -> float:
        return self.cp * T

    def temperature_from_enthalpy(self, h: float) -> float:
        return h / self.cp

    def sonic_velocity(self, p: float, T: float) -> float:
        return math.sqrt(self.gamma * self.s * self.R * T)

    def drho_dp(self, p: float, T: float) -> float:
        # rho = p / (s R T)  ->  d rho / d p = 1 / (s R T)
        return 1.0 / (self.s * self.R * T)


def helium() -> IdealGas:
    """Helium, the working fluid in the paper's benchmark cases."""
    return IdealGas(name="helium", R=2077.0, gamma=1.667)


def air() -> IdealGas:
    return IdealGas(name="air", R=287.0, gamma=1.4)
