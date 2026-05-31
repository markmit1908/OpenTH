"""Incompressible (liquid) fluid model.

For liquid flows the paper works with the total pressure p0 = p + rho V^2/2 + rho g z
and a constant density (eqs. 5-6, 10). This model provides that constant-density EOS.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from .base import FluidModel


@dataclass
class Incompressible(FluidModel):
    """Constant-density liquid with a finite bulk modulus for acoustic speed."""

    name: str
    rho: float          # density [kg/m^3]
    cp: float           # specific heat [J/(kg.K)]
    bulk_modulus: float = 2.2e9  # default ~ water [Pa]

    def density(self, p: float, T: float) -> float:
        return self.rho

    def enthalpy(self, T: float) -> float:
        return self.cp * T

    def temperature_from_enthalpy(self, h: float) -> float:
        return h / self.cp

    def sonic_velocity(self, p: float, T: float) -> float:
        # a = sqrt(K / rho), the liquid's acoustic (water-hammer) speed.
        return math.sqrt(self.bulk_modulus / self.rho)

    def drho_dp(self, p: float, T: float) -> float:
        # Constant density: a pressure change drives no density change in the storage term.
        # (A finite bulk modulus could be modelled as rho / bulk_modulus if liquid
        # compressibility / water-hammer storage is needed; left at 0 for now.)
        return 0.0


def water() -> Incompressible:
    return Incompressible(name="water", rho=998.0, cp=4182.0, bulk_modulus=2.2e9)
