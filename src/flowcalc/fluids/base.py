"""Fluid property models (equations of state).

The solver needs, at minimum, a relationship between pressure, density and temperature
plus thermodynamic helpers (enthalpy, sonic velocity, specific heats). Greyvenstein's
gas-flow variant uses the equation of state ``p = s * rho * R * T`` (eq. 12), where ``s``
is the compressibility factor. Liquid (incompressible) flows use a constant density.

Every concrete fluid implements :class:`FluidModel` so the solver stays agnostic to the
particular working fluid (helium, air, water, ...).
"""

from __future__ import annotations

from abc import ABC, abstractmethod


class FluidModel(ABC):
    """Abstract equation of state + thermophysical properties for a working fluid.

    All quantities are SI: pressure [Pa], density [kg/m^3], temperature [K],
    enthalpy [J/kg], velocity [m/s].
    """

    name: str

    @abstractmethod
    def density(self, p: float, T: float) -> float:
        """Density rho(p, T). For an ideal gas this inverts eq. (12)."""

    @abstractmethod
    def enthalpy(self, T: float) -> float:
        """Static specific enthalpy h(T)."""

    @abstractmethod
    def sonic_velocity(self, p: float, T: float) -> float:
        """Speed of sound a(p, T). Used to evaluate Mach number and total quantities."""

    def total_enthalpy(self, T: float, velocity: float) -> float:
        """Total (stagnation) enthalpy h0 = h + V^2 / 2  (see paper, after eq. 3)."""
        return self.enthalpy(T) + 0.5 * velocity * velocity
