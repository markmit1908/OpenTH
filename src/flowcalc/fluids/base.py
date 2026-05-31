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

    def drho_dp(self, p: float, T: float) -> float:
        """Isothermal compressibility (d rho / d p)_T [s^2/m^2].

        This is the coefficient that couples a pressure correction to a density change in
        the transient continuity storage term (paper eq. 17). Default is a central finite
        difference; closed-form overrides are preferred where available.
        """
        dp = max(1.0, 1e-6 * abs(p))
        return (self.density(p + dp, T) - self.density(p - dp, T)) / (2.0 * dp)

    def total_enthalpy(self, T: float, velocity: float) -> float:
        """Total (stagnation) enthalpy h0 = h + V^2 / 2  (see paper, after eq. 3)."""
        return self.enthalpy(T) + 0.5 * velocity * velocity

    @abstractmethod
    def temperature_from_enthalpy(self, h: float) -> float:
        """Invert the static enthalpy relation: return T such that enthalpy(T) = h.

        Used by the energy solver to recover temperature from the solved total enthalpy
        (after subtracting the kinetic part): T = inverse_enthalpy(h0 - V^2/2).
        """
