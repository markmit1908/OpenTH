"""Pipe element: the canonical face in the paper (Fig. 1, Section 3).

A pipe contributes Darcy friction plus a convective-acceleration (momentum-flux) term.

Friction resistance, from f rho V|V| / (2D) with V = mdot / (rho A):

    dp_friction = K * mdot * |mdot|,   K = f * dx / (2 * D * A^2 * rho_face)

Convective term, from the steady momentum-flux change rho V dV across the (constant-area)
element, with mdot constant along it:

    dp_conv = (mdot / A) * (V_down - V_up),   V = mdot / (rho A)

This is the term the paper takes pains to retain (Section 2): neglecting it causes up to
~13% error in the pressure ratio at outlet Mach 0.7 (Fig. 2). For incompressible flow
rho is constant so V_up = V_down and dp_conv vanishes.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from ..network.element import Element


@dataclass
class Pipe(Element):
    length: float = 1.0          # dx [m]
    diameter: float = 0.0        # D [m]
    friction_factor: float = 0.02  # Darcy f (constant in the paper's examples)
    angle: float = 0.0           # theta: centreline-to-vertical angle [rad] (gravity; unused yet)

    @property
    def area(self) -> float:
        """Cross-sectional flow area A [m^2]."""
        return 0.25 * math.pi * self.diameter * self.diameter

    def resistance(self, rho_face: float) -> float:
        A = self.area
        return self.friction_factor * self.length / (2.0 * self.diameter * A * A * rho_face)

    def convective_dp(self, rho_up: float, rho_down: float, mdot: float) -> float:
        A = self.area
        v_up = mdot / (rho_up * A)
        v_down = mdot / (rho_down * A)
        return (mdot / A) * (v_down - v_up)

    def inertance(self) -> float:
        return self.length / self.area
