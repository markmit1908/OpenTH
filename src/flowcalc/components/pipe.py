"""Pipe element: the canonical face in the paper (Fig. 1, Section 3).

A pipe contributes a momentum balance with Darcy friction, fluid inertia (the transient
``rho dV/dt`` term) and gravity. The link coefficients ``a_plus`` / ``a_minus`` follow
eqs. (21)-(22).

NOTE: :meth:`momentum_coeffs` is deliberately left unimplemented. The exact coefficient
algebra (eqs. 21-22, including the effective friction factor f~ of eq. 11) is the core
numerical work of this project and must be transcribed *carefully* from the source PDF
(``docs/papers/``) and validated against the steady-state benchmark in eqs. (35)-(37)
before it can be trusted. The geometry/friction helpers below are the verified pieces.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from ..network.element import Element, MomentumCoeffs


@dataclass
class Pipe(Element):
    length: float = 1.0          # dx [m]
    diameter: float = 0.0        # D [m]
    friction_factor: float = 0.02  # Darcy f (constant in the paper's examples)
    angle: float = 0.0           # theta: angle of centreline to vertical [rad]

    @property
    def area(self) -> float:
        """Cross-sectional flow area A_e [m^2]."""
        return 0.25 * math.pi * self.diameter * self.diameter

    def friction_resistance(self, rho_face: float) -> float:
        """Quadratic friction resistance coefficient, i.e. the factor K in dp_fric = K Q|Q|.

        From the friction term f rho V|V| / (2D) integrated over the face:
            K = f * dx * rho / (2 * D * A^2)
        """
        A = self.area
        return self.friction_factor * self.length * rho_face / (2.0 * self.diameter * A * A)

    def momentum_coeffs(self, dt: float, alpha: float) -> MomentumCoeffs:
        # TODO(core): implement eqs. (21)-(22). Needs face density rho_e, face total
        # pressure p0_e, the preliminary flow Q_bar (from the current iteration), the
        # effective friction factor f~ (eq. 11), and the gravity term rho g dx cos(theta).
        # Validate the resulting steady state against eqs. (35)-(37) before relying on it.
        raise NotImplementedError(
            "Pipe.momentum_coeffs: transcribe eqs. (21)-(22) from docs/papers/ and "
            "validate against the steady-state benchmark before use."
        )

    def mass_flow(self) -> float:
        if self.fluid is None:
            raise ValueError(f"pipe {self.id!r} has no fluid assigned")
        # Upwind face density per the paper's upstream treatment of convected quantities.
        upwind = self.upstream if self.Q >= 0 else self.downstream
        rho = self.fluid.density(upwind.state.p0, upwind.state.T)
        return rho * self.Q
