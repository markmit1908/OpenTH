"""Valve: a non-pipe component whose resistance varies in time.

The paper's transient examples (Section 5.2-5.3) are driven by *sudden valve closures*.
A valve is a variable flow resistance; as it closes the resistance grows without bound and
the flow is choked off. This is the simplest non-pipe component and the template for the
others (pump, compressor, turbine, orifice, heat exchanger), each of which overrides the
momentum closure (``resistance`` / ``convective_dp``) with its own characteristic.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass

from ..network.element import Element


@dataclass
class Valve(Element):
    """Variable-resistance valve.

    Parameters
    ----------
    k_open : float
        Loss coefficient when fully open, in the resistance ``K = k_open / opening^2``.
    opening : Callable[[float], float]
        Fractional opening in (0, 1] as a function of time; ``opening(t) -> 0`` shuts it.
        Default is always fully open. The steady solve uses ``opening(0.0)``.
    """

    k_open: float = 1.0
    opening: Callable[[float], float] = lambda t: 1.0
    _t: float = 0.0  # time at which to evaluate the opening (set by the solver)

    def resistance(self, rho_face: float) -> float:
        frac = self.opening(self._t)
        if frac <= 0.0:
            return math.inf  # shut: no flow
        # Density-scaled quadratic loss; K such that dp = K * mdot * |mdot|.
        return self.k_open / (frac * frac * rho_face)
