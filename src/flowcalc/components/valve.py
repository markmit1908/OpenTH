"""Valve: a non-pipe component whose resistance varies in time.

The paper's transient examples (Section 5.2-5.3) are driven by *sudden valve closures*.
A valve is modelled as a variable flow resistance; closing drives the resistance to
infinity (flow -> 0). This is the simplest non-pipe component and the template for the
others (pump, compressor, turbine, orifice, heat exchanger), each of which overrides the
momentum closure with its own characteristic.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

from ..network.element import Element, MomentumCoeffs


@dataclass
class Valve(Element):
    """Variable-resistance valve.

    Parameters
    ----------
    cv_open : float
        Flow coefficient when fully open.
    opening : Callable[[float], float]
        Fractional opening in [0, 1] as a function of time; ``opening(t)=0`` is shut.
        Default is always-open. A sudden closure at t0 is ``lambda t: 0.0 if t >= t0 else 1.0``.
    """

    cv_open: float = 1.0
    opening: Callable[[float], float] = lambda t: 1.0

    def momentum_coeffs(self, dt: float, alpha: float) -> MomentumCoeffs:
        # TODO(core): a valve imposes dp = (Q|Q|) / (cv * opening)^2 across the face;
        # assemble the linearised link coefficients consistently with the pipe's eqs.
        # (21)-(22) but with the valve resistance in place of pipe friction.
        raise NotImplementedError("Valve.momentum_coeffs not yet implemented")

    def mass_flow(self) -> float:
        if self.fluid is None:
            raise ValueError(f"valve {self.id!r} has no fluid assigned")
        upwind = self.upstream if self.Q >= 0 else self.downstream
        return self.fluid.density(upwind.state.p0, upwind.state.T) * self.Q
