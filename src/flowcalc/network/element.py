"""Element = a face / branch connecting two nodes (where flow is solved).

Each element carries a mass flow and supplies the discretized **momentum equation** to
the solver. In the segregated (SIMPLE / PCIM) scheme the momentum balance across a face
is written as

    p_up - p_down  =  K * mdot * |mdot|  +  C

where ``K`` is the (density-dependent) quadratic friction *resistance* and ``C`` is an
explicit convective-acceleration (momentum-flux) term. Greyvenstein eliminates the
convective term analytically by working in *total* pressure (paper eqs. 4-11); carrying it
as an explicit ``C`` here is mathematically equivalent and keeps each component's closure
self-contained.

Linearising the friction term gives the pressure-correction link coefficient used to
assemble the continuity system (paper eqs. 20-22):

    mdot' = d * (p'_up - p'_down),   d = 1 / (2 K |mdot|)

A *pipe* computes ``K`` from Darcy friction and ``C`` from the density change along it;
*non-pipe* components (valves, pumps, ...) override ``resistance``/``convective_dp`` with
their own characteristic. This is the single extension point that makes the method work
uniformly across component types.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field

from ..fluids import FluidModel
from .node import Node


@dataclass
class Element(ABC):
    """Abstract face/branch connecting an upstream node to a downstream node.

    Sign convention: positive ``mdot`` flows upstream -> downstream.
    """

    id: str
    upstream: Node
    downstream: Node
    mdot: float = 0.0           # mass flow rate [kg/s] on this face
    fluid: FluidModel | None = field(default=None)

    @abstractmethod
    def resistance(self, rho_face: float) -> float:
        """Quadratic friction resistance K, so that the friction pressure drop is K*mdot*|mdot|.

        ``rho_face`` is the representative density at the face [kg/m^3].
        """

    def convective_dp(self, rho_up: float, rho_down: float, mdot: float) -> float:
        """Static-pressure change from momentum flux (convective acceleration).

        Zero by default (e.g. incompressible constant-area flow). Pipes override it for
        the density change along the element. See paper Section 2 on why this term matters
        at higher Mach number (Fig. 2).
        """
        return 0.0

    def _require_fluid(self) -> FluidModel:
        if self.fluid is None:
            raise ValueError(f"element {self.id!r} has no fluid assigned")
        return self.fluid

    def density_at(self, node: Node) -> float:
        """Fluid density evaluated at one of this element's end nodes."""
        fluid = self._require_fluid()
        return fluid.density(node.state.p0, node.state.T)

    def face_density(self) -> float:
        """Representative (arithmetic-mean) density across the face."""
        return 0.5 * (self.density_at(self.upstream) + self.density_at(self.downstream))
