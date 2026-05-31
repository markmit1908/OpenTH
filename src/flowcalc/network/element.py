"""Element = a face / branch connecting two nodes (where flow is solved).

Each element carries a volumetric flow rate ``Q`` and supplies the discretized
**momentum equation** to the solver. The crux of the PCIM method is that, for each
element, the corrected flow can be written as a *linear* function of the pressure
corrections at its two end nodes (paper eq. 20):

    Q'_i = a_minus * p'_upstream  -  a_plus * p'_downstream

The ``a_plus`` / ``a_minus`` link coefficients (eqs. 21-22) depend on the element's
physics (friction, area, inertia, gravity). A *pipe* computes them from Darcy friction;
*non-pipe* components (valves, pumps, compressors, ...) override the momentum closure
with their own characteristic. This is the single extension point that makes the method
work uniformly across component types.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field

from ..fluids import FluidModel
from .node import Node


@dataclass
class MomentumCoeffs:
    """Link coefficients for one element's contribution to the pressure-correction system.

    See eqs. (20)-(22). ``a_plus`` multiplies the downstream node's pressure correction,
    ``a_minus`` the upstream node's.
    """

    a_plus: float
    a_minus: float


@dataclass
class Element(ABC):
    """Abstract face/branch connecting an upstream node to a downstream node.

    Parameters
    ----------
    id : str
        Unique label.
    upstream, downstream : Node
        The two control volumes this element connects. Sign convention: positive ``Q``
        flows upstream -> downstream.
    """

    id: str
    upstream: Node
    downstream: Node
    Q: float = 0.0          # volumetric flow rate [m^3/s] on this face
    fluid: FluidModel | None = field(default=None)

    @abstractmethod
    def momentum_coeffs(self, dt: float, alpha: float) -> MomentumCoeffs:
        """Assemble this element's momentum link coefficients (eqs. 21-22).

        Parameters
        ----------
        dt : float
            Time step. (Steady state is recovered as dt -> infinity / the transient
            terms dropped.)
        alpha : float
            Time-integration weighing factor in [0.5, 1]; see the solver. alpha=1 is
            fully implicit, alpha=0.5 Crank-Nicolson, alpha=0.6 a good compromise.
        """

    @abstractmethod
    def mass_flow(self) -> float:
        """Mass flow rate = rho_face * Q, using the upwind/face density."""
