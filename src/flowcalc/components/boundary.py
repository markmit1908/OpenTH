"""Boundary conditions.

The paper's examples use two kinds (Section 5): a fixed pressure (e.g. inlet at 700 kPa,
or outlet at 200 kPa) and a fixed mass flow (e.g. an outflow held at the steady-state
value, or zero when a valve shuts).

  * :class:`PressureBoundary` is a *Dirichlet* condition: it pins a node's pressure and
    temperature, and the solver removes that node from the pressure-correction unknowns.
  * :class:`MassFlowBoundary` is a *source*: the node stays an unknown (its pressure is
    solved) and the imposed mass flow enters the node's continuity balance.
"""

from __future__ import annotations

from dataclasses import dataclass

from ..network.node import Node


@dataclass
class PressureBoundary:
    """Hold a node's pressure (and temperature) fixed."""

    node: Node
    p: float
    T: float

    def apply(self) -> None:
        self.node.is_boundary = True
        self.node.state.p0 = self.p
        self.node.state.T = self.T
        self.node.fixed_temperature = True  # pressure boundaries also pin temperature


@dataclass
class MassFlowBoundary:
    """Impose a mass flow into (+) or out of (-) a node [kg/s].

    Used both for steady demand and for transient events such as a valve closing, where
    ``mdot`` is driven to zero. The node remains a pressure unknown. If ``T`` is given, the
    node's temperature is also pinned (appropriate for an inflow boundary when the energy
    equation is solved).
    """

    node: Node
    mdot: float
    T: float | None = None

    def apply(self) -> None:
        self.node.mass_source = self.mdot
        if self.T is not None:
            self.node.state.T = self.T
            self.node.fixed_temperature = True
