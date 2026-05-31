"""Boundary conditions.

The paper's examples use two kinds (Section 5): a fixed total pressure (e.g. inlet at
700 kPa) and a fixed mass flow (e.g. outflow held at the steady-state value, or zero when
a valve shuts). Boundaries are attached to nodes; a node marked ``is_boundary`` is removed
from the pressure-correction unknowns and instead has its value imposed here.
"""

from __future__ import annotations

from dataclasses import dataclass

from ..network.node import Node


@dataclass
class PressureBoundary:
    """Hold a node's total pressure (and temperature) fixed."""

    node: Node
    p0: float
    T: float

    def apply(self) -> None:
        self.node.is_boundary = True
        self.node.state.p0 = self.p0
        self.node.state.T = self.T


@dataclass
class MassFlowBoundary:
    """Impose a mass flow into/out of a node (negative = outflow).

    Used both for steady demand and for transient events such as a valve closing,
    where ``mdot`` is driven to zero (see :class:`flowcalc.components.valve.Valve`).
    """

    node: Node
    mdot: float

    def apply(self) -> None:
        self.node.is_boundary = True
