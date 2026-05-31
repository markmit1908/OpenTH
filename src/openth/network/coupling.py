"""Thermal couplings between nodes (heat transfer without flow).

A :class:`HeatExchanger` ties two nodes together thermally: heat ``UA*(T_hot - T_cold)``
flows from the hot node to the cold node. It is **not** a flow `Element` — it transfers heat
only, entering the energy equation (so it does nothing unless ``solve_energy`` is set). The
two nodes are usually on separate flow streams (the two sides of a recuperator/heat
exchanger), which in OpenTH's single-fluid model is a same-fluid exchanger such as a
gas-cycle recuperator.

This is a *lumped* (single node-pair) model; a distributed counter-flow exchanger is built by
coupling several node pairs along the two streams.
"""

from __future__ import annotations

from dataclasses import dataclass

from .node import Node


@dataclass
class HeatExchanger:
    """Lumped thermal coupling: heat ``UA*(T_hot - T_cold)`` flows hot -> cold."""

    id: str
    hot: Node
    cold: Node
    UA: float  # overall heat-transfer conductance [W/K]

    def duty(self) -> float:
        """Current heat duty [W] transferred from the hot node to the cold node."""
        return self.UA * (self.hot.state.T - self.cold.state.T)
