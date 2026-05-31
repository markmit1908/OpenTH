"""Network = the graph of nodes and elements, plus index bookkeeping for the solver.

The network owns the topology and assigns each *non-boundary* node a row index in the
pressure-correction system. For a single series pipeline the resulting matrix is
tridiagonal (solve with the Thomas algorithm); for branching networks it is a general
sparse matrix (paper, end of Section 4).
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .coupling import HeatExchanger
from .element import Element
from .node import Node


@dataclass
class Network:
    nodes: dict[str, Node] = field(default_factory=dict)
    elements: dict[str, Element] = field(default_factory=dict)
    heat_exchangers: dict[str, HeatExchanger] = field(default_factory=dict)

    def add_node(self, node: Node) -> Node:
        if node.id in self.nodes:
            raise ValueError(f"duplicate node id: {node.id!r}")
        self.nodes[node.id] = node
        return node

    def add_element(self, element: Element) -> Element:
        if element.id in self.elements:
            raise ValueError(f"duplicate element id: {element.id!r}")
        self.elements[element.id] = element
        return element

    def add_heat_exchanger(self, hx: HeatExchanger) -> HeatExchanger:
        """Register a thermal coupling (heat transfer only, used by the energy solve)."""
        if hx.id in self.heat_exchangers:
            raise ValueError(f"duplicate heat-exchanger id: {hx.id!r}")
        self.heat_exchangers[hx.id] = hx
        return hx

    def elements_at(self, node: Node) -> list[Element]:
        """All elements incident on ``node`` (used to assemble its continuity balance)."""
        return [
            e for e in self.elements.values()
            if e.upstream is node or e.downstream is node
        ]

    def solve_order(self) -> list[Node]:
        """Non-boundary nodes, in a stable order, indexed into the linear system."""
        return [n for n in self.nodes.values() if not n.is_boundary]

    def validate(self) -> None:
        """Cheap structural checks before a solve."""
        for e in self.elements.values():
            for end in (e.upstream, e.downstream):
                if self.nodes.get(end.id) is not end:
                    raise ValueError(
                        f"element {e.id!r} references node {end.id!r} not in this network"
                    )
        if not any(n.is_boundary for n in self.nodes.values()):
            raise ValueError("network has no boundary nodes; the system is singular")
