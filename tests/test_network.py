"""Tests for network topology assembly and validation."""

import pytest

from flowcalc.components import Pipe
from flowcalc.fluids import helium
from flowcalc.network import Network, Node


def build_series_pipeline(n_cells: int = 3) -> Network:
    net = Network()
    fluid = helium()
    nodes = [net.add_node(Node(id=f"n{i}", volume=1.0)) for i in range(n_cells + 1)]
    for i in range(n_cells):
        net.add_element(
            Pipe(
                id=f"p{i}",
                upstream=nodes[i],
                downstream=nodes[i + 1],
                fluid=fluid,
                length=5.0,
                diameter=0.5,
            )
        )
    nodes[0].is_boundary = True
    nodes[-1].is_boundary = True
    return net


def test_series_pipeline_topology():
    net = build_series_pipeline(3)
    assert len(net.nodes) == 4
    assert len(net.elements) == 3
    # interior nodes are the unknowns
    assert [n.id for n in net.solve_order()] == ["n1", "n2"]


def test_elements_at_interior_node_has_two():
    net = build_series_pipeline(3)
    assert len(net.elements_at(net.nodes["n1"])) == 2


def test_pipe_area():
    net = build_series_pipeline(1)
    pipe = net.elements["p0"]
    assert pipe.area == pytest.approx(3.141592653589793 * 0.25 * 0.5 * 0.5)


def test_duplicate_node_rejected():
    net = Network()
    net.add_node(Node(id="x"))
    with pytest.raises(ValueError, match="duplicate node"):
        net.add_node(Node(id="x"))


def test_validate_requires_boundary():
    net = Network()
    net.add_node(Node(id="a"))
    with pytest.raises(ValueError, match="no boundary nodes"):
        net.validate()
