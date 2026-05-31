"""Validation tests for the steady-state PCIM core.

Each test checks the solver against an independently-computable answer:
  * incompressible series -> exact resistance sum,
  * incompressible parallel -> flow split by sqrt(resistance ratio) + mass conservation,
  * compressible isothermal pipe -> the closed-form isothermal pipe-flow relation.
"""

import math

import pytest

from flowcalc.components import MassFlowBoundary, Pipe, PressureBoundary
from flowcalc.fluids import helium, water
from flowcalc.network import Network, Node
from flowcalc.solver import PCIMSolver, SolverConfig


def pipe_resistance(fluid, length, diameter, f, p, T):
    """Independent K such that dp = K * mdot * |mdot| (matches Pipe.resistance)."""
    area = 0.25 * math.pi * diameter**2
    rho = fluid.density(p, T)
    return f * length / (2.0 * diameter * area * area * rho)


def test_incompressible_series_exact():
    fluid = water()
    n_pipes, L, D, f = 4, 10.0, 0.2, 0.02
    p_out, mdot, T = 300e3, 50.0, 300.0

    net = Network()
    nodes = [net.add_node(Node(id=f"n{i}")) for i in range(n_pipes + 1)]
    for nd in nodes:
        nd.state.T = T
    for i in range(n_pipes):
        net.add_element(Pipe(id=f"p{i}", upstream=nodes[i], downstream=nodes[i + 1],
                             fluid=fluid, length=L, diameter=D, friction_factor=f))
    PressureBoundary(node=nodes[-1], p=p_out, T=T).apply()
    MassFlowBoundary(node=nodes[0], mdot=mdot).apply()

    result = PCIMSolver(net, SolverConfig(relaxation=0.7)).steady_state()
    assert result.converged

    k = pipe_resistance(fluid, L, D, f, p_out, T)  # density constant -> K constant
    dp_per_pipe = k * mdot * mdot
    # Node i (from inlet) sits i pipes upstream of the outlet.
    for i, nd in enumerate(nodes):
        expected = p_out + (n_pipes - i) * dp_per_pipe
        assert nd.state.p0 == pytest.approx(expected, rel=1e-6)
    # Every pipe carries the full mass flow.
    for e in net.elements.values():
        assert e.mdot == pytest.approx(mdot, rel=1e-6)


def test_incompressible_parallel_split_and_conservation():
    fluid = water()
    D, f, p_out, mdot, T = 0.2, 0.02, 300e3, 60.0, 300.0

    net = Network()
    a = net.add_node(Node(id="A"))
    b = net.add_node(Node(id="B"))
    a.state.T = b.state.T = T
    # Two parallel paths A->B with different lengths (hence different resistance).
    net.add_element(Pipe(id="p1", upstream=a, downstream=b, fluid=fluid,
                         length=10.0, diameter=D, friction_factor=f))
    net.add_element(Pipe(id="p2", upstream=a, downstream=b, fluid=fluid,
                         length=40.0, diameter=D, friction_factor=f))
    PressureBoundary(node=b, p=p_out, T=T).apply()
    MassFlowBoundary(node=a, mdot=mdot).apply()

    result = PCIMSolver(net).steady_state()
    assert result.converged

    m1 = net.elements["p1"].mdot
    m2 = net.elements["p2"].mdot
    # Mass conservation at A.
    assert m1 + m2 == pytest.approx(mdot, rel=1e-6)
    # Equal pressure drop -> K1 m1^2 = K2 m2^2 -> m1/m2 = sqrt(K2/K1) = sqrt(L2/L1).
    assert m1 / m2 == pytest.approx(math.sqrt(40.0 / 10.0), rel=1e-4)


def test_compressible_isothermal_matches_closed_form():
    fluid = helium()
    L, D, f, p_out, T = 100.0, 0.5, 0.02, 200e3, 300.0
    n_cells = 20
    area = 0.25 * math.pi * D**2
    mach = 0.3
    mdot = fluid.density(p_out, T) * (mach * fluid.sonic_velocity(p_out, T)) * area

    net = Network()
    nodes = [net.add_node(Node(id=f"n{i}")) for i in range(n_cells + 1)]
    for nd in nodes:
        nd.state.T = T
    dx = L / n_cells
    for i in range(n_cells):
        net.add_element(Pipe(id=f"p{i}", upstream=nodes[i], downstream=nodes[i + 1],
                             fluid=fluid, length=dx, diameter=D, friction_factor=f))
    PressureBoundary(node=nodes[-1], p=p_out, T=T).apply()
    MassFlowBoundary(node=nodes[0], mdot=mdot).apply()

    result = PCIMSolver(net, SolverConfig(relaxation=0.5)).steady_state()
    assert result.converged

    # Closed form: p1^2 - p2^2 = G^2 R T (fL/D + 2 ln(p1/p2)).
    g = mdot / area
    c = g * g * fluid.R * T
    p1_exact = p_out
    for _ in range(200):
        p1_exact = math.sqrt(p_out**2 + c * (f * L / D + 2.0 * math.log(p1_exact / p_out)))

    assert nodes[0].state.p0 == pytest.approx(p1_exact, rel=2e-3)


def test_mass_conserved_at_interior_nodes():
    fluid = helium()
    net = Network()
    nodes = [net.add_node(Node(id=f"n{i}")) for i in range(5)]
    for nd in nodes:
        nd.state.T = 300.0
    for i in range(4):
        net.add_element(Pipe(id=f"p{i}", upstream=nodes[i], downstream=nodes[i + 1],
                             fluid=fluid, length=5.0, diameter=0.3, friction_factor=0.02))
    PressureBoundary(node=nodes[-1], p=200e3, T=300.0).apply()
    MassFlowBoundary(node=nodes[0], mdot=5.0).apply()

    PCIMSolver(net).steady_state()
    flows = [e.mdot for e in net.elements.values()]
    # Series: every interior balance => identical mass flow through each pipe.
    assert max(flows) - min(flows) < 1e-6 * max(flows)
