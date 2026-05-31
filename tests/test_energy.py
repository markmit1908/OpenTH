"""Validation tests for the energy-equation coupling (non-isothermal flow).

  * adiabatic flow conserves total enthalpy h0 along the pipe, and temperature falls as the
    gas accelerates (expansion cooling) -- the kinetic part of h0;
  * a heat source raises the downstream total enthalpy by exactly Q / mdot (steady energy
    balance);
  * the transient energy solve is consistent with the steady one (marches to it);
  * with the energy equation off, temperatures are untouched (isothermal).

Conditions are kept at moderate Mach: the segregated pressure<->temperature coupling, like
the isothermal pressure-driven solve, is only robust below M ~ 0.4.
"""

import math

import pytest

from openth.components import Pipe, PressureBoundary
from openth.fluids import helium
from openth.network import Network, Node
from openth.solver import PCIMSolver, SolverConfig

T_IN = 300.0


def heated_pipeline(p_in: float = 230e3, heat_node: int | None = None,
                    heat: float = 0.0, n: int = 5) -> tuple[Network, list[Node]]:
    net = Network()
    fluid = helium()
    diameter, dx = 0.5, 5.0
    area = 0.25 * math.pi * diameter**2
    nodes = [net.add_node(Node(id=f"n{i}", volume=dx * area)) for i in range(n + 1)]
    for nd in nodes:
        nd.state.T = T_IN
    for i in range(n):
        net.add_element(Pipe(id=f"p{i}", upstream=nodes[i], downstream=nodes[i + 1],
                             fluid=fluid, length=dx, diameter=diameter, friction_factor=0.02))
    PressureBoundary(node=nodes[0], p=p_in, T=T_IN).apply()
    PressureBoundary(node=nodes[-1], p=200e3, T=T_IN).apply()
    if heat_node is not None:
        nodes[heat_node].heat_source = heat
    return net, nodes


def test_adiabatic_conserves_total_enthalpy_and_cools():
    net, nodes = heated_pipeline()
    result = PCIMSolver(net, SolverConfig(relaxation=0.4, solve_energy=True)).steady_state()
    assert result.converged

    h0 = [nd.state.h0 for nd in nodes]
    # Total enthalpy is uniform (adiabatic, no work extracted): h0_i == h0_inlet.
    spread = (max(h0[1:-1]) - min(h0[1:-1])) / h0[0]
    assert spread < 1e-6

    # Temperature falls monotonically downstream as the gas accelerates (expansion cooling).
    interior_T = [nd.state.T for nd in nodes[1:-1]]
    assert all(a > b for a, b in zip(interior_T, interior_T[1:], strict=False))
    assert interior_T[0] < T_IN


def test_heat_addition_raises_total_enthalpy_by_q_over_mdot():
    q = 5e5  # W
    net, nodes = heated_pipeline(heat_node=2, heat=q)
    result = PCIMSolver(net, SolverConfig(relaxation=0.4, solve_energy=True)).steady_state()
    assert result.converged

    mdot = net.elements["p0"].mdot
    # Steady energy balance: total enthalpy jumps by Q / mdot across the heated node.
    rise = nodes[4].state.h0 - nodes[1].state.h0
    assert rise == pytest.approx(q / mdot, rel=1e-3)


def test_transient_energy_matches_steady():
    # Steady reference with energy.
    net_s, nodes_s = heated_pipeline()
    PCIMSolver(net_s, SolverConfig(relaxation=0.4, solve_energy=True)).steady_state()
    steady_T = [nd.state.T for nd in nodes_s]

    # Transient with energy, marched to steady with constant boundaries.
    net_t, nodes_t = heated_pipeline()
    solver = PCIMSolver(net_t, SolverConfig(alpha=0.6, relaxation=0.4, solve_energy=True))
    for k in range(500):
        assert solver.step(dt=0.05, t=0.05 * (k + 1)).converged
    transient_T = [nd.state.T for nd in nodes_t]

    for ts, tt in zip(steady_T, transient_T, strict=True):
        assert tt == pytest.approx(ts, rel=1e-4)


def test_energy_off_leaves_temperature_isothermal():
    net, nodes = heated_pipeline()
    PCIMSolver(net, SolverConfig(relaxation=0.4)).steady_state()  # solve_energy defaults False
    assert all(nd.state.T == T_IN for nd in nodes)
