"""Validation tests for the Pump / compressor component.

  * a pump pushes flow against an adverse pressure gradient, settling on the operating
    point where its head curve meets the system (here, the imposed pressure difference);
  * with zero head the element is passive and flow follows the pressure gradient;
  * a compressor's shaft work raises the gas total enthalpy (and temperature) by exactly
    the specific work the energy equation adds.
"""

import math

import pytest

from openth.components import Pipe, PressureBoundary, Pump
from openth.fluids import helium, water
from openth.network import Network, Node
from openth.solver import PCIMSolver, SolverConfig


def solve(net: Network, energy: bool = False):
    cfg = SolverConfig(relaxation=0.4 if energy else 0.6, solve_energy=energy)
    return PCIMSolver(net, cfg).steady_state()


def test_pump_pushes_flow_against_adverse_gradient():
    fluid = water()
    head_shutoff, curve = 250e3, 500.0
    p_in, p_out = 200e3, 300e3  # outlet higher: flow would reverse without the pump

    net = Network()
    a = net.add_node(Node(id="A"))
    b = net.add_node(Node(id="B"))
    a.state.T = b.state.T = 300.0
    pump = net.add_element(Pump(id="pump", upstream=a, downstream=b, fluid=fluid,
                                head_shutoff=head_shutoff, curve=curve))
    PressureBoundary(node=a, p=p_in, T=300.0).apply()
    PressureBoundary(node=b, p=p_out, T=300.0).apply()

    result = solve(net)
    assert result.converged
    # Operating point: p_out - p_in = head_shutoff - curve * mdot^2.
    mdot_exact = math.sqrt((head_shutoff - (p_out - p_in)) / curve)
    assert pump.mdot > 0.0  # forced forward (uphill) by the pump
    assert pump.mdot == pytest.approx(mdot_exact, rel=1e-4)


def test_zero_head_pump_follows_pressure_gradient():
    fluid = water()
    net = Network()
    a = net.add_node(Node(id="A"))
    b = net.add_node(Node(id="B"))
    a.state.T = b.state.T = 300.0
    pump = net.add_element(Pump(id="pump", upstream=a, downstream=b, fluid=fluid,
                                head_shutoff=0.0, curve=500.0))
    PressureBoundary(node=a, p=300e3, T=300.0).apply()  # inlet higher
    PressureBoundary(node=b, p=200e3, T=300.0).apply()

    result = solve(net)
    assert result.converged
    assert pump.mdot == pytest.approx(math.sqrt(100e3 / 500.0), rel=1e-4)  # plain downhill flow


def test_compressor_raises_total_enthalpy_by_its_work():
    fluid = helium()
    T = 300.0
    area = 0.25 * math.pi * 0.3**2

    net = Network()
    inlet = net.add_node(Node(id="A"))
    mid = net.add_node(Node(id="C", volume=2.0 * area))
    outlet = net.add_node(Node(id="B"))
    for nd in (inlet, mid, outlet):
        nd.state.T = T
    comp = net.add_element(Pump(id="comp", upstream=inlet, downstream=mid, fluid=fluid,
                                head_shutoff=80e3, curve=200.0, efficiency=0.8))
    net.add_element(Pipe(id="pipe", upstream=mid, downstream=outlet, fluid=fluid,
                         length=5.0, diameter=0.3, friction_factor=0.02))
    PressureBoundary(node=inlet, p=200e3, T=T).apply()
    PressureBoundary(node=outlet, p=215e3, T=T).apply()

    result = solve(net, energy=True)
    assert result.converged

    work = comp.work_per_mass(comp.mdot, fluid.density(200e3, T))
    assert mid.state.h0 - inlet.state.h0 == pytest.approx(work, rel=1e-3)
    assert mid.state.T > T  # compression heats the gas
