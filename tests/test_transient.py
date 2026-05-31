"""Validation tests for the transient PCIM time step.

  * transient -> steady: with constant boundaries the transient must march to the same
    state the steady solver finds (the steady solution is the fixed point of the transient
    equations);
  * vessel fill: the compressibility storage term must conserve mass exactly in the
    scheme's theta-weighted sense, and the vessel must equilibrate to the supply pressure;
  * valve closure: a sudden closure must stop the flow through the valve and produce a
    water-hammer overpressure (which requires the momentum inertia term).
"""

import math

import pytest

from openth.components import MassFlowBoundary, Pipe, PressureBoundary, Valve
from openth.fluids import helium
from openth.network import Network, Node
from openth.solver import PCIMSolver, SolverConfig

T = 300.0


def series_pipeline(mdot_in: float, n: int = 5) -> tuple[Network, list[Node]]:
    net = Network()
    fluid = helium()
    diameter, dx = 0.5, 5.0
    area = 0.25 * math.pi * diameter**2
    nodes = [net.add_node(Node(id=f"n{i}", volume=dx * area)) for i in range(n + 1)]
    for nd in nodes:
        nd.state.T = T
    for i in range(n):
        net.add_element(Pipe(id=f"p{i}", upstream=nodes[i], downstream=nodes[i + 1],
                             fluid=fluid, length=dx, diameter=diameter, friction_factor=0.02))
    PressureBoundary(node=nodes[-1], p=200e3, T=T).apply()
    MassFlowBoundary(node=nodes[0], mdot=mdot_in).apply()
    return net, nodes


def test_transient_marches_to_steady_state():
    mdot_in = 30.0
    # Steady reference.
    net_s, nodes_s = series_pipeline(mdot_in)
    PCIMSolver(net_s, SolverConfig(relaxation=0.5)).steady_state()
    steady_p = [nd.state.p0 for nd in nodes_s]

    # Transient from rest, constant boundary, march out in time.
    net_t, nodes_t = series_pipeline(mdot_in)
    solver = PCIMSolver(net_t, SolverConfig(alpha=0.6, relaxation=0.5))
    for k in range(800):
        result = solver.step(dt=0.02, t=0.02 * (k + 1))
        assert result.converged
    transient_p = [nd.state.p0 for nd in nodes_t]

    for ps, pt in zip(steady_p, transient_p, strict=True):
        assert pt == pytest.approx(ps, rel=1e-5)
    assert net_t.elements["p0"].mdot == pytest.approx(mdot_in, rel=1e-5)


def test_vessel_fill_conserves_mass_and_equilibrates():
    fluid = helium()
    alpha = 0.6
    net = Network()
    supply = net.add_node(Node(id="S"))
    supply.state.T = T
    vessel = net.add_node(Node(id="V", volume=5.0))
    vessel.state.T = T
    vessel.state.p0 = 200e3
    net.add_element(Pipe(id="p", upstream=supply, downstream=vessel, fluid=fluid,
                         length=2.0, diameter=0.1, friction_factor=0.02))
    PressureBoundary(node=supply, p=300e3, T=T).apply()

    solver = PCIMSolver(net, SolverConfig(alpha=alpha, relaxation=0.6))
    rho_initial = fluid.density(vessel.state.p0, T)
    dt = 0.01
    influx = 0.0
    mdot_prev = 0.0
    for k in range(3000):
        solver.step(dt=dt, t=dt * (k + 1))
        mdot = net.elements["p"].mdot
        influx += dt * (alpha * mdot + (1.0 - alpha) * mdot_prev)  # theta-weighted integral
        mdot_prev = mdot

    # Equilibrium: vessel reaches the supply pressure, flow stops.
    assert vessel.state.p0 == pytest.approx(300e3, rel=1e-4)
    assert net.elements["p"].mdot == pytest.approx(0.0, abs=1e-4)

    # Conservation: mass gained by the vessel equals the integrated inflow.
    delta_mass = (fluid.density(vessel.state.p0, T) - rho_initial) * vessel.volume
    assert influx == pytest.approx(delta_mass, rel=1e-6)


def test_valve_closure_produces_water_hammer():
    fluid = helium()
    net = Network()
    inlet = net.add_node(Node(id="a"))
    mid = net.add_node(Node(id="b", volume=2.0 * 0.25 * math.pi * 0.5**2))
    outlet = net.add_node(Node(id="c"))
    for nd in (inlet, mid, outlet):
        nd.state.T = T
    mid.state.p0 = 700e3
    net.add_element(Pipe(id="pipe", upstream=inlet, downstream=mid, fluid=fluid,
                         length=20.0, diameter=0.5, friction_factor=0.02))
    valve = net.add_element(Valve(id="v", upstream=mid, downstream=outlet, fluid=fluid,
                                  k_open=1.0, opening=lambda t: 0.0 if t >= 0.05 else 1.0))
    PressureBoundary(node=inlet, p=700e3, T=T).apply()
    PressureBoundary(node=outlet, p=690e3, T=T).apply()

    solver = PCIMSolver(net, SolverConfig(alpha=0.6, relaxation=0.6))
    flow_before = 0.0
    peak_pressure = 0.0
    for k in range(400):
        t = 0.0009 * (k + 1)
        solver.step(dt=0.0009, t=t)
        if t < 0.05:
            flow_before = net.elements["pipe"].mdot
        else:
            peak_pressure = max(peak_pressure, mid.state.p0)
        assert math.isfinite(mid.state.p0)  # stays bounded

    assert flow_before > 0.0                      # flow was established before closure
    assert valve.mdot == pytest.approx(0.0)       # shut valve passes no flow
    assert peak_pressure > 700e3                   # water-hammer overpressure
