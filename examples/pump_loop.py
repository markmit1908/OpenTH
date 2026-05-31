"""Pump / compressor example.

Two demonstrations of the Pump component:

  1. A water pump driving flow *uphill* — from a low-pressure reservoir to a high-pressure
     one — settling on the operating point where the pump's head curve
     (head_shutoff - curve*mdot^2) balances the adverse pressure difference.
  2. A helium compressor raising both pressure and temperature: its shaft work shows up as
     a total-enthalpy (and temperature) rise in the energy equation.
"""

import math

from flowcalc.components import Pipe, PressureBoundary, Pump
from flowcalc.fluids import helium, water
from flowcalc.network import Network, Node
from flowcalc.solver import PCIMSolver, SolverConfig


def pump_uphill() -> None:
    fluid = water()
    net = Network()
    low = net.add_node(Node(id="low"))
    high = net.add_node(Node(id="high"))
    low.state.T = high.state.T = 300.0
    pump = net.add_element(Pump(id="pump", upstream=low, downstream=high, fluid=fluid,
                                head_shutoff=250e3, curve=500.0))
    PressureBoundary(node=low, p=200e3, T=300.0).apply()
    PressureBoundary(node=high, p=300e3, T=300.0).apply()

    PCIMSolver(net, SolverConfig(relaxation=0.6)).steady_state()
    print("1) Water pump, 200 kPa -> 300 kPa (100 kPa uphill):")
    print(f"   shutoff head 250 kPa -> delivers mdot = {pump.mdot:.2f} kg/s "
          f"(rise {pump.pressure_rise(pump.mdot)/1e3:.1f} kPa)\n")


def compressor() -> None:
    fluid = helium()
    T = 300.0
    area = 0.25 * math.pi * 0.3**2
    net = Network()
    inlet = net.add_node(Node(id="inlet"))
    mid = net.add_node(Node(id="mid", volume=2.0 * area))
    outlet = net.add_node(Node(id="outlet"))
    for nd in (inlet, mid, outlet):
        nd.state.T = T
    net.add_element(Pump(id="comp", upstream=inlet, downstream=mid, fluid=fluid,
                         head_shutoff=80e3, curve=200.0, efficiency=0.8))
    net.add_element(Pipe(id="pipe", upstream=mid, downstream=outlet, fluid=fluid,
                         length=5.0, diameter=0.3, friction_factor=0.02))
    PressureBoundary(node=inlet, p=200e3, T=T).apply()
    PressureBoundary(node=outlet, p=215e3, T=T).apply()

    PCIMSolver(net, SolverConfig(relaxation=0.4, solve_energy=True)).steady_state()
    print("2) Helium compressor (efficiency 0.8), inlet 200 kPa / 300 K:")
    print(f"   downstream T = {mid.state.T:.2f} K (rise {mid.state.T - T:+.2f} K from shaft work)")


def main() -> None:
    pump_uphill()
    compressor()


if __name__ == "__main__":
    main()
