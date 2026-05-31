"""Non-isothermal example: adiabatic expansion cooling vs. heat addition.

Helium flows through a 25 m pipe between two pressure boundaries (230 -> 200 kPa, inlet
300 K). With the energy equation enabled the solver tracks temperature:

  * adiabatic (no heat): total enthalpy h0 is conserved, so the gas *cools* as it
    accelerates down the pressure gradient (the kinetic part of h0 grows);
  * with a heat source: the total enthalpy rises by Q / mdot, warming the gas.

This is the machinery behind the paper's non-isothermal cases (Fig. 3, 6, 7).
"""

import math

from flowcalc.components import Pipe, PressureBoundary
from flowcalc.fluids import helium
from flowcalc.network import Network, Node
from flowcalc.solver import PCIMSolver, SolverConfig

T_IN, P_IN, P_OUT = 300.0, 230e3, 200e3
N, DX, D = 5, 5.0, 0.5


def build(heat: float = 0.0) -> tuple[Network, list[Node]]:
    net = Network()
    fluid = helium()
    area = 0.25 * math.pi * D**2
    nodes = [net.add_node(Node(id=f"n{i}", volume=DX * area)) for i in range(N + 1)]
    for nd in nodes:
        nd.state.T = T_IN
    for i in range(N):
        net.add_element(Pipe(id=f"p{i}", upstream=nodes[i], downstream=nodes[i + 1],
                             fluid=fluid, length=DX, diameter=D, friction_factor=0.02))
    PressureBoundary(node=nodes[0], p=P_IN, T=T_IN).apply()
    PressureBoundary(node=nodes[-1], p=P_OUT, T=T_IN).apply()
    if heat:
        nodes[2].heat_source = heat
    return net, nodes


def report(label: str, heat: float) -> None:
    net, nodes = build(heat)
    PCIMSolver(net, SolverConfig(relaxation=0.4, solve_energy=True)).steady_state()
    mdot = net.elements["p0"].mdot
    temps = " ".join(f"{nd.state.T:6.2f}" for nd in nodes)
    print(f"{label:18} mdot={mdot:6.3f}  T[K]: {temps}")


def main() -> None:
    print("Helium, 230 -> 200 kPa, inlet 300 K (energy equation on)\n")
    report("adiabatic", 0.0)
    report("heat +0.5 MW @ n2", 5e5)


if __name__ == "__main__":
    main()
