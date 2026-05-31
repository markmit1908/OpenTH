"""Paper benchmark setup: steady isothermal helium flow through a 100 m pipeline.

Reproduces the configuration from Greyvenstein (2002), Section 5.1:
  - helium, 100 m long, 0.5 m diameter pipe, f = 0.02
  - inlet T = 300 K, outlet total pressure = 200 kPa, 20 length increments
  - initial pressures (except inlet) = outlet pressure, initial velocities = 0

This script *builds* the network today; running the solve is gated on PCIMSolver being
implemented (it currently raises NotImplementedError). Once the core lands, this becomes
a validation case against Figure 2 in the paper.
"""

from flowcalc.components import Pipe, PressureBoundary
from flowcalc.fluids import helium
from flowcalc.network import Network, Node
from flowcalc.solver import PCIMSolver, SolverConfig

N_CELLS = 20
LENGTH = 100.0
DIAMETER = 0.5
FRICTION = 0.02
T_IN = 300.0
P_OUT = 200e3


def build() -> Network:
    net = Network()
    fluid = helium()
    dx = LENGTH / N_CELLS
    cell_volume = dx * 0.25 * 3.141592653589793 * DIAMETER**2

    nodes = [net.add_node(Node(id=f"n{i}", volume=cell_volume)) for i in range(N_CELLS + 1)]
    for i in range(N_CELLS):
        net.add_element(
            Pipe(
                id=f"pipe{i}",
                upstream=nodes[i],
                downstream=nodes[i + 1],
                fluid=fluid,
                length=dx,
                diameter=DIAMETER,
                friction_factor=FRICTION,
            )
        )

    # Outlet total pressure fixed; inlet temperature fixed. (Inlet pressure is the result.)
    PressureBoundary(node=nodes[-1], p0=P_OUT, T=T_IN).apply()
    nodes[0].state.T = T_IN
    nodes[0].is_boundary = True
    return net


def main() -> None:
    net = build()
    print(f"Built {len(net.nodes)} nodes / {len(net.elements)} elements (helium pipeline).")
    solver = PCIMSolver(net, SolverConfig(alpha=0.6))
    try:
        result = solver.steady_state()
        print(result)
    except NotImplementedError as exc:
        print(f"[pending] solver core not implemented yet: {exc}")


if __name__ == "__main__":
    main()
