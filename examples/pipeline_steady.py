"""Paper benchmark: steady isothermal helium flow through a 100 m pipeline (Section 5.1).

Configuration from Greyvenstein (2002), Fig. 2:
  - helium, 100 m long, 0.5 m diameter pipe, f = 0.02, 20 length increments
  - isothermal at 300 K, outlet pressure 200 kPa

We impose the outlet pressure and an inlet mass flow sized to hit a target outlet Mach
number, solve the steady state, and compare the resulting pressure ratio P1/P2 against the
closed-form isothermal compressible pipe-flow relation

    p1^2 - p2^2 = G^2 R T [ f L / D + 2 ln(p1 / p2) ],   G = mdot / A

which includes the convective-acceleration ("2 ln") term the paper retains.
"""

import math

from flowcalc.components import MassFlowBoundary, Pipe, PressureBoundary
from flowcalc.fluids import helium
from flowcalc.network import Network, Node
from flowcalc.solver import PCIMSolver, SolverConfig

N_CELLS = 20
LENGTH = 100.0
DIAMETER = 0.5
FRICTION = 0.02
T = 300.0
P_OUT = 200e3


def build(mdot: float) -> Network:
    net = Network()
    fluid = helium()
    dx = LENGTH / N_CELLS
    area = 0.25 * math.pi * DIAMETER**2
    cell_volume = dx * area

    nodes = [net.add_node(Node(id=f"n{i}", volume=cell_volume)) for i in range(N_CELLS + 1)]
    for node in nodes:
        node.state.T = T  # isothermal
    for i in range(N_CELLS):
        net.add_element(
            Pipe(id=f"pipe{i}", upstream=nodes[i], downstream=nodes[i + 1], fluid=fluid,
                 length=dx, diameter=DIAMETER, friction_factor=FRICTION)
        )

    PressureBoundary(node=nodes[-1], p=P_OUT, T=T).apply()
    MassFlowBoundary(node=nodes[0], mdot=mdot).apply()
    return net


def analytical_inlet_pressure(mdot: float) -> float:
    """Solve p1^2 - p2^2 = G^2 R T (fL/D + 2 ln(p1/p2)) for p1 by fixed-point iteration."""
    fluid = helium()
    area = 0.25 * math.pi * DIAMETER**2
    g = mdot / area
    c = g * g * fluid.R * T
    p1 = P_OUT
    for _ in range(200):
        p1 = math.sqrt(P_OUT**2 + c * (FRICTION * LENGTH / DIAMETER + 2.0 * math.log(p1 / P_OUT)))
    return p1


def main() -> None:
    fluid = helium()
    area = 0.25 * math.pi * DIAMETER**2
    a_out = fluid.sonic_velocity(P_OUT, T)
    rho_out = fluid.density(P_OUT, T)

    print(f"{'M_out':>6} {'P1/P2 (PCIM)':>14} {'P1/P2 (exact)':>15} {'err %':>8} {'iters':>6}")
    for mach in (0.1, 0.2, 0.3, 0.4, 0.5):
        mdot = rho_out * (mach * a_out) * area
        net = build(mdot)
        solver = PCIMSolver(net, SolverConfig(relaxation=0.5))
        result = solver.steady_state()

        p1 = net.nodes["n0"].state.p0
        ratio = p1 / P_OUT
        exact = analytical_inlet_pressure(mdot) / P_OUT
        err = 100.0 * (ratio - exact) / exact
        print(f"{mach:6.2f} {ratio:14.4f} {exact:15.4f} {err:8.2f} {result.iterations:6d}")


if __name__ == "__main__":
    main()
