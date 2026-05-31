"""Paper benchmark: pressure-vessel blow-down (Section 5.4, Fig. 10).

Helium at 300 K flows from an upstream vessel through a 10 m x 0.1 m pipe (f = 0.02) into a
downstream vessel held at 650 kPa. The upstream pressure decays as

    p1(t) = 650 + 50 * exp(-0.004 t)  [kPa]

so the driving pressure difference (and hence the mass flow) decays to zero over ~30 min.
This is the slow transient where the implicit method is dramatically faster than explicit
schemes (the paper used a 10 s time step).

Because the transient is slow, the flow is essentially quasi-steady: at each instant the
mass flow should match the steady-state flow for the instantaneous pressures. We exploit
that as a check, comparing the transient flow against a frozen-pressure steady solve.
"""

import math

from flowcalc.components import Pipe, PressureBoundary
from flowcalc.fluids import helium
from flowcalc.network import Network, Node
from flowcalc.solver import PCIMSolver, SolverConfig

T = 300.0
P2 = 650e3
LENGTH, DIAMETER, FRICTION = 10.0, 0.1, 0.02


def p1_of_t(t: float) -> float:
    return (650.0 + 50.0 * math.exp(-0.004 * t)) * 1e3


def build() -> tuple[Network, Node, Node]:
    net = Network()
    fluid = helium()
    up = net.add_node(Node(id="up"))
    down = net.add_node(Node(id="down"))
    up.state.T = down.state.T = T
    net.add_element(Pipe(id="pipe", upstream=up, downstream=down, fluid=fluid,
                         length=LENGTH, diameter=DIAMETER, friction_factor=FRICTION))
    PressureBoundary(node=up, p=p1_of_t(0.0), T=T).apply()
    PressureBoundary(node=down, p=P2, T=T).apply()
    return net, up, down


def quasi_steady_flow(p1: float) -> float:
    net, up, down = build()
    up.state.p0 = p1
    PCIMSolver(net, SolverConfig(relaxation=0.5)).steady_state()
    return net.elements["pipe"].mdot


def main() -> None:
    net, up, _ = build()
    solver = PCIMSolver(net, SolverConfig(alpha=0.6, relaxation=0.6))
    dt = 10.0
    sample_times = {0, 250, 500, 1000, 1750}

    print(f"{'t [s]':>6} {'p1 [kPa]':>9} {'mdot [kg/s]':>12} {'quasi-steady':>13} {'err %':>7}")
    t = 0.0
    while t < 1750.0:
        t += dt
        up.state.p0 = p1_of_t(t)  # time-varying upstream pressure boundary
        solver.step(dt=dt, t=t)
        if int(t) in sample_times:
            mdot = net.elements["pipe"].mdot
            qs = quasi_steady_flow(p1_of_t(t))
            err = 100.0 * (mdot - qs) / qs if qs else 0.0
            print(f"{t:6.0f} {p1_of_t(t)/1e3:9.2f} {mdot:12.5f} {qs:13.5f} {err:7.2f}")


if __name__ == "__main__":
    main()
