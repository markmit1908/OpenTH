"""Port / connection API example (openth.circuit.Circuit).

Builds a small closed helium loop the vision's way — components with `.inlet`/`.outlet`
ports, joined with `connect` — instead of naming nodes by hand:

    blower -> hot pipe -> (heat added) -> return pipe -> back to blower

The blower (a Pump) drives the flow; one node is pinned as the pressure reference (a closed
loop's absolute pressure level is otherwise undefined). Compare with `pump_loop.py`, which
builds a similar case through the lower-level name-based `FlowModel` API.
"""

import openth as th


def main() -> None:
    loop = th.Circuit(fluid=th.Fluid("helium"))

    blower = loop.add(th.Pump(head_shutoff=120e3, curve=300.0, name="blower"))
    hot = loop.add(th.Pipe(length=15.0, diameter=0.3, n_cells=6, name="hot_leg"))
    cold = loop.add(th.Pipe(length=15.0, diameter=0.3, n_cells=6, name="cold_leg"))

    # blower -> hot -> cold -> blower
    loop.connect(blower.outlet, hot.inlet)
    loop.connect(hot.outlet, cold.inlet)
    loop.connect(cold.outlet, blower.inlet)

    loop.pressure_boundary(blower.inlet, p=300e3, T=300.0)  # reference pressure

    result = loop.solve_steady_state(relaxation=0.6)
    print("Closed helium loop driven by a blower:")
    print(f"  converged: {result.converged}")
    print(f"  loop mass flow: {loop.flow(blower):.3f} kg/s")
    print(f"  blower:  in {loop.pressure(blower.inlet)/1e3:.1f} kPa "
          f"-> out {loop.pressure(blower.outlet)/1e3:.1f} kPa "
          f"(rise {(loop.pressure(blower.outlet) - loop.pressure(blower.inlet))/1e3:+.1f} kPa)")
    print(f"  compiled to {len(loop.model.network.nodes)} nodes / "
          f"{len(loop.model.network.elements)} elements")


if __name__ == "__main__":
    main()
