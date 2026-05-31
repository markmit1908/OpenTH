"""Natural circulation: buoyancy drives flow around a closed loop with no pump.

A square gas loop is heated at the bottom of the right-hand leg and held at a reference
temperature at the top-left corner (the "cooler"). The heated gas is lighter, so it rises in
the right leg while the cooler, denser gas sinks in the left leg — the difference in
hydrostatic head (g*rho*z, with rho set by temperature) drives a steady circulation.

This needs the gravity term *and* the energy equation, and it bootstraps transiently: with
no temperature difference there is no buoyancy, so we march in time and watch the flow spin
up from rest. Switch gravity off (`gravity=0.0`) and the circulation disappears.
"""

import openth as th

H = 3.0  # loop height [m]


def build(heat_power: float) -> th.Model:
    m = th.Model(fluid=th.air())
    m.add_pipe("BL", "BR", length=3.0, diameter=0.1, n_cells=2, name="bottom")
    m.add_pipe("BR", "TR", length=H, diameter=0.1, n_cells=2, delta_elevation=+H, name="riser")
    m.add_pipe("TR", "TL", length=3.0, diameter=0.1, n_cells=2, name="top")
    m.add_pipe("TL", "BL", length=H, diameter=0.1, n_cells=2, delta_elevation=-H, name="downcomer")
    m.heat_source("BR", power=heat_power)            # heater at the base of the riser
    m.pressure_boundary("TL", p=200e3, T=300.0)      # pressure + cold reference
    return m


def main() -> None:
    m = build(3000.0)
    hist = m.run(dt=0.1, duration=20.0, record=("flow:riser", "T:BR"),
                 alpha=0.6, relaxation=0.4, solve_energy=True, steady_init=False)
    t, flow = hist["t"], hist["flow:riser"]

    print("Natural circulation in a heated closed loop (3 kW, no pump):\n")
    print(f"{'t [s]':>6} {'circulation [kg/s]':>20} {'riser-base T [K]':>18}")
    for i in range(0, len(t), len(t) // 8):
        print(f"{t[i]:6.1f} {flow[i]:20.5f} {hist['T:BR'][i]:18.1f}")
    print(f"\nSteady circulation: {flow[-1]:.5f} kg/s up the heated riser (buoyancy-driven).")


if __name__ == "__main__":
    main()
