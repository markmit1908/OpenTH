"""Heat exchanger / recuperator: two helium streams coupled by a lumped UA.

A hot stream (500 K) and a cold stream (300 K) run side by side, thermally coupled by a
`HeatExchanger` (heat `UA*(T_hot - T_cold)` flows hot -> cold). With the energy equation on,
the hot stream cools and the cold warms; the heat lost by one equals the heat gained by the
other. The two streams are hydraulically independent — they only exchange heat — which is the
recuperator arrangement of a closed-cycle gas reactor.
"""

import openth as th

UA = 2000.0  # W/K
MDOT = 0.5   # kg/s in each stream


def stream(model: th.Model, side: str, T_in: float) -> None:
    model.mass_flow_boundary(f"{side}_in", mdot=MDOT, T=T_in)
    model.add_pipe(f"{side}_in", f"hx_{side}", length=5, diameter=0.1, n_cells=1, name=f"{side}_a")
    model.add_pipe(f"hx_{side}", f"{side}_out", length=5, diameter=0.1, n_cells=1, name=f"{side}_b")
    model.pressure_boundary(f"{side}_out", p=300e3, T=T_in)


def main() -> None:
    he = th.helium()
    cp = he.enthalpy(1.0)
    m = th.Model(fluid=he)
    stream(m, "h", 500.0)
    stream(m, "c", 300.0)
    m.add_heat_exchanger("hx_h", "hx_c", UA=UA)

    m.steady_state(relaxation=0.5, max_outer_iterations=2000, solve_energy=True)
    t_hot, t_cold = m.temperature("hx_h"), m.temperature("hx_c")

    print(f"Heat exchanger (UA = {UA:.0f} W/K), {MDOT} kg/s helium each side:\n")
    print(f"  hot stream : 500.0 K -> {t_hot:6.1f} K")
    print(f"  cold stream: 300.0 K -> {t_cold:6.1f} K")
    print(f"  duty       : {MDOT * cp * (500.0 - t_hot) / 1e3:6.1f} kW lost by the hot stream")
    print(f"               {MDOT * cp * (t_cold - 300.0) / 1e3:6.1f} kW gained by the cold stream")


if __name__ == "__main__":
    main()
