"""Validation test for the HeatExchanger component.

Two helium streams (hot and cold) run side by side, coupled by a lumped UA. The hot stream
cools, the cold warms, and energy is conserved: the heat lost by the hot stream equals the
heat gained by the cold stream equals UA*(T_hot - T_cold).
"""

import pytest

import openth as th


def _stream(model, side, T_in, mdot):
    model.mass_flow_boundary(f"{side}_in", mdot=mdot, T=T_in)
    model.add_pipe(f"{side}_in", f"hx_{side}", length=5, diameter=0.1, n_cells=1, name=f"{side}_a")
    model.add_pipe(f"hx_{side}", f"{side}_out", length=5, diameter=0.1, n_cells=1, name=f"{side}_b")
    model.pressure_boundary(f"{side}_out", p=300e3, T=T_in)


def test_heat_exchanger_conserves_energy():
    he = th.helium()
    cp = he.enthalpy(1.0)
    mdot, UA = 0.5, 2000.0
    T_hot_in, T_cold_in = 500.0, 300.0

    m = th.Model(fluid=he)
    _stream(m, "h", T_hot_in, mdot)
    _stream(m, "c", T_cold_in, mdot)
    m.add_heat_exchanger("hx_h", "hx_c", UA=UA)

    result = m.steady_state(relaxation=0.5, max_outer_iterations=2000, solve_energy=True)
    assert result.converged

    t_hot, t_cold = m.temperature("hx_h"), m.temperature("hx_c")
    assert t_hot < T_hot_in           # hot stream cooled
    assert t_cold > T_cold_in         # cold stream warmed

    q_hot = mdot * cp * (T_hot_in - t_hot)     # heat lost by hot
    q_cold = mdot * cp * (t_cold - T_cold_in)  # heat gained by cold
    q_hx = UA * (t_hot - t_cold)               # exchanger duty
    assert q_hot == pytest.approx(q_cold, rel=2e-3)   # energy conserved
    assert q_hx == pytest.approx(q_hot, rel=2e-3)      # matches UA*dT
