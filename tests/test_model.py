"""Tests for the FlowModel high-level facade."""

import math

import pytest

from flowcalc.fluids import helium
from flowcalc.model import FlowModel


def test_add_pipe_subdivides_into_cells():
    model = FlowModel(fluid=helium())
    model.add_pipe("inlet", "outlet", length=100.0, diameter=0.5, n_cells=20, name="pipe")
    # 20 segments (faces) and 21 nodes (2 named ends + 19 internal).
    assert len(model.network.elements) == 20
    assert len(model.network.nodes) == 21
    # Each interior node owns one full cell volume (two half-cells); ends own a half-cell.
    area = 0.25 * math.pi * 0.5**2
    dx = 100.0 / 20
    assert model.node("pipe#c0").volume == pytest.approx(dx * area)
    assert model.node("inlet").volume == pytest.approx(0.5 * dx * area)


def test_steady_pipeline_matches_closed_form():
    fluid = helium()
    T, p_out, L, D, f = 300.0, 200e3, 100.0, 0.5, 0.02
    area = 0.25 * math.pi * D**2
    mdot = fluid.density(p_out, T) * 0.3 * fluid.sonic_velocity(p_out, T) * area

    model = FlowModel(fluid=fluid)
    model.add_pipe("inlet", "outlet", length=L, diameter=D, friction_factor=f,
                   n_cells=20, name="pipe")
    model.pressure_boundary("outlet", p=p_out, T=T)
    model.mass_flow_boundary("inlet", mdot=mdot, T=T)

    result = model.steady_state(relaxation=0.5, max_outer_iterations=600)
    assert result.converged

    g = mdot / area
    c = g * g * fluid.R * T
    p1 = p_out
    for _ in range(200):
        p1 = math.sqrt(p_out**2 + c * (f * L / D + 2.0 * math.log(p1 / p_out)))
    assert model.pressure("inlet") == pytest.approx(p1, rel=2e-3)
    assert model.flow_through("pipe") == pytest.approx(mdot, rel=1e-4)


def test_time_varying_pressure_boundary_is_applied_during_run():
    # A vessel fed through a pipe from a ramping supply pressure should track it (slow fill).
    model = FlowModel(fluid=helium())
    model.add_pipe("supply", "vessel", length=2.0, diameter=0.1, n_cells=1, name="pipe")
    model.node("vessel").volume = 5.0

    def supply(t: float) -> float:
        return 200e3 + 1e3 * t  # 1 kPa/s ramp

    model.pressure_boundary("supply", p=supply, T=300.0)

    hist = model.run(dt=0.1, duration=2.0, record=("p:supply",),
                     alpha=0.6, relaxation=0.6, steady_init=False)
    # The supply boundary was driven by the callable: last value = supply(2.0).
    assert hist["p:supply"][-1] == pytest.approx(supply(2.0), rel=1e-9)
    # The vessel pressure rose toward the (higher) supply.
    assert model.pressure("vessel") > 200e3
