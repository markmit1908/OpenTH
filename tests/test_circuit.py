"""Tests for the port/connection layer (openth.circuit.Circuit)."""

import math

import pytest

import openth as th


def test_connections_merge_nodes_and_match_flowmodel():
    # Two 50 m pipes connected in series must equal one 100 m pipe in the name-based model.
    c = th.Circuit(fluid=th.Fluid("helium"))
    p1 = c.add(th.Pipe(length=50, diameter=0.5, n_cells=10))
    p2 = c.add(th.Pipe(length=50, diameter=0.5, n_cells=10))
    c.connect(p1.outlet, p2.inlet)
    c.pressure_boundary(p1.inlet, p=250e3, T=300)
    c.pressure_boundary(p2.outlet, p=200e3, T=300)
    assert c.solve_steady_state(relaxation=0.5, max_outer_iterations=800).converged

    ref = th.Model(fluid=th.Fluid("helium"))
    ref.add_pipe("a", "b", length=100, diameter=0.5, n_cells=20)
    ref.pressure_boundary("a", p=250e3, T=300)
    ref.pressure_boundary("b", p=200e3, T=300)
    ref.steady_state(relaxation=0.5, max_outer_iterations=800)

    assert c.flow(p1) == pytest.approx(c.flow(p2))                 # continuity across the join
    assert c.flow(p1) == pytest.approx(ref.flow_through("a->b"), rel=1e-6)


def test_pump_drives_closed_loop():
    c = th.Circuit(fluid=th.Fluid("water"))
    pump = c.add(th.Pump(head_shutoff=150e3, curve=400.0))
    pipe = c.add(th.Pipe(length=20, diameter=0.3, n_cells=5))
    c.connect(pump.outlet, pipe.inlet)
    c.connect(pipe.outlet, pump.inlet)                 # closed loop
    c.pressure_boundary(pump.inlet, p=200e3, T=300)    # pressure reference

    assert c.solve_steady_state(relaxation=0.6).converged
    # Operating point: pump head balances the loop friction -> head_shutoff = curve*mdot^2.
    assert c.flow(pump) == pytest.approx(math.sqrt(150e3 / 400.0), rel=1e-3)
    assert c.flow(pump) > 0.0


def test_pressure_accessor_at_ports():
    c = th.Circuit(fluid=th.Fluid("helium"))
    pipe = c.add(th.Pipe(length=10, diameter=0.5, n_cells=4))
    c.pressure_boundary(pipe.inlet, p=300e3, T=300)
    c.pressure_boundary(pipe.outlet, p=200e3, T=300)
    c.solve_steady_state(relaxation=0.5, max_outer_iterations=800)
    assert c.pressure(pipe.inlet) == pytest.approx(300e3)
    assert c.pressure(pipe.outlet) == pytest.approx(200e3)


def test_missing_fluid_raises():
    c = th.Circuit()
    pipe = c.add(th.Pipe(length=10, diameter=0.5))
    c.pressure_boundary(pipe.inlet, p=300e3, T=300)
    c.pressure_boundary(pipe.outlet, p=200e3, T=300)
    with pytest.raises(ValueError, match="no fluid"):
        c.solve_steady_state()


def test_multi_fluid_not_supported():
    c = th.Circuit()
    c.add(th.Pipe(length=10, diameter=0.5, fluid=th.helium()))
    c.add(th.Pipe(length=10, diameter=0.5, fluid=th.water()))
    with pytest.raises(NotImplementedError):
        c.compile()
