"""Tests for JSON model (de)serialization."""

import pytest

import openth as th


def _pipeline() -> th.Model:
    m = th.Model(fluid=th.Fluid("helium"))
    m.add_pipe("inlet", "outlet", length=100, diameter=0.5, n_cells=20, name="pipe")
    m.pressure_boundary("outlet", p=200e3, T=300)
    m.mass_flow_boundary("inlet", mdot=30, T=300)
    return m


def test_roundtrip_solves_identically():
    original = _pipeline()
    restored = th.Model.from_dict(original.to_dict())
    original.steady_state(relaxation=0.5, max_outer_iterations=600)
    restored.steady_state(relaxation=0.5, max_outer_iterations=600)
    assert restored.pressure("inlet") == pytest.approx(original.pressure("inlet"), rel=1e-9)


def test_save_load_file_roundtrip(tmp_path):
    path = str(tmp_path / "model.json")
    _pipeline().save(path)
    restored = th.Model.load(path)
    restored.steady_state(relaxation=0.5, max_outer_iterations=600)
    assert restored.pressure("inlet") > 200e3   # rebuilt and solvable


def test_roundtrip_preserves_full_structure():
    m = th.Model(fluid=th.IdealGas(name="co2", R=189.0, gamma=1.29))
    m.add_pipe("a", "b", length=5, diameter=0.3, n_cells=2, delta_elevation=2.0)
    m.add_pump("b", "c", head_shutoff=50e3, curve=100.0, name="pmp")
    m.set_volume("c", 4.0)
    m.set_initial("c", p=300e3, T=350)
    m.add_heat_exchanger("b", "c", UA=1000.0, name="hx")
    m.pressure_boundary("a", p=200e3, T=300)

    r = th.Model.from_dict(m.to_dict())
    assert sorted(r.network.nodes) == sorted(m.network.nodes)
    assert sorted(r.network.elements) == sorted(m.network.elements)
    assert sorted(r.network.heat_exchangers) == sorted(m.network.heat_exchangers)
    assert r.node("c").volume == 4.0
    assert r.node("c").state.p0 == 300e3
    assert r.node("b").elevation == pytest.approx(2.0)
    assert isinstance(r.fluid, th.IdealGas) and r.fluid.R == 189.0


def test_builtin_fluid_shorthand():
    data = {
        "schema": 1, "fluid": {"kind": "water"}, "default_temperature": 300.0,
        "build": [
            {"op": "add_pipe", "upstream": "a", "downstream": "b",
             "length": 10, "diameter": 0.2, "n_cells": 1, "name": "p"},
            {"op": "pressure_boundary", "node": "a", "p": 3e5, "T": 300},
            {"op": "pressure_boundary", "node": "b", "p": 2e5, "T": 300},
        ],
    }
    model = th.Model.from_dict(data)
    assert model.fluid.name == "water"
    model.steady_state()
    assert model.flow_through("p") > 0.0


def test_callable_boundary_is_not_serializable():
    m = th.Model(fluid=th.helium())
    m.add_pipe("x", "y", length=1, diameter=0.1)
    m.pressure_boundary("x", p=lambda t: 1e5 + t, T=300)
    with pytest.raises(ValueError, match="callable"):
        m.to_dict()


def test_unknown_schema_rejected():
    with pytest.raises(ValueError, match="schema"):
        th.Model.from_dict({"schema": 999, "fluid": {"kind": "helium"}, "build": []})
