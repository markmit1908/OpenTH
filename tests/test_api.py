"""The top-level ``import openth as th`` convenience API."""

import pytest

import openth as th


def test_top_level_names_are_exported():
    for name in ("Model", "FlowModel", "Fluid", "Pipe", "Valve", "Pump",
                 "PressureBoundary", "MassFlowBoundary", "PCIMSolver", "SolverConfig",
                 "helium", "air", "water", "IdealGas", "Incompressible"):
        assert hasattr(th, name), f"openth.{name} should be exported"
    assert th.Model is th.FlowModel


def test_fluid_factory():
    assert th.Fluid("helium").name == "helium"
    assert th.Fluid("HE ").name == "helium"     # case/space-insensitive alias
    assert th.Fluid("air").R == th.air().R
    with pytest.raises(ValueError):
        th.Fluid("unobtanium")


def test_build_and_solve_via_top_level_api():
    model = th.Model(fluid=th.Fluid("helium"))
    model.add_pipe("inlet", "outlet", length=100, diameter=0.5, n_cells=20)
    model.pressure_boundary("outlet", p=200e3, T=300)
    model.mass_flow_boundary("inlet", mdot=30, T=300)

    result = model.steady_state(relaxation=0.5, max_outer_iterations=600)
    assert result.converged
    assert model.pressure("inlet") > 200e3   # inlet pressure built up to drive the flow
