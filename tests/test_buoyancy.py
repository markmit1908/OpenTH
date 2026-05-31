"""Validation tests for the gravity / buoyancy term.

  * a static liquid column reproduces the hydrostatic pressure p = p_top + rho*g*H;
  * a heated closed loop develops buoyancy-driven natural circulation (up the heated leg),
    which vanishes when gravity is switched off.
"""

import warnings

import pytest

import openth as th

G = 9.80665


def test_hydrostatic_column():
    """Vertical incompressible column, dead-ended at the bottom: no flow, and the bottom
    pressure is the top pressure plus the hydrostatic head rho*g*H."""
    water = th.water()
    height = 10.0
    m = th.Model(fluid=water)
    m.add_pipe("top", "bot", length=height, diameter=0.1, n_cells=4,
               delta_elevation=-height, name="col")
    m.node("bot").volume = 0.0
    m.pressure_boundary("top", p=100e3, T=300.0)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # zero-flow case won't hit the mass-residual tol
        m.steady_state(relaxation=0.7)

    assert abs(m.flow_through("col")) < 1e-4                       # no flow
    assert m.pressure("bot") == pytest.approx(100e3 + water.rho * G * height, rel=1e-6)


def _loop(heat: float):
    """Square gas loop BL(0)->BR(0)->TR(H)->TL(H)->BL, heated at the bottom of the riser,
    with TL pinned as the pressure + cold reference."""
    h = 3.0
    m = th.Model(fluid=th.air())
    m.add_pipe("BL", "BR", length=3.0, diameter=0.1, n_cells=2, name="bottom")
    m.add_pipe("BR", "TR", length=h, diameter=0.1, n_cells=2, delta_elevation=+h, name="riser")
    m.add_pipe("TR", "TL", length=3.0, diameter=0.1, n_cells=2, name="top")
    m.add_pipe("TL", "BL", length=h, diameter=0.1, n_cells=2, delta_elevation=-h, name="downcomer")
    m.heat_source("BR", power=heat)
    m.pressure_boundary("TL", p=200e3, T=300.0)
    return m


def _circulation(heat: float, gravity: float = G) -> float:
    m = _loop(heat)
    hist = m.run(dt=0.2, duration=8.0, record=("flow:riser",), alpha=0.6, relaxation=0.4,
                 solve_energy=True, steady_init=False, gravity=gravity)
    return hist["flow:riser"][-1]


def test_natural_circulation_is_buoyancy_driven():
    with_gravity = _circulation(3000.0, gravity=G)
    without_gravity = _circulation(3000.0, gravity=0.0)
    assert with_gravity > 0.005            # circulates, up the heated riser (positive)
    assert without_gravity < 0.1 * with_gravity   # no buoyancy -> essentially no circulation


def test_circulation_increases_with_heating():
    assert _circulation(6000.0) > _circulation(1500.0) > 0.0
