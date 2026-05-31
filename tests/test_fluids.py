"""Tests for the fluid property models (the verified, implemented pieces)."""

import math

from openth.fluids import helium, water


def test_ideal_gas_equation_of_state():
    he = helium()
    p, T = 700e3, 300.0
    rho = he.density(p, T)
    # eq. (12): p = s rho R T  ->  round-trips back to p.
    assert math.isclose(rho * he.s * he.R * T, p, rel_tol=1e-12)


def test_helium_specific_heats_consistent():
    he = helium()
    # cp - cv = R, and gamma = cp/cv.
    assert math.isclose(he.cp - he.cv, he.R, rel_tol=1e-12)
    assert math.isclose(he.cp / he.cv, he.gamma, rel_tol=1e-12)


def test_helium_sonic_velocity_order_of_magnitude():
    # Helium at 300 K: ~ 1000 m/s.
    a = helium().sonic_velocity(700e3, 300.0)
    assert 900.0 < a < 1100.0


def test_total_enthalpy_adds_kinetic_energy():
    he = helium()
    T, V = 300.0, 50.0
    assert math.isclose(he.total_enthalpy(T, V), he.enthalpy(T) + 0.5 * V * V)


def test_incompressible_density_constant():
    w = water()
    assert w.density(1e5, 300.0) == w.density(5e5, 350.0) == w.rho
