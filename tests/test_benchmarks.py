"""Tests for the paper's Section 5 benchmark builders/runners.

Transient cases use short durations here to keep the suite fast (the qualitative behaviour
appears within the first cycle); the default durations in :mod:`flowcalc.benchmarks`
reproduce the paper's full runs.
"""

from flowcalc import benchmarks


def test_registry_covers_all_four_sections():
    assert set(benchmarks.BENCHMARKS) == {
        "steady_pipeline", "valve_closure", "branching_network", "blowdown",
    }


def test_steady_pipeline_matches_closed_form():
    res = benchmarks.run_steady_pipeline(mach=0.3)
    assert res["converged"]
    assert res["P1/P2 (PCIM)"] == res["P1/P2 (exact)"]


def test_valve_closure_produces_water_hammer():
    res = benchmarks.run_valve_closure(t_close=0.01, dt=0.0009, duration=0.06)
    assert res["water hammer"] is True
    assert res["peak p_valve [kPa]"] > res["pre-closure p_valve [kPa]"]


def test_branching_network_runs():
    res = benchmarks.run_branching_network(t_close=0.005, dt=0.0009, duration=0.03)
    assert res["steps"] > 0
    assert res["p_A peak [kPa]"] >= res["p_A min [kPa]"]


def test_blowdown_decays_monotonically():
    res = benchmarks.run_blowdown(dt=10.0, duration=300.0)
    assert res["monotonic decay"] is True
    assert res["final mdot [kg/s]"] < res["initial mdot [kg/s]"]


def test_run_unknown_benchmark_raises():
    import pytest
    with pytest.raises(KeyError):
        benchmarks.run("nonexistent")
