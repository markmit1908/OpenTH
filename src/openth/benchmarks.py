"""The paper's Section 5 test cases, built on :class:`~openth.model.FlowModel`.

Each case has a ``build_*`` function (returns a ready-to-solve model) and a ``run_*``
function (solves it and returns a concise result summary). :data:`BENCHMARKS` maps short
names to ``(description, run_fn)`` and :func:`run` dispatches by name -- used by the
``openth benchmark`` CLI command.

References are to G. P. Greyvenstein, Int. J. Numer. Meth. Engng 2002; 53:1127-1143.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from typing import Any

from .fluids import helium
from .model import FlowModel

T0 = 300.0  # all paper cases use helium at 300 K


def _area(diameter: float) -> float:
    return 0.25 * math.pi * diameter * diameter


def _closes_at(t_close: float) -> Callable[[float], float]:
    """Valve opening schedule: fully open until ``t_close``, then shut."""
    def opening(t: float) -> float:
        return 0.0 if t >= t_close else 1.0
    return opening


# --- Section 5.1: steady isothermal pipeline -------------------------------------------

def build_steady_pipeline(mach: float = 0.3) -> FlowModel:
    """100 m x 0.5 m helium pipeline, f=0.02, 20 increments; outlet 200 kPa, inlet flow
    sized to a target outlet Mach number (Fig. 2)."""
    fluid = helium()
    p_out = 200e3
    mdot = fluid.density(p_out, T0) * mach * fluid.sonic_velocity(p_out, T0) * _area(0.5)
    model = FlowModel(fluid=fluid)
    model.add_pipe("inlet", "outlet", length=100.0, diameter=0.5,
                   friction_factor=0.02, n_cells=20, name="pipe")
    model.pressure_boundary("outlet", p=p_out, T=T0)
    model.mass_flow_boundary("inlet", mdot=mdot, T=T0)
    return model


def run_steady_pipeline(mach: float = 0.3) -> dict[str, Any]:
    model = build_steady_pipeline(mach)
    result = model.steady_state(relaxation=0.5, max_outer_iterations=600)
    p_out = 200e3
    # Closed-form isothermal pipe law for comparison.
    g = model.flow_through("pipe") / _area(0.5)
    c = g * g * helium().R * T0
    p1 = p_out
    for _ in range(200):
        p1 = math.sqrt(p_out**2 + c * (0.02 * 100.0 / 0.5 + 2.0 * math.log(p1 / p_out)))
    return {
        "converged": result.converged,
        "mdot [kg/s]": round(model.flow_through("pipe"), 4),
        "P1/P2 (PCIM)": round(model.pressure("inlet") / p_out, 4),
        "P1/P2 (exact)": round(p1 / p_out, 4),
    }


# --- Section 5.2: sudden valve closure in a single pipe --------------------------------

def build_valve_closure(t_close: float = 0.01) -> FlowModel:
    """20 m x 0.5 m helium pipe, f=0.02, 20 increments; inlet 700 kPa, a valve at the
    downstream end shuts at ``t_close`` -> water-hammer pressure waves (Fig. 4)."""
    model = FlowModel(fluid=helium())
    model.add_pipe("inlet", "valve_in", length=20.0, diameter=0.5,
                   friction_factor=0.02, n_cells=20, name="pipe")
    model.add_valve("valve_in", "outlet", k_open=300.0, opening=_closes_at(t_close),
                    name="valve")
    model.pressure_boundary("inlet", p=700e3, T=T0)
    model.pressure_boundary("outlet", p=650e3, T=T0)
    return model


def run_valve_closure(t_close: float = 0.01, dt: float = 0.0009,
                      duration: float = 0.5) -> dict[str, Any]:
    model = build_valve_closure(t_close)
    hist = model.run(dt, duration, record=("p:valve_in",),
                     alpha=0.6, relaxation=0.6)
    pressures = hist["p:valve_in"]
    times = hist["t"]
    pre = [p for p, t in zip(pressures, times, strict=True) if t < t_close]
    post = [p for p, t in zip(pressures, times, strict=True) if t >= t_close]
    return {
        "pre-closure p_valve [kPa]": round((pre[-1] if pre else pressures[0]) / 1e3, 1),
        "peak p_valve [kPa]": round(max(post) / 1e3, 1),
        "min p_valve [kPa]": round(min(post) / 1e3, 1),
        "water hammer": max(post) > (pre[-1] if pre else pressures[0]),
    }


# --- Section 5.3: sudden valve closures in a branching network -------------------------

def build_branching_network(t_close: float = 0.01) -> FlowModel:
    """Five 10 m x 0.5 m pipes (10 increments each), inlet 700 kPa; branches 1 and 2 have
    valves that shut at ``t_close`` while branch 3 keeps a constant outflow (Fig. 8)."""
    model = FlowModel(fluid=helium())

    def pipe(upstream: str, downstream: str, name: str) -> None:
        model.add_pipe(upstream, downstream, length=10.0, diameter=0.5,
                       friction_factor=0.02, n_cells=10, name=name)

    pipe("inlet", "A", "pipe_in")   # inlet -> junction A
    pipe("A", "B", "pipe_AB")       # A -> junction B
    pipe("B", "out3", "pipe_3")     # B -> branch-3 outlet
    pipe("A", "v1_in", "pipe_1")    # A -> branch-1 (valved)
    pipe("B", "v2_in", "pipe_2")    # B -> branch-2 (valved)
    model.add_valve("v1_in", "out1", k_open=50.0, opening=_closes_at(t_close), name="valve1")
    model.add_valve("v2_in", "out2", k_open=50.0, opening=_closes_at(t_close), name="valve2")
    model.pressure_boundary("inlet", p=700e3, T=T0)
    model.pressure_boundary("out1", p=650e3, T=T0)
    model.pressure_boundary("out2", p=650e3, T=T0)
    model.mass_flow_boundary("out3", mdot=-11.61, T=T0)  # constant outflow from branch 3
    return model


def run_branching_network(t_close: float = 0.01, dt: float = 0.0009,
                          duration: float = 0.3) -> dict[str, Any]:
    model = build_branching_network(t_close)
    hist = model.run(dt, duration, record=("p:A",), alpha=0.6, relaxation=0.5)
    return {
        "p_A initial [kPa]": round(hist["p:A"][0] / 1e3, 1),
        "p_A peak [kPa]": round(max(hist["p:A"]) / 1e3, 1),
        "p_A min [kPa]": round(min(hist["p:A"]) / 1e3, 1),
        "steps": len(hist["t"]),
    }


# --- Section 5.4: blow-down of a pressure vessel through a pipe ------------------------

def build_blowdown() -> FlowModel:
    """10 m x 0.1 m helium pipe (40 increments) between two vessels; the upstream pressure
    decays as p1(t) = 650 + 50 exp(-0.004 t) kPa, downstream held at 650 kPa (Fig. 10)."""
    model = FlowModel(fluid=helium())
    model.add_pipe("up", "down", length=10.0, diameter=0.1,
                   friction_factor=0.02, n_cells=40, name="pipe")
    model.pressure_boundary("up", p=lambda t: (650.0 + 50.0 * math.exp(-0.004 * t)) * 1e3, T=T0)
    model.pressure_boundary("down", p=650e3, T=T0)
    return model


def run_blowdown(dt: float = 10.0, duration: float = 1750.0) -> dict[str, Any]:
    model = build_blowdown()
    hist = model.run(dt, duration, record=("flow:pipe",), alpha=0.6, relaxation=0.6)
    flows = hist["flow:pipe"]
    return {
        "initial mdot [kg/s]": round(flows[0], 4),
        "final mdot [kg/s]": round(flows[-1], 4),
        "monotonic decay": all(a >= b - 1e-9 for a, b in zip(flows, flows[1:], strict=False)),
        "steps": len(flows),
    }


BENCHMARKS: dict[str, tuple[str, Callable[[], dict[str, Any]]]] = {
    "steady_pipeline": ("Section 5.1 - steady isothermal pipeline", run_steady_pipeline),
    "valve_closure": ("Section 5.2 - sudden valve closure (water hammer)", run_valve_closure),
    "branching_network": ("Section 5.3 - branching network valve closures", run_branching_network),
    "blowdown": ("Section 5.4 - pressure-vessel blow-down", run_blowdown),
}


def run(name: str) -> dict[str, Any]:
    """Run a benchmark by name and return its result summary."""
    if name not in BENCHMARKS:
        raise KeyError(f"unknown benchmark {name!r}; choose from {sorted(BENCHMARKS)}")
    return BENCHMARKS[name][1]()
