# OpenTH User Guide

OpenTH solves **fluid flow networks** — steady and transient, liquid and gas, isothermal
and non-isothermal — using Greyvenstein's implicit Pressure-Correction (PCIM) method. This
guide covers building and running models through the high-level `FlowModel` API, and
generating/running the benchmark cases from the source paper.

For the underlying method and how it maps to the paper's equations, see
[`theory.md`](theory.md). For project status and limitations, see
[`project-status.md`](project-status.md).

## Contents

1. [Install](#1-install)
2. [Core concepts](#2-core-concepts)
3. [Building a model](#3-building-a-model)
4. [Boundary conditions](#4-boundary-conditions)
5. [Solving: steady state and transient](#5-solving-steady-state-and-transient)
6. [Reading results](#6-reading-results)
7. [Fluids](#7-fluids)
8. [Non-isothermal flow (the energy equation)](#8-non-isothermal-flow-the-energy-equation)
9. [The paper's test models (benchmarks)](#9-the-papers-test-models-benchmarks)
10. [Worked example: build a paper case from scratch](#10-worked-example-build-a-paper-case-from-scratch)
11. [Limitations](#11-limitations)

---

## 1. Install

Requires Python ≥ 3.10.

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"          # add ",llm" for the optional LLM interface
```

Quick check:

```bash
openth benchmark               # lists the paper's test cases
python -m pytest                 # runs the test suite
```

## 2. Core concepts

OpenTH uses a finite-volume discretization (the paper's Fig. 1):

- **Nodes** are control-volume cell centres. They carry the scalar state — total pressure
  `p₀`, density `ρ`, temperature `T`, total enthalpy `h₀` — and a control volume.
- **Elements** are faces/branches connecting two nodes. They carry the **mass flow** `ṁ`
  and define the momentum closure (a pipe's friction, a pump's head, …).
- **Boundaries** fix conditions: a `PressureBoundary` pins a node's pressure (and
  temperature); a `MassFlowBoundary` imposes a mass flow and leaves the pressure to be
  solved.

**Units are SI throughout**: pressure Pa, density kg/m³, temperature K, enthalpy J/kg,
length m, area m², mass flow kg/s, time s, power W.

A real pipe is *subdivided* into several cells; `FlowModel.add_pipe(..., n_cells=N)` does
this for you, creating the internal nodes and assigning control volumes.

## 3. Building a model

Everything goes through the high-level model, exposed at the top level as `th.Model` (an
alias of `openth.model.FlowModel`). Nodes are created automatically the first time you name
them.

```python
import openth as th

model = th.Model(fluid=th.Fluid("helium"))   # one working fluid per model
```

`th.Fluid("helium")` is a convenience factory for the built-ins; you can equally pass
`th.helium()` or a custom `th.IdealGas(...)` / `th.Incompressible(...)` (see
[Fluids](#7-fluids)). The lower-level imports (`from openth.model import FlowModel`,
`from openth.fluids import helium`) still work — `import openth as th` just collects the
common names in one place.

### Pipes

```python
model.add_pipe(
    "inlet", "outlet",       # upstream node, downstream node (created on first use)
    length=100.0,            # m
    diameter=0.5,            # m
    friction_factor=0.02,    # Darcy f (default 0.02)
    n_cells=20,              # subdivide into 20 finite-volume cells (default 1)
    name="pipe",             # optional; defaults to "inlet->outlet"
)
```

`n_cells=20` builds 20 pipe segments and 19 internal nodes between `inlet` and `outlet`
(21 nodes, 20 elements total). More cells = finer spatial resolution.

### Valves

A valve is a variable resistance; an `opening(t)` schedule in `(0, 1]` lets it close in time.

```python
model.add_valve(
    "pipe_end", "outlet",
    k_open=300.0,                                   # loss coefficient when fully open
    opening=lambda t: 0.0 if t >= 0.05 else 1.0,    # shut at t = 0.05 s
    name="valve",
)
```

### Pumps and compressors

A pump/compressor adds a pressure rise following a quadratic curve
`rise = head_shutoff - curve·ṁ²`. For a compressor the shaft work also heats the gas (only
relevant when the energy equation is on).

```python
model.add_pump(
    "low", "high",
    head_shutoff=250e3,      # Pa at zero flow
    curve=500.0,             # Pa/(kg/s)^2 fall-off
    efficiency=0.8,          # affects the compressor temperature rise (default 1.0)
    name="pump",
)
```

### Heat input

```python
model.heat_source("mid", power=5e5)   # 5e5 W into this node (energy equation only)
```

## 4. Boundary conditions

A network needs at least one `pressure_boundary` (otherwise the pressure level is undefined).

```python
# Fixed pressure (Dirichlet) — also fixes temperature at that node:
model.pressure_boundary("outlet", p=200e3, T=300.0)

# Imposed mass flow (+ into the node, - out of it); pressure stays an unknown.
# Pass T to also fix the inlet temperature (for the energy equation):
model.mass_flow_boundary("inlet", mdot=30.0, T=300.0)

# Time-varying pressure: pass a callable of time t [s].
import math
model.pressure_boundary("up", p=lambda t: (650 + 50 * math.exp(-0.004 * t)) * 1e3, T=300.0)
```

Time-varying pressure boundaries are re-evaluated automatically at each step of `run()`.

## 5. Solving: steady state and transient

### Steady state

```python
result = model.steady_state(relaxation=0.5)   # keyword args go to SolverConfig
print(result.converged, result.iterations, result.residual)
```

### Transient

`run(dt, duration, ...)` solves the steady state first (as the initial condition, unless
`steady_init=False`), then marches in time. `record` is a tuple of history keys; the return
value maps each key (plus `"t"`) to a list of values.

```python
history = model.run(
    dt=0.0009, duration=0.5,
    record=("p:valve_in", "flow:pipe"),   # see "Reading results" for key formats
    alpha=0.6, relaxation=0.6,            # SolverConfig options
)
times = history["t"]
pressures = history["p:valve_in"]
```

### `SolverConfig` options

Pass these as keyword arguments to `steady_state()` / `run()`:

| Option | Default | Meaning |
|--------|---------|---------|
| `alpha` | `0.6` | Time-integration weighing in `[0.5, 1]` (transient only). 1 = fully implicit; 0.5 = Crank–Nicolson; 0.6 is the recommended compromise. |
| `relaxation` | `0.7` | Under-relaxation in `(0, 1]`. Lower = more stable, slower. Use ~0.4–0.5 for higher-Mach or non-isothermal cases. |
| `solve_energy` | `False` | Solve the energy equation for temperature (non-isothermal). |
| `max_outer_iterations` | `200` | Iteration cap per solve/step. Raise for stiff/high-Mach cases. |
| `tol` | `1e-8` | Convergence tolerance on the (relative) mass-imbalance residual. |

## 6. Reading results

After a solve, query the committed state by node/element name:

```python
model.pressure("inlet")        # total pressure [Pa] at a node
model.temperature("inlet")     # temperature [K]
model.flow("valve")            # mass flow [kg/s] through a single element
model.flow_through("pipe")     # mean mass flow through a (subdivided) pipe by its name
```

For `run()`, the `record` keys are `"<kind>:<id>"`:

| Key | Meaning |
|-----|---------|
| `p:NODE` | total pressure at a node |
| `T:NODE` | temperature at a node |
| `mdot:ELEMENT` | mass flow through one element (e.g. a valve, or one pipe segment) |
| `flow:PIPENAME` | mean mass flow through a whole subdivided pipe |

Internal cell nodes of a subdivided pipe are named `"<pipe>#c0"`, `"<pipe>#c1"`, …, and the
segments `"<pipe>#s0"`, …; record those if you need a specific interior location. Named
endpoints and junctions are usually what you want.

## 7. Fluids

```python
import openth as th

th.Fluid("helium")   # factory for built-ins: "helium" (or "he"), "air", "water"
th.helium()          # ideal gas, R=2077, gamma=1.667  (same as Fluid("helium"))
th.air()             # ideal gas, R=287, gamma=1.4
th.water()           # incompressible, rho=998

# Custom gas / liquid:
th.IdealGas(name="co2", R=189.0, gamma=1.29)
th.Incompressible(name="oil", rho=850.0, cp=1900.0)
```

One fluid per model.

## 8. Non-isothermal flow (the energy equation)

By default the flow is isothermal (temperatures stay at their boundary/initial values).
Enable the energy equation with `solve_energy=True`:

```python
model.steady_state(solve_energy=True, relaxation=0.4)
```

With it on:

- **Adiabatic** flow conserves total enthalpy `h₀`, so a gas *cools* as it accelerates down
  a pressure gradient (expansion cooling).
- A **`heat_source`** raises downstream total enthalpy by `Q/ṁ`.
- A **compressor** raises temperature by its shaft work.

Inflow boundaries should fix temperature (`pressure_boundary` does; for a
`mass_flow_boundary` inlet, pass `T=`). Use `relaxation≈0.4` for non-isothermal cases.

## 9. The paper's test models (benchmarks)

`openth.benchmarks` generates and runs the four Section 5 cases from Greyvenstein (2002).
Each has a `build_*` function (returns a ready-to-solve `FlowModel`) and a `run_*` function
(solves it and returns a result summary), collected in the `BENCHMARKS` registry.

| Name | Section | What it shows |
|------|---------|---------------|
| `steady_pipeline` | 5.1 | Steady isothermal 100 m helium pipeline; matches the closed-form pipe law. |
| `valve_closure` | 5.2 | Sudden valve closure in a 20 m pipe → water-hammer overpressure. |
| `branching_network` | 5.3 | Valve closures in a five-pipe branching network. |
| `blowdown` | 5.4 | Pressure-vessel blow-down; mass flow decays to zero (slow transient). |

### From the command line

```bash
openth benchmark                 # list all four
openth benchmark steady_pipeline # generate + run one, print the summary
openth benchmark blowdown
```

### From Python

```python
from openth import benchmarks

# Run by name (uses the paper's parameters):
summary = benchmarks.run("valve_closure")
print(summary)   # {'pre-closure p_valve [kPa]': 698.4, 'peak p_valve [kPa]': 766.2, ...}

# Or build the model and drive it yourself for full control over outputs:
model = benchmarks.build_blowdown()
history = model.run(dt=10.0, duration=1750.0, record=("flow:pipe",),
                    alpha=0.6, relaxation=0.6)
mdot_t = history["flow:pipe"]      # mass-flow decay over the 30-minute blow-down

# Run functions also take parameters (e.g. a target Mach, or a shorter transient):
benchmarks.run_steady_pipeline(mach=0.5)
benchmarks.run_valve_closure(t_close=0.01, dt=0.0009, duration=0.5)
```

The matching `examples/` scripts (`pipeline_steady.py`, `blowdown_transient.py`,
`heated_pipe.py`, `pump_loop.py`) show the same cases built by hand with extra commentary
and comparisons to analytical solutions.

## 10. Worked example: build a paper case from scratch

Reproducing the Section 5.4 blow-down with the high-level API:

```python
import math
from openth.model import FlowModel
from openth.fluids import helium

model = FlowModel(fluid=helium(), default_temperature=300.0)

# 10 m x 0.1 m pipe, f = 0.02, 40 finite-volume cells, between two vessels.
model.add_pipe("up", "down", length=10.0, diameter=0.1,
               friction_factor=0.02, n_cells=40, name="pipe")

# Upstream pressure decays; downstream held constant.
model.pressure_boundary("up", p=lambda t: (650 + 50 * math.exp(-0.004 * t)) * 1e3, T=300.0)
model.pressure_boundary("down", p=650e3, T=300.0)

# Slow transient: 10 s steps over 30 minutes.
history = model.run(dt=10.0, duration=1750.0, record=("flow:pipe",),
                    alpha=0.6, relaxation=0.6)

for t, mdot in zip(history["t"], history["flow:pipe"]):
    if int(t) % 250 == 0:
        print(f"t = {t:6.0f} s   mdot = {mdot:.4f} kg/s")
```

## 11. Limitations

- **Near-choking flow**: pressure-driven cases are robust up to ~Mach 0.74 (the isothermal
  choking limit). Above that no steady subsonic solution exists; near it the solve can be
  mesh-sensitive — prefer coarser meshes / lower `relaxation`, and sanity-check that
  `result.converged` is `True` and the flow is non-zero.
- **Isothermal by default**: set `solve_energy=True` for temperature transport. Wall heat
  transfer is not modelled (`heat_source` is a fixed power, not temperature-dependent).
- **Component set**: pipes, valves, pumps/compressors. Turbines (negative-head pump),
  orifices (pure resistance), and heat exchangers are not yet provided as named components
  but follow the same element interface.
- **Performance**: pure-Python prototype; large networks or long fast transients will be
  slow until the planned C/C++ kernels land.
