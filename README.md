# OpenTH

**An open-source, Python-first, reactor-grade 1D thermal-hydraulic system code.**

OpenTH aims to be *the OpenMC of reactor thermal-hydraulics*: a transparent, scriptable, testable,
and extensible thermal-fluid network solver for steady-state and transient reactor systems —
filling the gap between Flownex-style system modelling and RELAP-class reactor transient analysis,
with advanced reactors (gas-cooled, molten-salt) as the first target. See the
**[Vision Statement](docs/vision-statement.md)** for the concept and positioning, the
**[Roadmap & Product Strategy](docs/roadmap.md)** for the phased development plan and status,
and the **[Backlog](docs/backlog.md)** for the requested component/capability catalog.

## Vision (in brief)

- **Open and developer-grade.** Permissive license, public repo, version-controlled models,
  CI-tested public benchmarks, transparent equations — vs. controlled/licensed codes like RELAP.
- **Python is the moat.** A clean model-building API (à la OpenMC) over a high-performance solver
  core; coupling to OpenMC, MOOSE, Cardinal, NekRS for neutronics/CFD.
- **Phased physics.** v0.1 incompressible → v0.2 compressible gas → real-fluid properties → heat
  structures → controls/trips → point kinetics → (eventually) two-phase. *Avoid full two-fluid
  two-phase complexity early — "that is where codes go to die."*
- **MVP target.** A helium-cooled reactor loop: core heat input, blower, recuperator/heat
  exchanger, heat structures, controls, and a transient pump coastdown.
- **Validation first.** A public benchmark suite (analytic friction → Edwards blowdown → natural
  circulation → pump coastdown → PWR/gas/MSR loops → coupled neutronics-T/H) run on every release.

## Current implementation

The repository today is the **single-phase PCIM flow-network kernel** — roughly the v0.1–v0.2
foundation of the vision above. It is built on Greyvenstein's implicit **Pressure Correction
(PC / PCIM)** method: a segregated (SIMPLE-style) finite-volume solver for the continuity, momentum
and energy equations, handling **liquid and gas**, **isothermal and non-isothermal**, and **fast and
slow transients**, extensible to non-pipe components (valves, pumps/compressors, …).

> Reference: G. P. Greyvenstein, *"An implicit method for the analysis of transient flows in
> pipe networks"*, Int. J. Numer. Meth. Engng 2002; 53:1127–1143
> (see [`docs/papers/`](docs/papers/), [`docs/theory.md`](docs/theory.md)).

## Status

Early prototype. The **steady-state and transient solvers work and are validated**, including
**non-isothermal flow**: steady reproduces the closed-form isothermal compressible pipe-flow relation
to ~0% error (`examples/pipeline_steady.py`); the transient marches to that steady state, conserves
mass, and captures water-hammer and the blow-down decay (`examples/blowdown_transient.py`); the energy
equation conserves total enthalpy in adiabatic flow and matches `Q/ṁ` for heat addition
(`examples/heated_pipe.py`, `tests/`). A compressible pressure-correction term keeps the solve robust
up to ~Mach 0.74 (the isothermal choking limit). Components so far: pipes, valves, and
pumps/compressors (`examples/pump_loop.py`); more are next.
Python-first for clarity and rapid iteration; performance-critical kernels will later be reimplemented
in C/C++ behind the same interfaces — see [`native/`](native/). A two-way LLM interface (build/query
networks in natural language) is planned in [`src/openth/llm/`](src/openth/llm/).

## Quick start

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e ".[dev,llm]"

pytest                              # run tests
openth benchmark                  # list the paper's Section 5 test cases
openth benchmark steady_pipeline  # generate + run one
```

> Requires Python ≥ 3.10. (The system Python on this machine is 3.9.6 — create a venv with a newer
> interpreter, e.g. via `pyenv` or `python3.12 -m venv .venv`.)

**In a notebook:** `pip install -e ".[notebook]"` (adds JupyterLab + matplotlib), then
`jupyter lab` and open [`examples/quickstart.ipynb`](examples/quickstart.ipynb) — a steady
solve, a transient plot, the `Circuit` API, and the benchmarks.

Build a model with the high-level API (`th.Model` is the builder; add elements by name,
set boundaries, then solve):

```python
import openth as th

model = th.Model(fluid=th.Fluid("helium"))
model.add_pipe("inlet", "outlet", length=100, diameter=0.5, n_cells=20)  # 20 FV cells
model.pressure_boundary("outlet", p=200e3, T=300)
model.mass_flow_boundary("inlet", mdot=30, T=300)

model.steady_state()                       # or: model.run(dt, duration, record=...)
print(model.pressure("inlet"), model.flow_through("inlet->outlet"))
```

Or build it the **component/connection** way with `th.Circuit` — add components and join
their `.inlet`/`.outlet` ports (nice for loops):

```python
import openth as th

c = th.Circuit(fluid=th.Fluid("helium"))
pump = c.add(th.Pump(head_shutoff=120e3, curve=300.0))
hot  = c.add(th.Pipe(length=15, diameter=0.3, n_cells=6))
cold = c.add(th.Pipe(length=15, diameter=0.3, n_cells=6))
c.connect(pump.outlet, hot.inlet); c.connect(hot.outlet, cold.inlet); c.connect(cold.outlet, pump.inlet)
c.pressure_boundary(pump.inlet, p=300e3, T=300)   # closed loop: pin a reference pressure

c.solve_steady_state()
print(c.flow(pump))
```

The paper's Section 5 cases live in `openth.benchmarks` (and the `openth benchmark` CLI).

**Full instructions are in the [User Guide](docs/user-guide.md)** — building models, steady
vs. transient runs, fluids, the energy equation, and generating/running the paper's tests.

## Layout

| Path | Purpose |
|------|---------|
| `src/openth/network/` | Topology: `Node` (control volume / cell centre) + `Element` (face / branch) + `Network` graph |
| `src/openth/components/` | Concrete elements: `Pipe`, `Valve`, `Pump`/compressor, boundaries; non-pipe components added here |
| `src/openth/fluids/` | Equations of state / fluid property models (`IdealGas`, incompressible) |
| `src/openth/solver/` | The PCIM solver, time integration, linear solves |
| `src/openth/model.py` | `FlowModel` — high-level model-building facade |
| `src/openth/benchmarks.py` | The paper's Section 5 test cases (`openth benchmark`) |
| `src/openth/io/` | Network (de)serialization |
| `src/openth/llm/` | Two-way LLM interface (optional) |
| `examples/` | Runnable cases, including paper benchmarks |
| `docs/` | Theory notes and the source papers |
| `native/` | C/C++ acceleration (stub for now) |

## License

MIT — see [LICENSE](LICENSE).
