# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

**OpenTH** (renamed from the prototype "FlowCalc") aims to be an open-source, Python-first,
reactor-grade 1D thermal-hydraulic system code — the OpenMC of reactor thermal hydraulics.
See `docs/vision-statement.md` for the founding concept and roadmap.

The repository today is the **single-phase flow-network kernel** (the v0.1–v0.2 foundation
of that vision): a solver implementing Greyvenstein's implicit **Pressure-Correction
(PC / PCIM)** method (Int. J. Numer. Meth. Engng 2002; 53:1127–1143; PDF in `docs/papers/`,
distilled in `docs/theory.md`). Python-first prototype; performance kernels are earmarked
for later C/C++ reimplementation (`native/`), and a two-way LLM interface is planned
(`src/openth/llm/`). The Python package imports as `openth`; end-user docs (the `FlowModel`
API, running the benchmarks) are in `docs/user-guide.md`.

**Project status.** The **steady-state and transient solvers are implemented and
validated**, including **non-isothermal flow** (energy-equation coupling). Steady matches
the closed-form isothermal compressible pipe-flow law to ~0% error; the transient marches
to the steady fixed point, conserves mass exactly, and reproduces water-hammer and the
blow-down decay; the energy solve conserves total enthalpy in adiabatic flow (~1e-9) and
matches `Q/ṁ` for heat addition. See `examples/` and `tests/`. With the compressible
pressure-correction term (eqs. 25-27) the solver is robust up to ~Mach 0.74 — the isothermal
choking limit `1/√γ`, beyond which no steady subsonic solution exists. More non-pipe
component models are the main remaining work — see "Implementing the solver core".

## Commands

```bash
# Environment (needs Python >= 3.10; system python here is 3.9.6 — use pyenv/python3.12)
python3.12 -m venv .venv && source .venv/bin/activate
pip install -e ".[dev,llm]"

python -m pytest                        # all tests (prefer `python -m pytest`, see note)
python -m pytest tests/test_steady.py   # one file
python -m pytest tests/test_steady.py::test_incompressible_series_exact   # one test
python -m pytest --cov=openth         # with coverage

ruff check . && ruff format .           # lint + format
mypy                                    # type-check (config in pyproject.toml)

python examples/pipeline_steady.py      # paper §5.1 steady benchmark: PCIM vs analytical
python examples/blowdown_transient.py   # paper §5.4 transient blow-down vs quasi-steady
python examples/heated_pipe.py          # non-isothermal: expansion cooling vs heat addition
python examples/pump_loop.py            # pump pushing flow uphill + compressor temperature rise
openth benchmark                      # list the paper's Section 5 test cases
openth benchmark blowdown             # generate + run one of them
openth --version                      # console entry point (openth.cli)
```

A repo-root `conftest.py` puts `src/` on the path so `python -m pytest` works even without
the editable install. If a bare `import openth` / `openth`/`pytest` console script
can't find the package (editable `.pth` not picked up), prefix with `PYTHONPATH=src`
(and `MYPYPATH=src` for mypy) or re-run `pip install -e .`.

## Architecture

The code mirrors the method's finite-volume structure (`docs/theory.md` maps every module
to the paper's equations — read it before touching the numerics):

- **`network/`** — topology. `Node` is a finite-volume **control volume / cell centre**
  carrying the scalar state (`p₀, ρ, T, h₀`) plus an optional `mass_source`; `Element` is a
  **face / branch** carrying mass flow `mdot`. `Network` owns the graph. The
  cell-centre-vs-face split is the paper's discretization (Fig. 1) and is load-bearing.
- **`components/`** — concrete `Element`s. `Pipe` is canonical; non-pipe components
  (`Valve`, `Pump`/compressor, and future turbine/orifice) plug in by overriding the closure
  hooks: `resistance` / `convective_dp` / `inertance` / `flow_area`, plus `head` (a
  pump/compressor pressure rise added to the momentum drive) and `work_per_mass` (shaft work
  raising total enthalpy in the energy solve). `PressureBoundary` pins a node's pressure (and
  temperature) and marks it `is_boundary` (removed from unknowns); `MassFlowBoundary` instead
  sets the node's `mass_source` (and optionally a fixed inlet `T`) and leaves the pressure an
  unknown. A `HeatExchanger` (in `network/coupling.py`, `Network.heat_exchangers`) is **not**
  an element — it's a thermal coupling between two nodes (heat `UA·ΔT`, hot→cold), applied in
  the energy solve only.
- **`fluids/`** — equations of state behind the `FluidModel` ABC. `IdealGas` implements
  `p = sρRT` (the paper's gas variant; `helium()` is the benchmark fluid); `Incompressible`
  the liquid variant. Each also provides `drho_dp` (compressibility) and
  `temperature_from_enthalpy` for the transient/energy solves.
- **`solver/`** — `PCIMSolver` (segregated SIMPLE-style loops), `SolverConfig` (notably
  `alpha ∈ [0.5,1]` for the transient, `relaxation ∈ (0,1]`, `solve_energy` for the
  non-isothermal energy solve, and `gravity` [m/s²] for the buoyancy term), and `linear.py`
  (Thomas for series pipelines, scipy sparse for networks — the solver currently uses sparse).
  Gravity enters the momentum drive as `g·ρ_face·(z_up−z_down)` via `Node.elevation` (set by
  `FlowModel.add_pipe(delta_elevation=…)`); inert when elevations are 0. With the energy
  solve on, temperature-driven density differences make this the **buoyancy** that drives
  natural circulation (which bootstraps transiently, not from a cold steady start).
- **`model.py`** — `FlowModel`, the high-level builder/facade most user code and the
  benchmarks go through: name-based nodes (auto-created), `add_pipe(..., n_cells=N)` that
  subdivides into finite-volume cells and assigns control volumes, callable (time-varying)
  pressure boundaries, and `steady_state()` / `run(dt, duration, record=...)`.
- **`circuit.py`** — the **port/connection layer** (`Circuit`): component objects with
  `.inlet`/`.outlet` ports (`circuit.Pipe`/`Valve`/`Pump`), `connect()` to join ports, and
  boundaries on ports. `Circuit.compile()` resolves connected ports into nodes (union-find)
  and emits into a `FlowModel`, so all physics is the validated `FlowModel`/solver path.
  Single-fluid for now. This is the vision's `loop = th.Model()` connection style.
- **`__init__.py`** — re-exports the common API for `import openth as th`: `Model`
  (alias of name-based `FlowModel`), `Circuit` (port-based), `Fluid(name)` factory,
  `helium`/`air`/`water`, the **port-based** `Pipe`/`Valve`/`Pump` (from `circuit`; the
  Element-level ones stay in `openth.components`), and `PCIMSolver`/`SolverConfig`.
  `__init__` must not import `llm` (optional `[llm]` extra).
- **`benchmarks.py`** — the paper's Section 5 cases built on `FlowModel` (`build_*`/`run_*`
  + a `BENCHMARKS` registry); exposed via `openth benchmark`.
- **`io/`** — JSON model save/load. `FlowModel` records its build calls as **directives**
  (`_directives`); `io.serialize` writes them + the fluid (high-level, hand-editable, not the
  expanded mesh) and replays them on load. `model.save`/`Model.load`/`to_dict`/`from_dict`,
  and `openth run model.json`. This declarative payload is also the LLM exchange format.
  Callables (time-varying boundaries/valve schedules) aren't serializable.
- **`llm/`** — optional two-way LLM interface. **The core must never import this**; the
  `anthropic` SDK lives behind the `[llm]` extra.

### The central abstraction

Each element exposes its momentum closure as **`resistance(rho_face) -> K`** (friction, so
that `Δp = K·mdot·|mdot|`) plus optional **`convective_dp(...)`** (momentum-flux term) and
**`head()`** (a pump/compressor pressure rise added to the momentum drive). The solver
linearizes friction into the pressure-correction conductance `d = 1/(2K|mdot|)` (paper
eq. 20); `head` is a constant source so it doesn't change `d`. These hooks are the single
interface that lets the same pressure-correction machinery handle pipes and arbitrary
non-pipe components uniformly — preserve them when adding components. (Carrying `convective_dp` explicitly is equivalent to
the paper's elimination of convective acceleration via total pressure; it's what makes the
high-Mach pressure ratio correct — see `docs/theory.md`.)

## Implementing the solver core

`PCIMSolver.steady_state` and `PCIMSolver.step` are both implemented (segregated SIMPLE;
steady is the dt→∞ limit of Section 4, transient adds the storage/inertia terms and the
`alpha` time-weighting), and the **energy equation** (eqs. 29-34) is coupled when
`SolverConfig.solve_energy` is set. All three are validated (`tests/`). The remaining work:

- **More non-pipe components** (`Pipe`/`Valve`/`Pump`/`HeatExchanger` exist; turbine/orifice
  and the rest of the [backlog](docs/backlog.md) catalog to come), each via the closure hooks
  (`resistance`/`convective_dp`/`inertance`/`flow_area`/`head`/`work_per_mass`). A turbine is
  a `Pump` with negative head / negative work; an orifice is a `Pipe`-like pure resistance.
- Wall heat transfer (temperature-dependent `q̇`) and heat structures; currently
  `Node.heat_source` is constant and `HeatExchanger` is a lumped UA·ΔT node coupling.
- **Transonic / near-choking robustness**: the solver is solid up to ~Mach 0.74 (the
  isothermal choking limit), but near choking it is mesh-sensitive — fine meshes can still
  collapse to the trivial zero-flow state (which then mass-balances to a false "converged").
  Genuine transonic handling (choked-flow boundary treatment) is future work.

Implementation notes worth preserving:
- The continuity equation carries a **compressible upwind convection of `p'`** — the density
  change `ρ'=∂ρ/∂p·p'` of the convected mass (eqs. 25-27; coefficient `κ=(∂ρ/∂p)/ρ=1/p` for
  an ideal gas). It scales with the mass flow, so it vanishes at low Mach (recovering
  incompressible SIMPLE) and dominates at high Mach, where it is what makes the solve
  converge. Without it the pressure-driven solve diverges above ~Mach 0.4.
- The flow update is **under-relaxed** (`SolverConfig.relaxation`): the convective term
  couples flow to itself and diverges at high Mach in pressure-driven problems.
- A network with **no interior pressure unknowns** (e.g. one element between two pressure
  boundaries — the blow-down case) is handled by `_converge_flows_only`: momentum alone
  fixes the flows.
- The transient is **mass-conservative by construction** (the storage term equals the
  θ-weighted flux integral); `tests/test_transient.py` checks this to ~1e-9.
- **Energy is solved on a converged (mass-conserving) flow field** — the energy outer loop
  alternates a *full* pressure/flow solve with one `h₀` update. Solving it on a
  non-conserving field makes the convective diagonal (nodal outflow) mismatch the inflow
  RHS and `h₀` blows up. The `h₀→T` step (via `½V²`) is also change-limited per iteration
  because the T↔ρ↔V kinetic coupling can otherwise run away (low ρ → high V → negative T).
- The pressure update is **floored positive** (`_apply_pressure_correction` caps the drop to
  ≥10% of the current `p`): a mass-flow inlet into a low-resistance pipe can otherwise
  overshoot to negative `p` → ρ≤0 → singular resistance → NaN. Needed for stiff and
  **hydraulically-disconnected** networks (e.g. the two independent streams of a heat
  exchanger, coupled only thermally).
- A **`HeatExchanger`** adds an implicit conductive link `g=UA/cp` between its two nodes in
  the energy matrix (to the RHS for a fixed-T side) — stable for any UA, energy-conserving.
- **Choking guard**: past the choking limit (`M=1/√γ` isothermal) there is no steady subsonic
  solution and the solve collapses to `ṁ=0`, which *trivially* mass-balances → a false
  "converged". `_guard_steady` downgrades such a result to `converged=False` by checking the
  flows actually satisfy the momentum balance (`_momentum_residual`) and that pressures are
  finite; a genuine no-flow state (hydrostatic, drive balanced by gravity) passes. `max_mach()`
  / the `"mach"` record key report proximity to choking. Simulating choked flow is future work.

Transcribe equations directly from `docs/papers/` (do **not** rely on memory or OCR), and
**validate against the paper's benchmarks before trusting results**. Add each validated
benchmark as a test (e.g. `tests/test_steady.py`, `test_energy.py`, `test_heat_exchanger.py`).

## Conventions

- SI units throughout (Pa, kg/m³, K, J/kg, m, m³/s, s). State this when adding quantities.
- Work in **total** pressure `p₀` and **total** enthalpy `h₀` (how the paper eliminates the
  convective-acceleration term) — not static `p`/`h`.
- Keep a pure-Python reference even after native kernels land; select backend at runtime so
  tests can exercise both and cross-check them.
