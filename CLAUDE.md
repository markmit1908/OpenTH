# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

FlowCalc is a **fluid flow-network solver** implementing Greyvenstein's implicit
**Pressure-Correction (PC / PCIM)** method (Int. J. Numer. Meth. Engng 2002; 53:1127–1143;
PDF in `docs/papers/`, distilled in `docs/theory.md`). Python-first prototype; performance
kernels are earmarked for later C/C++ reimplementation (`native/`), and a two-way LLM
interface is planned (`src/flowcalc/llm/`).

**Project status.** The **steady-state and transient solvers are implemented and
validated**, including **non-isothermal flow** (energy-equation coupling). Steady matches
the closed-form isothermal compressible pipe-flow law to ~0% error; the transient marches
to the steady fixed point, conserves mass exactly, and reproduces water-hammer and the
blow-down decay; the energy solve conserves total enthalpy in adiabatic flow (~1e-9) and
matches `Q/ṁ` for heat addition. See `examples/` and `tests/`. The segregated
pressure↔temperature coupling is robust below ~Mach 0.4 for pressure-driven flow (as with
the isothermal pressure-driven solve); higher-Mach robustness and more non-pipe components
are the remaining work — see "Implementing the solver core".

## Commands

```bash
# Environment (needs Python >= 3.10; system python here is 3.9.6 — use pyenv/python3.12)
python3.12 -m venv .venv && source .venv/bin/activate
pip install -e ".[dev,llm]"

python -m pytest                        # all tests (prefer `python -m pytest`, see note)
python -m pytest tests/test_steady.py   # one file
python -m pytest tests/test_steady.py::test_incompressible_series_exact   # one test
python -m pytest --cov=flowcalc         # with coverage

ruff check . && ruff format .           # lint + format
mypy                                    # type-check (config in pyproject.toml)

python examples/pipeline_steady.py      # paper §5.1 steady benchmark: PCIM vs analytical
python examples/blowdown_transient.py   # paper §5.4 transient blow-down vs quasi-steady
python examples/heated_pipe.py          # non-isothermal: expansion cooling vs heat addition
flowcalc --version                      # console entry point (flowcalc.cli)
```

A repo-root `conftest.py` puts `src/` on the path so `python -m pytest` works even without
the editable install. If a bare `import flowcalc` / `flowcalc`/`pytest` console script
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
  (`Valve`, and future pump/compressor/turbine/orifice/heat-exchanger) plug in by
  overriding `resistance` / `convective_dp` / `inertance` / `flow_area`. `PressureBoundary`
  pins a node's pressure (and temperature) and marks it `is_boundary` (removed from
  unknowns); `MassFlowBoundary` instead sets the node's `mass_source` (and optionally a
  fixed inlet `T`) and leaves the pressure an unknown.
- **`fluids/`** — equations of state behind the `FluidModel` ABC. `IdealGas` implements
  `p = sρRT` (the paper's gas variant; `helium()` is the benchmark fluid); `Incompressible`
  the liquid variant. Each also provides `drho_dp` (compressibility) and
  `temperature_from_enthalpy` for the transient/energy solves.
- **`solver/`** — `PCIMSolver` (segregated SIMPLE-style loops), `SolverConfig` (notably
  `alpha ∈ [0.5,1]` for the transient, `relaxation ∈ (0,1]`, and `solve_energy` to enable
  the non-isothermal energy solve), and `linear.py` (Thomas for series pipelines, scipy
  sparse for networks — the solver currently uses sparse).
- **`io/`** — declarative dict/JSON (de)serialization; also the LLM exchange payload.
- **`llm/`** — optional two-way LLM interface. **The core must never import this**; the
  `anthropic` SDK lives behind the `[llm]` extra.

### The central abstraction

Each element exposes its momentum closure as **`resistance(rho_face) -> K`** (friction, so
that `Δp = K·mdot·|mdot|`) plus an optional **`convective_dp(...)`** (momentum-flux term).
The solver linearizes friction into the pressure-correction conductance `d = 1/(2K|mdot|)`
(paper eq. 20). This pair of methods is the single interface that lets the same
pressure-correction machinery handle pipes and arbitrary non-pipe components uniformly —
preserve it when adding components. (Carrying `convective_dp` explicitly is equivalent to
the paper's elimination of convective acceleration via total pressure; it's what makes the
high-Mach pressure ratio correct — see `docs/theory.md`.)

## Implementing the solver core

`PCIMSolver.steady_state` and `PCIMSolver.step` are both implemented (segregated SIMPLE;
steady is the dt→∞ limit of Section 4, transient adds the storage/inertia terms and the
`alpha` time-weighting), and the **energy equation** (eqs. 29-34) is coupled when
`SolverConfig.solve_energy` is set. All three are validated (`tests/`). The remaining work:

- **Higher-Mach / harder networks**: the segregated pressure↔temperature coupling is robust
  below ~Mach 0.4 for pressure-driven flow; above that the convective feedback strains the
  simple under-relaxation. A more robust coupling (e.g. coupled solve, line search, or the
  paper's exact total-pressure linearization) would extend the range.
- **More non-pipe components** (`Valve`/`Pipe` exist; pump/compressor/turbine/orifice/
  heat-exchanger to come), each via `resistance`/`convective_dp`/`inertance`/`flow_area`.
- Wall heat transfer (temperature-dependent `q̇`); currently `Node.heat_source` is constant.

Implementation notes worth preserving:
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

Transcribe equations directly from `docs/papers/` (do **not** rely on memory or OCR), and
**validate against the paper's benchmarks before trusting results**. Add each validated
benchmark as a test (see `tests/test_steady.py`, `test_transient.py`, `test_energy.py`).

## Conventions

- SI units throughout (Pa, kg/m³, K, J/kg, m, m³/s, s). State this when adding quantities.
- Work in **total** pressure `p₀` and **total** enthalpy `h₀` (how the paper eliminates the
  convective-acceleration term) — not static `p`/`h`.
- Keep a pure-Python reference even after native kernels land; select backend at runtime so
  tests can exercise both and cross-check them.
