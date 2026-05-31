# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

FlowCalc is a **fluid flow-network solver** implementing Greyvenstein's implicit
**Pressure-Correction (PC / PCIM)** method (Int. J. Numer. Meth. Engng 2002; 53:1127–1143;
PDF in `docs/papers/`, distilled in `docs/theory.md`). Python-first prototype; performance
kernels are earmarked for later C/C++ reimplementation (`native/`), and a two-way LLM
interface is planned (`src/flowcalc/llm/`).

**Project status.** The **steady-state and transient solvers are implemented and
validated** (isothermal). Steady matches the closed-form isothermal compressible
pipe-flow law to ~0% error; the transient marches to the steady solution as its fixed
point, conserves mass exactly in the storage term, and reproduces water-hammer
overpressure and the blow-down decay. See `examples/` and `tests/`. **Non-isothermal flow
(energy-equation coupling) is the main remaining gap** — see "Implementing the solver
core" below.

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
  overriding `resistance` / `convective_dp`. `PressureBoundary` pins a node's pressure and
  marks it `is_boundary` (removed from unknowns); `MassFlowBoundary` instead sets the
  node's `mass_source` and leaves it an unknown.
- **`fluids/`** — equations of state behind the `FluidModel` ABC. `IdealGas` implements
  `p = sρRT` (the paper's gas variant; `helium()` is the benchmark fluid); `Incompressible`
  the liquid variant.
- **`solver/`** — `PCIMSolver` (segregated SIMPLE-style loop), `SolverConfig` (notably
  `alpha ∈ [0.5,1]` for the transient, and `relaxation ∈ (0,1]` for the steady solve), and
  `linear.py` (Thomas for series pipelines, scipy sparse for networks — the solver
  currently uses sparse).
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
`alpha` time-weighting) and validated — **but only for isothermal flow**: temperature is
held fixed and the energy equation is not yet solved. The remaining work:

- **Energy equation / non-isothermal flow** (eqs. 29-34): solve for total enthalpy `h₀`
  (upwinded, tridiagonal/sparse) alternately with the pressure-correction loop, and update
  temperatures. This unlocks the paper's adiabatic and heat-transfer cases (Fig. 3, 6, 7).
- Non-pipe component closures (`Valve.resistance`/`Pipe` exist; pump/compressor/turbine/
  orifice/heat-exchanger to come), each via `resistance`/`convective_dp`/`inertance`.

Implementation notes worth preserving:
- The flow update is **under-relaxed** (`SolverConfig.relaxation`): the convective term
  couples flow to itself and diverges at high Mach in pressure-driven steady problems.
- A network with **no interior pressure unknowns** (e.g. one element between two pressure
  boundaries — the blow-down case) is handled by `_converge_flows_only`: momentum alone
  fixes the flows.
- The transient is **mass-conservative by construction** (the storage term equals the
  θ-weighted flux integral); `tests/test_transient.py` checks this to ~1e-9.

Transcribe equations directly from `docs/papers/` (do **not** rely on memory or OCR), and
**validate against the paper's benchmarks before trusting results**. Add each validated
benchmark as a test (see `tests/test_steady.py`, `tests/test_transient.py`).

## Conventions

- SI units throughout (Pa, kg/m³, K, J/kg, m, m³/s, s). State this when adding quantities.
- Work in **total** pressure `p₀` and **total** enthalpy `h₀` (how the paper eliminates the
  convective-acceleration term) — not static `p`/`h`.
- Keep a pure-Python reference even after native kernels land; select backend at runtime so
  tests can exercise both and cross-check them.
