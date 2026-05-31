# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

FlowCalc is a **fluid flow-network solver** implementing Greyvenstein's implicit
**Pressure-Correction (PC / PCIM)** method (Int. J. Numer. Meth. Engng 2002; 53:1127–1143;
PDF in `docs/papers/`, distilled in `docs/theory.md`). Python-first prototype; performance
kernels are earmarked for later C/C++ reimplementation (`native/`), and a two-way LLM
interface is planned (`src/flowcalc/llm/`).

**Project status: scaffold.** The control flow, interfaces, fluid models, topology and
linear-solve back ends exist and are tested. The *numerical core* of the solver is
deliberately stubbed with `NotImplementedError` — see "Implementing the solver core" below.

## Commands

```bash
# Environment (needs Python >= 3.10; system python here is 3.9.6 — use pyenv/python3.12)
python3.12 -m venv .venv && source .venv/bin/activate
pip install -e ".[dev,llm]"

pytest                                  # all tests
pytest tests/test_fluids.py             # one file
pytest tests/test_network.py::test_pipe_area   # one test
pytest --cov=flowcalc                   # with coverage

ruff check . && ruff format .           # lint + format
mypy                                    # type-check (config in pyproject.toml)

python examples/pipeline_steady.py      # paper §5.1 benchmark setup (solve is gated)
flowcalc --version                      # console entry point (flowcalc.cli)
```

## Architecture

The code mirrors the method's finite-volume structure (`docs/theory.md` maps every module
to the paper's equations — read it before touching the numerics):

- **`network/`** — topology. `Node` is a finite-volume **control volume / cell centre**
  carrying the scalar state (`p₀, ρ, T, h₀`); `Element` is a **face / branch** carrying
  volumetric flow `Q`. `Network` owns the graph and assigns linear-system row indices.
  The cell-centre-vs-face split is the paper's discretization (Fig. 1) and is load-bearing.
- **`components/`** — concrete `Element`s. `Pipe` is canonical; non-pipe components
  (`Valve`, and future pump/compressor/turbine/orifice/heat-exchanger) plug in by
  overriding `momentum_coeffs`. Boundaries (`PressureBoundary`, `MassFlowBoundary`) mark a
  node `is_boundary`, removing it from the unknowns.
- **`fluids/`** — equations of state behind the `FluidModel` ABC. `IdealGas` implements
  `p = sρRT` (the paper's gas variant; `helium()` is the benchmark fluid); `Incompressible`
  the liquid variant.
- **`solver/`** — `PCIMSolver` (segregated SIMPLE-style loop), `SolverConfig` (notably
  `alpha ∈ [0.5,1]`, default 0.6), and `linear.py` (Thomas for series pipelines, scipy
  sparse for networks).
- **`io/`** — declarative dict/JSON (de)serialization; also the LLM exchange payload.
- **`llm/`** — optional two-way LLM interface. **The core must never import this**; the
  `anthropic` SDK lives behind the `[llm]` extra.

### The central abstraction

Each element linearizes its corrected flow as a function of its two end-node pressure
corrections — `Q'ᵢ = a⁻·p'ᵢ − a⁺·p'ᵢ₊₁` (paper eq. 20), returned by
`Element.momentum_coeffs() -> MomentumCoeffs`. This single interface is what lets the same
pressure-correction machinery handle pipes and arbitrary non-pipe components uniformly.
Preserve it when adding components.

## Implementing the solver core

The stubbed numerics (search for `NotImplementedError` / `TODO(core)`) are the project's
real work, concentrated in:

- `components/pipe.py::Pipe.momentum_coeffs` — eqs. (21)-(22), incl. effective friction `f̃`.
- `solver/pcim.py::PCIMSolver.steady_state` / `.step` — the segregated loop, eqs. (12)-(34).

Transcribe equations directly from `docs/papers/` (do **not** rely on memory or OCR), and
**validate against the paper's benchmarks before trusting results**: steady pipeline vs.
the ODE solution of eqs. (35)-(37) (`examples/pipeline_steady.py`, Fig. 2/3), then the
transient cases (Fig. 4-10). Add each validated benchmark as a test.

## Conventions

- SI units throughout (Pa, kg/m³, K, J/kg, m, m³/s, s). State this when adding quantities.
- Work in **total** pressure `p₀` and **total** enthalpy `h₀` (how the paper eliminates the
  convective-acceleration term) — not static `p`/`h`.
- Keep a pure-Python reference even after native kernels land; select backend at runtime so
  tests can exercise both and cross-check them.
