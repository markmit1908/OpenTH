# Project status & scaffold summary

Snapshot of the initial FlowCalc scaffold (2026-05-30). Living document — update as the
solver core lands.

## What was built

A Python-first scaffold for **FlowCalc**, structured directly around Greyvenstein's PC/PCIM
method (the segregated, SIMPLE-style implicit pressure-correction scheme described in
[`papers/`](papers/Greyvenstein-2001-implicit-transient-pipe-networks.pdf), distilled in
[`theory.md`](theory.md)).

```
FlowCalc/
├── CLAUDE.md              ← guidance for future Claude Code sessions
├── README.md  LICENSE(MIT)  pyproject.toml  .gitignore
├── src/flowcalc/
│   ├── network/      Node (control volume) · Element (face) · Network graph
│   ├── components/   Pipe · Valve · pressure/mass-flow boundaries
│   ├── fluids/       FluidModel ABC · IdealGas (p=sρRT, helium) · Incompressible
│   ├── solver/       PCIMSolver · SolverConfig (α∈[0.5,1]) · Thomas + sparse
│   ├── io/           declarative dict/JSON (de)serialization
│   ├── llm/          two-way LLM interface (optional [llm] extra; core never imports it)
│   └── cli.py
├── tests/            fluids · network · linear  (11 passing)
├── examples/         pipeline_steady.py — the paper's §5.1 helium benchmark
├── docs/
│   ├── papers/       Greyvenstein-2001-...pdf
│   ├── theory.md     maps every module to the paper's equations
│   └── project-status.md   ← this file
└── native/           stub for later C/C++ kernels
```

## Key design decision

The architecture mirrors the paper's discretization: **scalars (p₀, ρ, T, h₀) at cell
centres (`Node`); flow `Q` at faces (`Element`)**. The load-bearing abstraction is
`Element.momentum_coeffs() → MomentumCoeffs(a_plus, a_minus)` — the linearized eq. (20).
That single interface is what lets the same pressure-correction machinery handle pipes
*and* future non-pipe components (pumps, compressors, turbines, valves) uniformly, exactly
as the paper describes.

## Verified

- `pytest` → **11 passed** (fluid EOS, topology, Thomas solver)
- `ruff check` → clean · `mypy` → clean (`py.typed` marker present)
- `examples/pipeline_steady.py` builds the 21-node helium pipeline and cleanly reports the
  solver core as pending
- `flowcalc --version` works

## Deliberately left as the real work

The numerical core is stubbed with `NotImplementedError` + `TODO(core)` rather than guessed
at — physics transcribed from OCR should not be shipped as if correct:

- `Pipe.momentum_coeffs` — eqs. (21)–(22)
- `PCIMSolver.steady_state` / `.step` — the segregated loop, eqs. (12)–(34)

Both [`theory.md`](theory.md) and [`../CLAUDE.md`](../CLAUDE.md) point to these and stress
validating against the paper's benchmarks (Figs. 2–10) before trusting output.

## Environment note

The system Python is 3.9.6, but the project targets **≥ 3.10**. The `.venv` was created with
Homebrew **Python 3.12** (`/opt/homebrew/bin/python3.12`):

```bash
/opt/homebrew/bin/python3.12 -m venv .venv && source .venv/bin/activate
pip install -e ".[dev,llm]"
```

## Suggested next step

Implement the steady-state core, starting with `Pipe.momentum_coeffs`, and validate against
the eq. (35)–(37) ODE benchmark in `examples/pipeline_steady.py`.
