# Project status & scaffold summary

FlowCalc status snapshot (started 2026-05-30). Living document — update as the solver
evolves. The steady, transient, and energy (non-isothermal) solves are now implemented and
validated; the layout/decisions below trace how it got here.

## What was built

A Python-first implementation of **FlowCalc**, structured directly around Greyvenstein's
PC/PCIM method (the segregated, SIMPLE-style implicit pressure-correction scheme described
in [`papers/`](papers/Greyvenstein-2001-implicit-transient-pipe-networks.pdf), distilled in
[`theory.md`](theory.md)).

```
FlowCalc/
├── CLAUDE.md              ← guidance for future Claude Code sessions
├── README.md  LICENSE(MIT)  pyproject.toml  .gitignore
├── src/flowcalc/
│   ├── network/      Node (control volume) · Element (face) · Network graph
│   ├── components/   Pipe · Valve · Pump/compressor · pressure/mass-flow boundaries
│   ├── fluids/       FluidModel ABC · IdealGas (p=sρRT, helium) · Incompressible
│   ├── solver/       PCIMSolver · SolverConfig (α∈[0.5,1]) · Thomas + sparse
│   ├── model.py      FlowModel — high-level builder/facade
│   ├── benchmarks.py paper Section 5 test cases (build_* / run_*, `flowcalc benchmark`)
│   ├── io/           declarative dict/JSON (de)serialization
│   ├── llm/          two-way LLM interface (optional [llm] extra; core never imports it)
│   └── cli.py
├── tests/            fluids·network·linear·steady·transient·energy·pump·model·benchmarks (35 passing)
├── examples/         pipeline_steady · blowdown_transient · heated_pipe · pump_loop
├── docs/
│   ├── papers/       Greyvenstein-2001-...pdf
│   ├── theory.md     maps every module to the paper's equations
│   ├── user-guide.md end-user manual (FlowModel API + running the benchmarks)
│   └── project-status.md   ← this file
└── native/           stub for later C/C++ kernels
```

## Key design decision

The architecture mirrors the paper's discretization: **scalars (p₀, ρ, T, h₀) at cell
centres (`Node`); mass flow `ṁ` at faces (`Element`)**. The load-bearing abstraction is the
per-element momentum closure — `Element.resistance(ρ)` (friction `Δp = K·ṁ|ṁ|`) plus
optional `convective_dp` / `inertance` / `flow_area`. The solver linearizes it into the
pressure-correction sensitivity `sₑ` (eq. 20). That single interface is what lets the same
machinery handle pipes *and* non-pipe components (valves now; pumps, compressors, turbines
to come) uniformly, exactly as the paper describes.

## Verified

- `python -m pytest` → **35 passed** (fluid EOS, topology, Thomas solver, steady — incl.
  high-Mach pressure-driven — transient — march-to-steady, mass conservation, water-hammer
  — energy: adiabatic h₀ conservation, heat addition, transient↔steady consistency — the
  pump/compressor: uphill flow, operating point, temperature rise — and the `FlowModel`
  facade + the four Section 5 benchmarks)
- `ruff check` → clean · `mypy` → clean (`py.typed` marker present)
- `examples/pipeline_steady.py` reproduces the analytical pressure ratio to **0.00% up to
  Mach 0.5**; `examples/blowdown_transient.py` tracks quasi-steady to <1%;
  `examples/heated_pipe.py` shows expansion cooling and heat addition
- `flowcalc --version` works

## Steady + transient + energy core — done ✅

`PCIMSolver.steady_state` (dt→∞ limit), `PCIMSolver.step` (transient), and the **energy
equation** (`SolverConfig.solve_energy`) are segregated SIMPLE iterations with Picard
density updates. Per-face momentum is `Δp = K·mdot·|mdot| + C` (+ inertia `(Δx/A) dṁ/dt` in
the transient); continuity assembles the pressure-correction system (with a `V·∂ρ/∂p/Δt`
storage diagonal in the transient). The energy solve transports total enthalpy `h₀` by
upwind convection on the converged flow field, alternating with the pressure loop.

- Steady validated vs. `p1²−p2² = G²RT(fL/D + 2ln(p1/p2))`.
- Transient validated: marches to the steady fixed point (~1e-9), conserves mass (~1e-9),
  produces water-hammer, reproduces the blow-down decay.
- Energy validated: adiabatic h₀ conserved (~1e-9, with expansion cooling); heat addition
  matches `Q/ṁ`; transient energy matches the steady energy solve.
- Robustness: a **compressible upwind `p'`-convection term** in the continuity equation
  (the density correction `ρ'=∂ρ/∂p·p'` of the convected mass, eqs. 25-27) extends the
  pressure-driven range to ~Mach 0.74 (the choking limit); flow under-relaxation tames the
  convective feedback; `_converge_flows_only` handles no-interior-unknown networks; the
  energy `h₀→T` step is change-limited to keep the kinetic coupling stable.

## User model & benchmarks — done ✅

`flowcalc.model.FlowModel` is the high-level facade: name-based nodes (auto-created),
`add_pipe(..., n_cells=N)` that subdivides a pipe into N finite-volume cells and assigns
each node its control volume, `add_valve`/`add_pump`, callable (time-varying) pressure
boundaries, and `steady_state()` / `run(dt, duration, record=...)`. It reproduces the
hand-built pipeline result to ~0%.

`flowcalc.benchmarks` builds and runs the paper's four Section 5 cases on `FlowModel`
(steady pipeline, valve closure / water hammer, branching-network valve closures,
pressure-vessel blow-down), exposed via `flowcalc benchmark [name]`.

## Still the real work

- **Transonic / near-choking robustness**: solid to ~Mach 0.74 (isothermal choking limit),
  but near choking the solve is mesh-sensitive and can collapse to a false zero-flow state.
- More non-pipe components (`Pipe`/`Valve`/`Pump` exist; turbine = negative-head pump,
  orifice = pure resistance, heat exchanger to come) via the closure hooks.
- Wall heat transfer (temperature-dependent `q̇`); `Node.heat_source` is currently constant.

Both [`theory.md`](theory.md) and [`../CLAUDE.md`](../CLAUDE.md) stress transcribing
equations directly from the PDF and validating against the paper's benchmarks before
trusting output.

## Environment note

The system Python is 3.9.6, but the project targets **≥ 3.10**. The `.venv` was created with
Homebrew **Python 3.12** (`/opt/homebrew/bin/python3.12`):

```bash
/opt/homebrew/bin/python3.12 -m venv .venv && source .venv/bin/activate
pip install -e ".[dev,llm]"
```

## Suggested next step

Add a turbine (negative-head `Pump`) and an orifice (pure-resistance `Pipe`) via the same
closure hooks, and/or harden the near-choking regime (choked-flow boundary treatment) to
remove the mesh-sensitive zero-flow collapse.
