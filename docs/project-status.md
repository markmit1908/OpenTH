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

- `python -m pytest` → **18 passed** (fluid EOS, topology, Thomas solver, steady-state
  validation, and transient validation: march-to-steady, mass conservation, water-hammer)
- `ruff check` → clean · `mypy` → clean (`py.typed` marker present)
- `examples/pipeline_steady.py` reproduces the analytical isothermal pressure ratio to
  **0.00% error up to Mach 0.5**; `examples/blowdown_transient.py` tracks quasi-steady to <1%
- `flowcalc --version` works

## Steady + transient core — done ✅ (isothermal)

Both `PCIMSolver.steady_state` (dt→∞ limit) and `PCIMSolver.step` (transient) are
segregated SIMPLE iterations with Picard density updates and a fixed temperature field.
Per-face momentum is `Δp = K·mdot·|mdot| + C` (+ inertia `(Δx/A) dṁ/dt` in the transient);
continuity assembles the pressure-correction system (with a `V·∂ρ/∂p/Δt` storage diagonal
in the transient), solved sparsely.

- Steady validated vs. `p1²−p2² = G²RT(fL/D + 2ln(p1/p2))`.
- Transient validated: marches to the steady fixed point (~1e-9), conserves mass in the
  θ-weighted sense (~1e-9), produces water-hammer overpressure, and reproduces the
  blow-down mass-flow decay.
- Robustness: flow under-relaxation (`SolverConfig.relaxation`) tames the high-Mach
  convective feedback; `_converge_flows_only` handles networks with no interior pressure
  unknowns (single element between two pressure boundaries).

## Still the real work

- **Energy equation / non-isothermal flow** (eqs. 29–34) — solve for total enthalpy `h₀`
  alternately with the pressure loop and update temperatures. Unlocks Fig. 3, 6, 7.
- Non-pipe component closures beyond `Valve` (pump/compressor/turbine/orifice/heat-exchanger).

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

Couple the **energy equation** (non-isothermal flow): solve total enthalpy `h₀` alternately
with the pressure-correction loop and update temperatures, enabling the adiabatic and
heat-transfer cases (Figs. 3, 6, 7).
