# FlowCalc

A fluid **flow network calculator** built on Greyvenstein's implicit **Pressure Correction (PC / PCIM)**
method for transient flow in pipe networks.

The solver is segregated (SIMPLE-style): it discretizes the continuity, momentum and energy equations
with a finite-volume scheme and solves them iteratively via a **pressure-correction equation** at each
time step. The same machinery handles **liquid and gas**, **isothermal and non-isothermal**, and **fast
and slow transients**, and extends naturally to **non-pipe components** (valves, pumps, compressors,
turbines, orifices, heat exchangers).

> Reference: G. P. Greyvenstein, *"An implicit method for the analysis of transient flows in pipe
> networks"*, Int. J. Numer. Meth. Engng 2002; 53:1127–1143.
> See [`docs/papers/`](docs/papers/) and [`docs/theory.md`](docs/theory.md).

## Status

Early prototype. The **steady-state and transient solvers work and are validated**, including
**non-isothermal flow**: steady reproduces the closed-form isothermal compressible pipe-flow relation
to ~0% error (`examples/pipeline_steady.py`); the transient marches to that steady state, conserves
mass, and captures water-hammer and the blow-down decay (`examples/blowdown_transient.py`); the energy
equation conserves total enthalpy in adiabatic flow and matches `Q/ṁ` for heat addition
(`examples/heated_pipe.py`, `tests/`). A compressible pressure-correction term keeps the solve robust
up to ~Mach 0.74 (the isothermal choking limit); more component models are next.
Python-first for clarity and rapid iteration; performance-critical kernels will later be reimplemented
in C/C++ behind the same interfaces — see [`native/`](native/). A two-way LLM interface (build/query
networks in natural language) is planned in [`src/flowcalc/llm/`](src/flowcalc/llm/).

## Quick start

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e ".[dev,llm]"

pytest                       # run tests
python examples/pipeline_steady.py
```

> Requires Python ≥ 3.10. (The system Python on this machine is 3.9.6 — create a venv with a newer
> interpreter, e.g. via `pyenv` or `python3.12 -m venv .venv`.)

## Layout

| Path | Purpose |
|------|---------|
| `src/flowcalc/network/` | Topology: `Node` (control volume / cell centre) + `Element` (face / branch) + `Network` graph |
| `src/flowcalc/components/` | Concrete elements: `Pipe`, boundaries, `Valve`; non-pipe components added here |
| `src/flowcalc/fluids/` | Equations of state / fluid property models (`IdealGas`, incompressible) |
| `src/flowcalc/solver/` | The PCIM solver, time integration, linear solves |
| `src/flowcalc/io/` | Network (de)serialization |
| `src/flowcalc/llm/` | Two-way LLM interface (optional) |
| `examples/` | Runnable cases, including paper benchmarks |
| `docs/` | Theory notes and the source papers |
| `native/` | C/C++ acceleration (stub for now) |

## License

MIT — see [LICENSE](LICENSE).
