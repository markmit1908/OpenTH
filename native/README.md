# native/ — C/C++ acceleration (stub)

This directory is reserved for performance-critical kernels reimplemented in C/C++ and
exposed to Python. Nothing is built here yet; the pure-Python implementation in
`src/openth/` is the reference.

## What will live here

The hot paths in the PCIM solver, behind the *same* Python interfaces so the switch is
transparent:

- **Linear solves** — the tridiagonal (Thomas) and sparse pressure-correction / energy
  systems (`openth.solver.linear`).
- **Coefficient assembly** — per-element momentum link coefficients (eqs. 21-22) and the
  continuity/energy coefficients (eqs. 24-34), looped over every face/cell each iteration.

## Intended approach (decide when we get here)

- Bindings via **pybind11** or **nanobind**, built with **scikit-build-core** / CMake.
- Keep a pure-Python fallback; select the backend at runtime so tests can run either.
- Validate the native path against the Python reference on the paper benchmarks before
  making it the default.
