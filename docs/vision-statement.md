---
jdate: 2026-05-31
---

# OpenTH — Vision Statement

> The founding concept for OpenTH. The repository today implements the single-phase PCIM
> flow-network kernel (the v0.1–v0.2 foundation described below); this document is the north
> star for where it is headed. For the sequenced plan see [`roadmap.md`](roadmap.md); for the
> implementation snapshot see [`project-status.md`](project-status.md).

Concept: build OpenTH — an open-source, Python-first, reactor-grade 1D thermal-hydraulic system code.

OpenMC’s winning pattern is: open core engine, Python model-building API, modern data formats, strong docs, verification culture, and community trust. OpenMC is community-developed, open source, supports CSG/CAD geometry, HDF5 nuclear data, MPI/OpenMP, and Python workflows.  

Your analog would target the gap between Flownex-style system modeling and RELAP-style reactor transient analysis. Flownex is positioned around system-level flow/heat-transfer simulation and design optimization, while RELAP5-3D is used for reactor safety analysis, design, training, and advanced reactor transient/accident analysis but remains controlled/licensed through INL.  

Core thesis

Build an open, extensible thermal-hydraulics platform for advanced reactors:

A modern, Python-native, open-source 1D/0D thermal-fluid network solver for steady-state and transient reactor systems, with optional two-phase, compressible, heat-structure, controls, and reactor-coupling capability.

Not “replace RELAP on day one.” Instead:

Be the OpenMC of reactor thermal hydraulics: transparent, scriptable, testable, extensible, and trusted by advanced reactor teams.

Product positioning

The code should sit here:

|   |   |   |
|---|---|---|
|Need|Existing tool|New code opportunity|
|Fast plant/system design|Flownex|Open, scriptable, reactor-focused|
|Licensing-grade LWR safety|RELAP/TRACE|Not first target|
|Advanced reactor R&D|Mix of custom scripts, RELAP, Modelica, Flownex|Strong target|
|Coupled neutronics/T-H|Bespoke workflows|Python-native coupling|
|Open validation benchmark|Fragmented|Major opportunity|

Initial scope

Start with a single-phase compressible/incompressible network solver.

Minimum physical objects:

|   |   |
|---|---|
|Component|Purpose|
|Volume / junction|Control-volume state|
|Pipe|Friction, inertia, heat transfer|
|Pump / blower|Head-flow relation|
|Valve / orifice|Loss coefficient, choking later|
|Heat structure|Wall thermal mass, conduction|
|Heat exchanger|Primary/secondary coupling|
|Boundary condition|Pressure, mass flow, temperature|
|Point kinetics source|Optional reactor power transient|
|Control block|PID, trips, logic|

Then expand to:

1. gas-cooled reactors
2. molten-salt loops
3. water/steam two-phase
4. sodium / lead / FLiBe / helium properties
5. reactor transient benchmarks
6. uncertainty/sensitivity workflows

Architecture

Use a two-layer design:

1. High-performance solver core

Written in C++ or Rust, exposed to Python.

Responsibilities:

- nonlinear residual assembly
- sparse matrix solve
- time integration
- component equations
- fluid property calls
- event handling
- Jacobians, ideally automatic or semi-automatic

2. Python modeling layer

This is where you follow OpenMC.

Example style:

import openth as th

  

helium = th.Fluid("helium")

  

loop = th.Model()

core = loop.add(th.HeatExchanger("core", power=5e6))

pipe1 = loop.add(th.Pipe(length=4.0, diameter=0.08, fluid=helium))

blower = loop.add(th.Pump(curve="blower_curve.csv"))

  

loop.connect(core.outlet, pipe1.inlet)

loop.connect(pipe1.outlet, blower.inlet)

loop.connect(blower.outlet, core.inlet)

  

loop.solve_steady_state()

loop.transient(t_end=500)

The Python API is the moat.

Numerical approach

Use a modern implicit formulation:

- finite-volume 1D network formulation
- conservation of mass, momentum, energy
- pressure-based or fully coupled Newton solve
- implicit time stepping for stiffness
- sparse linear algebra
- event handling for trips, valve actions, control logic
- steady-state solver as “transient to equilibrium” or direct nonlinear solve

For the first version, avoid full RELAP-like two-fluid two-phase complexity. That is where codes go to die.

A good phased model stack:

|   |   |
|---|---|
|Phase|Physics|
|v0.1|incompressible single-phase|
|v0.2|compressible ideal gas|
|v0.3|real-fluid properties|
|v0.4|heat structures and conjugate heat transfer|
|v0.5|controls and trips|
|v0.6|reactor point kinetics|
|v0.7|homogeneous equilibrium two-phase|
|v1.0|drift-flux two-phase|
|later|two-fluid model|

Differentiator vs Flownex

Do not compete on GUI first.

Compete on:

- open source
- scriptability
- version-controlled models
- reproducible simulations
- CI-tested benchmarks
- Python ecosystem integration
- transparent equations
- reactor-specific examples
- coupling to OpenMC, MOOSE, Cardinal, NekRS, etc.

Flownex has a polished commercial system-simulation environment. Your wedge is developer-grade reactor simulation infrastructure.

Differentiator vs RELAP

RELAP’s strength is validation history and safety pedigree. Your advantage is openness and extensibility.

Make the pitch:

|   |   |
|---|---|
|RELAP|New code|
|controlled distribution|open development|
|legacy input decks|Python models|
|hard to modify|modular components|
|licensing pedigree|transparency and CI validation|
|LWR heritage|advanced-reactor-first|

INL says RELAP5-3D is used for reactor safety analysis, design, training, and advanced reactor analysis, but it is also a controlled 10 CFR 810 code with license/export review requirements. That friction is exactly the opening.  

Validation strategy

This is the most important part.

Create a public benchmark suite from day one:

1. analytic pipe friction tests
2. Edwards pipe blowdown, later
3. natural circulation loop
4. heat exchanger transient
5. pump coastdown
6. simple PWR loop
7. gas-cooled reactor loop
8. molten salt drain/freeze transient
9. point-kinetics reactivity insertion
10. coupled neutronics/T-H demonstration with OpenMC

Every PR should run regression tests. Every release should publish benchmark reports.

Governance model

Use the OpenMC pattern:

- permissive license: BSD-3 or MIT
- public GitHub repo
- strong documentation
- example gallery
- discussion forum
- contributor guide
- code of conduct
- formal verification test suite
- citation paper
- annual workshop/user meeting eventually

Avoid GPL if you want industrial adoption.

Suggested name

Options:

- OpenTH
- OpenFlowSys
- ThermoNet
- ReTH
- OpenREL
- HeliosTH
- AstraTH
- NexusTH

My favorite: OpenTH. Clear, boring, credible.

MVP target

The first credible demo should be:

A helium-cooled reactor loop with core heat input, blower, recuperator/heat exchanger, heat structures, controls, and transient pump coastdown.

That is advanced-reactor-relevant, avoids water two-phase initially, and shows value fast.

The strategic insight

Do not start by building “open RELAP.” Start by building:

A modern open thermal-fluid simulation kernel that can grow into RELAP-class capability.

The order matters. Start with clean architecture, testable physics, excellent Python UX, and public benchmarks. The nuclear-grade credibility comes from disciplined verification and validation over time.