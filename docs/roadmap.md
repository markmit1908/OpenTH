# OpenTH Roadmap & Product Strategy

This document turns the [vision statement](vision-statement.md) into a sequenced development
plan: where OpenTH is today, the phased path forward, and the validation and governance that
make it credible. For the current implementation detail see
[project-status.md](project-status.md); for the method see [theory.md](theory.md).

*Status legend:* ✅ done · ◑ partial · ☐ planned. Timing is by horizon — **Now** (current
focus), **Next** (following milestones), **Later** (post-MVP / long-range) — not hard dates.

---

## 1. Strategy in one page

**What OpenTH is.** An open-source, Python-first, reactor-grade 1D/0D thermal-hydraulic
system code — *"the OpenMC of reactor thermal hydraulics"*: transparent, scriptable,
testable, extensible.

**The gap we target.** Between Flownex-style system modelling and RELAP-class reactor
transient analysis, **advanced reactors first** (gas-cooled, molten-salt), where today's
options are a fragmented mix of custom scripts, Modelica, controlled/licensed codes, and
spreadsheets.

**The wedge.**
- *vs Flownex* — open source, scriptable, version-controlled models, CI-tested benchmarks,
  Python-ecosystem integration (not a GUI-first commercial environment).
- *vs RELAP/TRACE* — open development and modular, extensible physics instead of controlled
  distribution and legacy input decks; advanced-reactor-first instead of LWR heritage.

**Operating principles** (the moat is discipline, not any single feature):
1. **Clean two-layer architecture** — a Python modelling API over a solver core; the API is
   the product.
2. **Validation-first** — a public, CI-run benchmark suite from day one; nuclear-grade
   credibility is earned through verification over time.
3. **Physics depth in the right order** — a clean single-phase kernel → reactor features →
   (much later) two-phase. *"Do not start by building open RELAP."*
4. **Permissive, open governance** — MIT/BSD, public repo, strong docs, citable.

---

## 2. Where we are today (mid-2026)

The repository is a **validated single-phase flow-network kernel** — effectively the
**v0.1–v0.2** foundation complete, with parts of v0.3.

- ✅ **Solver core** — Greyvenstein implicit Pressure-Correction (PCIM), segregated
  SIMPLE-style: **steady**, **transient** (θ-weighted, water-hammer-capable), and the
  **energy equation** (non-isothermal). Compressible *and* incompressible.
- ✅ **Robustness** to ~Mach 0.74 (isothermal choking limit) via the compressible
  pressure-correction term.
- ✅ **Components** — pipe (friction + inertia), valve, pump/compressor (head curve + shaft
  work), pressure / mass-flow boundaries, constant nodal heat source.
- ✅ **Two modelling APIs** — name-based `FlowModel` and the port/connection `Circuit`
  (`import openth as th`), plus a declarative model dict (partial).
- ✅ **Validation** — the paper's four Section 5 benchmarks (steady pipeline vs. closed-form,
  water-hammer, branching network, blow-down), the `openth benchmark` CLI, unit tests, ruff +
  mypy clean.
- ✅ **Docs & onboarding** — theory notes, a user guide, a runnable `quickstart.ipynb`,
  VS Code config, MIT license, public GitHub repo.

---

## 3. Phased physics roadmap

The vision's version stack, with status and the gate that "closes" each phase:

| Version | Theme | Status | Closing validation gate |
|---------|-------|--------|--------------------------|
| v0.1 | Incompressible single-phase | ✅ | analytic friction (series/parallel) |
| v0.2 | Compressible ideal gas (steady + transient + energy) | ✅ | isothermal pipe law; water-hammer; blow-down |
| v0.3 | Real-fluid properties + buoyancy/gravity | ◑ Now | natural-circulation loop |
| v0.4 | Heat structures & conjugate heat transfer | ☐ Next | heat-exchanger transient |
| v0.5 | Controls, trips & rotating-machine dynamics | ☐ Next | pump coastdown |
| v0.6 | Reactor point kinetics | ☐ Next | reactivity-insertion transient |
| **MVP** | **Gas-cooled reactor loop demo** | ☐ Next | helium loop (see §7) |
| v0.7 | Homogeneous-equilibrium two-phase | ☐ Later | Edwards pipe blow-down |
| v1.0 | Drift-flux two-phase | ☐ Later | boiling-channel / loop benchmarks |
| — | Two-fluid model | ☐ Later (research) | — |

**v0.3 — real fluids + buoyancy** *(Now).* Add an `equation-of-state` layer beyond ideal
gas / constant density: helium real properties, then sodium / lead / FLiBe / water-steam
property tables behind the existing `FluidModel` ABC. Activate the **gravity/buoyancy term**
in the momentum balance (the `Pipe.angle`/`Node.elevation` fields exist but are not yet used)
— the prerequisite for natural circulation.

**v0.4 — heat structures.** A `HeatStructure` (wall thermal mass + 1D conduction) and a
`HeatExchanger` (primary/secondary coupling), plus wall heat-transfer correlations so
`q̇` becomes temperature-dependent (today `Node.heat_source` is a fixed power).

**v0.5 — controls & machines.** A `ControlBlock` layer (PID, trips, logic, time-/event-driven
actions) and rotating-machine dynamics (pump speed, inertia, **coastdown**); generalise the
current static valve schedule into the same event framework.

**v0.6 — point kinetics.** A `PointKinetics` power source (reactivity feedback hooks) so a
reactor transient can drive the thermal-hydraulics.

**v0.7+ — two-phase.** Deliberately last: homogeneous-equilibrium first, then drift-flux,
then (research) two-fluid. *"That is where codes go to die"* — enter only on a proven kernel.

---

## 4. Component roadmap

The vision's minimum object set, mapped to the code:

| Component | Purpose | Status |
|-----------|---------|--------|
| Volume / junction | control-volume state | ✅ `Node` |
| Pipe | friction, inertia, (wall) heat transfer | ◑ friction + inertia; wall HT ☐ |
| Pump / blower / compressor | head–flow relation, shaft work | ✅ `Pump` (speed dynamics ☐) |
| Valve / orifice | loss coefficient; choking later | ◑ `Valve`; orifice via pure-resistance pipe; choked-flow ☐ |
| Turbine | negative-head/work machine | ☐ (a `Pump` with negative head/work) |
| Heat structure | wall thermal mass, conduction | ☐ |
| Heat exchanger | primary/secondary coupling | ☐ |
| Boundary condition | pressure, mass flow, temperature | ✅ |
| Point-kinetics source | reactor power transient | ☐ |
| Control block | PID, trips, logic | ☐ |

New components plug in through the existing per-element closure hooks
(`resistance` / `convective_dp` / `inertance` / `flow_area` / `head` / `work_per_mass`) and
node hooks (`heat_source`), so most additions don't touch the solver.

---

## 5. Validation & benchmark roadmap

A public, CI-run benchmark suite is the credibility engine. Every PR runs regression tests;
every release publishes benchmark reports.

| # | Benchmark | Exercises | Status |
|---|-----------|-----------|--------|
| 1 | Analytic pipe friction | steady momentum/continuity | ✅ |
| 2 | Pressure-vessel blow-down | compressible transient | ✅ (single-phase) |
| 3 | Natural-circulation loop | buoyancy + energy | ☐ (needs v0.3 gravity) |
| 4 | Heat-exchanger transient | conjugate heat transfer | ☐ (v0.4) |
| 5 | Pump coastdown | inertia + machine dynamics | ☐ (v0.5) |
| 6 | Simple PWR loop | integrated single-phase loop | ☐ |
| 7 | Gas-cooled reactor loop | the **MVP** demo | ☐ |
| 8 | Molten-salt drain / freeze | property + phase handling | ☐ (Later) |
| 9 | Point-kinetics reactivity insertion | neutronics coupling | ☐ (v0.6) |
| 10 | Coupled neutronics/T-H (OpenMC) | external coupling | ☐ (Later) |
| — | Edwards pipe blow-down | two-phase critical flow | ☐ (v0.7) |

Already in-repo today: the Greyvenstein Section 5 set (`openth benchmark`) covers #1 and the
single-phase part of #2.

---

## 6. Cross-cutting workstreams

- **Numerics / robustness** — transonic & choked-flow boundary treatment (today solid to
  ~Mach 0.74, mesh-sensitive near choking); stronger coupled/Newton option for stiff cases.
- **Performance (native core)** — reimplement the hot kernels (linear solves, coefficient
  assembly) in C/C++/Rust behind the same Python interfaces (`native/`), keeping the
  pure-Python reference for cross-checking. Sparse/Jacobian, MPI/OpenMP later.
- **APIs & UX** — `FlowModel` + `Circuit` (done); finish declarative JSON/YAML
  (de)serialization (`io/`), results/post-processing helpers, an example gallery.
- **External coupling** — Python-native coupling to OpenMC (neutronics), and to
  MOOSE / Cardinal / NekRS for higher-fidelity sub-models.
- **LLM interface** — the optional two-way `openth.llm` layer (build/query models in natural
  language) behind the `[llm]` extra; never imported by the core.
- **Docs, CI & governance** — set up CI (tests + benchmark regression on every PR);
  contributor guide + code of conduct; a citation/methods paper; an eventual user workshop.

---

## 7. MVP milestone — gas-cooled reactor loop

The first credible demonstrator (vision §"MVP target"):

> A helium-cooled reactor loop with **core heat input, blower, recuperator/heat exchanger,
> heat structures, controls, and a transient pump coastdown.**

Advanced-reactor-relevant, avoids water two-phase, and shows value fast. Gap to close:

| Needed | Status |
|--------|--------|
| Helium loop topology (pipes/junctions) | ✅ |
| Blower (pump) | ✅ |
| Core heat input | ◑ (constant `heat_source`; wall/structure ☐) |
| Recuperator / heat exchanger | ☐ (v0.4) |
| Heat structures | ☐ (v0.4) |
| Controls / trips | ☐ (v0.5) |
| Pump coastdown (machine dynamics) | ☐ (v0.5) |

So the MVP is gated on **v0.4 + v0.5**; it then becomes benchmark #7 and the flagship example.

---

## 8. Near-term next steps (Now → Next)

1. **Activate gravity/buoyancy** in the momentum balance (use `angle`/`elevation`) and add a
   **natural-circulation** benchmark (#3). *(smallest step with the biggest reactor payoff.)*
2. **Real-fluid property layer** behind `FluidModel` (helium real props first), with a clean
   table/correlation interface.
3. **Heat structures + heat exchanger** components and wall heat-transfer correlations (v0.4).
4. **Pump coastdown / machine dynamics** and a first **control block** (v0.5).
5. **CI** running the test suite + benchmark regressions on every PR.
6. Opportunistic: turbine/orifice as named components; finish `io` JSON round-trip;
   near-choking robustness.

---

## 9. Release & governance

- **License** MIT (permissive — industrial adoption; avoid GPL). Public GitHub repo. ✅
- **Versioning** — semantic-ish, tracking the phase milestones above (currently 0.0.x).
- **Quality gates** — ruff + mypy + pytest clean on every change; benchmark regression in CI;
  validated cases become tests.
- **Community** (as adoption grows) — contributor guide, code of conduct, example gallery,
  discussion forum, a citation/methods paper, and eventually a user meeting.

The order matters: **clean architecture, testable physics, excellent Python UX, and public
benchmarks first** — the nuclear-grade credibility compounds from there.
