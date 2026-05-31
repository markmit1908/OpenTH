# OpenTH Backlog

A capture of requested capabilities not yet built, to be sequenced into the
[roadmap](roadmap.md). Grouped by theme. Status: ✅ done · ◑ partial · ☐ planned.

---

## 1. Component catalog

A target set of element/component types (codes from established system codes such as
Flownex). OpenTH adds these through the per-element closure hooks
(`resistance` / `convective_dp` / `inertance` / `flow_area` / `head` / `work_per_mass`) and
node hooks (`heat_source`), so most are data + a small class, not solver changes.

### Ducts / pipes / resistances

| Code | Component | Status / notes |
|------|-----------|----------------|
| DW | Darcy–Weisbach pipe | ✅ `Pipe` (constant `f`; add Colebrook/Moody `f(Re, ε/D)`) |
| CP | Compressible pipe | ✅ (compressible PCIM is the core) |
| DR | Resistance duct | ◑ `Pipe`-like pure resistance; add as a named component |
| HW | Hazen–Williams pipe | ☐ (liquid head-loss correlation) |
| DG | Duct with area change | ☐ (variable `flow_area`; diffuser/nozzle, convective term) |

### Orifices / restrictors

| Code | Component | Status / notes |
|------|-----------|----------------|
| RL | Restrictor with loss coefficient | ◑ `Valve` (loss-coefficient resistance) |
| RD | Restrictor with discharge coefficient | ☐ (Cd-based) |
| OR | British-standard orifice | ☐ (BS-1042 correlation) |
| SO | Special orifice (Lucas correlation) | ☐ |
| CO | Cooling orifice | ☐ |

### Rotating machines

| Code | Component | Status / notes |
|------|-----------|----------------|
| FA | Fan / pump | ✅ `Pump` (quadratic head curve) |
| CM | Compressor | ✅ `Pump` (head + shaft work raises h₀) |
| TU | Turbine | ☐ (`Pump` with negative head / negative work) |
| PDC | Positive-displacement compressor | ☐ |
| CA | Cold-air unit | ☐ |

### Heat exchangers

| Code | Component | Status / notes |
|------|-----------|----------------|
| HX | Heat exchanger | ✅ `HeatExchanger` (lumped UA·ΔT coupling) |
| RX | Recuperator | ◑ via `HeatExchanger` (same-fluid); add effectiveness/NTU + distributed |
| UX | Duct heat exchanger | ☐ |
| EV | Evaporator | ☐ (needs two-phase) |
| CX | Gas/liquid complex heat exchanger | ☐ (needs multi-fluid) |

### Controls & relief

| Code | Component | Status / notes |
|------|-----------|----------------|
| CNT | PID controller | ☐ (control-block layer, v0.5) |
| PR | Pressure-relief valve | ☐ (event/threshold valve) |
| PG | Pressure-regulating valve | ☐ (setpoint-controlled valve) |

### Reactor & generic

| Code | Component | Status / notes |
|------|-----------|----------------|
| PB / PBR | Pebble-bed reactor core | ☐ (porous-media friction + heat + point kinetics) |
| PMR | Prismatic-fuel reactor core | ☐ |
| GE | General empirical relationship | ☐ (user-supplied closure) |
| US | User-specified curve | ☐ (table-driven head/loss) |
| SE | Special element | ☐ (extension hook) |

---

## 2. Heat-transport elements

First-class heat-transfer paths (today heat enters only via a constant `Node.heat_source`):

- **Convection** — wall↔fluid with `h·A·(T_wall − T_fluid)` correlations (Dittus–Boelter, etc.).
- **Conduction** — through heat structures (wall thermal mass + 1D conduction); ties to v0.4.
- **Radiation** — `σεA(T₁⁴ − T₂⁴)` between surfaces (high-temperature gas reactors).

These generalise the lumped `HeatExchanger` coupling into a proper heat-structure / surface
network coupled to the flow energy equation.

---

## 3. Fluids: mixtures → multiphase

- **Gas mixtures** — multi-species gas with mixture properties (mole/mass fractions, mixing
  rules for `R`, `cp`, `γ`, transport); species transport in the energy/continuity solve.
- **Multiphase mixtures** — then liquid+gas: homogeneous-equilibrium → drift-flux → two-fluid
  (the v0.7+ physics track). Needs quality/void handling and phase-change EOS.

Prerequisite: lift the current single-fluid-per-network restriction (needed for gas/liquid
heat exchangers, `CX`, and any mixed loop).

---

## 4. Fluid-property databases to develop

Behind the `FluidModel` ABC, a property layer (ideally with a table/correlation backend and
optional CoolProp/NIST REFPROP bridge):

- **Ideal gas** — ✅ (`IdealGas`, with compressibility factor `s`).
- **Incompressible liquid** — ✅ (`Incompressible`; add thermal expansion for liquid buoyancy).
- **Real gases** — helium, CO₂ (sCO₂ cycles), nitrogen, air, argon (real-gas EOS / tables).
- **Gas mixtures** — air as a mixture, combustion products, He–Xe, etc. (see §3).
- **Water / steam** — IAPWS-IF97 (subcooled, saturation, superheat; the two-phase workhorse).
- **Liquid metals** — sodium, lead, lead–bismuth eutectic (LBE) — for fast reactors.
- **Molten salts** — FLiBe, FLiNaK, chloride salts (MSR coolants/fuels).
- **Refrigerants / working fluids** — for cold-air units and bottoming cycles.
- **Transport properties** — viscosity, conductivity, Prandtl (for friction `f(Re)` and
  heat-transfer correlations), not just `ρ`/`cp`.

---

## 5. Tooling: web editor, solver-in-browser, doc export

- **Web-based network editor + file input** — a browser UI to lay out nodes/components and
  edit the declarative model (the `io/` JSON/YAML schema; finish round-trip serialization
  first). Build/query also drives the planned two-way LLM interface.
- **Online editor that runs the solver** — for **beginner users**: a hosted editor that
  invokes OpenTH and returns results/plots in the browser (no install). Pairs with the
  notebook/quickstart as the on-ramp; likely a small web service wrapping `FlowModel`/`Circuit`.
- **draw.io export** — export a network design to **draw.io** (diagrams.net) format for
  documentation, so models can be embedded in reports/specs. (Import from draw.io is a natural
  complement.)

---

These items feed the roadmap; the near-term sequence there is real-fluid properties → heat
structures/conjugate heat transfer → controls → point kinetics → the gas-cooled-loop MVP.
