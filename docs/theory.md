# Theory: the implicit Pressure-Correction (PC / PCIM) method

Working notes mapping Greyvenstein (2002) onto the FlowCalc code. The source paper is in
[`papers/`](papers/Greyvenstein-2001-implicit-transient-pipe-networks.pdf). Equation
numbers below refer to that paper.

## The idea in one paragraph

Transient compressible flow in a network is governed by **continuity**, **momentum** and
**energy** PDEs (eqs. 1–3). Discretized with a finite-volume scheme (Fig. 1), they become
a coupled algebraic system for pressure, flow, density and enthalpy at the new time level.
Rather than solve that monolithically, the method is **segregated** (SIMPLE-style): from a
guessed pressure field it computes preliminary flows and densities, then derives a
**pressure-correction equation** (from continuity) whose solution nudges the fields toward
satisfying all the equations simultaneously. Iterating to convergence at each time step
gives an implicit scheme that is stable for large time steps — its key advantage over
explicit methods (MOC, Lax–Wendroff) on slow transients.

## Discretization layout (Fig. 1) → code

| Quantity | Lives at | Code |
|----------|----------|------|
| pressure `p₀`, density `ρ`, temperature `T`, enthalpy `h₀` | **cell centre** | `network.Node` / `NodeState` |
| volumetric flow `Q`, velocity `V` | **cell face** | `network.Element` (e.g. `components.Pipe`) |

Subscripts `e`/`w` = the east/west (downstream/upstream) faces of a cell; superscript `o`
= previous time level; a prime `'` = a correction; an overbar = a preliminary value.

## Governing equations

- **Continuity** (eq. 1), integrated over a control volume → eq. (13).
- **Momentum** (eq. 2). The convective-acceleration term is *retained* but eliminated by
  recasting in terms of **total pressure** `p₀` (eqs. 4–10): for liquids via eq. (6), for
  gases via eq. (10) with the **effective friction factor** `f̃` (eq. 11). Integrated over
  a face-centred control volume → eq. (14).
- **Energy** (eq. 3), integrated → eq. (15).
- **Equation of state** for gas: `p = s ρ R T` (eq. 12), with compressibility factor `s`.
  → `fluids.IdealGas`. Liquids use constant density → `fluids.Incompressible`.

## Solution algorithm (Section 4) → `solver.PCIMSolver`

1. Preliminary pressures `p̄₀` (previous time level / iteration).
2. Preliminary flows `Q̄` from momentum (eq. 14); preliminary densities `ρ̄` from eq. (12).
3. The corrected flow on a face is **linear in the end-node pressure corrections**:
   `Q'ᵢ = a⁻·p'ᵢ − a⁺·p'ᵢ₊₁` (eq. 20), with link coefficients `a⁺`, `a⁻` (eqs. 21–22).
   → `Element.momentum_coeffs` returning `MomentumCoeffs(a_plus, a_minus)`.
4. Substituting into continuity gives the **pressure-correction equation**
   `cP·p'ᵢ = cE·p'ᵢ₊₁ + cW·p'ᵢ₋₁ + bᵢ` (eqs. 24–28).
5. Solve it: **Thomas algorithm** for a series pipeline, **sparse** solve for networks
   (`solver.linear.thomas` / `sparse_solve`).
6. Update `p₀, Q, ρ` with the corrections (eqs. 16–20); iterate (inner loop).
7. Solve the **energy equation** for `h₀` (eqs. 29–34, upwinded), tridiagonal/sparse.
8. Iterate momentum/continuity ↔ energy to convergence; advance to the next step.

## Time integration

Weighing factor **α ∈ [0.5, 1]** (`SolverConfig.alpha`):

- α = 1 → fully implicit, 1st-order in time (more numerical damping).
- α = 0.5 → Crank–Nicolson, 2nd-order, but the scheme becomes **unstable near 0.5**.
- **α = 0.6** → the paper's recommended accuracy/stability compromise (default in code).

The scheme conserves mass (important for closed loops like the PBMR Brayton cycle that
motivated the method).

## Implementation status

The **steady-state solve** (`PCIMSolver.steady_state`, the dt→∞ limit), the **transient
time step** (`PCIMSolver.step`), and the **energy equation** (`solve_energy=True`) are all
implemented as segregated SIMPLE iterations with Picard density updates. The transient adds
the finite-volume storage term `V(ρ−ρᵒ)/Δt`, the momentum inertia `(Δx/A)(ṁ−ṁᵒ)/Δt`, and
the θ-weighting `α`; its pressure-correction coefficients are `s = α/(I/Δt + 2αK|ṁ|)`,
`cP = V(∂ρ/∂p)/Δt + α·Σsₑ`, `cnb = α·sₑ`.

The **energy solve** transports total enthalpy `h₀` by upwind convection on the
(already-converged, mass-conserving) flow field (eqs. 29–34): per node
`kₚ h₀ᵢ = Σ kⱼ h₀ⱼ + rᵢ` with `k = (α or 1)·max(±ṁ, 0)`, a storage diagonal `Vρ/Δt`, and a
RHS carrying heat `Q̇`, pressure work `V(p−pᵒ)/Δt`, and old-time fluxes. Temperature comes
back as `T = (h₀ − ½V²)/cₚ`. The non-isothermal outer loop alternates a full pressure solve
with one energy update until both converge.

Validated (`tests/`): steady matches the closed-form isothermal pipe law to ~0%; the
transient marches to the steady fixed point (~1e-9), conserves mass exactly in the
θ-weighted sense (~1e-9), and produces water-hammer + the blow-down decay; the energy solve
conserves h₀ in adiabatic flow (~1e-9, with expansion cooling) and matches `Q/ṁ` for heat
addition. The pressure↔temperature coupling is robust below ~Mach 0.4 for pressure-driven
flow (higher-Mach robustness is future work).

## Benchmarks to validate against (Section 5)

1. **Steady isothermal/non-isothermal pipeline** — compare pressure ratio vs. outlet Mach
   number to the ODE benchmark (eqs. 35–37). Implemented + validated:
   [`examples/pipeline_steady.py`](../examples/pipeline_steady.py). (Fig. 2/3.) ✅
2. **Sudden valve closure** in a 20 m pipe — pressure-wave amplitude/frequency vs. MOC and
   Lax–Wendroff (Fig. 4–7). Water-hammer overpressure validated in
   `tests/test_transient.py` (isothermal; quantitative match to the figures pending). ◑
3. **Branching network** valve closures (Fig. 8–9).
4. **Pressure-vessel blow-down** — slow transient where the implicit method is ~70–1000×
   faster than explicit methods (Fig. 10). Implemented + checked vs. quasi-steady:
   [`examples/blowdown_transient.py`](../examples/blowdown_transient.py). ✅
5. **Non-isothermal flow** (Fig. 3, 6, 7) — energy equation. Adiabatic h₀ conservation
   (expansion cooling) and heat addition (`Q/ṁ`) validated in `tests/test_energy.py`;
   demo in [`examples/heated_pipe.py`](../examples/heated_pipe.py). ◑ (robust below ~M 0.4)

## Beyond pipes

The momentum closure (step 3) is per-element, so **non-pipe components** — valves, pumps,
compressors, turbines, orifices, heat exchangers — slot in by supplying their own
`resistance` / `convective_dp` / `inertance`. `components.Valve` is the first such template.
