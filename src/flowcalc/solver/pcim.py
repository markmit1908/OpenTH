"""PCIM: the implicit Pressure-Correction solver (Greyvenstein 2002).

This is the heart of FlowCalc. The algorithm is segregated (SIMPLE-style); one time step
proceeds as:

    1. Take preliminary pressures p0_bar (previous time level, or previous iteration).
    2. Solve preliminary flows Q_bar from the momentum equation (eq. 14) and preliminary
       densities rho_bar from the equation of state (eq. 12).
    3. Assemble the pressure-correction equation from continuity (eqs. 24-28):
           cP * p'_i = cE * p'_{i+1} + cW * p'_{i-1} + b_i
       Solve it (Thomas for a series pipeline; sparse for networks).
    4. Update p0, Q, rho with the corrections (eqs. 16-20).
    5. Repeat 2-4 (inner pressure iterations) until the correction is small.
    6. Solve the energy equation for h0 (eqs. 29-34), again tridiagonal/sparse.
    7. Repeat 2-6 (outer iterations) until convergence.
    8. Commit to the new time level and advance.

Steady state is the same loop with the transient (dt) terms switched off.

The numerical core (steps 2-6) is intentionally unimplemented: it must be transcribed
carefully from ``docs/papers/`` and validated against the paper's benchmarks. The control
flow, configuration and linear-solve back ends around it are in place.
"""

from __future__ import annotations

from .base import Solver, StepResult


class PCIMSolver(Solver):
    """Segregated implicit pressure-correction solver."""

    def steady_state(self) -> StepResult:
        self.network.validate()
        # TODO(core): run the segregated loop with transient terms disabled.
        raise NotImplementedError(
            "PCIMSolver.steady_state: implement the pressure-correction loop "
            "(eqs. 12-28) and validate against eqs. (35)-(37)."
        )

    def step(self, dt: float, t: float) -> StepResult:
        self.network.validate()
        # TODO(core): one transient step -- momentum/continuity pressure-correction inner
        # loop (eqs. 13-28) alternating with the energy solve (eqs. 29-34) until converged.
        raise NotImplementedError(
            "PCIMSolver.step: implement the transient pressure-correction time step "
            "(eqs. 13-34)."
        )
