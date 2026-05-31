"""PCIM: the implicit Pressure-Correction solver (Greyvenstein 2002).

This module implements the **steady-state** core. The transient time-marching step is
still pending (see :meth:`PCIMSolver.step`); the steady solve below is also the natural
initial condition for a future transient run, exactly as in the paper (Section 5.1).

Steady-state algorithm (segregated / SIMPLE, the dt -> infinity limit of Section 4):

    repeat until the mass imbalance is below tolerance:
      1. Evaluate face densities from the current node pressures/temperatures (eq. 12).
      2. For each face, satisfy the momentum balance for the current pressures:
             p_up - p_down = K * mdot|mdot| + C
         -> mdot* = sign(dp - C) * sqrt(|dp - C| / K)
         where K is the friction resistance and C the convective term (Element).
      3. Linearise: mdot' = d * (p'_up - p'_down), d = 1 / (2 K |mdot*|)  (eqs. 20-22).
      4. Assemble the pressure-correction equation from the nodal mass balance
         (eqs. 24-28, transient terms dropped):
             (sum_e d_e) p'_i  -  sum_nb d_e p'_nb  =  -O_i
         where O_i is the net mass outflow from node i (the continuity residual).
      5. Solve the linear system for p' (sparse; Thomas for a pure series chain).
      6. Update pressures p_i += relaxation * p'_i and the face flows.

Densities are refreshed each iteration (Picard), so at convergence the equation of state,
momentum and continuity are all satisfied simultaneously. The steady solve assumes a
fixed temperature field (isothermal / energy equation not yet coupled here).
"""

from __future__ import annotations

import math

import numpy as np
from scipy.sparse import csr_matrix

from .base import Solver, StepResult
from .linear import sparse_solve


class PCIMSolver(Solver):
    """Segregated implicit pressure-correction solver."""

    # Relative floor on |mdot| when forming the linearised conductance, to avoid a
    # singular system on the first iterations when flows are near zero.
    _MDOT_FLOOR = 1e-9

    def steady_state(self) -> StepResult:
        net = self.network
        net.validate()
        self._initialise_state()

        unknowns = net.solve_order()  # nodes whose pressure is solved (not Dirichlet)
        index = {node.id: i for i, node in enumerate(unknowns)}
        n = len(unknowns)
        if n == 0:
            return StepResult(time=0.0, iterations=0, residual=0.0, converged=True)

        cfg = self.config
        residual = math.inf
        scale = self._mass_scale()

        for it in range(1, cfg.max_outer_iterations + 1):
            self._update_flows()

            # Assemble the pressure-correction system: A p' = b.
            rows: list[int] = []
            cols: list[int] = []
            data: list[float] = []
            b = np.zeros(n)

            for node in unknowns:
                i = index[node.id]
                diag = 0.0
                # Continuity residual O_i = net mass outflow - imposed source.
                outflow = -node.mass_source
                for e in net.elements_at(node):
                    d = self._conductance(e)
                    diag += d
                    other = e.downstream if e.upstream is node else e.upstream
                    sign = 1.0 if e.upstream is node else -1.0  # +mdot leaves an upstream node
                    outflow += sign * e.mdot
                    if not other.is_boundary:
                        rows.append(i)
                        cols.append(index[other.id])
                        data.append(-d)
                rows.append(i)
                cols.append(i)
                data.append(diag)
                b[i] = -outflow

            residual = float(np.max(np.abs(b))) / scale
            if residual < cfg.tol:
                return StepResult(time=0.0, iterations=it, residual=residual, converged=True)

            a = csr_matrix((data, (rows, cols)), shape=(n, n))
            p_corr = sparse_solve(a, b)

            for node in unknowns:
                node.state.p0 += cfg.relaxation * float(p_corr[index[node.id]])

        return StepResult(time=0.0, iterations=cfg.max_outer_iterations,
                          residual=residual, converged=False)

    def step(self, dt: float, t: float) -> StepResult:
        self.network.validate()
        # TODO(core): transient time step -- add the storage terms (V/dt, previous-time
        # level) to continuity/momentum and couple the energy equation (eqs. 13-34).
        raise NotImplementedError(
            "PCIMSolver.step: transient time-marching not implemented yet; steady_state() works."
        )

    # ---- helpers -------------------------------------------------------------------

    def _initialise_state(self) -> None:
        """Fill in any unset pressures/temperatures from the fixed boundaries."""
        nodes = list(self.network.nodes.values())
        fixed_p = [n.state.p0 for n in nodes if n.is_boundary and n.state.p0 > 0.0]
        temps = [n.state.T for n in nodes if n.state.T > 0.0]
        if not fixed_p:
            raise ValueError("steady_state needs at least one PressureBoundary")
        if not temps:
            raise ValueError("steady_state needs a temperature on at least one node")
        p_guess = sum(fixed_p) / len(fixed_p)
        t_guess = sum(temps) / len(temps)
        for node in nodes:
            if node.state.p0 <= 0.0:
                node.state.p0 = p_guess
            if node.state.T <= 0.0:
                node.state.T = t_guess  # isothermal fill

    def _update_flows(self) -> None:
        """Recompute each face's mass flow from the current pressures (step 2 above)."""
        for e in self.network.elements.values():
            rho_up = e.density_at(e.upstream)
            rho_down = e.density_at(e.downstream)
            rho_face = 0.5 * (rho_up + rho_down)
            k = e.resistance(rho_face)
            if not math.isfinite(k) or k <= 0.0:  # e.g. a shut valve
                e.mdot = 0.0
                continue
            conv = e.convective_dp(rho_up, rho_down, e.mdot)
            drive = (e.upstream.state.p0 - e.downstream.state.p0) - conv
            e.mdot = math.copysign(math.sqrt(abs(drive) / k), drive)

    def _conductance(self, e) -> float:
        """Linearised friction conductance d = 1 / (2 K |mdot|) for face ``e`` (eq. 20)."""
        rho_face = e.face_density()
        k = e.resistance(rho_face)
        if not math.isfinite(k) or k <= 0.0:
            return 0.0
        mdot_floor = self._MDOT_FLOOR * self._mass_scale()
        return 1.0 / (2.0 * k * max(abs(e.mdot), mdot_floor))

    def _mass_scale(self) -> float:
        """A representative mass-flow magnitude for relative tolerances / flooring."""
        sources = [abs(n.mass_source) for n in self.network.nodes.values() if n.mass_source]
        flows = [abs(e.mdot) for e in self.network.elements.values()]
        scale = max(sources + flows, default=0.0)
        return scale if scale > 0.0 else 1.0
