"""PCIM: the implicit Pressure-Correction solver (Greyvenstein 2002).

Implements both the **steady-state** solve and the **transient** time step. Both share the
same segregated (SIMPLE) inner loop; the transient adds the finite-volume storage/inertia
terms and the paper's time-integration weighing factor ``alpha`` (Section 3).

Per face, the momentum balance at the new time level is (paper eq. 14, in mass-flow form):

    I (mdot - mdot_o)/dt + alpha*[ K mdot|mdot| + C - (p_up - p_down) ]
                         + (1-alpha)*[ K mdot|mdot| + C - (p_up - p_down) ]_o  = 0

where I = dx/A is the inertance, K the friction resistance, C the convective term. Solving
this for mdot given the current pressures, and substituting into the nodal continuity
balance (eq. 13)

    V (rho - rho_o)/dt + alpha*O + (1-alpha)*O_o = 0,   O = net mass outflow,

gives the pressure-correction equation (eqs. 24-28). Its coefficients are

    s   = alpha / (I/dt + 2 alpha K |mdot|)         (flow sensitivity to dp, eq. 20)
    cP  = V (drho/dp)/dt + alpha * sum_e s_e         (storage + link diagonal)
    cnb = alpha * s_e
    b   = -(continuity residual)

The steady solve is the dt -> infinity limit (storage and inertia drop out, s -> 1/(2K|m|)).
Densities are refreshed each iteration (Picard) so the equation of state, momentum and
continuity are all satisfied at convergence. Temperature is held fixed (isothermal); the
energy equation (eqs. 29-34) is not yet coupled.
"""

from __future__ import annotations

import math

import numpy as np
from scipy.sparse import csr_matrix

from ..network.element import Element
from ..network.node import Node
from .base import Solver, StepResult
from .linear import sparse_solve


class PCIMSolver(Solver):
    """Segregated implicit pressure-correction solver."""

    # Floor on |mdot| (relative to the characteristic mass flow) when forming the
    # linearised conductance d = 1/(2K|mdot|). Quadratic resistance makes d blow up as
    # mdot -> 0, which makes near-zero-flow branches ill-conditioned; this keeps the system
    # well-conditioned without affecting the converged result (where flows are finite).
    _MDOT_FLOOR = 1e-3

    # ---- steady state --------------------------------------------------------------

    def steady_state(self) -> StepResult:
        net = self.network
        net.validate()
        self._initialise_state()

        unknowns = net.solve_order()  # nodes whose pressure is solved (not Dirichlet)
        index = {node.id: i for i, node in enumerate(unknowns)}
        n = len(unknowns)
        cfg = self.config
        if n == 0:
            # No pressure unknowns (e.g. one element between two pressure boundaries):
            # the flows are fixed by the boundary pressures; just converge them.
            return self._converge_flows_only(self._update_flows_steady, time=0.0)

        residual = math.inf
        scale = self._mass_scale()

        for it in range(1, cfg.max_outer_iterations + 1):
            self._update_flows_steady()

            rows: list[int] = []
            cols: list[int] = []
            data: list[float] = []
            b = np.zeros(n)

            for node in unknowns:
                i = index[node.id]
                diag = 0.0
                outflow = -node.mass_source
                for e in net.elements_at(node):
                    d = self._conductance_steady(e)
                    diag += d
                    other = e.downstream if e.upstream is node else e.upstream
                    sign = 1.0 if e.upstream is node else -1.0
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

    # ---- transient time step -------------------------------------------------------

    def step(self, dt: float, t: float) -> StepResult:
        net = self.network
        net.validate()
        self._initialise_state()
        for e in net.elements.values():
            e.set_time(t)  # update time-dependent elements (e.g. a closing valve)

        cfg = self.config
        alpha = cfg.alpha

        # Snapshot previous-time-level (superscript o) quantities from the committed state.
        rho_o = {node.id: self._density(node) for node in net.nodes.values()}
        mdot_o = {e.id: e.mdot for e in net.elements.values()}
        resid_o = {e.id: self._old_momentum_residual(e, rho_o, mdot_o, alpha)
                   for e in net.elements.values()}
        outflow_o = {node.id: self._net_outflow(node, mdot_o) for node in net.nodes.values()}

        unknowns = net.solve_order()
        index = {node.id: i for i, node in enumerate(unknowns)}
        n = len(unknowns)
        if n == 0:
            return self._converge_flows_only(
                lambda: self._update_flows_transient(dt, alpha, mdot_o, resid_o), time=t)

        residual = math.inf
        scale = self._mass_scale()

        for it in range(1, cfg.max_outer_iterations + 1):
            self._update_flows_transient(dt, alpha, mdot_o, resid_o)

            rows: list[int] = []
            cols: list[int] = []
            data: list[float] = []
            b = np.zeros(n)

            for node in unknowns:
                i = index[node.id]
                # Storage diagonal: V (drho/dp) / dt  (compressibility coupling, eq. 17).
                diag = node.volume * self._drho_dp(node) / dt
                outflow = self._net_outflow(node)
                for e in net.elements_at(node):
                    s = self._sensitivity_transient(e, dt, alpha)
                    diag += alpha * s
                    other = e.downstream if e.upstream is node else e.upstream
                    if not other.is_boundary:
                        rows.append(i)
                        cols.append(index[other.id])
                        data.append(-alpha * s)
                rows.append(i)
                cols.append(i)
                data.append(diag)
                # Transient continuity residual G_i (eq. 13).
                storage = node.volume * (self._density(node) - rho_o[node.id]) / dt
                b[i] = -(storage + alpha * outflow + (1.0 - alpha) * outflow_o[node.id])

            residual = float(np.max(np.abs(b))) / scale
            if residual < cfg.tol:
                return StepResult(time=t, iterations=it, residual=residual, converged=True)

            a = csr_matrix((data, (rows, cols)), shape=(n, n))
            p_corr = sparse_solve(a, b)
            for node in unknowns:
                node.state.p0 += cfg.relaxation * float(p_corr[index[node.id]])

        return StepResult(time=t, iterations=cfg.max_outer_iterations,
                          residual=residual, converged=False)

    def _converge_flows_only(self, update_flows, time: float) -> StepResult:
        """Iterate the face-flow update to convergence when there are no pressure unknowns.

        Applies when every node has a fixed pressure (e.g. a single element between two
        pressure boundaries): momentum alone determines the flows, so we iterate the
        (under-relaxed, convective-lagged) flow update until it stops changing.
        """
        cfg = self.config
        change = math.inf
        for it in range(1, cfg.max_outer_iterations + 1):
            previous = {e.id: e.mdot for e in self.network.elements.values()}
            update_flows()
            delta = max((abs(e.mdot - previous[e.id]) for e in self.network.elements.values()),
                        default=0.0)
            change = delta / self._mass_scale()  # scale reflects the now-updated flows
            if change < cfg.tol:
                return StepResult(time=time, iterations=it, residual=change, converged=True)
        return StepResult(time=time, iterations=cfg.max_outer_iterations,
                          residual=change, converged=False)

    # ---- flow updates --------------------------------------------------------------

    def _update_flows_steady(self) -> None:
        """Recompute each face's mass flow from the current pressures (steady momentum).

        The flow is under-relaxed: the convective term couples flow to itself, which is a
        positive feedback that can diverge at high Mach in pressure-driven problems.
        """
        relax = self.config.relaxation
        for e in self.network.elements.values():
            rho_up = e.density_at(e.upstream)
            rho_down = e.density_at(e.downstream)
            k = e.resistance(0.5 * (rho_up + rho_down))
            if not math.isfinite(k) or k <= 0.0:
                e.mdot = 0.0
                continue
            conv = e.convective_dp(rho_up, rho_down, e.mdot)
            drive = (e.upstream.state.p0 - e.downstream.state.p0) - conv
            target = math.copysign(math.sqrt(abs(drive) / k), drive)
            e.mdot += relax * (target - e.mdot)

    def _update_flows_transient(self, dt: float, alpha: float, mdot_o: dict[str, float],
                                resid_o: dict[str, float]) -> None:
        """Solve each face's transient momentum balance for the mass flow."""
        for e in self.network.elements.values():
            rho_up = e.density_at(e.upstream)
            rho_down = e.density_at(e.downstream)
            k = e.resistance(0.5 * (rho_up + rho_down))
            if not math.isfinite(k) or k <= 0.0:
                e.mdot = 0.0
                continue
            a = e.inertance() / dt          # inertia coefficient
            quad = alpha * k                # quadratic friction coefficient
            conv = e.convective_dp(rho_up, rho_down, e.mdot)
            dp = e.upstream.state.p0 - e.downstream.state.p0
            # a*mdot + quad*mdot|mdot| = a*mdot_o + alpha*dp - alpha*conv - resid_o
            rhs = a * mdot_o[e.id] + alpha * dp - alpha * conv - resid_o[e.id]
            e.mdot = self._solve_momentum(a, quad, rhs)

    # ---- coefficient helpers -------------------------------------------------------

    def _conductance_steady(self, e: Element) -> float:
        """Linearised friction conductance d = 1 / (2 K |mdot|) for face ``e`` (eq. 20)."""
        k = e.resistance(e.face_density())
        if not math.isfinite(k) or k <= 0.0:
            return 0.0
        return 1.0 / (2.0 * k * max(abs(e.mdot), self._MDOT_FLOOR * self._mass_scale()))

    def _sensitivity_transient(self, e: Element, dt: float, alpha: float) -> float:
        """Flow sensitivity s = alpha / (I/dt + 2 alpha K |mdot|) for face ``e``."""
        k = e.resistance(e.face_density())
        if not math.isfinite(k) or k <= 0.0:
            return 0.0
        a = e.inertance() / dt
        floor = self._MDOT_FLOOR * self._mass_scale()
        return alpha / (a + 2.0 * alpha * k * max(abs(e.mdot), floor))

    def _old_momentum_residual(self, e: Element, rho_o: dict[str, float],
                               mdot_o: dict[str, float], alpha: float) -> float:
        """The (1-alpha)-weighted old-time momentum term carried into this step."""
        rho_up = rho_o[e.upstream.id]
        rho_down = rho_o[e.downstream.id]
        k = e.resistance(0.5 * (rho_up + rho_down))
        m_o = mdot_o[e.id]
        if not math.isfinite(k) or k <= 0.0:
            friction = 0.0
            conv = 0.0
        else:
            friction = k * m_o * abs(m_o)
            conv = e.convective_dp(rho_up, rho_down, m_o)
        dp_o = e.upstream.state.p0 - e.downstream.state.p0
        return (1.0 - alpha) * (friction + conv - dp_o)

    @staticmethod
    def _solve_momentum(a: float, b: float, rhs: float) -> float:
        """Solve the monotonic equation ``a*m + b*m*|m| = rhs`` for m (a, b >= 0)."""
        if b <= 0.0:
            return rhs / a if a > 0.0 else 0.0
        if a <= 0.0:
            return math.copysign(math.sqrt(abs(rhs) / b), rhs)
        mag = (-a + math.sqrt(a * a + 4.0 * b * abs(rhs))) / (2.0 * b)
        return math.copysign(mag, rhs)

    # ---- state helpers -------------------------------------------------------------

    def _net_outflow(self, node: Node, mdot: dict[str, float] | None = None) -> float:
        """Net mass outflow O_i from a node: signed sum of face flows minus its source.

        Uses ``mdot[e.id]`` when a map is given (e.g. the previous-time-level flows),
        otherwise the elements' current ``mdot``.
        """
        out = -node.mass_source
        for e in self.network.elements_at(node):
            flow = e.mdot if mdot is None else mdot[e.id]
            out += (1.0 if e.upstream is node else -1.0) * flow
        return out

    def _node_fluid(self, node: Node):
        for e in self.network.elements_at(node):
            if e.fluid is not None:
                return e.fluid
        raise ValueError(f"node {node.id!r} has no adjacent element with a fluid")

    def _density(self, node: Node) -> float:
        return self._node_fluid(node).density(node.state.p0, node.state.T)

    def _drho_dp(self, node: Node) -> float:
        return self._node_fluid(node).drho_dp(node.state.p0, node.state.T)

    def _initialise_state(self) -> None:
        """Fill in any unset pressures/temperatures from the fixed boundaries."""
        nodes = list(self.network.nodes.values())
        fixed_p = [n.state.p0 for n in nodes if n.is_boundary and n.state.p0 > 0.0]
        temps = [n.state.T for n in nodes if n.state.T > 0.0]
        if not fixed_p:
            raise ValueError("solver needs at least one PressureBoundary")
        if not temps:
            raise ValueError("solver needs a temperature on at least one node")
        p_guess = sum(fixed_p) / len(fixed_p)
        t_guess = sum(temps) / len(temps)
        for node in nodes:
            if node.state.p0 <= 0.0:
                node.state.p0 = p_guess
            if node.state.T <= 0.0:
                node.state.T = t_guess  # isothermal fill

    def _mass_scale(self) -> float:
        """A representative mass-flow magnitude for relative tolerances / flooring."""
        sources = [abs(n.mass_source) for n in self.network.nodes.values() if n.mass_source]
        flows = [abs(e.mdot) for e in self.network.elements.values()]
        scale = max(sources + flows, default=0.0)
        return scale if scale > 0.0 else 1.0
