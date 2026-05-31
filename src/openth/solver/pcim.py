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

        if not cfg.solve_energy:
            residual, iters = self._pressure_loop_steady(unknowns, index, n)
            return StepResult(time=0.0, iterations=iters, residual=residual,
                              converged=residual < cfg.tol)

        # Non-isothermal: alternate a full pressure/flow solve with one energy update,
        # until both the mass imbalance and the temperature change are below tolerance.
        # Energy is solved on a mass-conserving flow field, which keeps it well-posed.
        residual, total = math.inf, 0
        for _ in range(cfg.max_outer_iterations):
            residual, iters = self._pressure_loop_steady(unknowns, index, n)
            total += iters
            energy_change = self._solve_energy_steady()
            if residual < cfg.tol and energy_change < cfg.tol:
                return StepResult(time=0.0, iterations=total, residual=residual,
                                  converged=True)
        return StepResult(time=0.0, iterations=total, residual=residual, converged=False)

    def _pressure_loop_steady(self, unknowns: list[Node], index: dict[str, int],
                              n: int) -> tuple[float, int]:
        """Iterate flows + the pressure-correction equation to convergence (fixed T)."""
        cfg = self.config
        scale = self._mass_scale()
        residual = math.inf
        for it in range(1, cfg.max_outer_iterations + 1):
            self._update_flows_steady()
            rows: list[int] = []
            cols: list[int] = []
            data: list[float] = []
            b = np.zeros(n)
            kappa = {node.id: self._kappa(node) for node in unknowns}
            for node in unknowns:
                i = index[node.id]
                diag = 0.0
                outflow = -node.mass_source
                for e in self.network.elements_at(node):
                    d = self._conductance_steady(e)
                    other = e.downstream if e.upstream is node else e.upstream
                    mass_out = (1.0 if e.upstream is node else -1.0) * e.mdot
                    outflow += mass_out
                    # Diffusion (flow change s*dp'); compressible upwind convection of p'
                    # (density change rho'=drho/dp*p' of the convected mass -- eqs. 25-27),
                    # which scales with the mass flow and stabilises high-Mach flow.
                    diag += d
                    off = -d
                    if mass_out >= 0.0:
                        diag += mass_out * kappa[node.id]        # outflow: upwind is node i
                    else:
                        off += mass_out * self._kappa(other)     # inflow: upwind is neighbour
                    if not other.is_boundary:
                        rows.append(i)
                        cols.append(index[other.id])
                        data.append(off)
                rows.append(i)
                cols.append(i)
                data.append(diag)
                b[i] = -outflow
            residual = float(np.max(np.abs(b))) / scale
            if residual < cfg.tol:
                return residual, it
            a = csr_matrix((data, (rows, cols)), shape=(n, n))
            p_corr = sparse_solve(a, b)
            for node in unknowns:
                dp = cfg.relaxation * float(p_corr[index[node.id]])
                self._apply_pressure_correction(node, dp)
        return residual, cfg.max_outer_iterations

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
        p_o = {node.id: node.state.p0 for node in net.nodes.values()}
        h0_o = {node.id: node.state.h0 for node in net.nodes.values()}

        unknowns = net.solve_order()
        index = {node.id: i for i, node in enumerate(unknowns)}
        n = len(unknowns)
        if n == 0:
            return self._converge_flows_only(
                lambda: self._update_flows_transient(dt, alpha, mdot_o, resid_o), time=t)

        loop = (unknowns, index, n, dt, alpha, mdot_o, resid_o, outflow_o, rho_o)
        if not cfg.solve_energy:
            residual, iters = self._pressure_loop_transient(*loop)
            return StepResult(time=t, iterations=iters, residual=residual,
                              converged=residual < cfg.tol)

        residual, total = math.inf, 0
        for _ in range(cfg.max_outer_iterations):
            residual, iters = self._pressure_loop_transient(*loop)
            total += iters
            energy_change = self._solve_energy_transient(dt, alpha, rho_o, mdot_o, h0_o, p_o)
            if residual < cfg.tol and energy_change < cfg.tol:
                return StepResult(time=t, iterations=total, residual=residual, converged=True)
        return StepResult(time=t, iterations=total, residual=residual, converged=False)

    def _pressure_loop_transient(self, unknowns: list[Node], index: dict[str, int], n: int,
                                 dt: float, alpha: float, mdot_o: dict[str, float],
                                 resid_o: dict[str, float], outflow_o: dict[str, float],
                                 rho_o: dict[str, float]) -> tuple[float, int]:
        """Iterate flows + the transient pressure-correction equation to convergence."""
        cfg = self.config
        scale = self._mass_scale()
        residual = math.inf
        for it in range(1, cfg.max_outer_iterations + 1):
            self._update_flows_transient(dt, alpha, mdot_o, resid_o)
            rows: list[int] = []
            cols: list[int] = []
            data: list[float] = []
            b = np.zeros(n)
            kappa = {node.id: self._kappa(node) for node in unknowns}
            for node in unknowns:
                i = index[node.id]
                # Storage diagonal: V (drho/dp) / dt  (compressibility coupling, eq. 17).
                diag = node.volume * self._drho_dp(node) / dt
                outflow = self._net_outflow(node)
                for e in self.network.elements_at(node):
                    s = self._sensitivity_transient(e, dt, alpha)
                    other = e.downstream if e.upstream is node else e.upstream
                    mass_out = (1.0 if e.upstream is node else -1.0) * e.mdot
                    # Diffusion (flow change) + compressible upwind p' convection (density
                    # change of the convected mass), as in the steady solve.
                    diag += alpha * s
                    off = -alpha * s
                    if mass_out >= 0.0:
                        diag += alpha * mass_out * kappa[node.id]
                    else:
                        off += alpha * mass_out * self._kappa(other)
                    if not other.is_boundary:
                        rows.append(i)
                        cols.append(index[other.id])
                        data.append(off)
                rows.append(i)
                cols.append(i)
                data.append(diag)
                # Transient continuity residual G_i (eq. 13).
                storage = node.volume * (self._density(node) - rho_o[node.id]) / dt
                b[i] = -(storage + alpha * outflow + (1.0 - alpha) * outflow_o[node.id])
            residual = float(np.max(np.abs(b))) / scale
            if residual < cfg.tol:
                return residual, it
            a = csr_matrix((data, (rows, cols)), shape=(n, n))
            p_corr = sparse_solve(a, b)
            for node in unknowns:
                dp = cfg.relaxation * float(p_corr[index[node.id]])
                self._apply_pressure_correction(node, dp)
        return residual, cfg.max_outer_iterations

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
            # e.head() is a pump/compressor pressure rise; _gravity_head is the hydrostatic
            # term (buoyancy when densities differ) — both assist forward flow.
            drive = ((e.upstream.state.p0 - e.downstream.state.p0)
                     + e.head() + self._gravity_head(e, 0.5 * (rho_up + rho_down)) - conv)
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
            dp = (e.upstream.state.p0 - e.downstream.state.p0
                  + e.head() + self._gravity_head(e, 0.5 * (rho_up + rho_down)))
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
        dp_o = (e.upstream.state.p0 - e.downstream.state.p0
                + e.head() + self._gravity_head(e, 0.5 * (rho_up + rho_down)))
        return (1.0 - alpha) * (friction + conv - dp_o)

    def _gravity_head(self, e: Element, rho_face: float) -> float:
        """Hydrostatic pressure that assists forward flow: g * rho_face * (z_up - z_down).

        Density-weighted, so when two legs of a loop carry fluid at different temperatures
        (hence densities) their heads differ and the imbalance drives natural circulation.
        Zero when the gravity is off or the two ends are at the same elevation.
        """
        dz = e.upstream.elevation - e.downstream.elevation
        if dz == 0.0:
            return 0.0
        return self.config.gravity * rho_face * dz

    @staticmethod
    def _solve_momentum(a: float, b: float, rhs: float) -> float:
        """Solve the monotonic equation ``a*m + b*m*|m| = rhs`` for m (a, b >= 0)."""
        if b <= 0.0:
            return rhs / a if a > 0.0 else 0.0
        if a <= 0.0:
            return math.copysign(math.sqrt(abs(rhs) / b), rhs)
        mag = (-a + math.sqrt(a * a + 4.0 * b * abs(rhs))) / (2.0 * b)
        return math.copysign(mag, rhs)

    # ---- energy equation -----------------------------------------------------------

    def _solve_energy_steady(self) -> float:
        return self._solve_energy(dt=1.0, alpha=1.0, rho_o=None, mdot_o=None,
                                  h0_o=None, p_o=None)

    def _solve_energy_transient(self, dt: float, alpha: float, rho_o: dict[str, float],
                                mdot_o: dict[str, float], h0_o: dict[str, float],
                                p_o: dict[str, float]) -> float:
        return self._solve_energy(dt, alpha, rho_o, mdot_o, h0_o, p_o)

    def _solve_energy(self, dt: float, alpha: float, rho_o: dict[str, float] | None,
                      mdot_o: dict[str, float] | None, h0_o: dict[str, float] | None,
                      p_o: dict[str, float] | None) -> float:
        """Assemble and solve the total-enthalpy transport equation (eqs. 29-34).

        Upwind convection of h0 by the (already-solved) face mass flows, plus a storage
        term and heat input. ``rho_o is None`` selects the steady form (no time terms,
        convection fully weighted). Updates each unknown node's h0 and temperature, and
        returns the maximum relative temperature change (for convergence).
        """
        net = self.network
        steady = rho_o is None
        conv_w = 1.0 if steady else alpha

        # Temperature-Dirichlet nodes (pressure boundaries, fixed-T inlets): set h0 from
        # the fixed temperature plus the local kinetic energy.
        for node in net.nodes.values():
            if node.is_boundary or node.fixed_temperature:
                fluid = self._node_fluid(node)
                v = self._node_velocity(node)
                node.state.h0 = fluid.enthalpy(node.state.T) + 0.5 * v * v

        unknowns = [n for n in net.nodes.values()
                    if not (n.is_boundary or n.fixed_temperature)]
        if not unknowns:
            return 0.0
        index = {node.id: i for i, node in enumerate(unknowns)}
        m = len(unknowns)

        rows: list[int] = []
        cols: list[int] = []
        data: list[float] = []
        rhs = np.zeros(m)

        for node in unknowns:
            i = index[node.id]
            kp = (self._density(node) * node.volume / dt) if not steady else 0.0
            r = node.heat_source
            if not steady:
                assert rho_o is not None and h0_o is not None and p_o is not None
                r += node.volume * rho_o[node.id] * h0_o[node.id] / dt
                r += node.volume * (node.state.p0 - p_o[node.id]) / dt  # pressure work
            old_flux = 0.0
            for e in net.elements_at(node):
                sign = 1.0 if e.upstream is node else -1.0
                other = e.downstream if e.upstream is node else e.upstream
                mdot_out = sign * e.mdot
                if mdot_out >= 0.0:
                    kp += conv_w * mdot_out  # outflow: upwind enthalpy is this node's
                else:
                    k_in = conv_w * (-mdot_out)  # inflow from neighbour: upwind is other
                    if other.is_boundary or other.fixed_temperature:
                        r += k_in * other.state.h0
                    else:
                        rows.append(i)
                        cols.append(index[other.id])
                        data.append(-k_in)
                    # Pump/compressor work added to the fluid arriving through this face.
                    r += k_in * e.work_per_mass(e.mdot, self._density(other))
                if not steady:
                    assert mdot_o is not None and h0_o is not None and rho_o is not None
                    out_o = sign * mdot_o[e.id]
                    if out_o >= 0.0:
                        old_flux += out_o * h0_o[node.id]
                    else:
                        w_o = e.work_per_mass(mdot_o[e.id], rho_o[other.id])
                        old_flux += out_o * (h0_o[other.id] + w_o)
            if not steady:
                r -= (1.0 - alpha) * old_flux
            if kp < self._MDOT_FLOOR * self._mass_scale():
                # Near-stagnant node with negligible storage: the h0 balance is degenerate,
                # so leave its enthalpy unchanged (identity row) rather than divide by ~0.
                rows.append(i)
                cols.append(i)
                data.append(1.0)
                rhs[i] = node.state.h0
            else:
                rows.append(i)
                cols.append(i)
                data.append(kp)
                rhs[i] = r

        # Heat-exchanger couplings: a conductive link g = UA/cp on the total enthalpy that
        # transfers heat hot -> cold (implicit, energy-conserving; ignores the small kinetic
        # part of h0 at the low-velocity exchanger nodes).
        for hx in net.heat_exchangers.values():
            g = hx.UA / self._cp(hx.hot)
            for node, other in ((hx.cold, hx.hot), (hx.hot, hx.cold)):
                if node.id not in index:
                    continue  # this side's temperature is pinned (boundary/fixed-T)
                i = index[node.id]
                rows.append(i)
                cols.append(i)
                data.append(g)
                if other.id in index:
                    rows.append(i)
                    cols.append(index[other.id])
                    data.append(-g)
                else:
                    rhs[i] += g * other.state.h0

        matrix = csr_matrix((data, (rows, cols)), shape=(m, m))
        h0_new = sparse_solve(matrix, rhs)

        max_change = 0.0
        relax = self.config.relaxation
        for node in unknowns:
            fluid = self._node_fluid(node)
            v = self._node_velocity(node)
            updated = node.state.h0 + relax * (float(h0_new[index[node.id]]) - node.state.h0)
            node.state.h0 = updated
            t_old = node.state.T
            t_new = fluid.temperature_from_enthalpy(updated - 0.5 * v * v)
            # Limit the temperature change per outer iteration: the kinetic term couples
            # T <-> rho <-> V, which can run away (low rho -> high V -> negative T) far from
            # the solution. Bounding the step keeps T strictly positive and the coupling
            # stable; it does not affect the converged result.
            t_new = min(max(t_new, 0.7 * t_old), 1.3 * t_old)
            node.state.T = t_new
            max_change = max(max_change, abs(t_new - t_old) / t_old)
        return max_change

    def _node_velocity(self, node: Node) -> float:
        """Representative speed at a node for the kinetic part of h0 (max face velocity)."""
        rho = self._density(node)
        speed = 0.0
        for e in self.network.elements_at(node):
            area = e.flow_area()
            if math.isfinite(area) and area > 0.0:
                speed = max(speed, abs(e.mdot) / (rho * area))
        return speed

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

    def _cp(self, node: Node) -> float:
        """Specific heat cp = dh/dT [J/(kg.K)] of the node's fluid (h is linear in T here)."""
        fluid = self._node_fluid(node)
        return fluid.enthalpy(1.0) - fluid.enthalpy(0.0)

    @staticmethod
    def _apply_pressure_correction(node: Node, dp: float) -> None:
        """Apply a pressure correction, but keep the pressure positive: a too-large negative
        correction would drive density <= 0 (and the resistance singular). Capping the drop
        per iteration keeps the solve stable for stiff / disconnected networks; it does not
        change the converged result (where the correction -> 0)."""
        node.state.p0 = max(node.state.p0 + dp, 0.1 * node.state.p0)

    def _drho_dp(self, node: Node) -> float:
        return self._node_fluid(node).drho_dp(node.state.p0, node.state.T)

    def _kappa(self, node: Node) -> float:
        """Relative compressibility (drho/dp)/rho at a node [1/Pa].

        The coefficient of the upwind pressure-correction convection term that makes the
        continuity equation compressible (eqs. 25-27); for an ideal gas it is 1/p.
        """
        rho = self._density(node)
        return self._drho_dp(node) / rho if rho > 0.0 else 0.0

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
        if self.config.solve_energy:
            # Initialise total enthalpy from temperature (kinetic part added once flows are
            # known); only fill unset values so a committed transient state is preserved.
            for node in nodes:
                if node.state.h0 <= 0.0:
                    node.state.h0 = self._node_fluid(node).enthalpy(node.state.T)

    def _mass_scale(self) -> float:
        """A representative mass-flow magnitude for relative tolerances / flooring."""
        sources = [abs(n.mass_source) for n in self.network.nodes.values() if n.mass_source]
        flows = [abs(e.mdot) for e in self.network.elements.values()]
        scale = max(sources + flows, default=0.0)
        return scale if scale > 0.0 else 1.0
