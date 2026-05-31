"""High-level model-building facade.

`FlowModel` wraps the low-level :class:`~flowcalc.network.Network` / `Element` / boundary
objects in a fluent, name-based API so a user can describe a network in a few lines instead
of wiring nodes and elements by hand:

    model = FlowModel(fluid=helium())
    model.add_pipe("inlet", "outlet", length=100, diameter=0.5, n_cells=20)
    model.pressure_boundary("outlet", p=200e3, T=300)
    model.mass_flow_boundary("inlet", mdot=30, T=300)
    result = model.steady_state()
    print(model.pressure("inlet"))

Key conveniences:

* **Nodes are created on first reference** by name; no explicit `Node` objects.
* **`add_pipe(..., n_cells=N)`** subdivides a pipe into N finite-volume cells (the paper's
  discretization, Fig. 1), inserting the internal nodes and assigning each node its
  control-volume size automatically.
* **Boundaries accept a callable** ``p=lambda t: ...`` for time-varying pressure (e.g. the
  blow-down case), updated automatically during :meth:`run`.
* **`run(dt, duration, record=...)`** solves the steady state as the initial condition, then
  time-marches, returning the recorded time histories.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass, field

from .components import MassFlowBoundary, Pipe, PressureBoundary, Pump, Valve
from .fluids import FluidModel
from .network import Network, Node
from .solver import PCIMSolver, SolverConfig, StepResult


@dataclass
class FlowModel:
    """A fluid-network model: build it with the ``add_*`` / ``*_boundary`` methods, then
    solve with :meth:`steady_state` or :meth:`run`."""

    fluid: FluidModel
    default_temperature: float = 300.0
    network: Network = field(default_factory=Network)
    _time_pressures: dict[str, Callable[[float], float]] = field(default_factory=dict, init=False)
    _solver: PCIMSolver | None = field(default=None, init=False)

    # ---- topology ------------------------------------------------------------------

    def node(self, name: str) -> Node:
        """Return the node ``name``, creating it (at the default temperature) if needed."""
        nd = self.network.nodes.get(name)
        if nd is None:
            nd = self.network.add_node(Node(id=name))
            nd.state.T = self.default_temperature
        return nd

    def add_pipe(self, upstream: str, downstream: str, *, length: float, diameter: float,
                 friction_factor: float = 0.02, n_cells: int = 1,
                 name: str | None = None) -> FlowModel:
        """Add a pipe from ``upstream`` to ``downstream``, subdivided into ``n_cells``
        finite-volume segments (with ``n_cells - 1`` internal nodes)."""
        if n_cells < 1:
            raise ValueError("n_cells must be >= 1")
        name = name or f"{upstream}->{downstream}"
        area = 0.25 * math.pi * diameter * diameter
        dx = length / n_cells
        half_cell = 0.5 * dx * area  # each segment gives half its cell volume to each end node
        prev = self.node(upstream)
        for k in range(n_cells):
            nxt = self.node(downstream) if k == n_cells - 1 else self.node(f"{name}#c{k}")
            prev.volume += half_cell
            nxt.volume += half_cell
            seg_id = name if n_cells == 1 else f"{name}#s{k}"
            self.network.add_element(Pipe(id=seg_id, upstream=prev, downstream=nxt,
                                          fluid=self.fluid, length=dx, diameter=diameter,
                                          friction_factor=friction_factor))
            prev = nxt
        return self

    def add_valve(self, upstream: str, downstream: str, *, k_open: float = 1.0,
                  opening: Callable[[float], float] | None = None,
                  name: str | None = None) -> FlowModel:
        name = name or f"valve_{upstream}->{downstream}"
        valve = Valve(id=name, upstream=self.node(upstream), downstream=self.node(downstream),
                      fluid=self.fluid, k_open=k_open)
        if opening is not None:
            valve.opening = opening
        self.network.add_element(valve)
        return self

    def add_pump(self, upstream: str, downstream: str, *, head_shutoff: float, curve: float,
                 efficiency: float = 1.0, name: str | None = None) -> FlowModel:
        name = name or f"pump_{upstream}->{downstream}"
        self.network.add_element(Pump(id=name, upstream=self.node(upstream),
                                      downstream=self.node(downstream), fluid=self.fluid,
                                      head_shutoff=head_shutoff, curve=curve,
                                      efficiency=efficiency))
        return self

    # ---- boundary conditions -------------------------------------------------------

    def pressure_boundary(self, node: str, *, p: float | Callable[[float], float],
                          T: float | None = None) -> FlowModel:
        """Fix a node's pressure (and temperature). ``p`` may be a callable of time."""
        nd = self.node(node)
        temperature = self.default_temperature if T is None else T
        if callable(p):
            self._time_pressures[nd.id] = p
            PressureBoundary(node=nd, p=float(p(0.0)), T=temperature).apply()
        else:
            PressureBoundary(node=nd, p=p, T=temperature).apply()
        return self

    def mass_flow_boundary(self, node: str, *, mdot: float, T: float | None = None) -> FlowModel:
        MassFlowBoundary(node=self.node(node), mdot=mdot, T=T).apply()
        return self

    def heat_source(self, node: str, *, power: float) -> FlowModel:
        """Impose a heat input rate [W] at a node (for the energy equation)."""
        self.node(node).heat_source = power
        return self

    # ---- solving -------------------------------------------------------------------

    def solver(self, **config: object) -> PCIMSolver:
        self._solver = PCIMSolver(self.network, SolverConfig(**config))  # type: ignore[arg-type]
        return self._solver

    def steady_state(self, **config: object) -> StepResult:
        return self.solver(**config).steady_state()

    def run(self, dt: float, duration: float, *, record: tuple[str, ...] = (),
            steady_init: bool = True, **config: object) -> dict[str, list[float]]:
        """Time-march for ``duration`` in steps of ``dt``.

        If ``steady_init`` the steady state is solved first as the initial condition.
        ``record`` is a tuple of ``"<kind>:<id>"`` keys (``p``/``T`` for a node, ``mdot`` for
        an element, ``flow`` for a whole subdivided pipe); the returned dict maps each key
        (plus ``"t"``) to its time history.
        """
        solver = self.solver(**config)
        if steady_init:
            solver.steady_state()
        history: dict[str, list[float]] = {"t": []}
        for key in record:
            history[key] = []
        t = 0.0
        for _ in range(int(round(duration / dt))):
            t += dt
            for name, fn in self._time_pressures.items():
                self.network.nodes[name].state.p0 = float(fn(t))
            solver.step(dt, t)
            history["t"].append(t)
            for key in record:
                history[key].append(self._query(key))
        return history

    # ---- results -------------------------------------------------------------------

    def pressure(self, node: str) -> float:
        return self.network.nodes[node].state.p0

    def temperature(self, node: str) -> float:
        return self.network.nodes[node].state.T

    def flow(self, element: str) -> float:
        return self.network.elements[element].mdot

    def flow_through(self, pipe_name: str) -> float:
        """Mean mass flow through a (possibly subdivided) pipe added under ``pipe_name``."""
        segs = [e for eid, e in self.network.elements.items()
                if eid == pipe_name or eid.startswith(f"{pipe_name}#s")]
        if not segs:
            raise KeyError(f"no pipe named {pipe_name!r}")
        return sum(e.mdot for e in segs) / len(segs)

    def _query(self, key: str) -> float:
        kind, _, ident = key.partition(":")
        if kind == "p":
            return self.pressure(ident)
        if kind == "T":
            return self.temperature(ident)
        if kind in ("m", "mdot"):
            return self.flow(ident)
        if kind == "flow":
            return self.flow_through(ident)
        raise ValueError(f"unknown record key {key!r} (use p:/T:/mdot:/flow:)")
