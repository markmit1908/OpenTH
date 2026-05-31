"""Port / connection modelling layer — the vision's component-and-connection style.

Instead of naming nodes by hand (the :class:`~openth.model.FlowModel` way), you build a
circuit from **component objects with `.inlet` / `.outlet` ports** and `connect` them::

    import openth as th

    c = th.Circuit(fluid=th.Fluid("helium"))
    pipe  = c.add(th.Pipe(length=50, diameter=0.5, n_cells=10))
    pump  = c.add(th.Pump(head_shutoff=80e3, curve=200.0))
    c.connect(pump.outlet, pipe.inlet)
    c.connect(pipe.outlet, pump.inlet)        # a closed loop
    c.pressure_boundary(pump.inlet, p=200e3, T=300)   # reference pressure
    c.solve_steady_state()
    print(c.flow(pipe))

A `Circuit` *compiles* to a `FlowModel`: each group of connected ports becomes one network
node, and each component emits its element(s) between its two port-nodes (pipes still
subdivide via ``n_cells``). All the solver/physics is the validated `FlowModel` machinery.

Single working fluid per circuit (multi-fluid coupling, e.g. a two-sided heat exchanger, is
future work). A network needs at least one pressure boundary — for a closed loop, pin one
node as the pressure reference.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from typing import cast

from .fluids import FluidModel
from .model import FlowModel
from .solver import StepResult


class Port:
    """A connection point on a component (its ``inlet`` or ``outlet``)."""

    def __init__(self, component: Component, role: str) -> None:
        self.component = component
        self.role = role

    def __repr__(self) -> str:
        return f"Port({self.component.kind}.{self.role})"


class Component:
    """Base class for port-based components. Concrete components are dataclasses that add
    their physical parameters and implement :meth:`_emit`."""

    inlet: Port
    outlet: Port
    name: str | None

    def _init_ports(self) -> None:
        self.inlet = Port(self, "inlet")
        self.outlet = Port(self, "outlet")

    @property
    def kind(self) -> str:
        return type(self).__name__

    def _emit(self, model: FlowModel, inlet: str, outlet: str, label: str) -> None:
        """Add this component's element(s) to ``model`` between the two resolved nodes."""
        raise NotImplementedError


@dataclass
class Pipe(Component):
    """A pipe (subdivided into ``n_cells`` finite-volume cells)."""

    length: float
    diameter: float
    fluid: FluidModel | None = None
    friction_factor: float = 0.02
    n_cells: int = 1
    name: str | None = None

    def __post_init__(self) -> None:
        self._init_ports()

    def _emit(self, model: FlowModel, inlet: str, outlet: str, label: str) -> None:
        model.add_pipe(inlet, outlet, length=self.length, diameter=self.diameter,
                       friction_factor=self.friction_factor, n_cells=self.n_cells, name=label)


@dataclass
class Valve(Component):
    """A valve / orifice (variable resistance; ``opening(t)`` can close it in time)."""

    k_open: float = 1.0
    opening: Callable[[float], float] | None = None
    fluid: FluidModel | None = None
    name: str | None = None

    def __post_init__(self) -> None:
        self._init_ports()

    def _emit(self, model: FlowModel, inlet: str, outlet: str, label: str) -> None:
        model.add_valve(inlet, outlet, k_open=self.k_open, opening=self.opening, name=label)


@dataclass
class Pump(Component):
    """A pump / compressor with a quadratic head curve ``head_shutoff - curve*mdot|mdot|``."""

    head_shutoff: float
    curve: float
    efficiency: float = 1.0
    fluid: FluidModel | None = None
    name: str | None = None

    def __post_init__(self) -> None:
        self._init_ports()

    def _emit(self, model: FlowModel, inlet: str, outlet: str, label: str) -> None:
        model.add_pump(inlet, outlet, head_shutoff=self.head_shutoff, curve=self.curve,
                       efficiency=self.efficiency, name=label)


class Circuit:
    """A component-and-connection model. Add components, connect their ports, set
    boundaries, then solve. Compiles to a :class:`~openth.model.FlowModel`."""

    def __init__(self, fluid: FluidModel | None = None, default_temperature: float = 300.0):
        self.fluid = fluid
        self.default_temperature = default_temperature
        self._components: list[Component] = []
        self._connections: list[tuple[Port, Port]] = []
        self._pressure_bcs: list[tuple[Port, float | Callable[[float], float], float | None]] = []
        self._mass_bcs: list[tuple[Port, float, float | None]] = []
        self._model: FlowModel | None = None
        self._port_node: dict[int, str] = {}
        self._comp_label: dict[int, str] = {}

    # ---- building ------------------------------------------------------------------

    def add(self, component: Component) -> Component:
        """Register a component and return it (so you can reference its ports)."""
        self._components.append(component)
        self._model = None
        return component

    def connect(self, port_a: Port, port_b: Port) -> Circuit:
        """Join two ports — they resolve to the same network node."""
        self._connections.append((port_a, port_b))
        self._model = None
        return self

    def pressure_boundary(self, port: Port, *, p: float | Callable[[float], float],
                          T: float | None = None) -> Circuit:
        self._pressure_bcs.append((port, p, T))
        self._model = None
        return self

    def mass_flow_boundary(self, port: Port, *, mdot: float, T: float | None = None) -> Circuit:
        self._mass_bcs.append((port, mdot, T))
        self._model = None
        return self

    # ---- compilation ---------------------------------------------------------------

    def compile(self) -> FlowModel:
        """Resolve connected ports into nodes and build the underlying ``FlowModel``."""
        fluid = self._resolve_fluid()
        node_for_port = self._resolve_nodes()
        self._port_node = node_for_port

        model = FlowModel(fluid=fluid, default_temperature=self.default_temperature)
        self._comp_label = {}
        for i, comp in enumerate(self._components):
            label = comp.name or f"{comp.kind.lower()}{i}"
            self._comp_label[id(comp)] = label
            comp._emit(model, node_for_port[id(comp.inlet)],
                       node_for_port[id(comp.outlet)], label)

        for port, p, temperature in self._pressure_bcs:
            model.pressure_boundary(node_for_port[id(port)], p=p, T=temperature)
        for port, mdot, temperature in self._mass_bcs:
            model.mass_flow_boundary(node_for_port[id(port)], mdot=mdot, T=temperature)

        self._model = model
        return model

    def _resolve_fluid(self) -> FluidModel:
        fluids = [self.fluid] if self.fluid is not None else []
        fluids += [f for c in self._components if (f := getattr(c, "fluid", None)) is not None]
        distinct = {id(f): f for f in fluids}
        if not distinct:
            raise ValueError(
                "circuit has no fluid; pass Circuit(fluid=...) or set it on a component")
        if len(distinct) > 1:
            raise NotImplementedError("multi-fluid circuits are not supported yet")
        return next(iter(distinct.values()))

    def _resolve_nodes(self) -> dict[int, str]:
        # Union-find over ports: connected ports share a node.
        parent: dict[int, int] = {}

        def find(x: int) -> int:
            parent.setdefault(x, x)
            root = x
            while parent[root] != root:
                root = parent[root]
            while parent[x] != root:
                parent[x], x = root, parent[x]
            return root

        def union(a: int, b: int) -> None:
            parent[find(a)] = find(b)

        for comp in self._components:
            find(id(comp.inlet))
            find(id(comp.outlet))
        for a, b in self._connections:
            union(id(a), id(b))

        group_node: dict[int, str] = {}
        node_for_port: dict[int, str] = {}
        counter = 0
        for comp in self._components:
            for port in (comp.inlet, comp.outlet):
                root = find(id(port))
                if root not in group_node:
                    group_node[root] = f"n{counter}"
                    counter += 1
                node_for_port[id(port)] = group_node[root]
        return node_for_port

    def _ensure_compiled(self) -> FlowModel:
        if self._model is None:
            self.compile()
        assert self._model is not None
        return self._model

    # ---- solving -------------------------------------------------------------------

    def solve_steady_state(self, **config: object) -> StepResult:
        return self._ensure_compiled().steady_state(**config)

    def transient(self, dt: float, t_end: float, *,
                  record: tuple[tuple[str, object], ...] = (),
                  steady_init: bool = True, **config: object) -> dict[str, list[float]]:
        """Time-march to ``t_end`` in steps ``dt``. ``record`` is a tuple of
        ``(kind, target)`` pairs: kind ``"p"``/``"T"`` with a Port, or ``"flow"`` with a
        Component. Returns a dict keyed by readable labels."""
        model = self._ensure_compiled()
        model_keys: list[str] = []
        labels: list[str] = []
        for kind, target in record:
            model_keys.append(self._record_key(kind, target))
            labels.append(self._record_label(kind, target))
        history = model.run(dt, t_end, record=tuple(model_keys),
                            steady_init=steady_init, **config)
        out: dict[str, list[float]] = {"t": history["t"]}
        for mk, label in zip(model_keys, labels, strict=True):
            out[label] = history[mk]
        return out

    # ---- results -------------------------------------------------------------------

    def pressure(self, port: Port) -> float:
        return self._ensure_compiled().pressure(self._port_node[id(port)])

    def temperature(self, port: Port) -> float:
        return self._ensure_compiled().temperature(self._port_node[id(port)])

    def flow(self, component: Component) -> float:
        """Mean mass flow through a component."""
        self._ensure_compiled()
        return self._model.flow_through(self._comp_label[id(component)])  # type: ignore[union-attr]

    @property
    def model(self) -> FlowModel:
        """The compiled underlying ``FlowModel`` (for advanced/low-level access)."""
        return self._ensure_compiled()

    # ---- record helpers ------------------------------------------------------------

    def _record_key(self, kind: str, target: object) -> str:
        if kind in ("p", "T"):
            return f"{kind}:{self._port_node[id(target)]}"
        if kind in ("flow", "mdot"):
            return f"flow:{self._comp_label[id(target)]}"
        raise ValueError(f"unknown record kind {kind!r} (use 'p', 'T', or 'flow')")

    def _record_label(self, kind: str, target: object) -> str:
        if kind in ("p", "T"):
            port = cast("Port", target)
            return f"{kind}:{self._comp_label[id(port.component)]}.{port.role}"
        return f"flow:{self._comp_label[id(target)]}"
