"""OpenTH: implicit pressure-correction solver for transient flow in fluid networks.

After G. P. Greyvenstein, "An implicit method for the analysis of transient flows in
pipe networks", Int. J. Numer. Meth. Engng 2002; 53:1127-1143.

The common entry points are re-exported here so you can ``import openth as th``::

    import openth as th

    model = th.Model(fluid=th.Fluid("helium"))
    model.add_pipe("inlet", "outlet", length=100, diameter=0.5, n_cells=20)
    model.pressure_boundary("outlet", p=200e3, T=300)
    model.mass_flow_boundary("inlet", mdot=30, T=300)
    model.steady_state()
    print(model.pressure("inlet"))

``th.Model`` is the name-based builder (alias of :class:`~openth.model.FlowModel`); add
elements with ``add_pipe`` / ``add_valve`` / ``add_pump`` and set ``pressure_boundary`` /
``mass_flow_boundary``.

For the **component-and-connection** style (the vision's ``loop = th.Model()`` sketch), use
:class:`~openth.circuit.Circuit` with port-based components::

    c = th.Circuit(fluid=th.Fluid("helium"))
    pump = c.add(th.Pump(head_shutoff=80e3, curve=200.0))
    pipe = c.add(th.Pipe(length=50, diameter=0.5, n_cells=10))
    c.connect(pump.outlet, pipe.inlet)
    c.connect(pipe.outlet, pump.inlet)
    c.pressure_boundary(pump.inlet, p=200e3, T=300)   # reference pressure
    c.solve_steady_state()

Here ``th.Pipe`` / ``th.Valve`` / ``th.Pump`` are the *port-based component specs* (the
Element-level classes used internally live in ``openth.components``). The optional LLM
interface lives in ``openth.llm`` and is **not** imported here (it needs the ``[llm]`` extra).
"""

from .circuit import Circuit, Component, Pipe, Port, Pump, Valve
from .components import MassFlowBoundary, PressureBoundary
from .fluids import Fluid, FluidModel, IdealGas, Incompressible, air, helium, water
from .model import FlowModel
from .network import Element, HeatExchanger, Network, Node
from .solver import PCIMSolver, SolverConfig

__version__ = "0.0.1"

#: Vision-style alias: ``th.Model(...)`` is the high-level :class:`~openth.model.FlowModel`.
Model = FlowModel

__all__ = [
    # high-level builders
    "Model", "FlowModel", "Circuit",
    # fluids
    "Fluid", "FluidModel", "IdealGas", "Incompressible", "helium", "air", "water",
    # port-based components (for Circuit) + connection primitives
    "Pipe", "Valve", "Pump", "Component", "Port", "HeatExchanger",
    # boundary helpers (low-level) + topology / solver
    "PressureBoundary", "MassFlowBoundary",
    "Network", "Node", "Element", "PCIMSolver", "SolverConfig",
    "__version__",
]
