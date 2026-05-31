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

``th.Model`` is the high-level builder (alias of :class:`~openth.model.FlowModel`); add
elements with ``add_pipe`` / ``add_valve`` / ``add_pump`` and set ``pressure_boundary`` /
``mass_flow_boundary``. A port/connection-style API (``model.connect(a.outlet, b.inlet)``,
as sketched in ``docs/vision-statement.md``) is future work — use the name-based ``add_*``
methods for now. The optional LLM interface lives in ``openth.llm`` and is **not** imported
here (it needs the ``[llm]`` extra).
"""

from .components import MassFlowBoundary, Pipe, PressureBoundary, Pump, Valve
from .fluids import Fluid, FluidModel, IdealGas, Incompressible, air, helium, water
from .model import FlowModel
from .network import Element, Network, Node
from .solver import PCIMSolver, SolverConfig

__version__ = "0.0.1"

#: Vision-style alias: ``th.Model(...)`` is the high-level :class:`~openth.model.FlowModel`.
Model = FlowModel

__all__ = [
    # high-level
    "Model", "FlowModel",
    # fluids
    "Fluid", "FluidModel", "IdealGas", "Incompressible", "helium", "air", "water",
    # components
    "Pipe", "Valve", "Pump", "PressureBoundary", "MassFlowBoundary",
    # low-level topology / solver
    "Network", "Node", "Element", "PCIMSolver", "SolverConfig",
    "__version__",
]
