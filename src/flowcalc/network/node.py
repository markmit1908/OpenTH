"""Node = a finite-volume control volume (cell centre).

In Greyvenstein's scheme (Fig. 1), the *cell-centred* quantities are pressure, density
and temperature; volumetric flow rate and velocity live on the *faces* (see
:mod:`flowcalc.network.element`). A :class:`Node` therefore carries the scalar state and
the control-volume geometry, and is where the **continuity** (eq. 13) and **energy**
(eq. 15) equations are integrated.
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class NodeState:
    """Primary unknowns stored at a node, at the *current* time/iteration level.

    The solver also keeps the previous-time-level values (superscript ``o`` in the
    paper); see :class:`flowcalc.solver.base.SolverState`.
    """

    p0: float = 0.0      # total pressure [Pa]      (primary continuity unknown)
    T: float = 0.0       # temperature [K]
    rho: float = 0.0     # density [kg/m^3]          (from equation of state)
    h0: float = 0.0      # total enthalpy [J/kg]     (primary energy unknown)


@dataclass
class Node:
    """A control volume in the network.

    Parameters
    ----------
    id : str
        Unique label.
    volume : float
        Control-volume size V [m^3] (the ``V`` in eqs. 13, 15). Boundary/junction
        nodes may have zero volume (no storage).
    elevation : float
        Height z [m] for the gravitational term rho g z.
    """

    id: str
    volume: float = 0.0
    elevation: float = 0.0
    state: NodeState = field(default_factory=NodeState)

    # Set True for nodes whose *pressure* is imposed (Dirichlet); the solver excludes
    # these rows from the pressure-correction system. A node with an imposed mass flow
    # (see MassFlowBoundary) stays an unknown and instead carries a non-zero mass_source.
    is_boundary: bool = False

    # Imposed mass flow into the node [kg/s] (+ = inflow, - = outflow); a source term in
    # the continuity balance. Zero for ordinary interior nodes.
    mass_source: float = 0.0

    # Imposed heat input rate into the control volume [W] for the energy equation
    # (the V*q_dot term of eq. 15). Zero = adiabatic.
    heat_source: float = 0.0

    # If True, the temperature (total enthalpy) is held fixed in the energy solve
    # (Dirichlet), e.g. at a flow inlet. Pressure-fixed boundary nodes are treated as
    # temperature-fixed too. Set automatically by the boundary helpers.
    fixed_temperature: bool = False
