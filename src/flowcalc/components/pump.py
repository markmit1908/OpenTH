"""Pump / compressor: a non-pipe component that adds pressure (and, for a compressor,
energy) to the flow.

The device is modelled by a quadratic characteristic curve relating the pressure rise to
the mass flow:

    Δp_rise(mdot) = head_shutoff - curve * mdot * |mdot|

i.e. maximum head at zero flow, falling as the flow increases. This maps cleanly onto the
solver's existing momentum closure: the constant ``head_shutoff`` is reported via
:meth:`head` (added to the momentum driving pressure, so the pump can push flow against an
adverse pressure gradient), and the flow-dependent slope ``curve`` is the friction
``resistance`` (so the pressure-correction sensitivity ``1/(2 K |mdot|)`` is unchanged).

For a **compressor**, the shaft work also raises the gas total enthalpy; :meth:`work_per_mass`
returns the specific work (≈ Δp_rise / (ρ·efficiency)) which the energy equation adds to the
downstream node, producing the temperature rise that distinguishes a compressor from an
isothermal pump.
"""

from __future__ import annotations

from dataclasses import dataclass

from ..network.element import Element

# Floor on the curve slope used as the momentum resistance: a pump with no flow-dependence
# (curve = 0) would give a singular momentum solve (division by zero), so keep it positive.
_MIN_CURVE = 1e-9


@dataclass
class Pump(Element):
    """Pump or compressor with a quadratic head curve.

    Parameters
    ----------
    head_shutoff : float
        Pressure rise at zero flow [Pa] (the top of the pump curve).
    curve : float
        Quadratic fall-off [Pa/(kg/s)^2]: ``rise = head_shutoff - curve * mdot|mdot|``.
    efficiency : float
        Isentropic-ish efficiency in (0, 1]; only affects the enthalpy rise (a compressor
        with lower efficiency heats the gas more for the same pressure rise).
    """

    head_shutoff: float = 0.0
    curve: float = _MIN_CURVE
    efficiency: float = 1.0

    def resistance(self, rho_face: float) -> float:
        return max(self.curve, _MIN_CURVE)

    def head(self) -> float:
        return self.head_shutoff

    def pressure_rise(self, mdot: float) -> float:
        """Actual pressure rise at the current mass flow (the operating point)."""
        return self.head_shutoff - self.curve * mdot * abs(mdot)

    def work_per_mass(self, mdot: float, rho: float) -> float:
        if mdot <= 0.0 or rho <= 0.0:
            return 0.0  # no work added in the reverse direction
        rise = max(self.pressure_rise(mdot), 0.0)
        return rise / (rho * self.efficiency)
