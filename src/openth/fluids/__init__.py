"""Fluid property models."""

from .base import FluidModel
from .ideal_gas import IdealGas, air, helium
from .incompressible import Incompressible, water

_BUILTIN = {"helium": helium, "he": helium, "air": air, "water": water}


def Fluid(name: str) -> FluidModel:
    """Return a built-in fluid by name, e.g. ``Fluid("helium")``.

    Known names: ``helium`` (``he``), ``air``, ``water``. For anything else, construct an
    :class:`IdealGas` or :class:`Incompressible` directly with its properties.
    """
    try:
        return _BUILTIN[name.strip().lower()]()
    except KeyError:
        raise ValueError(
            f"unknown fluid {name!r}; built-ins are {sorted(set(_BUILTIN))} "
            "(or build an IdealGas / Incompressible directly)"
        ) from None


__all__ = ["Fluid", "FluidModel", "IdealGas", "Incompressible", "air", "helium", "water"]
