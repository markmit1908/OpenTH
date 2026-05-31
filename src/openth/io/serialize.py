"""JSON (de)serialization for :class:`~openth.model.FlowModel`.

A model is saved as its **high-level build directives** (the ``add_pipe`` / ``add_valve`` /
boundary / ... calls, with their arguments) plus the fluid — not the expanded finite-volume
network — so the JSON is compact and hand-editable. Loading replays the directives onto a
fresh model, reproducing it exactly. This declarative payload is also what the planned
two-way LLM interface (:mod:`openth.llm`) exchanges.

Schema (v1)::

    {
      "schema": 1,
      "fluid": {"kind": "ideal_gas", "name": "helium", "R": 2077.0, "gamma": 1.667, "s": 1.0},
      "default_temperature": 300.0,
      "build": [
        {"op": "add_pipe", "upstream": "inlet", "downstream": "outlet",
         "length": 100.0, "diameter": 0.5, "friction_factor": 0.02,
         "n_cells": 20, "delta_elevation": 0.0, "name": "pipe"},
        {"op": "pressure_boundary", "node": "outlet", "p": 200000.0, "T": 300.0},
        {"op": "mass_flow_boundary", "node": "inlet", "mdot": 30.0, "T": 300.0}
      ]
    }

``fluid.kind`` may also be a built-in shorthand — ``"helium"``, ``"air"``, ``"water"`` — for
hand-authoring. Callable arguments (time-varying boundaries, valve schedules) cannot be
written to JSON; serializing such a model raises ``ValueError``.
"""

from __future__ import annotations

import json
from typing import TYPE_CHECKING, Any

from ..fluids import FluidModel, IdealGas, Incompressible, air, helium, water

if TYPE_CHECKING:
    from ..model import FlowModel

SCHEMA = 1
_BUILTIN_FLUIDS = {"helium": helium, "air": air, "water": water}


# --- fluids ----------------------------------------------------------------------------

def fluid_to_dict(fluid: FluidModel) -> dict[str, Any]:
    if isinstance(fluid, IdealGas):
        return {"kind": "ideal_gas", "name": fluid.name, "R": fluid.R,
                "gamma": fluid.gamma, "s": fluid.s}
    if isinstance(fluid, Incompressible):
        return {"kind": "incompressible", "name": fluid.name, "rho": fluid.rho,
                "cp": fluid.cp, "bulk_modulus": fluid.bulk_modulus}
    raise ValueError(f"cannot serialize fluid of type {type(fluid).__name__}")


def fluid_from_dict(data: dict[str, Any]) -> FluidModel:
    kind = data["kind"]
    if kind in _BUILTIN_FLUIDS:
        return _BUILTIN_FLUIDS[kind]()
    if kind == "ideal_gas":
        return IdealGas(name=data.get("name", "gas"), R=data["R"], gamma=data["gamma"],
                        s=data.get("s", 1.0))
    if kind == "incompressible":
        return Incompressible(name=data.get("name", "liquid"), rho=data["rho"], cp=data["cp"],
                              bulk_modulus=data.get("bulk_modulus", 2.2e9))
    raise ValueError(f"unknown fluid kind {kind!r}")


# --- model -----------------------------------------------------------------------------

def model_to_dict(model: FlowModel) -> dict[str, Any]:
    """Serialize a :class:`~openth.model.FlowModel` to a JSON-ready dict."""
    build: list[dict[str, Any]] = []
    for op, kwargs in model._directives:
        for key, value in kwargs.items():
            if callable(value):
                raise ValueError(
                    f"cannot serialize {op!r}: argument {key!r} is a callable "
                    "(e.g. a time-varying boundary or valve schedule); remove it to save."
                )
        build.append({"op": op, **kwargs})
    return {
        "schema": SCHEMA,
        "fluid": fluid_to_dict(model.fluid),
        "default_temperature": model.default_temperature,
        "build": build,
    }


def model_from_dict(data: dict[str, Any]) -> FlowModel:
    """Rebuild a :class:`~openth.model.FlowModel` from :func:`model_to_dict` output."""
    from ..model import FlowModel  # local import avoids a cycle (model imports io lazily)

    schema = data.get("schema", SCHEMA)
    if schema != SCHEMA:
        raise ValueError(f"unsupported model schema {schema} (this build reads {SCHEMA})")
    model = FlowModel(fluid=fluid_from_dict(data["fluid"]),
                      default_temperature=data.get("default_temperature", 300.0))
    for raw in data["build"]:
        step = dict(raw)
        op = step.pop("op")
        getattr(model, op)(**step)
    return model


def save_model(model: FlowModel, path: str) -> None:
    """Write ``model`` to ``path`` as JSON."""
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(model_to_dict(model), fh, indent=2)
        fh.write("\n")


def load_model(path: str) -> FlowModel:
    """Load a :class:`~openth.model.FlowModel` from a JSON file."""
    with open(path, encoding="utf-8") as fh:
        return model_from_dict(json.load(fh))
