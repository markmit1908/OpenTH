"""(De)serialization of networks to/from plain dicts (JSON/YAML-friendly).

A declarative network description is also the natural payload exchanged with the LLM
interface (see :mod:`openth.llm`): the model emits a dict, this module turns it into a
:class:`~openth.network.Network`, and results serialize back the same way.

Schema (v0, subject to change):

    {
      "fluid": {"kind": "ideal_gas", "name": "helium", "R": 2077.0, "gamma": 1.667},
      "nodes": [{"id": "n1", "volume": 0.0, "elevation": 0.0}, ...],
      "elements": [
        {"id": "p1", "kind": "pipe", "upstream": "n1", "downstream": "n2",
         "length": 5.0, "diameter": 0.5, "friction_factor": 0.02}, ...
      ],
      "boundaries": [{"kind": "pressure", "node": "n1", "p0": 700e3, "T": 300.0}, ...]
    }
"""

from __future__ import annotations

from typing import Any

from ..network import Network, Node

# TODO: implement the element/fluid/boundary registries and round-trip (to_dict/from_dict).


def from_dict(data: dict[str, Any]) -> Network:
    """Build a Network from a declarative description (see module docstring)."""
    net = Network()
    for nd in data.get("nodes", []):
        net.add_node(
            Node(
                id=nd["id"],
                volume=nd.get("volume", 0.0),
                elevation=nd.get("elevation", 0.0),
            )
        )
    # TODO: instantiate elements (pipe/valve/...) and apply boundaries via a kind registry.
    return net


def to_dict(network: Network) -> dict[str, Any]:
    """Serialize a Network back to a declarative description."""
    raise NotImplementedError("to_dict not yet implemented")
