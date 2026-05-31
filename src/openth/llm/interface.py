"""Two-way interface between OpenTH and an LLM.

Goal (planned): let a model both *drive* the solver (build a network, run a case, read
results back) and be *driven by* it (have the solver hand structured state to a model for
explanation, design suggestions, or anomaly detection). The two directions:

  * inbound  -- LLM -> OpenTH: a natural-language or structured request is turned into a
    network description (see :mod:`openth.io.serialize`) and a simulation request.
  * outbound -- OpenTH -> LLM: solver state/results are summarised into a payload the
    model can reason over (tool results, design feedback).

Kept dependency-light: the core solver must never import this module. The ``anthropic``
SDK lives behind the optional ``[llm]`` extra. When wiring this up, prefer the latest
Claude model and enable prompt caching for the (large, static) tool/schema preamble.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from ..network import Network


@dataclass
class ToolSpec:
    """Description of a OpenTH capability exposed to the model as a tool."""

    name: str
    description: str
    input_schema: dict[str, Any]


def network_tools() -> list[ToolSpec]:
    """Tool specs the model can call to build and query networks.

    TODO: define build_network / add_element / run_steady_state / run_transient /
    get_results tools backed by :mod:`openth.io` and :mod:`openth.solver`.
    """
    raise NotImplementedError("network_tools not yet defined")


def summarize(network: Network) -> dict[str, Any]:
    """Produce an LLM-friendly summary of a network and (eventually) its results."""
    return {
        "n_nodes": len(network.nodes),
        "n_elements": len(network.elements),
        # TODO: include solved state (pressures, flows, temperatures) once the solver runs.
    }
