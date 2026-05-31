"""Optional two-way LLM interface. Requires the ``[llm]`` extra; never imported by core."""

from .interface import ToolSpec, network_tools, summarize

__all__ = ["ToolSpec", "network_tools", "summarize"]
