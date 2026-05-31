"""Solvers: the PCIM pressure-correction method and its linear-algebra back ends."""

from .base import Solver, SolverConfig, StepResult
from .linear import sparse_solve, thomas
from .pcim import PCIMSolver

__all__ = ["PCIMSolver", "Solver", "SolverConfig", "StepResult", "sparse_solve", "thomas"]
