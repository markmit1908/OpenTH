"""Solver interfaces and configuration shared by all solver back ends."""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass

from ..network import Network


@dataclass
class SolverConfig:
    """Tunables for the segregated pressure-correction solve.

    Attributes
    ----------
    alpha : float
        Time-integration weighing factor in [0.5, 1] (paper Section 3). alpha=1 is fully
        implicit (1st-order), alpha=0.5 is Crank-Nicolson (2nd-order but prone to
        instability); alpha=0.6 is the recommended accuracy/stability compromise.
        (Unused by the steady-state solve.)
    relaxation : float
        Pressure under-relaxation factor in (0, 1] for the segregated (SIMPLE) iteration;
        smaller is more stable but slower. 0.7 is a typical default.
    max_outer_iterations : int
        Maximum segregated iterations (per time step, or for the steady solve).
    max_pressure_iterations : int
        Inner pressure-correction iterations per outer iteration (transient).
    tol : float
        Convergence tolerance on the (relative) mass-imbalance residual.
    """

    alpha: float = 0.6
    relaxation: float = 0.7
    max_outer_iterations: int = 200
    max_pressure_iterations: int = 20
    tol: float = 1e-8
    # Solve the energy equation for temperature (non-isothermal). Default off: the flow is
    # isothermal and temperatures stay at their initial/boundary values.
    solve_energy: bool = False
    # Gravitational acceleration [m/s^2] for the hydrostatic/buoyancy term. Acts through the
    # node elevations, so it is inert when all elevations are zero (flat networks). Set to 0
    # to disable gravity entirely.
    gravity: float = 9.80665

    def __post_init__(self) -> None:
        if not 0.5 <= self.alpha <= 1.0:
            raise ValueError("alpha must lie in [0.5, 1.0] for a stable scheme (paper Section 3)")
        if not 0.0 < self.relaxation <= 1.0:
            raise ValueError("relaxation must lie in (0, 1]")


@dataclass
class StepResult:
    """Outcome of advancing one time step (or one steady-state solve)."""

    time: float
    iterations: int
    residual: float
    converged: bool


class Solver(ABC):
    """Base class for flow-network solvers."""

    def __init__(self, network: Network, config: SolverConfig | None = None) -> None:
        self.network = network
        self.config = config or SolverConfig()

    @abstractmethod
    def steady_state(self) -> StepResult:
        """Solve for the steady state (used standalone and as a transient initial condition)."""

    @abstractmethod
    def step(self, dt: float, t: float) -> StepResult:
        """Advance the solution by one time step of size ``dt`` ending at time ``t``."""

    def run(self, dt: float, duration: float) -> list[StepResult]:
        """March from the current state over ``duration`` in steps of ``dt``."""
        results: list[StepResult] = []
        t = 0.0
        while t < duration - 1e-12:
            t += dt
            results.append(self.step(dt, t))
        return results
