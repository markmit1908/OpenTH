"""Linear-system back ends for the pressure-correction and energy equations.

The discretized equations (24) and (29) form, per node i:

    cP * x_i = cE * x_{i+1} + cW * x_{i-1} + b_i

For a single series pipeline this is tridiagonal and is solved directly with the Thomas
algorithm (paper, Section 4). For branching networks the same coefficients assemble into
a general sparse matrix, solved with scipy. These are exactly the kernels earmarked for a
later C/C++ reimplementation (see ``native/``).
"""

from __future__ import annotations

import numpy as np


def thomas(a: np.ndarray, b: np.ndarray, c: np.ndarray, d: np.ndarray) -> np.ndarray:
    """Solve a tridiagonal system by the Thomas algorithm.

    Parameters
    ----------
    a : sub-diagonal (a[0] unused), b : diagonal, c : super-diagonal (c[-1] unused),
    d : right-hand side. All length n. Returns the solution vector.
    """
    n = len(b)
    cp = np.empty(n)
    dp = np.empty(n)
    cp[0] = c[0] / b[0]
    dp[0] = d[0] / b[0]
    for i in range(1, n):
        m = b[i] - a[i] * cp[i - 1]
        cp[i] = c[i] / m
        dp[i] = (d[i] - a[i] * dp[i - 1]) / m
    x = np.empty(n)
    x[-1] = dp[-1]
    for i in range(n - 2, -1, -1):
        x[i] = dp[i] - cp[i] * x[i + 1]
    return x


def sparse_solve(matrix, rhs: np.ndarray) -> np.ndarray:
    """Solve a general sparse linear system (scipy CSR/CSC) for network topologies."""
    from scipy.sparse.linalg import spsolve  # local import keeps scipy optional at import time

    return spsolve(matrix.tocsr(), rhs)
