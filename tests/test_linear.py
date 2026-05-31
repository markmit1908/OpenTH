"""Tests for the linear-solve back ends (Thomas algorithm)."""

import numpy as np

from flowcalc.solver import thomas


def test_thomas_matches_dense_solve():
    # Build a random diagonally-dominant tridiagonal system and compare to numpy.
    rng = np.random.default_rng(0)
    n = 6
    a = np.concatenate([[0.0], rng.uniform(-1, 1, n - 1)])      # sub-diagonal
    c = np.concatenate([rng.uniform(-1, 1, n - 1), [0.0]])      # super-diagonal
    b = np.abs(a) + np.abs(c) + rng.uniform(1, 2, n)            # diagonally dominant
    x_true = rng.uniform(-5, 5, n)

    M = np.diag(b) + np.diag(a[1:], -1) + np.diag(c[:-1], 1)
    d = M @ x_true

    x = thomas(a, b, c, d)
    np.testing.assert_allclose(x, x_true, rtol=1e-10, atol=1e-10)
