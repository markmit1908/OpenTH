"""Minimal command-line entry point.

Wired up via the ``flowcalc`` console script (see pyproject.toml). For now it just reports
the version and is the place to add ``run <network.json>`` once the solver lands.
"""

from __future__ import annotations

import argparse

from . import __version__


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="flowcalc", description=__doc__)
    parser.add_argument("--version", action="version", version=f"flowcalc {__version__}")
    sub = parser.add_subparsers(dest="command")
    sub.add_parser("run", help="(planned) load a network description and solve it")

    args = parser.parse_args(argv)
    if args.command is None:
        parser.print_help()
        return 0
    if args.command == "run":
        raise SystemExit("`flowcalc run` is not implemented yet; see flowcalc.solver.PCIMSolver")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
