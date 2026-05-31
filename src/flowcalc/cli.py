"""Command-line entry point (the ``flowcalc`` console script; see pyproject.toml).

Currently supports ``flowcalc benchmark`` to generate and run the paper's Section 5 test
cases (see :mod:`flowcalc.benchmarks`).
"""

from __future__ import annotations

import argparse

from . import __version__, benchmarks


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="flowcalc", description=__doc__)
    parser.add_argument("--version", action="version", version=f"flowcalc {__version__}")
    sub = parser.add_subparsers(dest="command")

    bench = sub.add_parser("benchmark", help="run a paper Section 5 test case")
    bench.add_argument("name", nargs="?",
                       help="benchmark to run; omit to list the available ones")

    args = parser.parse_args(argv)
    if args.command is None:
        parser.print_help()
        return 0
    if args.command == "benchmark":
        return _run_benchmark(args.name)
    return 0


def _run_benchmark(name: str | None) -> int:
    if name is None:
        print("Available benchmarks (from Greyvenstein 2002, Section 5):\n")
        for key, (description, _) in benchmarks.BENCHMARKS.items():
            print(f"  {key:20} {description}")
        print("\nRun one with:  flowcalc benchmark <name>")
        return 0
    if name not in benchmarks.BENCHMARKS:
        print(f"unknown benchmark {name!r}; choose from: {', '.join(benchmarks.BENCHMARKS)}")
        return 1
    description, _ = benchmarks.BENCHMARKS[name]
    print(f"{name}: {description}\n")
    for key, value in benchmarks.run(name).items():
        print(f"  {key}: {value}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
