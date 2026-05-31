"""Command-line entry point (the ``openth`` console script; see pyproject.toml).

  * ``openth benchmark [name]`` — run a paper Section 5 test case (:mod:`openth.benchmarks`).
  * ``openth run <model.json>`` — load a model from JSON and solve its steady state.
"""

from __future__ import annotations

import argparse

from . import __version__, benchmarks


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="openth", description=__doc__)
    parser.add_argument("--version", action="version", version=f"openth {__version__}")
    sub = parser.add_subparsers(dest="command")

    bench = sub.add_parser("benchmark", help="run a paper Section 5 test case")
    bench.add_argument("name", nargs="?",
                       help="benchmark to run; omit to list the available ones")

    run = sub.add_parser("run", help="load a model .json and solve its steady state")
    run.add_argument("path", help="path to a model JSON file (see openth.io)")
    run.add_argument("--energy", action="store_true", help="solve the energy equation")
    run.add_argument("--relaxation", type=float, default=0.5, help="under-relaxation (0,1]")

    args = parser.parse_args(argv)
    if args.command is None:
        parser.print_help()
        return 0
    if args.command == "benchmark":
        return _run_benchmark(args.name)
    if args.command == "run":
        return _run_model(args.path, energy=args.energy, relaxation=args.relaxation)
    return 0


def _run_benchmark(name: str | None) -> int:
    if name is None:
        print("Available benchmarks (from Greyvenstein 2002, Section 5):\n")
        for key, (description, _) in benchmarks.BENCHMARKS.items():
            print(f"  {key:20} {description}")
        print("\nRun one with:  openth benchmark <name>")
        return 0
    if name not in benchmarks.BENCHMARKS:
        print(f"unknown benchmark {name!r}; choose from: {', '.join(benchmarks.BENCHMARKS)}")
        return 1
    description, _ = benchmarks.BENCHMARKS[name]
    print(f"{name}: {description}\n")
    for key, value in benchmarks.run(name).items():
        print(f"  {key}: {value}")
    return 0


def _run_model(path: str, *, energy: bool, relaxation: float) -> int:
    from .model import FlowModel

    model = FlowModel.load(path)
    result = model.steady_state(relaxation=relaxation, solve_energy=energy,
                                max_outer_iterations=2000)
    print(f"{path}: converged={result.converged} in {result.iterations} iterations\n")
    print(f"{'node':>16} {'p [kPa]':>12}" + ("  T [K]" if energy else ""))
    for name, node in model.network.nodes.items():
        line = f"{name:>16} {node.state.p0 / 1e3:12.2f}"
        if energy:
            line += f" {node.state.T:7.1f}"
        print(line)
    print(f"\n{'element':>16} {'mdot [kg/s]':>12}")
    for name, element in model.network.elements.items():
        print(f"{name:>16} {element.mdot:12.4f}")
    return 0 if result.converged else 1


if __name__ == "__main__":
    raise SystemExit(main())
