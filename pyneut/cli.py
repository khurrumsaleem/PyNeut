"""The `pyneut` command.

    pyneut check model.toml      validate the input, data, source and geometry
    pyneut run   model.toml      solve and write <model>_results/results.json

The input format is described in docs/input-format.md.
"""
import argparse
import logging
import os
import sys

from . import __version__
from .errors import InputError


def _load(args):
    from .model import Model
    model = Model.from_toml(args.input)
    if getattr(args, "data", None):
        model.run.data = os.path.abspath(args.data)
    return model


def _print_problems(title, items, stream):
    if items:
        print(f"{title}:", file=stream)
        for item in items:
            print(f"  - {item}", file=stream)


def cmd_check(args):
    try:
        model = _load(args)
    except InputError as e:
        _print_problems(f"{args.input}: {len(e.problems)} problem(s)", e.problems, sys.stderr)
        return 1
    errors, warnings = model.check(samples=args.samples)
    print(model.describe())
    _print_problems("warnings", warnings, sys.stdout)
    if errors:
        _print_problems(f"{len(errors)} error(s)", errors, sys.stderr)
        return 1
    print("OK")
    return 0


def cmd_run(args):
    try:
        model = _load(args)
        for key in ("engine", "particles", "seed", "workers"):
            value = getattr(args, key)
            if value is not None:
                setattr(model.run, key, value)
        model.build()
    except InputError as e:
        _print_problems(f"{args.input}: {len(e.problems)} problem(s)", e.problems, sys.stderr)
        return 1
    if model.run.engine == "vector" and any(t.mesh or t.type == "detector"
                                            for t in model.tallies.values()):
        print("error: mesh and detector tallies need the scalar engine", file=sys.stderr)
        return 1
    results = model.solve()
    print(results.summary())
    out = args.output or os.path.splitext(os.path.basename(args.input))[0] + "_results"
    path = results.save(out)
    print(f"results written to {path}")
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="pyneut",
        description="Continuous-energy Monte Carlo neutron transport.")
    parser.add_argument("--version", action="version", version=f"pyneut {__version__}")
    parser.add_argument("-v", "--verbose", action="store_true", help="debug logging")
    sub = parser.add_subparsers(dest="command")

    p = sub.add_parser("check", help="validate a model without running it")
    p.add_argument("input", help="TOML model file")
    p.add_argument("--data", help="nuclear-data directory (overrides run.data)")
    p.add_argument("--samples", type=int, default=20_000,
                   help="random points for the overlap check (default 20000)")
    p.set_defaults(func=cmd_check)

    p = sub.add_parser("run", help="solve a model and write its results")
    p.add_argument("input", help="TOML model file")
    p.add_argument("-o", "--output", help="results folder (default <model>_results)")
    p.add_argument("--data", help="nuclear-data directory (overrides run.data)")
    p.add_argument("--engine", choices=("vector", "scalar"))
    p.add_argument("--particles", type=int)
    p.add_argument("--seed", type=int)
    p.add_argument("--workers", type=int, help="processes (0 = all cores)")
    p.set_defaults(func=cmd_run)

    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                        format="%(levelname)s %(name)s: %(message)s")
    if not getattr(args, "func", None):
        parser.print_help()
        return 0
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
