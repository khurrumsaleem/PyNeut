"""The `pyneut` command.

Subcommands (`run`, `check`, `plot`) arrive with the TOML input format; see
docs/design/input-format.md.
"""
import argparse

from . import __version__


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="pyneut",
        description="Continuous-energy Monte Carlo neutron transport.")
    parser.add_argument("--version", action="version",
                        version=f"pyneut {__version__}")
    parser.parse_args(argv)
    parser.print_help()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
