"""Unit handling for model inputs.

Lengths are always plain numbers in cm. Energies may be a plain number in eV
or a string with a unit: "2 MeV", "25 keV", "0.0253 eV".
"""
import re

ENERGY_UNITS = {"eV": 1.0, "keV": 1.0e3, "MeV": 1.0e6}

_ENERGY_RE = re.compile(r"^\s*([-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)\s*([a-zA-Z]+)\s*$")


def parse_energy(value):
    """Return an energy in eV from a number (eV) or a string with a unit.

    Raises ValueError with a readable message for anything else.
    """
    if isinstance(value, bool):
        raise ValueError(f"expected an energy, got {value!r}")
    if isinstance(value, (int, float)):
        return float(value)
    if isinstance(value, str):
        m = _ENERGY_RE.match(value)
        if m and m.group(2) in ENERGY_UNITS:
            return float(m.group(1)) * ENERGY_UNITS[m.group(2)]
        units = ", ".join(ENERGY_UNITS)
        raise ValueError(f"cannot read energy {value!r}; write a number in eV "
                         f"or a number with a unit ({units}), e.g. \"2 MeV\"")
    raise ValueError(f"expected an energy, got {value!r}")
