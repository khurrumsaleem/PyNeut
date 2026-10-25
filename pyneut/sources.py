"""Source distributions for the input format.

A source is three independent distributions: where a neutron starts
(`space`), which way it goes (`direction`) and its energy (`energy`). Each is
built from its input table by `build_space`, `build_direction` and
`build_energy`, and every one samples `n` values at once from a
`numpy.random.Generator`.

`Source.sample(n, rng)` returns the vector engine's source dict (arrays x, y,
z, u, v, w, energy); `to_states` converts it to the scalar engine's
per-particle dicts.
"""
import numpy as np

from .physics import WATT_PARAMS
from .units import parse_energy


def _vec3(value, what):
    if (not isinstance(value, (list, tuple)) or len(value) != 3
            or not all(isinstance(v, (int, float)) and not isinstance(v, bool) for v in value)):
        raise ValueError(f"{what} must be a list of three numbers, got {value!r}")
    return np.array(value, dtype=float)


def _number(value, what, positive=False):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{what} must be a number, got {value!r}")
    if positive and value <= 0:
        raise ValueError(f"{what} must be positive, got {value!r}")
    return float(value)


def _isotropic(n, rng):
    mu = 1.0 - 2.0 * rng.random(n)
    phi = 2.0 * np.pi * rng.random(n)
    s = np.sqrt(np.maximum(0.0, 1.0 - mu * mu))
    return s * np.cos(phi), s * np.sin(phi), mu


# --------------------------------------------------------------------- space
class PointSpace:
    keys = {"xyz"}

    def __init__(self, xyz):
        self.xyz = _vec3(xyz, "xyz")

    def sample(self, n, rng):
        return tuple(np.full(n, c) for c in self.xyz)


class BoxSpace:
    keys = {"min", "max"}

    def __init__(self, min, max):
        self.lo, self.hi = _vec3(min, "min"), _vec3(max, "max")
        if np.any(self.hi <= self.lo):
            raise ValueError("box max must exceed min in every coordinate")

    def sample(self, n, rng):
        return tuple(self.lo[i] + (self.hi[i] - self.lo[i]) * rng.random(n) for i in range(3))


class SphereSpace:
    keys = {"center", "radius"}

    def __init__(self, center, radius):
        self.center = _vec3(center, "center")
        self.radius = _number(radius, "radius", positive=True)

    def sample(self, n, rng):
        r = self.radius * np.cbrt(rng.random(n))
        u, v, w = _isotropic(n, rng)
        return (self.center[0] + r * u, self.center[1] + r * v, self.center[2] + r * w)


class CylinderSpace:
    """Uniform in a finite cylinder: `axis`, `center`, `radius`, and the
    coordinate `range = [lo, hi]` along the axis."""
    keys = {"axis", "center", "radius", "range"}

    def __init__(self, axis, center, radius, range):
        if axis not in ("x", "y", "z"):
            raise ValueError(f"axis must be 'x', 'y' or 'z', got {axis!r}")
        self.axis = "xyz".index(axis)
        self.center = _vec3(center, "center")
        self.radius = _number(radius, "radius", positive=True)
        if (not isinstance(range, (list, tuple)) or len(range) != 2
                or not range[1] > range[0]):
            raise ValueError(f"range must be [lo, hi] with hi > lo, got {range!r}")
        self.range = (float(range[0]), float(range[1]))

    def sample(self, n, rng):
        r = self.radius * np.sqrt(rng.random(n))
        ang = 2.0 * np.pi * rng.random(n)
        along = self.range[0] + (self.range[1] - self.range[0]) * rng.random(n)
        a, b = [i for i in range(3) if i != self.axis]
        out = [None, None, None]
        out[a] = self.center[a] + r * np.cos(ang)
        out[b] = self.center[b] + r * np.sin(ang)
        out[self.axis] = along
        return tuple(out)


class CellSpace:
    """Uniform inside a cell, by rejection inside `bounds = {min, max}`.

    `contains(x, y, z)` is attached by the model once cells are built.
    """
    keys = {"cell", "bounds"}

    def __init__(self, cell, bounds):
        if not isinstance(cell, str):
            raise ValueError(f"cell must be a cell name, got {cell!r}")
        if not isinstance(bounds, dict) or set(bounds) != {"min", "max"}:
            raise ValueError("bounds must be a table with 'min' and 'max'")
        self.cell = cell
        self.box = BoxSpace(bounds["min"], bounds["max"])
        self.contains = None

    def sample(self, n, rng, max_rounds=1000):
        xs, ys, zs, have = [], [], [], 0
        for _ in range(max_rounds):
            x, y, z = self.box.sample(max(2 * (n - have), 64), rng)
            keep = np.array([self.contains(a, b, c) for a, b, c in zip(x, y, z)])
            xs.append(x[keep]); ys.append(y[keep]); zs.append(z[keep])
            have += int(keep.sum())
            if have >= n:
                break
        else:
            raise ValueError(f"could not place source points in cell {self.cell!r}; "
                             "check that its bounds overlap the cell")
        return tuple(np.concatenate(a)[:n] for a in (xs, ys, zs))


SPACES = {"point": PointSpace, "box": BoxSpace, "sphere": SphereSpace,
          "cylinder": CylinderSpace, "cell": CellSpace}


# ----------------------------------------------------------------- direction
class Isotropic:
    def sample(self, n, rng):
        return _isotropic(n, rng)


class Beam:
    keys = {"direction"}

    def __init__(self, direction):
        d = _vec3(direction, "direction")
        norm = np.linalg.norm(d)
        if norm == 0:
            raise ValueError("beam direction cannot be zero")
        self.d = d / norm

    def sample(self, n, rng):
        return tuple(np.full(n, c) for c in self.d)


# -------------------------------------------------------------------- energy
class Monoenergetic:
    def __init__(self, energy):
        self.energy = parse_energy(energy)
        if self.energy <= 0:
            raise ValueError("energy must be positive")

    def sample(self, n, rng):
        return np.full(n, self.energy)


class Discrete:
    keys = {"energies", "probabilities"}

    def __init__(self, energies, probabilities):
        if not isinstance(energies, list) or not energies:
            raise ValueError("energies must be a non-empty list")
        self.energies = np.array([parse_energy(e) for e in energies])
        p = np.array([_number(q, "probability") for q in probabilities])
        if p.size != self.energies.size or np.any(p < 0) or p.sum() <= 0:
            raise ValueError("probabilities must be non-negative, one per energy")
        self.p = p / p.sum()

    def sample(self, n, rng):
        return self.energies[rng.choice(self.energies.size, size=n, p=self.p)]


class Uniform:
    keys = {"min", "max"}

    def __init__(self, min, max):
        self.lo, self.hi = parse_energy(min), parse_energy(max)
        if not 0 < self.lo < self.hi:
            raise ValueError("uniform energy needs 0 < min < max")

    def sample(self, n, rng):
        return self.lo + (self.hi - self.lo) * rng.random(n)


def _maxwell(n, rng, kT):
    # f(E) ~ sqrt(E) exp(-E/kT) (algorithm C64, as physics.sample_maxwellian)
    r1, r2, r3 = (1.0 - rng.random(n) for _ in range(3))
    c = np.cos(0.5 * np.pi * r3)
    return -kT * (np.log(r1) + np.log(r2) * c * c)


class Maxwell:
    keys = {"kT"}

    def __init__(self, kT):
        self.kT = parse_energy(kT)
        if self.kT <= 0:
            raise ValueError("kT must be positive")

    def sample(self, n, rng):
        return _maxwell(n, rng, self.kT)


class Watt:
    """chi(E) ~ exp(-E/a) sinh(sqrt(b E)); `a` in eV, `b` in 1/eV, or the
    tabulated parameters of `nuclide` (U235, U233, Pu239)."""
    keys = {"nuclide", "a", "b"}

    def __init__(self, nuclide=None, a=None, b=None):
        if nuclide is not None:
            if a is not None or b is not None:
                raise ValueError("give either nuclide or a and b, not both")
            if nuclide not in WATT_PARAMS:
                known = ", ".join(WATT_PARAMS)
                raise ValueError(f"no Watt parameters for {nuclide!r} (known: {known}); give a and b")
            self.a, self.b = WATT_PARAMS[nuclide]
        else:
            if a is None or b is None:
                raise ValueError("watt needs nuclide, or both a and b")
            self.a, self.b = parse_energy(a), _number(b, "b", positive=True)

    def sample(self, n, rng):
        w = _maxwell(n, rng, self.a)
        a2b = self.a * self.a * self.b
        return w + a2b / 4.0 + (2.0 * rng.random(n) - 1.0) * np.sqrt(a2b * w)


ENERGIES = {"discrete": Discrete, "uniform": Uniform, "maxwell": Maxwell, "watt": Watt}


# --------------------------------------------------------------------- build
def _build(table, kinds, what):
    if not isinstance(table, dict):
        raise ValueError(f"{what} must be a table with a 'type'")
    if "type" not in table:
        raise ValueError(f"{what} needs a 'type' ({', '.join(kinds)})")
    kind = table["type"]
    if kind not in kinds:
        raise ValueError(f"unknown {what} type {kind!r}; expected one of {', '.join(kinds)}")
    cls = kinds[kind]
    params = {k: v for k, v in table.items() if k != "type"}
    unknown = set(params) - cls.keys
    if unknown:
        raise ValueError(f"unknown key(s) {', '.join(sorted(unknown))} for {what} "
                         f"type {kind!r}; allowed: {', '.join(sorted(cls.keys))}")
    try:
        return cls(**params)
    except TypeError:
        missing = sorted(set(cls.__init__.__code__.co_varnames[1:cls.__init__.__code__.co_argcount])
                         - set(params))
        raise ValueError(f"{what} type {kind!r} is missing {', '.join(missing)}") from None


def build_space(table):
    return _build(table, SPACES, "space")


def build_direction(value):
    if value == "isotropic" or value == {"type": "isotropic"}:
        return Isotropic()
    return _build(value, {"beam": Beam}, "direction")


def build_energy(value):
    if isinstance(value, (int, float, str)) and not isinstance(value, bool):
        return Monoenergetic(value)
    return _build(value, ENERGIES, "energy")


class Source:
    def __init__(self, space, direction, energy):
        self.space, self.direction, self.energy = space, direction, energy

    def sample(self, n, rng):
        x, y, z = self.space.sample(n, rng)
        u, v, w = self.direction.sample(n, rng)
        e = self.energy.sample(n, rng)
        return {"x": x, "y": y, "z": z, "u": u, "v": v, "w": w, "energy": e}


def to_states(arrays):
    """Vector-engine source arrays -> scalar-engine particle dicts."""
    w = np.clip(arrays["w"], -1.0, 1.0)
    theta = np.arccos(w)
    phi = np.arctan2(arrays["v"], arrays["u"])
    return [{"x": float(arrays["x"][i]), "y": float(arrays["y"][i]),
             "z": float(arrays["z"][i]), "theta": float(theta[i]),
             "phi": float(phi[i]), "energy": float(arrays["energy"][i]),
             "weight": 1.0, "has_interacted": False}
            for i in range(len(theta))]

