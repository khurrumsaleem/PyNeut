"""A complete PyNeut calculation: read from TOML, checked, solved.

    model = Model.from_toml("pb208_sphere.toml")
    model.run.particles = 20_000        # anything can be overridden
    results = model.solve()
    results.save("out/")

The input format is documented in docs/input-format.md. `Model.from_dict`
takes the same structure as a Python dict.

Reproducibility: the master `run.seed` feeds a numpy SeedSequence that is
split into one stream for sampling the source and one for transport. The
scalar engine gives every history its own child stream, so a result does not
depend on how a process pool schedules histories.
"""
import difflib
import hashlib
import json
import logging
import os
import platform
import re
import subprocess
import time
from dataclasses import asdict, dataclass

import numpy as np

from . import __version__
from .data import (check_nuclide_name, data_manifest, nuclide_file, read_awr,
                   resolve_data_dir)
from .errors import InputError, ProblemList
from .material import Material, Nuclide
from .medium import Box, Cylinder, Plane, Region, Sphere
from .regions import compile_region, parse_region, surface_names
from .sources import (CellSpace, Source, build_direction, build_energy,
                      build_space, to_states)
from .statistics import batch_edges, mean_sem
from .units import parse_energy

try:
    import tomllib
except ModuleNotFoundError:  # Python 3.10
    import tomli as tomllib

logger = logging.getLogger(__name__)

SECTIONS = ("run", "materials", "surfaces", "cells", "source", "tallies")
_NAME_RE = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")


# ------------------------------------------------------------------ helpers
def _suggest(word, options):
    match = difflib.get_close_matches(word, list(options), n=1)
    return f" (did you mean {match[0]!r}?)" if match else ""


def _check_keys(table, allowed, where, problems, required=()):
    if not isinstance(table, dict):
        problems.add(where, f"expected a table, got {table!r}")
        return False
    for key in table:
        if key not in allowed:
            problems.add(f"{where}.{key}" if where else key,
                         "unknown key" + _suggest(key, allowed))
    for key in required:
        if key not in table:
            problems.add(where, f"missing required key {key!r}")
    return True


def _is_int(v):
    return isinstance(v, int) and not isinstance(v, bool)


def _is_num(v):
    return isinstance(v, (int, float)) and not isinstance(v, bool)


def _vec3(v):
    return isinstance(v, list) and len(v) == 3 and all(_is_num(c) for c in v)


# --------------------------------------------------------------------- run
@dataclass
class RunSettings:
    mode: str
    particles: int
    engine: str = "vector"
    batches: int = 10
    generations: int = None
    inactive: int = None
    seed: int = 1
    transport: str = None          # resolved: implicit-capture (fixed source) / analog
    weight_cutoff: float = 1e-4
    survival: float = 0.1
    workers: int = 1
    data: str = None
    temperature: str = "294K"

    @property
    def implicit_capture(self):
        return self.transport == "implicit-capture"

    def settings(self):
        """The low-level Settings object the engines take."""
        from .settings import Settings
        if self.mode == "eigenvalue" or not self.implicit_capture:
            s = Settings("criticality", self.particles, self.batches)
        else:
            s = Settings("shielding", self.particles, self.batches)
            s.weight_cutoff = self.weight_cutoff
            s.roulette_survival_prob = self.survival
        return s


_RUN_KEYS = {f for f in RunSettings.__dataclass_fields__}


def _parse_run(t, problems):
    if not _check_keys(t, _RUN_KEYS, "run", problems, required=("mode", "particles")):
        return None
    mode = t.get("mode")
    if mode not in ("fixed-source", "eigenvalue"):
        problems.add("run.mode", f"must be 'fixed-source' or 'eigenvalue', got {mode!r}")
    engine = t.get("engine", "vector")
    if engine not in ("vector", "scalar"):
        problems.add("run.engine", f"must be 'vector' or 'scalar', got {engine!r}")
    for key, minimum in (("particles", 1), ("batches", 2), ("seed", 0), ("workers", 0)):
        if key in t and (not _is_int(t[key]) or t[key] < minimum):
            problems.add(f"run.{key}", f"must be an integer >= {minimum}, got {t[key]!r}")
    transport = t.get("transport")
    if mode == "eigenvalue":
        for key in ("generations", "inactive"):
            if key not in t:
                problems.add("run", f"eigenvalue mode needs {key!r}")
            elif not _is_int(t[key]) or t[key] < (1 if key == "generations" else 0):
                problems.add(f"run.{key}", f"must be a non-negative integer, got {t[key]!r}")
        if (_is_int(t.get("generations")) and _is_int(t.get("inactive"))
                and t["inactive"] >= t["generations"]):
            problems.add("run.inactive", "must be smaller than generations")
        if transport not in (None, "analog"):
            problems.add("run.transport", "eigenvalue mode is always analog")
        transport = "analog"
    else:
        for key in ("generations", "inactive"):
            if key in t:
                problems.add(f"run.{key}", "only used in eigenvalue mode")
        transport = transport or "implicit-capture"
        if transport not in ("implicit-capture", "analog"):
            problems.add("run.transport",
                         f"must be 'implicit-capture' or 'analog', got {transport!r}")
    if "weight_cutoff" in t and (not _is_num(t["weight_cutoff"]) or t["weight_cutoff"] < 0):
        problems.add("run.weight_cutoff", "must be a non-negative number")
    if "survival" in t and (not _is_num(t["survival"]) or not 0 < t["survival"] <= 1):
        problems.add("run.survival", "must be in (0, 1]")
    if "data" in t and not isinstance(t["data"], str):
        problems.add("run.data", "must be a path string")
    temp = t.get("temperature", "294K")
    if not (isinstance(temp, str) and re.fullmatch(r"\d+(\.\d+)?K", temp)):
        problems.add("run.temperature", f"must look like \"294K\", got {temp!r}")
    kwargs = {k: v for k, v in t.items() if k in _RUN_KEYS}
    kwargs["transport"] = transport
    try:
        return RunSettings(**kwargs)
    except TypeError:
        return None


# --------------------------------------------------------------- materials
@dataclass
class MaterialSpec:
    name: str
    density: float = None             # g/cm3
    nuclides: dict = None             # name -> fraction
    fractions: str = "atom"
    atom_densities: dict = None       # name -> atoms/b-cm

    def nuclide_names(self):
        return list((self.nuclides or self.atom_densities or {}).keys())

    def build(self, awr, temperature):
        """list[Nuclide] given {nuclide: AWR} from the data files."""
        if self.atom_densities is not None:
            return [Nuclide(n, d * 1e24, awr[n], temperature=temperature)
                    for n, d in self.atom_densities.items()]
        if self.fractions == "atom":
            return Material.mixture(self.name, self.density,
                                    [(n, awr[n], f) for n, f in self.nuclides.items()],
                                    temperature=temperature)
        # weight fractions: N_i = rho * N_A * w_i / M_i
        from .material import AVOGADRO, NEUTRON_MASS_AMU
        total = sum(self.nuclides.values())
        return [Nuclide(n, self.density * AVOGADRO * (w / total) / (awr[n] * NEUTRON_MASS_AMU),
                        awr[n], temperature=temperature)
                for n, w in self.nuclides.items()]


def _parse_nuclide_table(table, where, problems):
    if not isinstance(table, dict) or not table:
        problems.add(where, "must be a non-empty table like { U235 = 1.0 }")
        return False
    ok = True
    for n, v in table.items():
        msg = check_nuclide_name(n)
        if msg:
            problems.add(f"{where}.{n}", msg)
            ok = False
        if not _is_num(v) or v <= 0:
            problems.add(f"{where}.{n}", f"must be a positive number, got {v!r}")
            ok = False
    return ok


def _parse_material(name, t, problems):
    where = f"materials.{name}"
    if name == "void":
        problems.add(where, "'void' is reserved for empty cells")
        return None
    if not _check_keys(t, {"density", "nuclides", "fractions", "atom_densities"}, where, problems):
        return None
    if "atom_densities" in t:
        for key in ("density", "nuclides", "fractions"):
            if key in t:
                problems.add(f"{where}.{key}", "not used with atom_densities")
        if _parse_nuclide_table(t["atom_densities"], f"{where}.atom_densities", problems):
            return MaterialSpec(name, atom_densities=dict(t["atom_densities"]))
        return None
    if "nuclides" not in t or "density" not in t:
        problems.add(where, "give 'density' and 'nuclides', or 'atom_densities'")
        return None
    if not _is_num(t["density"]) or t["density"] <= 0:
        problems.add(f"{where}.density", "must be a positive number (g/cm3)")
        return None
    fractions = t.get("fractions", "atom")
    if fractions not in ("atom", "weight"):
        problems.add(f"{where}.fractions", f"must be 'atom' or 'weight', got {fractions!r}")
        return None
    if _parse_nuclide_table(t["nuclides"], f"{where}.nuclides", problems):
        return MaterialSpec(name, density=float(t["density"]),
                            nuclides=dict(t["nuclides"]), fractions=fractions)
    return None


# ---------------------------------------------------------------- surfaces
_SURFACE_PARAMS = {
    "plane": ("normal", "offset"),
    "x-plane": ("at",), "y-plane": ("at",), "z-plane": ("at",),
    "sphere": ("center", "radius"),
    "cylinder": ("axis", "center", "radius"),
    "box": ("min", "max"),
}


def _parse_surface(name, t, problems):
    where = f"surfaces.{name}"
    if not _NAME_RE.match(name):
        problems.add(where, "surface names may use letters, digits and '_' only "
                            "(they appear in region expressions)")
        return None
    if not isinstance(t, dict) or "type" not in t:
        problems.add(where, f"needs a 'type' ({', '.join(_SURFACE_PARAMS)})")
        return None
    kind = t["type"]
    if kind not in _SURFACE_PARAMS:
        problems.add(f"{where}.type", f"unknown surface type {kind!r}"
                     + _suggest(str(kind), _SURFACE_PARAMS))
        return None
    params = _SURFACE_PARAMS[kind]
    if not _check_keys(t, {"type", "boundary", *params}, where, problems, required=params):
        return None
    boundary = t.get("boundary", "vacuum")
    if boundary not in ("vacuum", "reflective"):
        problems.add(f"{where}.boundary", f"must be 'vacuum' or 'reflective', got {boundary!r}")
        return None
    bt = "reflective" if boundary == "reflective" else "transmission"
    n_before = len(problems.items)
    for key in params:
        v = t.get(key)
        if key in ("center", "normal", "min", "max") and not _vec3(v):
            problems.add(f"{where}.{key}", f"must be [x, y, z], got {v!r}")
        elif key in ("radius",) and (not _is_num(v) or v <= 0):
            problems.add(f"{where}.{key}", f"must be a positive number, got {v!r}")
        elif key in ("offset", "at") and not _is_num(v):
            problems.add(f"{where}.{key}", f"must be a number, got {v!r}")
        elif key == "axis" and v not in ("x", "y", "z"):
            problems.add(f"{where}.axis", f"must be 'x', 'y' or 'z', got {v!r}")
    if kind == "plane" and _vec3(t.get("normal")) and not any(t["normal"]):
        problems.add(f"{where}.normal", "cannot be zero")
    if kind == "box" and _vec3(t.get("min")) and _vec3(t.get("max")) \
            and not all(b > a for a, b in zip(t["min"], t["max"])):
        problems.add(where, "box max must exceed min in every coordinate")
    if len(problems.items) > n_before:
        return None
    if kind == "plane":
        a, b, c = t["normal"]
        return Plane(a, b, c, t["offset"], boundary_type=bt)
    if kind.endswith("-plane"):
        normal = {"x": (1, 0, 0), "y": (0, 1, 0), "z": (0, 0, 1)}[kind[0]]
        return Plane(*normal, t["at"], boundary_type=bt)
    if kind == "sphere":
        return Sphere(tuple(t["center"]), t["radius"], boundary_type=bt)
    if kind == "cylinder":
        return Cylinder(t["axis"], t["radius"], tuple(t["center"]), boundary_type=bt)
    (x0, y0, z0), (x1, y1, z1) = t["min"], t["max"]
    return Box(x0, x1, y0, y1, z0, z1, boundary_type=bt)


# ------------------------------------------------------------------- cells
@dataclass
class CellSpec:
    name: str
    region_text: str
    material: str
    priority: int = None
    expr: object = None
    geometry: object = None           # compiled surface/Region


def _parse_cell(name, t, surfaces, materials, problems):
    where = f"cells.{name}"
    if not _check_keys(t, {"region", "material", "priority"}, where, problems,
                       required=("region", "material")):
        return None
    ok = True
    try:
        expr = parse_region(t.get("region"))
    except ValueError as e:
        problems.add(f"{where}.region", str(e))
        expr, ok = None, False
    if expr is not None:
        for s in sorted(surface_names(expr)):
            if s not in surfaces:
                problems.add(f"{where}.region", f"unknown surface {s!r}" + _suggest(s, surfaces))
                ok = False
            elif surfaces[s] is None:      # the surface itself was invalid (reported)
                ok = False
    mat = t.get("material")
    if mat != "void" and mat not in materials:
        problems.add(f"{where}.material", f"unknown material {mat!r}"
                     + _suggest(str(mat), [*materials, "void"]))
        ok = False
    pri = t.get("priority")
    if pri is not None and not _is_int(pri):
        problems.add(f"{where}.priority", f"must be an integer, got {pri!r}")
        ok = False
    if not ok:
        return None
    geometry = compile_region(expr, surfaces)
    return CellSpec(name, t["region"], mat, pri, expr, geometry)


# ------------------------------------------------------------------ source
def _parse_source(t, cells, problems):
    if not _check_keys(t, {"space", "direction", "energy"}, "source", problems,
                       required=("space", "energy")):
        return None
    parts = {}
    for key, build in (("space", build_space), ("direction", build_direction),
                       ("energy", build_energy)):
        try:
            parts[key] = build(t.get(key, "isotropic"))
        except ValueError as e:
            problems.add(f"source.{key}", str(e))
    if len(parts) < 3:
        return None
    if isinstance(parts["space"], CellSpace):
        cell = cells.get(parts["space"].cell)
        if cell is None:
            problems.add("source.space.cell", f"unknown cell {parts['space'].cell!r}"
                         + _suggest(parts["space"].cell, cells))
            return None
        parts["space"].contains = Region(surfaces=[cell.geometry]).contains
    return Source(parts["space"], parts["direction"], parts["energy"])


# ----------------------------------------------------------------- tallies
@dataclass
class TallySpec:
    name: str
    type: str
    energy_bins: np.ndarray = None
    mesh: dict = None
    box: dict = None
    cells: list = None
    estimator: str = None

    @property
    def scored(self):
        """True for tallies the event Scorer handles (not keff/detector)."""
        return self.type in ("leakage", "absorption", "flux")


_TALLY_KEYS = {
    "leakage": {"type", "energy_bins"},
    "absorption": {"type", "cells", "energy_bins", "mesh"},
    "flux": {"type", "cells", "energy_bins", "mesh", "estimator"},
    "detector": {"type", "box"},
    "keff": {"type"},
}
_PLANNED = {"current": "surface-current tallies"}


def _parse_bins(v):
    if isinstance(v, list):
        edges = np.array([parse_energy(e) for e in v])
    elif isinstance(v, dict) and len(v) == 1 and next(iter(v)) in ("log", "lin"):
        kind, spec = next(iter(v.items()))
        if not (isinstance(spec, list) and len(spec) == 3 and _is_int(spec[2]) and spec[2] > 0):
            raise ValueError(f"{{{kind} = [min, max, n_bins]}} expected, got {spec!r}")
        lo, hi = parse_energy(spec[0]), parse_energy(spec[1])
        if kind == "log":
            if lo <= 0:
                raise ValueError("log bins need min > 0")
            edges = np.logspace(np.log10(lo), np.log10(hi), spec[2] + 1)
        else:
            edges = np.linspace(lo, hi, spec[2] + 1)
    else:
        raise ValueError("energy_bins is a list of edges, {log = [min, max, n]} "
                         "or {lin = [min, max, n]}")
    if edges.size < 2 or np.any(np.diff(edges) <= 0) or edges[0] < 0:
        raise ValueError("energy bin edges must be non-negative and increasing")
    return edges


def _parse_mesh(m):
    if not (isinstance(m, dict) and set(m) == {"min", "max", "dims"} and _vec3(m["min"])
            and _vec3(m["max"]) and isinstance(m["dims"], list) and len(m["dims"]) == 3
            and all(_is_int(d) and d > 0 for d in m["dims"])
            and all(b > a for a, b in zip(m["min"], m["max"]))):
        raise ValueError("must be {min = [x,y,z], max = [x,y,z], dims = [nx,ny,nz]}")
    return m


def _parse_tally(name, t, run, cells, problems):
    where = f"tallies.{name}"
    if not isinstance(t, dict) or "type" not in t:
        problems.add(where, f"needs a 'type' ({', '.join(_TALLY_KEYS)})")
        return None
    kind = t["type"]
    if kind in _PLANNED:
        problems.add(f"{where}.type", f"{_PLANNED[kind]} are not available yet")
        return None
    if kind not in _TALLY_KEYS:
        problems.add(f"{where}.type", f"unknown tally type {kind!r}" + _suggest(str(kind), _TALLY_KEYS))
        return None
    if not _check_keys(t, _TALLY_KEYS[kind], where, problems):
        return None
    spec = TallySpec(name, kind)
    if run is not None:
        if run.mode == "eigenvalue" and kind == "detector":
            problems.add(where, "detector tallies are fixed-source only")
            return None
        if run.mode == "fixed-source" and kind == "keff":
            problems.add(where, "keff tallies need run.mode = \"eigenvalue\"")
            return None
    ok = True
    if "energy_bins" in t:
        try:
            spec.energy_bins = _parse_bins(t["energy_bins"])
        except ValueError as e:
            problems.add(f"{where}.energy_bins", str(e))
            ok = False
    if "mesh" in t:
        try:
            spec.mesh = _parse_mesh(t["mesh"])
        except ValueError as e:
            problems.add(f"{where}.mesh", str(e))
            ok = False
    if "cells" in t:
        names = t["cells"]
        if not (isinstance(names, list) and names and all(isinstance(c, str) for c in names)):
            problems.add(f"{where}.cells", "must be a non-empty list of cell names")
            ok = False
        else:
            for c in names:
                if c not in cells:
                    problems.add(f"{where}.cells", f"unknown cell {c!r}" + _suggest(c, cells))
                    ok = False
            if len(set(names)) != len(names):
                problems.add(f"{where}.cells", "lists a cell twice")
                ok = False
            spec.cells = list(names)
    if kind == "flux":
        est = t.get("estimator", "collision" if spec.mesh else "track-length")
        if est not in ("track-length", "collision"):
            problems.add(f"{where}.estimator", f"must be 'track-length' or 'collision', got {est!r}")
            ok = False
        elif est == "track-length" and spec.mesh:
            problems.add(f"{where}.estimator", "a mesh flux uses the collision estimator "
                         "(track-length mesh flux is not available yet)")
            ok = False
        spec.estimator = est
    if kind == "detector":
        b = t.get("box")
        if not (isinstance(b, dict) and set(b) == {"min", "max"} and _vec3(b["min"])
                and _vec3(b["max"]) and all(hi > lo for lo, hi in zip(b["min"], b["max"]))):
            problems.add(f"{where}.box", "must be {min = [x,y,z], max = [x,y,z]}")
            return None
        spec.box = b
        if run is not None and run.engine == "vector":
            problems.add(where, "detector tallies need run.engine = \"scalar\" for now")
            return None
    return spec if ok else None


# ------------------------------------------------------------------- model
class Model:
    """A parsed, validated model. Build with `from_toml` or `from_dict`."""

    def __init__(self, spec, base_dir=".", input_path=None, input_text=None):
        problems = ProblemList()
        if not isinstance(spec, dict):
            raise InputError("the input must be a table")
        for key in spec:
            if key not in SECTIONS:
                problems.add(key, "unknown section" + _suggest(key, SECTIONS))
        missing = [key for key in ("run", "materials", "cells", "source") if key not in spec]
        if missing:
            for key in missing:
                problems.add("", f"missing [{key}] section")
            problems.raise_if_any()

        self.base_dir = os.path.abspath(base_dir)
        self.input_path = input_path
        self.input_text = input_text
        self.spec = spec

        self.run = _parse_run(spec["run"], problems)

        self.materials = {}
        mat_tables = spec["materials"]
        if not isinstance(mat_tables, dict):
            problems.add("materials", "must be a table of named materials")
            mat_tables = {}
        for name, t in mat_tables.items():
            m = _parse_material(name, t, problems)
            if m is not None:
                self.materials[name] = m

        self.surfaces = {}
        surf_tables = spec.get("surfaces", {})
        if isinstance(surf_tables, dict):
            for name, t in surf_tables.items():
                self.surfaces[name] = _parse_surface(name, t, problems)
        else:
            problems.add("surfaces", "must be a table of named surfaces")

        self.cells = {}
        cell_tables = spec["cells"]
        if not isinstance(cell_tables, dict) or not cell_tables:
            problems.add("cells", "must define at least one cell")
        else:
            known_mats = set(mat_tables)   # an invalid material is reported once, above
            for name, t in cell_tables.items():
                c = _parse_cell(name, t, self.surfaces, known_mats, problems)
                if c is not None:
                    self.cells[name] = c

        self.source = _parse_source(spec["source"], self.cells, problems)

        self.tallies = {}
        tally_tables = spec.get("tallies", {})
        if not isinstance(tally_tables, dict):
            problems.add("tallies", "must be a table of named tallies")
            tally_tables = {}
        for name, t in tally_tables.items():
            tally = _parse_tally(name, t, self.run, self.cells, problems)
            if tally is not None:
                self.tallies[name] = tally
        if sum(1 for t in self.tallies.values() if t.type == "detector") > 1:
            problems.add("tallies", "only one detector tally is supported for now")
        if self.run is not None and self.run.mode == "eigenvalue" and \
                not any(t.type == "keff" for t in self.tallies.values()):
            self.tallies["keff"] = TallySpec("keff", "keff")
        if self.run is not None and self.run.mode == "fixed-source" and not self.tallies:
            self.tallies["leakage"] = TallySpec("leakage", "leakage")

        problems.raise_if_any()
        self._built = None

    # ---------------------------------------------------------- constructors
    @classmethod
    def from_toml(cls, path):
        with open(path, "rb") as f:
            raw = f.read()
        try:
            spec = tomllib.loads(raw.decode("utf-8"))
        except tomllib.TOMLDecodeError as e:
            raise InputError(f"{path}: not valid TOML: {e}") from None
        return cls(spec, base_dir=os.path.dirname(os.path.abspath(path)),
                   input_path=os.path.abspath(path), input_text=raw)

    @classmethod
    def from_dict(cls, spec, base_dir="."):
        return cls(spec, base_dir=base_dir)

    # ----------------------------------------------------------------- data
    @property
    def data_dir(self):
        return resolve_data_dir(self.run.data, self.base_dir)

    def nuclides(self):
        names = []
        for m in self.materials.values():
            for n in m.nuclide_names():
                if n not in names:
                    names.append(n)
        return names

    def build(self):
        """Resolve nuclear data and build the transport objects.

        Returns (mediums, compositions). Raises InputError listing every
        nuclide whose data file is missing.
        """
        if self._built is not None:
            return self._built
        data_dir = self.data_dir
        problems = ProblemList()
        if not os.path.isdir(os.path.join(data_dir, "neutron")):
            problems.add("run.data", f"no nuclear data at {data_dir} (expected "
                         f"{os.path.join(data_dir, 'neutron')}/<nuclide>.h5); set run.data, "
                         "$PYNEUT_DATA, or run from a folder containing endfb/")
            problems.raise_if_any()
        awr = {}
        for n in self.nuclides():
            try:
                awr[n] = read_awr(data_dir, n)
            except FileNotFoundError:
                problems.add(f"nuclide {n}", f"no data file {nuclide_file(data_dir, n)}")
            except (KeyError, OSError) as e:
                problems.add(f"nuclide {n}", f"cannot read its atomic-weight ratio: {e}")
        problems.raise_if_any()
        temperature = float(self.run.temperature[:-1])
        compositions = {name: m.build(awr, temperature) for name, m in self.materials.items()}
        mediums = []
        for c in self.cells.values():
            mediums.append(Region(surfaces=[c.geometry], operation="intersection",
                                  name=c.name, priority=c.priority or 0,
                                  is_void=(c.material == "void"),
                                  composition=None if c.material == "void"
                                  else compositions[c.material]))
        self._built = (mediums, compositions, awr)
        return self._built

    # ---------------------------------------------------------------- solve
    def solve(self):
        """Run the calculation and return a Results object."""
        mediums, _, awr = self.build()
        run = self.run
        ss = np.random.SeedSequence(run.seed)
        source_ss, transport_ss = ss.spawn(2)
        n = run.particles
        source = self.source.sample(n, np.random.default_rng(source_ss))
        from .cross_section_read import CrossSectionReader
        reader = CrossSectionReader(self.data_dir, run.temperature)
        workers = run.workers or os.cpu_count()

        t0 = time.perf_counter()
        if run.mode == "eigenvalue":
            out = self._solve_eigenvalue(reader, mediums, source, transport_ss, workers)
        else:
            out = self._solve_fixed_source(reader, mediums, source, transport_ss, workers)
        wall = time.perf_counter() - t0

        histories = n * (run.generations if run.mode == "eigenvalue" else 1)
        info = {"mode": run.mode, "engine": run.engine, "transport": run.transport,
                "particles": n, "batches": run.batches if run.mode == "fixed-source" else None,
                "generations": run.generations, "inactive": run.inactive,
                "seed": run.seed, "workers": workers, "wall_time_s": wall,
                "histories_per_s": histories / wall if wall > 0 else None,
                "counters": out.pop("_counters", {})}
        return Results(out, info, self._provenance(awr))

    def _tally_defs(self):
        """TallyDefs for the event Scorer, in the order of self.cells."""
        from .tallies import TallyDef
        names = list(self.cells)
        defs = []
        for t in self.tallies.values():
            if not t.scored:
                continue
            cell_map = None
            if t.cells:
                cell_map = np.full(len(names), -1, dtype=np.int64)
                for i, c in enumerate(t.cells):
                    cell_map[names.index(c)] = i
            defs.append(TallyDef(t.name, t.type, cell_map=cell_map, cell_names=t.cells,
                                 mesh=t.mesh, energy_bins=t.energy_bins,
                                 estimator=t.estimator))
        return defs

    def _solve_fixed_source(self, reader, mediums, source, transport_ss, workers):
        from .tallies import Scorer
        run = self.run
        n = run.particles
        settings = run.settings()
        edges = batch_edges(n, run.batches)
        n_batches = len(edges) - 1
        hist_batch = np.searchsorted(edges, np.arange(n), side="right") - 1
        scorer = Scorer(self._tally_defs(), n_batches)
        detector = next((t for t in self.tallies.values() if t.type == "detector"), None)

        if run.engine == "scalar":
            bounds = None
            if detector is not None:
                (x0, y0, z0), (x1, y1, z1) = detector.box["min"], detector.box["max"]
                bounds = (x0, x1, y0, y1, z0, z1)
            states, seeds = to_states(source), transport_ss.spawn(n)
            n_chunks = max(1, min(n, 8 * workers))
            cuts = np.linspace(0, n, n_chunks + 1).astype(int)
            tasks = [(states[lo:hi], seeds[lo:hi], hist_batch[lo:hi], reader, mediums,
                      settings, bounds, scorer.defs, n_batches)
                     for lo, hi in zip(cuts[:-1], cuts[1:]) if hi > lo]
            if workers > 1:
                from multiprocessing import Pool
                with Pool(workers) as pool:
                    parts = pool.map(_scalar_chunk, tasks)
            else:
                parts = [_scalar_chunk(task) for task in tasks]
            detected = np.concatenate([p["detected"] for p in parts])
            counters = {}
            for p in parts:
                scorer.merge(p["scorer"])
                for k, v in p["counters"].items():
                    counters[k] = counters.get(k, 0) + v
        else:
            from . import vector_transport as vt
            scorer.set_batches(hist_batch)
            if workers > 1:
                res = vt.run_transport_parallel(
                    self.data_dir, mediums, source, settings, n_workers=workers,
                    seed=int(transport_ss.generate_state(1)[0]),
                    temperature=run.temperature, scorer=scorer)
            else:
                res = vt.run_transport(reader, mediums, source,
                                       np.random.default_rng(transport_ss), settings,
                                       scorer=scorer)
            detected, counters = None, res["counters"]

        out = scorer.finalize(np.diff(edges))
        out["_counters"] = counters
        if detector is not None:
            frac = [detected[lo:hi].mean() for lo, hi in zip(edges[:-1], edges[1:])]
            m, sem = mean_sem(frac)
            out[detector.name] = {"type": "detector", "fraction": m, "fraction_sem": sem,
                                  "box": detector.box}
        return {name: out[name] for name in [*self.tallies, "_counters"]}

    def _solve_eigenvalue(self, reader, mediums, source, transport_ss, workers):
        from .tallies import Scorer
        run = self.run
        settings = run.settings()
        seed = int(transport_ss.generate_state(1)[0] % (2 ** 31))
        n_active = run.generations - run.inactive
        defs = self._tally_defs()
        scorer = Scorer(defs, n_active) if defs else None
        kwargs = dict(generations=run.generations, inactive=run.inactive, seed=seed,
                      scorer=scorer)
        if run.engine == "scalar":
            from .criticality import run_keff
            from .vt_calc import VelocitySampler
            states = to_states(source)
            if workers > 1:
                from multiprocessing import Pool
                with Pool(workers) as pool:
                    k = run_keff(reader, mediums, states, 1.0, 1.0, VelocitySampler(1.0),
                                 settings, pool=pool, **kwargs)
            else:
                k = run_keff(reader, mediums, states, 1.0, 1.0, VelocitySampler(1.0),
                             settings, **kwargs)
        else:
            from . import vector_transport as vt
            if workers > 1:
                logger.warning("the vector k-eigenvalue solver runs in one process; "
                               "ignoring run.workers = %d", workers)
            k = vt.run_keff_vector(reader, mediums, source, settings, **kwargs)
        keff = {"type": "keff", "k_eff": k["k_eff"], "k_eff_sem": k["k_sem"],
                "k_source": k["k_source"], "k_source_sem": k["k_source_sem"],
                "active_generations": len(k["active_k"]),
                "cycles": [{"generation": g, "k_collision": kc, "k_source": ks,
                            "active": a} for g, kc, ks, a in k["cycles"]]}
        if len(k["cycles"]) < run.generations:
            keff["note"] = (f"fission chain died out after {len(k['cycles'])} "
                            "generations")
        out = {"_counters": {}}
        done = len(k["active_k"])
        if scorer is not None and done > 0:
            # tallies are per source neutron of each active generation
            out.update(scorer.finalize(np.full(done, run.particles)))
        for t in self.tallies.values():
            if t.type == "keff":
                out[t.name] = keff
            elif t.name not in out:
                out[t.name] = {"type": t.type, "note": "no active generations completed"}
        return {name: out[name] for name in [*self.tallies, "_counters"]}

    # ----------------------------------------------------------- provenance
    def resolved_input(self):
        """The input with every default filled in, as recorded in results."""
        spec = {k: v for k, v in self.spec.items()}
        spec["run"] = {k: v for k, v in asdict(self.run).items()}
        spec["run"]["data"] = self.data_dir
        tallies = {}
        for name, t in self.tallies.items():
            entry = {"type": t.type}
            if t.cells:
                entry["cells"] = t.cells
            if t.energy_bins is not None:
                entry["energy_bins"] = t.energy_bins.tolist()
            if t.mesh:
                entry["mesh"] = t.mesh
            if t.estimator:
                entry["estimator"] = t.estimator
            if t.box:
                entry["box"] = t.box
            tallies[name] = entry
        spec["tallies"] = tallies
        return spec

    def _provenance(self, awr):
        text = self.input_text
        return {
            "input_file": self.input_path,
            "input_sha256": hashlib.sha256(text).hexdigest() if text else None,
            "resolved_input": self.resolved_input(),
            "pyneut_version": __version__,
            "git": _git_state(),
            "python": platform.python_version(),
            "numpy": np.__version__,
            "data": {"path": self.data_dir, "manifest": data_manifest(self.data_dir),
                     "atomic_weight_ratios": awr},
            "created": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        }

    # ---------------------------------------------------------------- check
    def check(self, samples=20_000, seed=0):
        """Checks that need the geometry and data, without transport.

        Returns (errors, warnings), lists of strings. Input errors raise
        InputError when the model is constructed, so they are not repeated.
        """
        errors, warnings = [], []
        try:
            mediums, _, _ = self.build()
        except InputError as e:
            return list(e.problems), warnings
        from .vector_transport import contains_many
        rng = np.random.default_rng(seed)

        # source points that start outside every cell leak immediately
        try:
            src = self.source.sample(min(samples, 5000), rng)
        except ValueError as e:
            errors.append(f"source: {e}")
        else:
            inside = np.zeros(src["x"].size, dtype=bool)
            for m in mediums:
                inside |= contains_many(m, src["x"], src["y"], src["z"])
            frac_out = 1.0 - inside.mean()
            if frac_out == 1.0:
                errors.append("source: every sampled source point is outside all cells")
            elif frac_out > 0:
                warnings.append(f"source: {100 * frac_out:.1f}% of source points are outside "
                                "all cells; those neutrons leak without transport")

        # overlapping cells must have distinct priorities
        bounds = self.bounding_box()
        if bounds is None:
            warnings.append("overlap check skipped: could not bound the geometry "
                            "(it needs finite spheres, boxes, cylinders or axis planes)")
        else:
            lo, hi = bounds
            pts = [lo[i] + (hi[i] - lo[i]) * rng.random(samples) for i in range(3)]
            masks = [contains_many(m, *pts) for m in mediums]
            names = list(self.cells)
            seen = set()
            for i in range(len(mediums)):
                for j in range(i + 1, len(mediums)):
                    both = masks[i] & masks[j]
                    if not both.any() or (names[i], names[j]) in seen:
                        continue
                    seen.add((names[i], names[j]))
                    # an unset priority is 0, as in the transport kernels
                    pi = self.cells[names[i]].priority or 0
                    pj = self.cells[names[j]].priority or 0
                    k = int(np.argmax(both))
                    at = f"({pts[0][k]:.4g}, {pts[1][k]:.4g}, {pts[2][k]:.4g})"
                    if pi == pj:
                        errors.append(f"cells {names[i]!r} and {names[j]!r} overlap, e.g. at {at}; "
                                      "make their regions disjoint or give them different "
                                      "priorities (the higher one wins; unset means 0)")
        return errors, warnings

    def bounding_box(self):
        """Axis-aligned box around every finite surface the cells use, or None."""
        lo = [np.inf] * 3
        hi = [-np.inf] * 3

        def grow(axis, a, b):
            lo[axis] = min(lo[axis], a)
            hi[axis] = max(hi[axis], b)

        used = set().union(*(surface_names(c.expr) for c in self.cells.values()))
        for name in used:
            s = self.surfaces[name]
            if isinstance(s, Sphere):
                for i, c in enumerate((s.x0, s.y0, s.z0)):
                    grow(i, c - s.radius, c + s.radius)
            elif isinstance(s, Box):
                for p in s.surfaces:
                    for i, comp in enumerate((p.A, p.B, p.C)):
                        if abs(comp) == 1.0:
                            grow(i, -p.D * comp, -p.D * comp)
            elif isinstance(s, Cylinder):
                axis = "xyz".index(s.axis)
                for i, c in enumerate((s.x0, s.y0, s.z0)):
                    if i != axis:
                        grow(i, c - s.radius, c + s.radius)
            elif isinstance(s, Plane):
                for i, comp in enumerate((s.A, s.B, s.C)):
                    if abs(comp) == 1.0:
                        grow(i, -s.D * comp, -s.D * comp)
        if not all(np.isfinite(lo) & np.isfinite(hi)) or any(h <= l for l, h in zip(lo, hi)):
            return None
        pad = [0.01 * (h - l) for l, h in zip(lo, hi)]
        return ([l - p for l, p in zip(lo, pad)], [h + p for h, p in zip(hi, pad)])

    def describe(self):
        """A short human-readable summary of the model."""
        r = self.run
        lines = [f"mode {r.mode}, engine {r.engine}, {r.particles} particles"
                 + (f" x {r.generations} generations ({r.inactive} inactive)"
                    if r.mode == "eigenvalue" else f", {r.batches} batches")
                 + f", seed {r.seed}",
                 f"data: {self.data_dir}"]
        for name, m in self.materials.items():
            lines.append(f"material {name}: " + ", ".join(m.nuclide_names()))
        for name, c in self.cells.items():
            pri = f" (priority {c.priority})" if c.priority is not None else ""
            lines.append(f"cell {name}: {c.region_text}  ->  {c.material}{pri}")
        lines.append("tallies: " + ", ".join(f"{n} ({t.type})" for n, t in self.tallies.items()))
        return "\n".join(lines)


# ------------------------------------------------------- scalar engine work
def _scalar_chunk(task):
    """Run a contiguous slice of histories on the scalar engine with its own
    Scorer; returns the scorer totals, detector flags and event counters."""
    from .random_number_generator import RNGHandler
    from .simulation import simulate_single_particle
    from .tallies import Scorer
    from .vt_calc import VelocitySampler
    states, seeds, batches, reader, mediums, settings, bounds, defs, n_batches = task
    scorer = Scorer(defs, n_batches)
    placeholder = VelocitySampler(1.0)
    detected = np.zeros(len(states), dtype=bool)
    counters = {}
    for i, (state, seed) in enumerate(zip(states, seeds)):
        scorer.batch = int(batches[i])
        r = simulate_single_particle((state, reader, mediums, 1.0, 1.0, placeholder,
                                      bounds, False, RNGHandler(seed), settings, scorer))
        detected[i] = r["region_detected"]
        for k, v in r["physics_counts"].items():
            counters[k] = counters.get(k, 0) + v
    return {"scorer": scorer.state(), "detected": detected, "counters": counters}


def _git_state():
    here = os.path.dirname(os.path.abspath(__file__))
    try:
        commit = subprocess.run(["git", "-C", here, "rev-parse", "HEAD"], capture_output=True,
                                text=True, timeout=5)
        if commit.returncode != 0:
            return None
        dirty = subprocess.run(["git", "-C", here, "status", "--porcelain", "--", here],
                               capture_output=True, text=True, timeout=5)
        return {"commit": commit.stdout.strip(), "dirty": bool(dirty.stdout.strip())}
    except (OSError, subprocess.SubprocessError):
        return None


# ----------------------------------------------------------------- results
class Results:
    """Tally results, run information and provenance of one solve."""

    def __init__(self, tallies, info, provenance):
        self.tallies = tallies
        self.info = info
        self.provenance = provenance

    def __getitem__(self, name):
        return self.tallies[name]

    def __iter__(self):
        return iter(self.tallies)

    def to_dict(self):
        """JSON-ready results. Mesh arrays are left to the files save() writes."""
        def clean(v):
            if isinstance(v, dict):
                return {k: clean(x) for k, x in v.items()}
            if isinstance(v, (list, tuple)):
                return [clean(x) for x in v]
            if isinstance(v, np.ndarray):
                return clean(v.tolist())
            if isinstance(v, (np.floating, float)):
                return None if not np.isfinite(v) else float(v)
            if isinstance(v, np.integer):
                return int(v)
            if isinstance(v, np.bool_):
                return bool(v)
            return v
        tallies = {}
        for name, t in self.tallies.items():
            t = dict(t)
            if "x" in t.get("axes", ()):
                t.pop("mean"), t.pop("sem")
                t["files"] = [f"{name}.vtk", f"{name}.npz"]
            tallies[name] = t
        return clean({"tallies": tallies, "run": self.info, "provenance": self.provenance})

    def save(self, directory):
        """Write results.json, plus <tally>.vtk and <tally>.npz per mesh tally.

        The NPZ holds the full mean and sem arrays (axes as in the tally's
        `axes`); the VTK holds the mean summed over any cell and energy axes.
        """
        os.makedirs(directory, exist_ok=True)
        for name, t in self.tallies.items():
            if "x" not in t.get("axes", ()):
                continue
            from .mesh import MeshTally
            m = t["mesh"]
            axes = t["axes"]
            grid = t["mean"]
            for ax in ("energy", "cell"):
                if ax in axes:
                    grid = grid.sum(axis=axes.index(ax))
            mt = MeshTally(*zip(m["min"], m["max"]), tuple(m["dims"]))
            mt.data = grid
            mt.write_vtk(os.path.join(directory, f"{name}.vtk"))
            extra = {"energy_bins": t["energy_bins"]} if "energy_bins" in t else {}
            np.savez_compressed(os.path.join(directory, f"{name}.npz"), mean=t["mean"],
                                sem=t["sem"], axes=np.array(axes), min=m["min"],
                                max=m["max"], dims=m["dims"], **extra)
        path = os.path.join(directory, "results.json")
        with open(path, "w") as f:
            json.dump(self.to_dict(), f, indent=2)
        return path

    def summary(self):
        lines = []
        i = self.info
        lines.append(f"{i['mode']} on the {i['engine']} engine: {i['particles']} particles"
                     + (f" x {i['generations']} generations" if i["mode"] == "eigenvalue" else "")
                     + f", {i['wall_time_s']:.1f} s")

        def pm(v, s, fmt):
            if v is None:
                return "n/a"
            return f"{v:{fmt}}" + (f" ± {s:{fmt}}" if s is not None and np.isfinite(s) else "")

        def per_cell(t, fmt):
            if t.get("axes") != ["cell"]:
                return ""
            return "; " + ", ".join(f"{c} {pm(m, e, fmt)}"
                                    for c, m, e in zip(t["cells"], t["mean"], t["sem"]))

        def binned(t):
            axes = [a for a in t.get("axes", []) if a != "cell" or t["axes"] != ["cell"]]
            return f" ({' × '.join(axes)} bins on save)" if axes else ""

        for name, t in self.tallies.items():
            kind = t["type"]
            if "note" in t and kind != "keff":
                lines.append(f"  {name}: {t['note']}")
            elif kind == "leakage":
                lines.append(f"  {name}: leakage {pm(t['total'], t['total_sem'], '.5f')}, "
                             f"mean escape energy {pm(t['mean_energy'], t['mean_energy_sem'], '.4e')} eV"
                             + binned(t))
            elif kind == "absorption":
                lines.append(f"  {name}: absorbed {pm(t['total'], t['total_sem'], '.5f')} per source neutron"
                             + per_cell(t, ".5f") + binned(t))
            elif kind == "flux":
                lines.append(f"  {name}: flux ({t['estimator']}) {pm(t['total'], t['total_sem'], '.4g')} "
                             "cm per source neutron" + per_cell(t, ".4g") + binned(t))
            elif kind == "detector":
                lines.append(f"  {name}: reached detector {pm(t['fraction'], t['fraction_sem'], '.5f')}")
            elif kind == "keff":
                lines.append(f"  {name}: k_eff {pm(t['k_eff'], t['k_eff_sem'], '.5f')} "
                             f"(generation estimator {pm(t['k_source'], t['k_source_sem'], '.5f')})"
                             + (f"; {t['note']}" if "note" in t else ""))
        return "\n".join(lines)
