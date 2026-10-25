"""TOML input format: parser, model building, every example, and equivalence
with the hand-built low-level API.

Runs without ENDF data: a synthetic library (tests/synthetic_data.py) stands
in, so these tests check the code paths and bookkeeping, not physics values.
"""
import glob
import hashlib
import json
import os

import numpy as np
import pytest

import synthetic_data
from pyneut.cli import main as cli_main
from pyneut.data import data_manifest
from pyneut.errors import InputError
from pyneut.material import NEUTRON_MASS_AMU
from pyneut.medium import Cylinder, Plane, Region, Sphere
from pyneut.model import Model
from pyneut.regions import compile_region, parse_region
from pyneut.sources import build_energy, build_space, to_states
from pyneut.units import parse_energy

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXAMPLES = sorted(glob.glob(os.path.join(ROOT, "examples", "**", "*.toml"), recursive=True))


@pytest.fixture(scope="session")
def data_dir(tmp_path_factory):
    return synthetic_data.write_library(str(tmp_path_factory.mktemp("synthetic_endf")))


@pytest.fixture(autouse=True)
def _use_synthetic_data(data_dir, monkeypatch):
    monkeypatch.setenv("PYNEUT_DATA", data_dir)


PB_SPHERE = """
[run]
mode = "fixed-source"
engine = "{engine}"
particles = 400
batches = 10
seed = 99

[materials.lead]
density = 11.35
nuclides = {{ Pb208 = 1.0 }}

[surfaces]
shell = {{ type = "sphere", center = [0, 0, 0], radius = 10.0 }}

[cells]
shield = {{ region = "-shell", material = "lead" }}

[source]
space = {{ type = "point", xyz = [0, 0, 0] }}
direction = "isotropic"
energy = "2 MeV"

[tallies.leak]
type = "leakage"
energy_bins = {{ log = ["1 keV", "2.1 MeV", 20] }}
"""


def write(tmp_path, text, name="model.toml"):
    path = tmp_path / name
    path.write_text(text)
    return str(path)


# ------------------------------------------------------------------ parsing
@pytest.mark.parametrize("text, inside, outside", [
    ("-s", [(0, 0, 0)], [(20, 0, 0)]),
    ("+s", [(20, 0, 0)], [(0, 0, 0)]),
    ("-s & -p", [(-1, 0, 0)], [(1, 0, 0), (20, 0, 0)]),
    ("-s | -p", [(-20, 0, 0), (1, 0, 0)], [(20, 0, 0)]),
    ("~(-s & -p)", [(1, 0, 0), (20, 0, 0)], [(-1, 0, 0)]),
    ("-s & ~(-c)", [(0, 5, 0)], [(0, 0, 0), (0, 20, 0)]),
    ("-s & -p | -c", [(-1, 0, 0), (5, 0, 0)], [(5, 5, 0)]),   # & binds tighter
])
def test_region_expressions(text, inside, outside):
    surfaces = {"s": Sphere((0, 0, 0), 10), "p": Plane(1, 0, 0, 0),
                "c": Cylinder("x", 1, (0, 0, 0))}
    region = Region(surfaces=[compile_region(parse_region(text), surfaces)])
    assert all(region.contains(*pt) for pt in inside)
    assert not any(region.contains(*pt) for pt in outside)


@pytest.mark.parametrize("text, message", [
    ("-a -b", "missing operator"),
    ("shell", "needs a sign"),
    ("-a &", "ends where a term"),
    ("(-a", r"missing '\)'"),
    ("-a)", "unbalanced"),
    ("", "empty"),
])
def test_region_errors(text, message):
    with pytest.raises(ValueError, match=message):
        parse_region(text)


def test_energy_units():
    assert parse_energy(2e6) == parse_energy("2 MeV") == parse_energy("2000keV") == 2e6
    assert parse_energy("0.0253 eV") == 0.0253
    for bad in ("2 Mev", "2 meV", "fast", True):
        with pytest.raises(ValueError):
            parse_energy(bad)


def test_source_distributions():
    rng = np.random.default_rng(0)
    x, y, z = build_space({"type": "sphere", "center": [0, 0, 0], "radius": 2}).sample(50_000, rng)
    assert np.sqrt(x * x + y * y + z * z).mean() == pytest.approx(1.5, rel=0.01)
    e = build_energy({"type": "watt", "nuclide": "U235"}).sample(50_000, rng)
    assert e.mean() == pytest.approx(2.03e6, rel=0.02)
    e = build_energy({"type": "discrete", "energies": ["1 MeV", "2 MeV"],
                      "probabilities": [3, 1]}).sample(50_000, rng)
    assert (e == 1e6).mean() == pytest.approx(0.75, abs=0.01)


def test_input_errors_are_all_reported(tmp_path):
    bad = PB_SPHERE.format(engine="vector").replace('engine = "vector"', 'engin = "vector"') \
        .replace('material = "lead"', 'material = "leed"') \
        .replace('region = "-shell"', 'region = "-shell -shell"') \
        .replace('energy = "2 MeV"', 'energy = "2 Mev"')
    with pytest.raises(InputError) as err:
        Model.from_toml(write(tmp_path, bad))
    text = str(err.value)
    for expected in ("run.engin: unknown key (did you mean 'engine'?)",
                     "unknown material 'leed' (did you mean 'lead'?)",
                     "missing operator",
                     "cannot read energy '2 Mev'"):
        assert expected in text
    assert len(err.value.problems) == 4


def test_element_names_are_rejected(tmp_path):
    with pytest.raises(InputError, match="natural elements are not expanded"):
        Model.from_toml(write(tmp_path, PB_SPHERE.format(engine="vector")
                              .replace("Pb208", "Pb")))


def test_missing_data_is_an_error(tmp_path):
    model = Model.from_toml(write(tmp_path, PB_SPHERE.format(engine="vector")
                                  .replace("Pb208", "Pb206")))
    with pytest.raises(InputError, match="no data file .*Pb206.h5"):
        model.build()


def test_vector_engine_rejects_mesh_tallies(tmp_path):
    text = PB_SPHERE.format(engine="vector") + \
        '\n[tallies.map]\ntype = "absorption"\nmesh = { min = [-1,-1,-1], max = [1,1,1], dims = [2,2,2] }\n'
    with pytest.raises(InputError, match="need run.engine = \"scalar\""):
        Model.from_toml(write(tmp_path, text))


# ----------------------------------------------------------------- building
def test_awr_and_number_density_come_from_the_data(tmp_path, data_dir):
    model = Model.from_toml(write(tmp_path, PB_SPHERE.format(engine="vector")))
    _, compositions, awr = model.build()
    expected_awr = synthetic_data.MASSES["Pb208"] / NEUTRON_MASS_AMU
    assert awr["Pb208"] == pytest.approx(expected_awr)
    (pb,) = compositions["lead"]
    assert pb.atomic_weight_ratio == pytest.approx(expected_awr)
    assert pb.number_density == pytest.approx(11.35 * 6.022e23 / (expected_awr * NEUTRON_MASS_AMU))


def test_check_finds_overlaps_and_outside_sources(tmp_path):
    text = PB_SPHERE.format(engine="vector").replace(
        'shield = { region = "-shell", material = "lead" }',
        'shield = { region = "-shell", material = "lead" }\ncore = { region = "-shell", material = "void" }')
    errors, _ = Model.from_toml(write(tmp_path, text)).check(samples=2000)
    assert any("overlap" in e for e in errors)

    with_priority = text.replace('material = "void" }', 'material = "void", priority = 1 }')
    errors, _ = Model.from_toml(write(tmp_path, with_priority, "p.toml")).check(samples=2000)
    assert not errors

    outside = PB_SPHERE.format(engine="vector").replace("xyz = [0, 0, 0]", "xyz = [50, 0, 0]")
    errors, _ = Model.from_toml(write(tmp_path, outside, "o.toml")).check(samples=2000)
    assert any("outside all cells" in e for e in errors)


# ------------------------------------------------------------------ examples
@pytest.mark.parametrize("path", EXAMPLES, ids=lambda p: os.path.relpath(p, ROOT))
def test_example_checks_and_runs(path):
    model = Model.from_toml(path)
    errors, _ = model.check(samples=4000)
    assert errors == []
    model.run.workers = 1
    if os.path.basename(path) == "pincell_beavrs.toml":
        # a reflective lattice on flat synthetic cross sections barely absorbs,
        # so histories are very long: keep this one minimal
        model.run.particles, model.run.generations, model.run.inactive = 8, 2, 1
    elif model.run.mode == "eigenvalue":
        model.run.particles, model.run.generations, model.run.inactive = 60, 3, 1
    else:
        model.run.particles = 60
    results = model.solve()
    assert set(results) == set(model.tallies)
    json.dumps(results.to_dict())                     # serialisable


def test_examples_cover_every_feature():
    text = "\n".join(open(p).read() for p in EXAMPLES)
    for feature in ('"eigenvalue"', '"fixed-source"', '"scalar"', '"analog"',
                    '"reflective"', '"weight"', "atom_densities", "priority",
                    '"plane"', '"z-plane"', '"sphere"',
                    '"cylinder"', '"box"', '"point"', '"cell"', '"beam"',
                    '"uniform"', '"discrete"', '"maxwell"', '"watt"',
                    '"leakage"', '"absorption"', '"detector"', "mesh", "|", "~"):
        assert feature in text, feature


# --------------------------------------------------------------- equivalence
def _hand_built_pb(data_dir):
    from pyneut.cross_section_read import CrossSectionReader
    from pyneut.material import Material
    awr = synthetic_data.MASSES["Pb208"] / NEUTRON_MASS_AMU
    pb = Material.mixture("lead", 11.35, [("Pb208", awr, 1.0)])
    regions = [Region(surfaces=[Sphere((0, 0, 0), 10.0)], name="shield", composition=pb)]
    return CrossSectionReader(data_dir), regions


def _seeded_source(n, seed):
    """The model's documented seeding: SeedSequence(seed) -> (source, transport)."""
    source_ss, transport_ss = np.random.SeedSequence(seed).spawn(2)
    rng = np.random.default_rng(source_ss)
    mu = 1.0 - 2.0 * rng.random(n)
    phi = 2.0 * np.pi * rng.random(n)
    s = np.sqrt(1.0 - mu * mu)
    src = {"x": np.zeros(n), "y": np.zeros(n), "z": np.zeros(n),
           "u": s * np.cos(phi), "v": s * np.sin(phi), "w": mu, "energy": np.full(n, 2e6)}
    return src, transport_ss


def test_scalar_engine_matches_hand_built_api(tmp_path, data_dir):
    from pyneut.random_number_generator import RNGHandler
    from pyneut.settings import Settings
    from pyneut.simulation import simulate_single_particle
    from pyneut.statistics import escape_statistics
    from pyneut.vt_calc import VelocitySampler

    results = Model.from_toml(write(tmp_path, PB_SPHERE.format(engine="scalar"))).solve()

    reader, regions = _hand_built_pb(data_dir)
    src, transport_ss = _seeded_source(400, 99)
    settings = Settings("shielding", 400, 10)
    per = [simulate_single_particle((st, reader, regions, 1.0, 1.0, VelocitySampler(1.0), None,
                                     False, RNGHandler(child), settings))
           for st, child in zip(to_states(src), transport_ss.spawn(400))]
    stats = escape_statistics(np.array([r["final_weight"] for r in per]),
                              np.array([r["final_energy"] for r in per]), 400, 10)
    assert results["leak"]["total"] == stats["leakage"]
    assert results["leak"]["mean_energy"] == stats["avg_energy"]


def test_vector_engine_matches_hand_built_api(tmp_path, data_dir):
    from pyneut import vector_transport as vt
    from pyneut.settings import Settings
    from pyneut.statistics import escape_statistics

    results = Model.from_toml(write(tmp_path, PB_SPHERE.format(engine="vector"))).solve()

    reader, regions = _hand_built_pb(data_dir)
    src, transport_ss = _seeded_source(400, 99)
    res = vt.run_transport(reader, regions, src, np.random.default_rng(transport_ss),
                           Settings("shielding", 400, 10))
    stats = escape_statistics(res["leaked_weight"], res["escape_energy"], 400, 10)
    assert results["leak"]["total"] == stats["leakage"]
    assert results["leak"]["mean_energy"] == stats["avg_energy"]


# ------------------------------------------------------------ results / CLI
def test_results_record_provenance(tmp_path, data_dir):
    path = write(tmp_path, PB_SPHERE.format(engine="vector"))
    results = Model.from_toml(path).solve()
    saved = json.load(open(results.save(str(tmp_path / "out"))))
    prov = saved["provenance"]
    assert prov["input_sha256"] == hashlib.sha256(open(path, "rb").read()).hexdigest()
    assert prov["data"]["manifest"] == data_manifest(data_dir)
    assert prov["resolved_input"]["run"]["batches"] == 10
    assert prov["pyneut_version"]
    leak = saved["tallies"]["leak"]
    assert len(leak["spectrum"]) == 20 and len(leak["energy_bins"]) == 21


def test_cli_check_and_run(tmp_path, data_dir, capsys):
    path = write(tmp_path, PB_SPHERE.format(engine="vector"))
    assert cli_main(["check", path]) == 0
    out = tmp_path / "res"
    assert cli_main(["run", path, "-o", str(out), "--particles", "100", "--engine", "scalar"]) == 0
    saved = json.load(open(out / "results.json"))
    assert saved["run"]["particles"] == 100 and saved["run"]["engine"] == "scalar"
    bad = write(tmp_path, "[run]\nmode = 'nope'\n", "bad.toml")
    assert cli_main(["check", bad]) == 1
    assert "missing [materials] section" in capsys.readouterr().err
