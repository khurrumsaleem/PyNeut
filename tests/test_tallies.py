"""The event-based tally interface (pyneut.tallies) and its engine hooks.

Runs on the synthetic data library: flat, energy-independent cross sections
make several answers exact or analytic, which is what these tests use.
"""
import numpy as np
import pytest

import synthetic_data
from pyneut import vector_transport as vt
from pyneut.cross_section_read import CrossSectionReader
from pyneut.material import Nuclide
from pyneut.medium import Box, Region, Sphere
from pyneut.model import Model
from pyneut.random_number_generator import RNGHandler
from pyneut.settings import Settings
from pyneut.simulation import simulate_single_particle
from pyneut.sources import to_states
from pyneut.statistics import batch_edges
from pyneut.tallies import Scorer, TallyDef
from pyneut.vt_calc import VelocitySampler

N_B10 = 1e22          # atoms/cm3; synthetic B10: sigma_el 4 b, sigma_capture 3 b


@pytest.fixture(scope="module")
def reader(tmp_path_factory):
    path = synthetic_data.write_library(str(tmp_path_factory.mktemp("endf")),
                                        ["B10", "Pb208", "U235"])
    return CrossSectionReader(path)


def point_source(n, seed, energy=2e6):
    rng = np.random.default_rng(seed)
    mu = 1.0 - 2.0 * rng.random(n)
    phi = 2.0 * np.pi * rng.random(n)
    s = np.sqrt(1.0 - mu * mu)
    return {"x": np.zeros(n), "y": np.zeros(n), "z": np.zeros(n),
            "u": s * np.cos(phi), "v": s * np.sin(phi), "w": mu,
            "energy": np.full(n, energy)}


def two_shell():
    inner = [Nuclide("B10", N_B10, 9.93)]
    outer = [Nuclide("Pb208", 3.3e22, 206.19)]
    return [Region(surfaces=[Sphere((0, 0, 0), 5)], name="inner", priority=2, composition=inner),
            Region(surfaces=[Sphere((0, 0, 0), 10)], name="outer", priority=1, composition=outer)]


def run_engine(engine, reader, mediums, defs, n=2000, batches=10, seed=3, settings=None):
    """Run `n` histories on one engine with a Scorer; return (scorer, engine totals)."""
    settings = settings or Settings("shielding", n, batches)
    edges = batch_edges(n, batches)
    hist_batch = np.searchsorted(edges, np.arange(n), side="right") - 1
    scorer = Scorer(defs, batches)
    src = point_source(n, seed)
    if engine == "vector":
        scorer.set_batches(hist_batch)
        res = vt.run_transport(reader, mediums, src, np.random.default_rng(seed + 1),
                               settings, scorer=scorer)
        leaked, absorbed = res["leaked_weight"].sum(), res["absorbed_weight"]
    else:
        leaked = absorbed = 0.0
        seeds = np.random.SeedSequence(seed + 1).spawn(n)
        for i, (state, s) in enumerate(zip(to_states(src), seeds)):
            scorer.batch = int(hist_batch[i])
            r = simulate_single_particle((state, reader, mediums, 1.0, 1.0, VelocitySampler(1.0),
                                          None, False, RNGHandler(s), settings, scorer))
            leaked += r["final_weight"]
            absorbed += r["absorbed_weight"]
    return scorer.finalize(np.diff(edges)), leaked / n, absorbed / n


# -------------------------------------------------------------- scorer unit
def test_leakage_is_scored_per_neutron():
    # one history leaking a 1 MeV and a 9 MeV neutron (as after (n,2n)) must
    # fill the 1 and 9 MeV bins; a per-history average would put it all at 5
    bins = np.array([0.0, 2e6, 4e6, 6e6, 8e6, 10e6])
    scorer = Scorer([TallyDef("leak", "leakage", energy_bins=bins)], 1)
    scorer.set_batches(np.array([0]))
    scorer.leak(np.array([0, 0]), [0, 0], [0, 0], [0, 0], np.array([1e6, 9e6]), np.array([1.0, 1.0]))
    out = scorer.finalize([1])["leak"]
    assert list(out["mean"]) == [1.0, 0.0, 0.0, 0.0, 1.0]
    assert out["total"] == 2.0 and out["mean_energy"] == 5e6


def test_filters_bin_events():
    defs = [TallyDef("a", "absorption", cell_map=[1, -1, 0], cell_names=["c2", "c0"],
                     mesh={"min": [0, 0, 0], "max": [2, 2, 2], "dims": [2, 1, 1]},
                     energy_bins=[0.0, 1.0, 10.0])]
    scorer = Scorer(defs, 2)
    scorer.batch = 1
    scorer.absorb(None, 0, 1.5, 0.5, 0.5, 5.0, 2.0)     # cell 0 -> bin 1, x bin 1, E bin 1
    scorer.absorb(None, 1, 0.5, 0.5, 0.5, 5.0, 7.0)     # cell 1 is not tallied
    scorer.absorb(None, 2, 9.0, 0.5, 0.5, 0.5, 3.0)     # outside the mesh
    v = scorer.values["a"]
    assert v.shape == (2, 2, 2, 1, 1, 2)
    assert v[1, 1, 1, 0, 0, 1] == 2.0 and v.sum() == 2.0
    assert scorer.totals["a"][1] == 12.0              # totals ignore the filters


def test_track_length_mesh_is_refused():
    with pytest.raises(ValueError, match="collision estimator"):
        TallyDef("f", "flux", mesh={"min": [0, 0, 0], "max": [1, 1, 1], "dims": [1, 1, 1]},
                 estimator="track-length")
    assert TallyDef("f", "flux", mesh={"min": [0, 0, 0], "max": [1, 1, 1],
                                       "dims": [1, 1, 1]}).estimator == "collision"


# ------------------------------------------------------------- engine hooks
@pytest.mark.parametrize("engine", ["vector", "scalar"])
def test_scorer_does_not_change_transport(engine, reader):
    mediums = two_shell()
    plain = run_engine(engine, reader, mediums, [], n=300)
    scored = run_engine(engine, reader, mediums,
                        [TallyDef("f", "flux"), TallyDef("l", "leakage")], n=300)
    assert plain[1] == scored[1] and plain[2] == scored[2]


@pytest.mark.parametrize("engine", ["vector", "scalar"])
def test_totals_match_engine_bookkeeping(engine, reader):
    defs = [TallyDef("leak", "leakage", energy_bins=np.logspace(-5, 7.4, 40)),
            TallyDef("abs", "absorption", cell_map=[0, 1], cell_names=["inner", "outer"])]
    out, leaked, absorbed = run_engine(engine, reader, two_shell(), defs)
    assert out["leak"]["total"] == pytest.approx(leaked, rel=1e-10)
    assert out["abs"]["total"] == pytest.approx(absorbed, rel=1e-10)
    assert out["leak"]["total"] + out["abs"]["total"] == pytest.approx(1.0, abs=1e-12)
    # filters that cover everything partition the total
    assert out["abs"]["mean"].sum() == pytest.approx(out["abs"]["total"], rel=1e-10)
    assert out["leak"]["mean"].sum() == pytest.approx(out["leak"]["total"], rel=1e-10)


@pytest.mark.parametrize("engine", ["vector", "scalar"])
def test_flux_estimators_agree(engine, reader):
    defs = [TallyDef("tl", "flux", cell_map=[0, 1], cell_names=["inner", "outer"]),
            TallyDef("col", "flux", cell_map=[0, 1], cell_names=["inner", "outer"],
                     estimator="collision"),
            TallyDef("mesh", "flux", mesh={"min": [-10, -10, -10], "max": [10, 10, 10],
                                           "dims": [3, 3, 3]})]
    out, _, _ = run_engine(engine, reader, two_shell(), defs, n=4000)
    tl, col = out["tl"], out["col"]
    z = np.abs(tl["mean"] - col["mean"]) / np.hypot(tl["sem"], col["sem"])
    assert np.all(z < 4), z
    # the mesh covers the whole model, so its collision flux is the model total
    assert out["mesh"]["mean"].sum() == pytest.approx(col["total"], rel=1e-10)


def test_infinite_medium_flux_is_one_over_sigma_a(reader):
    # in an infinite medium every source neutron is absorbed, so the flux
    # integrated over the medium is exactly 1 / Sigma_a per source neutron
    b10 = [Nuclide("B10", N_B10, 9.93)]
    cube = [Region(surfaces=[Box(-20, 20, -20, 20, -20, 20, boundary_type="reflective")],
                   name="cube", composition=b10)]
    sigma_a = N_B10 * 3e-24
    for engine in ("vector", "scalar"):
        out, leaked, absorbed = run_engine(
            engine, reader, cube,
            [TallyDef("tl", "flux"), TallyDef("col", "flux", estimator="collision")],
            n=2000, settings=Settings("criticality", 2000, 10))      # analog
        assert leaked == 0.0 and absorbed == pytest.approx(1.0)
        for name in ("tl", "col"):
            t = out[name]
            assert abs(t["total"] - 1.0 / sigma_a) < 4 * t["total_sem"], (engine, name, t)


def test_engines_agree(reader):
    defs = [TallyDef("tl", "flux", cell_map=[0, 1], cell_names=["inner", "outer"]),
            TallyDef("leak", "leakage")]
    a, _, _ = run_engine("vector", reader, two_shell(), defs, n=4000, seed=11)
    b, _, _ = run_engine("scalar", reader, two_shell(), defs, n=4000, seed=11)
    z = np.abs(a["tl"]["mean"] - b["tl"]["mean"]) / np.hypot(a["tl"]["sem"], b["tl"]["sem"])
    assert np.all(z < 4), z
    zl = abs(a["leak"]["total"] - b["leak"]["total"]) / np.hypot(a["leak"]["total_sem"],
                                                                 b["leak"]["total_sem"])
    assert zl < 4


def test_parallel_vector_merges_scores(reader):
    n, batches = 2000, 10
    edges = batch_edges(n, batches)
    scorer = Scorer([TallyDef("leak", "leakage"), TallyDef("abs", "absorption")], batches)
    scorer.set_batches(np.searchsorted(edges, np.arange(n), side="right") - 1)
    res = vt.run_transport_parallel(reader.base_path, two_shell(), point_source(n, 5),
                                    Settings("shielding", n, batches), n_workers=2, seed=9,
                                    scorer=scorer)
    out = scorer.finalize(np.diff(edges))
    assert out["leak"]["total"] == pytest.approx(res["leaked_weight"].sum() / n, rel=1e-10)
    assert out["abs"]["total"] == pytest.approx(res["absorbed_weight"] / n, rel=1e-10)


# ------------------------------------------------------------------- models
GODIVA_LIKE = """
[run]
mode = "eigenvalue"
engine = "{engine}"
particles = 150
generations = 5
inactive = 2
seed = 4

[materials.fuel]
density = 18.74
nuclides = {{ U235 = 1.0 }}

[materials.lead]
density = 11.35
nuclides = {{ Pb208 = 1.0 }}

[surfaces]
core = {{ type = "sphere", center = [0, 0, 0], radius = 8 }}
refl = {{ type = "sphere", center = [0, 0, 0], radius = 14 }}

[cells]
core = {{ region = "-core", material = "fuel" }}
reflector = {{ region = "+core & -refl", material = "lead" }}

[source]
space = {{ type = "sphere", center = [0, 0, 0], radius = 8 }}
energy = {{ type = "watt", nuclide = "U235" }}

[tallies.flux]
type = "flux"
cells = ["core", "reflector"]
energy_bins = {{ log = ["1 keV", "20 MeV", 10] }}

[tallies.leakage]
type = "leakage"

[tallies.map]
type = "absorption"
mesh = {{ min = [-14, -14, -14], max = [14, 14, 14], dims = [4, 4, 4] }}
"""


@pytest.mark.parametrize("engine", ["vector", "scalar"])
def test_eigenvalue_tallies(engine, reader, tmp_path, monkeypatch):
    monkeypatch.setenv("PYNEUT_DATA", reader.base_path)
    path = tmp_path / "model.toml"
    path.write_text(GODIVA_LIKE.format(engine=engine))
    results = Model.from_toml(str(path)).solve()
    flux = results["flux"]
    assert flux["axes"] == ["cell", "energy"] and flux["mean"].shape == (2, 10)
    assert flux["total"] > 0 and np.all(np.isfinite(flux["sem"]))
    assert results["keff"]["active_generations"] == 3
    assert results["map"]["mean"].shape == (4, 4, 4)
    saved = results.save(str(tmp_path / "out"))
    assert (tmp_path / "out" / "map.vtk").exists() and saved.endswith("results.json")
