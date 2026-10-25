"""Reaction classification: every non-redundant reaction counts in the total.

The reader sorts a nuclide's reactions by their products: one with a neutron
among its products is transported as an emitting channel (with its neutron
yield), one without is absorption, and a redundant one (a sum of others) is
skipped. A small hand-made file checks the bookkeeping, both engines check
that a channel emitting no neutron removes the particle, and the real data,
when present, checks the totals against OpenMC's.
"""
import os
import pickle

import h5py
import numpy as np
import pytest

import synthetic_data
from pyneut import vector_transport as vt
from pyneut.cross_section_read import CrossSectionReader
from pyneut.material import Nuclide
from pyneut.medium import Region, Sphere
from pyneut.random_number_generator import RNGHandler
from pyneut.settings import Settings
from pyneut.simulation import _neutron_multiplicity, simulate_single_particle
from pyneut.sources import to_states
from pyneut.vt_calc import VelocitySampler

ENDF = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "endfb")


def write_nuclide(directory, name, reactions):
    """Write <directory>/neutron/<name>.h5 with flat cross sections.

    `reactions` maps MT -> (xs in barns, products, redundant), where products
    is a list of (particle, yield) and a yield is ("poly", coeffs) or
    ("tab", energies, values).
    """
    out = os.path.join(directory, "neutron")
    os.makedirs(out, exist_ok=True)
    energy = np.logspace(-5, np.log10(2e7), 50)
    with h5py.File(os.path.join(out, f"{name}.h5"), "w") as f:
        f.attrs["filetype"] = np.bytes_("data_neutron")
        root = f.create_group(name)
        root.attrs["atomic_weight_ratio"] = (synthetic_data.MASSES[name]
                                             / synthetic_data.NEUTRON_MASS_AMU)
        root.create_dataset("energy/294K", data=energy)
        for mt, (xs, products, redundant) in reactions.items():
            grp = root.create_group(f"reactions/reaction_{mt:03}")
            grp.attrs["mt"] = mt
            grp.attrs["redundant"] = int(redundant)
            d = grp.create_dataset("294K/xs", data=np.full(energy.size, xs))
            d.attrs["threshold_idx"] = 0
            for i, (particle, yld) in enumerate(products):
                pr = grp.create_group(f"product_{i}")
                pr.attrs["particle"] = np.bytes_(particle)
                if yld[0] == "poly":
                    ds = pr.create_dataset("yield", data=np.asarray(yld[1], float))
                    ds.attrs["type"] = np.bytes_("Polynomial")
                else:
                    ds = pr.create_dataset("yield", data=np.vstack([yld[1], yld[2]]))
                    ds.attrs["type"] = np.bytes_("Tabulated1D")
    return directory


N = [("neutron", ("poly", [1.0]))]


@pytest.fixture(scope="module")
def reader(tmp_path_factory):
    path = str(tmp_path_factory.mktemp("reactions"))
    write_nuclide(path, "Al27", {
        2: (1.0, N, False),
        16: (0.125, [("neutron", ("poly", [2.0]))], False),
        5: (2.0, [("neutron", ("tab", [1e6, 1e7], [0.5, 1.5])), ("alpha", ("poly", [1.0]))], False),
        22: (0.25, N + [("alpha", ("poly", [1.0]))], False),
        103: (100.0, [("proton", ("poly", [1.0]))], True),     # redundant: skipped
        600: (0.5, [("proton", ("poly", [1.0]))], False),
        107: (0.7, [], False),
        102: (0.03, [("photon", ("poly", [1.0]))], False),
    })
    # every collision that is not elastic emits no neutron at all
    write_nuclide(path, "O16", {
        2: (1.0, N, False),
        5: (3.0, [("neutron", ("poly", [0.0]))], False),
    })
    synthetic_data.write_library(path, ["Pb208"])        # no product data
    return CrossSectionReader(path)


def test_reactions_are_classified_by_their_products(reader):
    reader._scan_inelastic_mts("Al27")
    assert reader._available_inelastic_mts["Al27"] == [16, 5, 22]
    assert reader._absorption_mts["Al27"] == [102, 107, 600]


def test_total_is_the_sum_of_non_redundant_reactions(reader):
    reader._build_fast_table("Al27")
    t = reader._fast_tables["Al27"]
    np.testing.assert_allclose(t["el"], 1.0)
    np.testing.assert_allclose(t["in"], 0.125 + 2.0 + 0.25)
    np.testing.assert_allclose(t["cap"], 0.03 + 0.7 + 0.5)


def test_files_without_products_keep_the_legacy_classification(reader):
    reader._scan_inelastic_mts("Pb208")
    assert reader._available_inelastic_mts["Pb208"] == []
    assert reader._absorption_mts["Pb208"] == [102]


def test_neutron_yields(reader):
    assert reader.neutron_yield("Al27", 16, 1e7) == 2.0
    assert reader.neutron_yield("Al27", 22, 1e7) == 1.0
    assert reader.neutron_yield("Al27", 5, 5.5e6) == pytest.approx(1.0)
    assert reader.neutron_yield("Al27", 5, 1.4e7) == pytest.approx(1.5)
    assert reader.neutron_yield("Al27", 51, 1e7) == 1.0     # unknown -> 1


def test_reader_with_yields_pickles(reader):
    # the parallel engines send the reader to worker processes
    copy = pickle.loads(pickle.dumps(reader))
    assert copy.neutron_yield("Al27", 5, 1.4e7) == pytest.approx(1.5)


def test_multiplicity_is_unbiased(reader):
    rng = np.random.default_rng(1)
    for energy, mean in ((2e6, 0.5 + 1.0 / 9.0), (8e6, 0.5 + 7.0 / 9.0)):
        n = np.array([_neutron_multiplicity(reader, "Al27", 5, energy, rng)
                      for _ in range(20000)])
        assert set(np.unique(n)) <= {0, 1, 2}
        assert n.mean() == pytest.approx(mean, abs=4 * n.std() / np.sqrt(n.size))


def _sphere_run(engine, reader, n=2000, seed=5):
    mediums = [Region(surfaces=[Sphere((0, 0, 0), 30.0)], name="ball",
                      composition=[Nuclide("O16", 1e22, 15.86)])]
    settings = Settings("shielding", n, 10)
    rng = np.random.default_rng(seed)
    mu = 1.0 - 2.0 * rng.random(n)
    phi = 2.0 * np.pi * rng.random(n)
    s = np.sqrt(1.0 - mu * mu)
    src = {"x": np.zeros(n), "y": np.zeros(n), "z": np.zeros(n),
           "u": s * np.cos(phi), "v": s * np.sin(phi), "w": mu,
           "energy": np.full(n, 1.4e7)}
    if engine == "vector":
        res = vt.run_transport(reader, mediums, src, np.random.default_rng(seed + 1), settings)
        return (res["leaked_weight"].sum() / n, res["absorbed_weight"] / n,
                res["counters"]["removed_weight"] / n)
    leaked = absorbed = removed = 0.0
    for state, ss in zip(to_states(src), np.random.SeedSequence(seed + 1).spawn(n)):
        r = simulate_single_particle((state, reader, mediums, 1.0, 1.0, VelocitySampler(1.0),
                                      None, False, RNGHandler(ss), settings))
        leaked += r["final_weight"]
        absorbed += r["absorbed_weight"]
        removed += r["physics_counts"].get("removed_weight", 0.0)
    return leaked / n, absorbed / n, removed / n


@pytest.mark.parametrize("engine", ["vector", "scalar"])
def test_a_channel_emitting_no_neutron_removes_it(engine, reader):
    # Sigma_t = 0.04 /cm and P(removal | collision) = 3/4, so the fraction
    # leaking from the 30 cm sphere is well above the uncollided 30 %, and
    # every neutron either leaks or is removed. As in OpenMC, removal is not
    # absorption, and O16 here has no absorption reaction.
    leaked, absorbed, removed = _sphere_run(engine, reader)
    assert absorbed == 0.0
    assert leaked + removed == pytest.approx(1.0, abs=1e-12)
    uncollided = np.exp(-0.04 * 30.0)
    assert uncollided < leaked < 0.6


def test_engines_agree_on_removal(reader):
    lv = _sphere_run("vector", reader)[0]
    ls = _sphere_run("scalar", reader)[0]
    sigma = np.sqrt(lv * (1 - lv) / 2000)
    assert abs(lv - ls) < 4 * np.sqrt(2) * sigma


# ----------------------------------------------------------------- real data
@pytest.mark.skipif(not os.path.exists(os.path.join(ENDF, "neutron", "Fe56.h5")),
                    reason="ENDF data not present")
@pytest.mark.parametrize("nuclide", ["Fe56", "Al27", "O16", "Ni58", "B10", "U235"])
def test_total_matches_openmc(nuclide):
    openmc_data = pytest.importorskip("openmc.data")
    path = os.path.join(ENDF, "neutron", f"{nuclide}.h5")
    if not os.path.exists(path):
        pytest.skip(f"{nuclide}.h5 not present")
    ref = openmc_data.IncidentNeutron.from_hdf5(path)
    energies = np.array([1e-2, 1e3, 1e6, 5e6, 1e7, 1.4e7, 1.9e7])
    expected = ref[1].xs["294K"](energies)
    reader = CrossSectionReader(ENDF)
    reader._build_fast_table(nuclide)
    t = reader._fast_tables[nuclide]
    total = sum(np.interp(energies, t["grid"], t[k]) for k in ("el", "in", "cap", "fis"))
    np.testing.assert_allclose(total, expected, rtol=1e-8)
