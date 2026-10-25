<p align="center">
  <img src="docs/assets/pyneut-logo.jpg" alt="PyNeut logo" width="420">
</p>

# PyNeut — Python Neutronics

A continuous-energy Monte Carlo neutron transport code in Python. It reads
ENDF/B-VIII.0 data directly (OpenMC HDF5 format) and is validated against
OpenMC, against analytic results and against ICSBEP critical benchmarks.
**For education and research prototyping, not for production or safety work.**

[![tests](https://github.com/ybadr16/PyNeut/actions/workflows/tests.yml/badge.svg)](https://github.com/ybadr16/PyNeut/actions/workflows/tests.yml)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

<p align="center">
  <img src="docs/images/validation_summary.png"
       alt="Every PyNeut validation case plotted against OpenMC as a z-score; twelve of thirteen fall inside the two-sigma agreement band"
       width="88%">
</p>

## What it does

- **Physics:** elastic scattering with tabulated ENDF angular distributions,
  plus free-gas thermal motion below 10 eV. Discrete inelastic levels, and
  continuum and (n,2n)/(n,3n) emission through ENDF Law 61 / Kalbach–Mann.
  Capture, unresolved-resonance probability tables, and fission with ν̄ and a
  Watt spectrum.
- **Problems:** fixed-source shielding (analog, or implicit capture with Russian
  roulette) and k-eigenvalue by fission-source power iteration.
- **Geometry:** constructive solid geometry from planes, spheres, axis-aligned
  cylinders and boxes. Regions can be nested, have priorities, be voids, or
  have reflective boundaries.
- **Materials:** one isotope or a multi-isotope mixture per region. The reacting
  isotope is sampled at each collision.
- **Two engines that share the same physics:**
  - `simulation.py` follows one neutron at a time. It is the readable
    reference implementation.
  - `vector_transport.py` advances the whole particle bank one event at a time
    with NumPy. It is 5–64× faster on a single core.
  Both parallelise across cores with `multiprocessing`.
- **Outputs:** leakage and escape-energy tallies with batch-means errors, a
  Cartesian mesh tally written as VTK for ParaView, and trajectory recording.

## Validation at a glance

| Check | Result |
|-------|--------|
| 13 fixed-source cases vs OpenMC 0.15 (Pb208, Fe56, Be9, Al27, C12, O16, B10, Cd112; thermal to 14 MeV) | **12 OK, 1 warning, 0 failures** (OK means within 2σ) |
| ICSBEP Godiva / Jezebel-23 / Jezebel | k_eff within **0.8σ** of the benchmark and **1.3σ** of OpenMC |
| Bare U235 spheres, R = 6 / 8.7 / 11 cm, vs OpenMC eigenvalue mode | z = 1.5 / 2.0 / 2.2 |
| BEAVRS pin cell, k∞ with matched physics (5 seeds per code) | **+72 ± 98 pcm**, z = 0.74 |
| B10 + C12 mixture vs OpenMC | leakage within 0.1σ, escape energy within 0.4σ |
| Analytic: Beer–Lambert attenuation, elastic [α,1] and ξ | 13 / 13 pass |
| Cross-section reader vs `openmc.data` (elastic, capture, (n,2n)) | identical (0.000 %) for every isotope |

Three spectrum shapes are POOR even though their integrals agree: the two
(n,2n) spheres and B10 at 1 keV. Fe56 at 14 MeV is the one integral warning.
The [validation report](docs/index.md) has the full tables, the method, and
what is and is not understood about these cases.

<p align="center">
  <img src="docs/images/icsbep_benchmarks.png"
       alt="PyNeut and OpenMC k_eff for the Godiva, Jezebel-23 and Jezebel benchmarks against the published unity value"
       width="72%">
</p>

## Install

```bash
git clone https://github.com/ybadr16/PyNeut.git
cd PyNeut
pip install -e .               # runtime: numpy, h5py
pip install -e ".[test]"       # + pytest, trimesh, rtree
pip install -e ".[validation]" # + openmc, matplotlib (OpenMC comparison suite)
```

The package is imported as `src` (`from src.simulation import ...`).

### Nuclear data

The ENDF/B-VIII.0 files are **not** in the repository. Download the OpenMC HDF5
library from [openmc.org/data-libraries](https://openmc.org/data-libraries/)
and put the `*.h5` files in `endfb/neutron/` (for example
`endfb/neutron/Pb208.h5`). The fixed-source suite needs Al27, B10, Be9, C12,
Cd112, Fe56, O16 and Pb208. The k-eigenvalue cases add U233/U234/U235/U238,
Pu239/Pu240/Pu241 and Ga, and the pin cell needs its 36 nuclides.

> A nuclide name with no matching file gets a zero cross section, and no error
> is raised. Check that the spelling matches the file name.

## Quick start

A 2 MeV isotropic point source at the centre of a 10 cm Pb208 sphere:

```python
import numpy as np
from multiprocessing import Pool

from src.cross_section_read import CrossSectionReader
from src.material import Material, NEUTRON_MASS_AMU
from src.medium import Region, Sphere
from src.random_number_generator import RNGHandler
from src.settings import Settings
from src.simulation import simulate_single_particle
from src.tally import Tally
from src.vt_calc import VelocitySampler

AWR = 206.190  # Pb208 atomic-weight ratio (target mass / neutron mass), from the data file

if __name__ == "__main__":
    reader = CrossSectionReader("./endfb")
    lead = Material("Pb208", density=11.35,
                    atomic_mass=AWR * NEUTRON_MASS_AMU, atomic_weight_ratio=AWR)
    sampler = VelocitySampler(lead.kg_mass)

    regions = [Region(surfaces=[Sphere((0, 0, 0), 10.0)],
                      name="Lead", priority=1, element="Pb208")]

    settings = Settings(mode="shielding", particles=10_000)   # implicit capture
    n = settings.num_particles
    rngs = [RNGHandler(seed=12345 + i) for i in range(n)]

    # theta is the polar angle from +z; arccos(1 - 2ξ) makes the source isotropic
    states = [{"x": 0.0, "y": 0.0, "z": 0.0,
               "theta": np.arccos(1 - 2 * r.random()),
               "phi": r.uniform(0, 2 * np.pi),
               "energy": 2.0e6,          # eV
               "weight": 1.0,
               "has_interacted": False} for r in rngs]

    # simulate_single_particle takes one 10-tuple per source neutron
    args = [(s, reader, regions, lead.atomic_weight_ratio, lead.number_density,
             sampler, None, False, r, settings) for s, r in zip(states, rngs)]

    with Pool() as pool:
        results = pool.map(simulate_single_particle, args)

    tally = Tally()
    for res in results:
        tally.merge_partial_results(res)
    tally.print_summary(n)          # leakage and mean escape energy, ± batch-means σ
```

The [Quick Start guide](docs/QUICK_START.md) covers slabs and detectors, mesh
tallies, material mixtures, reflective boundaries, k-eigenvalue and the vector
engine. `main.py` is a complete mesh-tally run.

## Repository layout

```
src/                       the code (imported as `src`)
  cross_section_read.py      ENDF HDF5 reader: XS, angular/secondary distributions, ν̄, URR
  angular_distribution.py    tabulated elastic (MT 2) cosine sampling
  material.py                Material, Nuclide, Material.mixture
  medium.py                  CSG: Region, Plane, Cylinder, Sphere, Box
  geometry.py                ray–surface crossing search, reflection
  physics.py                 elastic/inelastic kinematics, evaporation, Watt χ
  vt_calc.py                 free-gas target-velocity sampling
  simulation.py              one-neutron-at-a-time transport kernel (reference engine)
  vector_transport.py        event-based NumPy engine + multiprocessing wrappers
  criticality.py             k-eigenvalue power iteration (scalar engine)
  settings.py, tally.py, statistics.py, mesh.py, random_number_generator.py
tests/                     pytest suite (165 tests; 36 need ENDF data and skip without it)
Validation/
  analytic_benchmarks.py     exact analytic checks, no OpenMC needed
  OpenMC_Comparison/         live PyNeut-vs-OpenMC suite, see its README
tools/                     figure scripts (validation summary, geometry maps, spectra)
docs/                      validation report, quick start, API reference, images
Animation/                 Manim trajectory replay
main.py                    end-to-end example with a mesh tally
```

## Testing

```bash
pytest                                                  # unit tests
python Validation/analytic_benchmarks.py                # analytic benchmarks (needs data)
cd Validation/OpenMC_Comparison && python run_all.py    # vs OpenMC (needs data + openmc)
```

## Main limitations

- Thermal scattering is free-gas only. There is no S(α,β), so bound moderators
  such as water and graphite are approximate.
- Fission uses a Watt χ and total ν̄. There is no tabulated ENDF χ and no
  separate delayed neutrons.
- Neutrons only. Charged particles and photons are not transported.
- Sources are built by hand, as per-particle dicts or arrays.
- Only one temperature (294 K) is loaded.
- Geometry covers planes, spheres and axis-aligned cylinders and boxes, with no
  rotations.

The full list, with how each item shows up in the validation, is in
[docs/index.md](docs/index.md#limitations).

## Roadmap

1. **Animation pipeline.** Replace the bare-coordinate JSON export with a
   self-describing bundle: serialised geometry, plus per-vertex
   energy/event/weight tagged by `parent_id`. A generic 3-D Manim scene would
   then need no hand-editing.
2. **Criticality.** Read the tabulated ENDF χ, separate prompt and delayed ν̄,
   add reflected and thermal-solution benchmarks, and ship a small fissile data
   subset so the benchmarks run without the full library.
3. **Performance.** The vector engine is done. A Numba-compiled kernel is in
   progress on the `numba-experiment` branch.
4. **Physics.** Energy-correlated (n,xn) multiplicity, S(α,β), multi-temperature
   data, a track-length estimator, built-in source samplers and a
   natural-abundance database.
5. **Code hygiene.**
   - Replace the 10-element argument tuple with a dataclass.
   - Keep the element-alias map in one place.
   - Raise an error when a nuclide's data file is missing, instead of returning
     zero.
   - Give the package a real name instead of `src`.

## Documentation

- [Validation report](docs/index.md): method, every result table, physics
  models, limitations
- [Quick Start](docs/QUICK_START.md): worked examples
- [API reference](docs/api-reference.md)
- [Validation suite guide](Validation/OpenMC_Comparison/README.md)
- [Contributing](CONTRIBUTING.md)

## References

- ENDF/B-VIII.0: https://www.nndc.bnl.gov/endf/
- OpenMC theory manual: https://docs.openmc.org/en/stable/methods/index.html
- Duderstadt & Hamilton, *Nuclear Reactor Analysis*
- J. Wood, *Computational Methods in Reactor Shielding*

## License and contact

MIT, see [LICENSE](LICENSE). Youssef, Nuclear and Radiation Engineering,
Alexandria University. Bug reports and questions go to
[GitHub issues](https://github.com/ybadr16/PyNeut/issues).
