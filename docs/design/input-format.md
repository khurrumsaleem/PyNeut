# Design: TOML input files and a tally system

Status: **PRs 1–3 implemented.** PR 1 is the `pyneut` package and PR 2 the
input format (`docs/input-format.md` is the user reference). PR 3 is the
tally interface (#6): `pyneut/tallies.py` and the engine hooks. Addresses
#26 (configuration file), #6 (tallies) and the API-ergonomics roadmap item.

**Where PR 3 differs from this proposal:**
- Flux has a track-length estimator as well as the collision estimator, and
  it is the default. It works in voids and has lower variance. A mesh flux
  still uses the collision estimator.
- Leakage became a per-neutron event score. The old per-history escape-energy
  histogram distorted spectra in multiplying problems.
- Tallies work in eigenvalue mode too, over the active generations.
- Detector tallies still take a box and run on the scalar engine.

**Where PR 2 differs from this proposal:**
- Detector tallies take a `box`, not a `cell`, because the kernel's detector
  is a box. They need the scalar engine, as do mesh tallies, until PR 3.
- `Model.from_dict` is the Python constructor. Separate `Material`, `Cell`
  and `Surface` builder classes were not needed.
- `pyneut plot` is not in PR 2.
- The source gained a `cylinder` space, needed for the pin-cell pellet.
- The Maxwellian takes `kT` as an energy rather than a temperature.
- On the vector engine, eigenvalue runs use one process, and `workers` is
  ignored with a warning.
- Unit strings deliberately leave out `meV`, so a case typo can't shift an
  energy by 10⁹.

## Goals

1. **One file describes one run** — data, materials, geometry, source,
   tallies and run settings — so a calculation is reproduced by keeping
   that file under version control. No Python needed for a standard run.
2. **Readable and diffable.** Comments allowed, one fact per line, names
   instead of numbers.
3. **Remove the known foot-guns** by construction: atomic-weight ratios
   are read from the data file (no more 207.2 vs 206.190), a missing
   nuclide file is an error (not a silent zero cross section), directions
   and sources come from named distributions (no hand-built `theta`).
4. **Both engines, one model.** The same file runs on the vector or the
   scalar engine; the engine is a setting.
5. **The Python API stays first-class.** The file is parsed into the same
   objects a script would build, so the validation suite and existing
   scripts keep working, and anything the file can say, Python can too.

Non-goals for the first version: natural-element expansion, multiple
temperatures, S(α,β), lattices/universes, general quadrics.

## Format: TOML

TOML (decided). Parsed with `tomllib` (standard library, Python 3.11+);
on 3.10 the identical `tomli` package is added as a conditional
dependency (`tomli; python_version < "3.11"`), so the 3.10 floor stays.

Units are fixed and documented: **cm, g/cm³, eV, K**. Numbers are in those
units. As a convenience, energies may also be written as strings with a unit
(`"2 MeV"`, `"25 keV"`, `"0.0253 eV"`); see open question 3.

## Example: the Pb208 validation sphere

```toml
# pb208_sphere.toml — 2 MeV point source in a 10 cm Pb208 sphere
[run]
mode      = "fixed-source"     # or "eigenvalue"
engine    = "vector"           # or "scalar"
particles = 100_000
seed      = 12345
data      = "../endfb"         # relative to this file

[materials.lead]
density  = 11.35               # g/cm3
nuclides = { Pb208 = 1.0 }     # atom fractions; AWR comes from Pb208.h5

[surfaces]
shell = { type = "sphere", center = [0, 0, 0], radius = 10.0 }

[cells]
shield = { region = "-shell", material = "lead" }

[source]
space     = { type = "point", xyz = [0, 0, 0] }
direction = "isotropic"
energy    = 2.0e6              # eV  (or "2 MeV")

[tallies.leakage]
type        = "leakage"
energy_bins = { log = [1e3, 2.1e7, 60] }   # 60 log-spaced bins
```

Run it with `pyneut run pb208_sphere.toml`.

## The sections

### `[run]`

| Key | Default | Meaning |
|-----|---------|---------|
| `mode` | — | `"fixed-source"` or `"eigenvalue"` |
| `engine` | `"vector"` | `"vector"` or `"scalar"` |
| `particles` | — | histories (fixed source) or neutrons per generation |
| `batches` | 10 | batch-means groups (fixed source) |
| `generations`, `inactive` | — | eigenvalue only |
| `seed` | 1 | master seed; per-history/per-worker streams derive from it |
| `transport` | `"implicit-capture"` | or `"analog"`; eigenvalue is always analog |
| `weight_cutoff`, `survival` | 1e-4, 0.1 | roulette parameters |
| `workers` | 1 | processes (`0` = all cores) |
| `data`, `temperature` | `"./endfb"`, `"294K"` | nuclear data |

### `[materials.<name>]`

Either a mass density with fractions, or explicit atom densities:

```toml
[materials.b4c]
density   = 2.52
nuclides  = { B10 = 4, C12 = 1 }        # atom fractions, need not sum to 1
# fractions = "weight"                  # optional; default "atom"

[materials.godiva]
atom_densities = { U234 = 4.8271e-4, U235 = 4.4994e-2, U238 = 2.4984e-3 }  # atoms/b-cm
```

Each nuclide's AWR is read from its HDF5 file (`atomic_weight_ratio`
attribute); molar masses are AWR × mₙ. A nuclide whose file is missing is a
load-time error naming the expected path. `void` is a built-in material name.

### `[surfaces]`

Each surface is named; `boundary` is `"vacuum"` (default, the neutron may
leave) or `"reflective"`.

| `type` | Parameters | `-name` (inside) means |
|--------|------------|------------------------|
| `plane` | `normal = [a,b,c]`, `offset = d` | a·x + b·y + c·z ≤ d |
| `x-plane`, `y-plane`, `z-plane` | `at = x0` | x ≤ x0 (etc.) |
| `sphere` | `center`, `radius` | inside the sphere |
| `cylinder` | `axis = "x"/"y"/"z"`, `center`, `radius` | inside the infinite cylinder |
| `box` | `min = [x,y,z]`, `max = [x,y,z]` | inside the box (a 6-plane macro-body) |

### `[cells]` and the region syntax

```toml
[surfaces]
outer  = { type = "cylinder", axis = "z", center = [0,0,0], radius = 10 }
top    = { type = "z-plane", at =  10 }
bottom = { type = "z-plane", at = -10 }
world  = { type = "box", min = [-20,-20,-20], max = [20,20,20] }

[cells]
shield = { region = "-outer & -top & +bottom", material = "lead" }
air    = { region = "-world & ~(-outer & -top & +bottom)", material = "void" }
```

- `-s` is the inside of surface `s`, `+s` the outside.
- `&` intersection, `|` union, `~` complement, parentheses group. Precedence:
  `~` > `&` > `|`.
- Any point in no cell is outside the model; a neutron that reaches it has
  leaked (the kernel's current behaviour).
- Overlaps: see open question 2.

Each expression compiles to the existing `Region` tree (intersection /
union / complement), so the transport kernels are untouched.

### `[source]`

Three independent distributions:

| Key | Options |
|-----|---------|
| `space` | `point {xyz}`, `box {min, max}`, `sphere {center, radius}` (uniform in volume), `cell {cell}` (uniform in a cell, by rejection in its bounding box) |
| `direction` | `"isotropic"`, `beam {direction = [u,v,w]}` |
| `energy` | a number, `discrete {energies, probabilities}`, `uniform {min, max}`, `watt {nuclide = "U235"}` or `watt {a, b}`, `maxwell {temperature}` |

The loader samples the source in whichever form the engine needs (dicts for
scalar, arrays for vector) — users never build per-particle states again.
For eigenvalue mode this is the first-generation source.

### `[tallies.<name>]`

This is where #6 lands. A tally is a **score** plus optional **filters**:

```toml
[tallies.leak_spectrum]
type        = "leakage"
energy_bins = { log = [1e3, 2.1e7, 60] }

[tallies.dose_map]
type  = "absorption"
mesh  = { min = [-15,-15,-15], max = [15,15,15], dims = [40,40,40] }

[tallies.shield_flux]
type        = "flux"              # collision estimator, phase 2
cells       = ["shield"]
energy_bins = [1e-5, 0.625, 1e5, 2e7]
```

| `type` | Scores | Filters |
|--------|--------|---------|
| `leakage` | weight leaving the model | `energy_bins` |
| `absorption` | weight removed by capture | `cells`, `mesh`, `energy_bins` |
| `flux` | collision estimator w/Σ_t (phase 2) | `cells`, `mesh`, `energy_bins` |
| `detector` | histories reaching a cell (existing binary count) | `cell` |
| `keff` | implicit in eigenvalue mode: collision and generation estimators | — |

Every tally reports mean and batch-means standard error per bin.

**What changes inside the engines.** Today the engines return fixed fields
(per-history leaked weight and mean escape energy; the scalar engine also
returns absorption positions, the vector engine only a total). Filters by
cell or mesh, and flux, need scoring *at the event*. Phase 2 therefore adds a
small scoring interface both engines call — `on_collision`, `on_absorb`,
`on_leak` with position, energy, weight and cell — that tallies subscribe to.
The vector engine calls it once per event pass with arrays, so the cost
stays vectorised. Validation: the new leakage tally must reproduce
`validation_results.csv` exactly at the same seed (scalar engine) before the
old fields are retired.

### Results

`pyneut run model.toml -o out/` writes:

- `out/results.json` — every tally (bins, mean, sem), k_eff when relevant,
  run time, histories/s, the engine, plus **provenance**: the resolved
  input (after defaults), its SHA-256, the PyNeut version and git commit,
  and the nuclear-data manifest (the `data_manifest` fingerprint the
  validation suite already uses). A result is then traceable to the exact
  input, code and data that produced it.
- `out/<tally>.vtk` for mesh tallies.
- A printed summary.

## Python API

```python
import pyneut

model = pyneut.Model.from_toml("pb208_sphere.toml")
model.run.particles = 20_000          # override anything before solving
results = model.solve()

leak = results["leak_spectrum"]
print(leak.total, leak.total_sem)
results.save("out/")
```

`Model` can also be built directly from `pyneut.Material`, `pyneut.Cell`,
`pyneut.Surface` … objects; `from_toml` is just one constructor. The current
low-level functions (`simulate_single_particle`, `run_transport`, `run_keff`)
remain and are what `Model.solve()` calls.

## CLI

```
pyneut run   model.toml [-o out/] [--engine scalar] [--particles N] [--seed S] [--workers N]
pyneut check model.toml          # parse, validate, resolve data, report; no transport
pyneut plot  model.toml [--slice z=0]   # material map via tools/plot_geometry
```

`pyneut check` catches, before any transport: unknown keys (with a
"did you mean" suggestion), undefined surface/material/cell names, missing
data files, a source outside every cell, and (by random point sampling)
overlapping cells and gaps.

## Package rename: `src` → `pyneut`

Decided. `src/` moves to `pyneut/`; every import in tests, validation
scripts, tools and docs changes in the same PR. `setup.py` installs the
`pyneut` package and the `pyneut` console script.

## Testing

- **Equivalence:** for each example file, the TOML path and a hand-built
  Python model produce *identical* scalar-engine results at the same seed
  (the scalar kernel is deterministic per seed). The validation cases become
  example files, and this check pins them to the current validated numbers.
- Parser tests: every section, every error message.
- Region-expression tests against `Region.contains` on random points.
- Source-sampler tests: moments of each distribution.

## Delivery plan

| PR | Contents |
|----|----------|
| 1 | Rename `src` → `pyneut`; console-script skeleton. Mechanical, no behaviour change. |
| 2 | `Model`, TOML loader, region parser, sources, `pyneut check/run`, existing tally types (leakage + spectrum, k_eff, scalar mesh, detector), `results.json` with provenance, `examples/` (Pb sphere, slab + detector, B4C mixture, Godiva, pin cell), equivalence tests. Closes #26. |
| 3 | Engine scoring interface; cell/mesh/energy filters in both engines; collision-estimator flux; vector-engine mesh tallies. Closes #6. |
| 4 (optional) | Validation suite reads its cases from `examples/` instead of Python dicts. |

## Decisions

All five open questions were settled with the proposed answers:

1. Intersection is written with an explicit `&`; a space is an error.
2. Cells must not overlap unless they have different `priority` values
   (unset means 0). `pyneut check` reports overlaps.
3. Energies accept unit strings (`"2 MeV"`); lengths are plain cm.
4. The default engine is the vector engine.
5. The rename to `pyneut` was a clean break, with no `src` shim.
