# Input format reference

A PyNeut model is one [TOML](https://toml.io) file. It describes the
materials, geometry, source, tallies and run settings of a calculation, so
keeping the file under version control is enough to reproduce a result.

```bash
pyneut check model.toml     # validate the input, data, source and geometry
pyneut run   model.toml     # solve; writes model_results/results.json
```

The [`examples/`](../examples/README.md) folder has an input for every
validated case and every feature. They are the quickest way to learn the
format.

**Units:** lengths in **cm**, densities in **g/cm³**, number densities in
**atoms/b-cm**, energies in **eV**. An energy can also be a string with a
unit: `"2 MeV"`, `"25 keV"`, `"0.0253 eV"`. The unit is case-sensitive, and
there is deliberately no `meV`, so a typo can't be off by a factor of 10⁹.

**Contents:**
[Example](#a-complete-example) · [`[run]`](#run) ·
[`[materials]`](#materials) · [`[surfaces]`](#surfaces) ·
[`[cells]`](#cells-and-region-expressions) · [`[source]`](#source) ·
[`[tallies]`](#tallies) · [Results](#results) · [CLI](#command-line) ·
[Python](#from-python) · [Errors and checks](#errors-and-checks)

## A complete example

```toml
# A 2 MeV point source at the centre of a 10 cm Pb208 sphere.
[run]
mode      = "fixed-source"
particles = 10_000
seed      = 12345

[materials.lead]
density  = 11.35                 # g/cm3
nuclides = { Pb208 = 1.0 }       # atom fractions

[surfaces]
shell = { type = "sphere", center = [0, 0, 0], radius = 10.0 }

[cells]
shield = { region = "-shell", material = "lead" }

[source]
space     = { type = "point", xyz = [0, 0, 0] }
direction = "isotropic"
energy    = "2 MeV"

[tallies.leakage]
type        = "leakage"
energy_bins = { log = ["1 keV", "2.1 MeV", 60] }
```

TOML notes:
- A `{ ... }` inline table must fit on one line. For anything longer, use a
  sub-table such as `[materials.core.atom_densities]` or `[source.space]`.
  The examples show both forms.
- Underscores are allowed in numbers: `100_000`.

## `[run]`

| Key | Default | Meaning |
|-----|---------|---------|
| `mode` | required | `"fixed-source"` or `"eigenvalue"` |
| `particles` | required | source histories (fixed source) or neutrons per generation (eigenvalue) |
| `engine` | `"vector"` | `"vector"` (fast) or `"scalar"` (one history at a time; the reference implementation) |
| `batches` | `10` | fixed source: number of batches for the batch-means standard errors |
| `generations`, `inactive` | required in eigenvalue mode | total generations, and the settling generations discarded before averaging |
| `seed` | `1` | master random seed |
| `transport` | `"implicit-capture"` | fixed source: or `"analog"`. Eigenvalue mode is always analog. |
| `weight_cutoff`, `survival` | `1e-4`, `0.1` | Russian-roulette threshold and survival probability (implicit capture) |
| `workers` | `1` | processes; `0` means every core |
| `data` | see below | nuclear-data directory, relative to the input file |
| `temperature` | `"294K"` | data temperature to load |

**Where the data comes from:** `run.data` if given, otherwise the
`PYNEUT_DATA` environment variable, otherwise `./endfb` in the working
directory. The directory must contain `neutron/<nuclide>.h5` files in OpenMC
format.

**Reproducibility:** `seed` feeds a NumPy `SeedSequence`, split into one
stream for sampling the source and one for transport. On the scalar engine,
each history gets its own child stream, so results don't depend on how
`workers` divides the work.

## `[materials]`

Each material is a named table. Give either a mass density with fractions:

```toml
[materials.b4c]
density   = 2.52
nuclides  = { B10 = 4, C12 = 1 }   # need not sum to 1
fractions = "atom"                 # default; or "weight"
```

or explicit atom densities in atoms/b-cm:

```toml
[materials.core.atom_densities]
U234 = 4.8271e-4
U235 = 4.4994e-2
U238 = 2.4984e-3
```

- Nuclide names match the data files: `U235` means `neutron/U235.h5`.
- Each nuclide's **atomic-weight ratio is read from its data file**, and its
  molar mass is AWR × mₙ. You never type an AWR.
- A nuclide with no data file is an error.
- Natural elements (`Fe`, `Pb`) are not expanded; list the isotopes.
- `void` is a reserved material name for empty cells.

## `[surfaces]`

Each surface has a name, a `type` with its parameters, and an optional
`boundary`: `"vacuum"` (the default; neutrons may leave) or `"reflective"`.
Names may use letters, digits and `_`.

| `type` | Parameters | `-name` (the inside) is |
|--------|------------|-------------------------|
| `sphere` | `center`, `radius` | inside the sphere |
| `cylinder` | `axis` (`"x"`/`"y"`/`"z"`), `center`, `radius` | inside the infinite cylinder |
| `box` | `min = [x,y,z]`, `max = [x,y,z]` | inside the box |
| `x-plane`, `y-plane`, `z-plane` | `at` | x ≤ at (or y, z) |
| `plane` | `normal = [a,b,c]`, `offset = d` | a·x + b·y + c·z ≤ d (the normal need not be a unit vector) |

```toml
[surfaces]
side   = { type = "cylinder", axis = "z", center = [0, 0, 0], radius = 10 }
top    = { type = "z-plane", at = 10 }
bottom = { type = "z-plane", at = -10 }

[surfaces.cube]
type     = "box"
min      = [-10, -10, -10]
max      = [10, 10, 10]
boundary = "reflective"
```

## `[cells]` and region expressions

A cell is a region of space filled with one material:

```toml
[cells]
shield = { region = "-side & -top & +bottom", material = "lead" }
air    = { region = "-world & ~(-side & -top & +bottom)", material = "void" }
```

| Syntax | Meaning |
|--------|---------|
| `-s` | the inside of surface `s` |
| `+s` | the outside of surface `s` |
| `A & B` | intersection |
| `A \| B` | union |
| `~A` | complement |
| `( ... )` | grouping |

- Precedence: `~` binds tightest, then `&`, then `|`.
- **Intersection must be written with `&`.** `"-a -b"` is an error rather
  than an implicit intersection, so a missing operator can't go unnoticed.
- A point inside no cell is outside the model. A neutron that gets there has
  leaked, so the model doesn't need an outer "graveyard" cell.

**Overlaps.** Cells should not overlap. `pyneut check` samples the geometry
and reports any overlap. If overlapping on purpose is more convenient (for
example nesting shells instead of subtracting them), give the cells different
`priority` values: the higher priority wins where they overlap, and an unset
priority counts as 0. `examples/features/priority_nesting.toml` and
`layered_shield.toml` build the same model both ways.

## `[source]`

A source has three independent parts:

```toml
[source]
space     = { type = "sphere", center = [0, 0, 0], radius = 5 }
direction = "isotropic"
energy    = { type = "watt", nuclide = "U235" }
```

| `space` type | Parameters | Neutrons start |
|--------------|------------|----------------|
| `point` | `xyz` | at a point |
| `box` | `min`, `max` | uniformly in a box |
| `sphere` | `center`, `radius` | uniformly in a sphere |
| `cylinder` | `axis`, `center`, `radius`, `range = [lo, hi]` (along the axis) | uniformly in a finite cylinder |
| `cell` | `cell`, `bounds = {min, max}` | uniformly in a cell (rejection sampling inside `bounds`, which must enclose the cell) |

| `direction` | Meaning |
|-------------|---------|
| `"isotropic"` | uniform over all directions (the default) |
| `{ type = "beam", direction = [u, v, w] }` | a fixed direction (normalised for you) |

| `energy` | Meaning |
|----------|---------|
| a number or `"2 MeV"` | monoenergetic |
| `{ type = "discrete", energies = [...], probabilities = [...] }` | lines (probabilities are normalised) |
| `{ type = "uniform", min = ..., max = ... }` | uniform between two energies |
| `{ type = "maxwell", kT = "1.3 MeV" }` | √E·exp(−E/kT) |
| `{ type = "watt", nuclide = "U235" }` | Watt fission spectrum (U235, U233, Pu239) |
| `{ type = "watt", a = ..., b = ... }` | Watt with explicit a (eV) and b (1/eV) |

In eigenvalue mode, the source is the first-generation guess.

## `[tallies]`

Each tally is named and has a `type`: what it scores. Optional **filters**
split the score into bins:
- `cells`: a list of cell names;
- `energy_bins`;
- `mesh`: a regular 3-D grid.

Filters combine, so a flux tally with `cells` and `energy_bins` gives a
spectrum per cell. Results are per source neutron, with batch-means standard
errors. In eigenvalue mode, they're per fission-source neutron, averaged over
the active generations.

| `type` | Scores | Filters | Engines |
|--------|--------|---------|---------|
| `flux` | neutron flux | `cells`, `energy_bins`, `mesh`, `estimator` | both |
| `absorption` | weight removed by capture and fission | `cells`, `energy_bins`, `mesh` | both |
| `leakage` | weight leaving the model, scored per neutron at the energy it leaves with | `energy_bins` | both |
| `detector` | fraction of histories that reach a box | `box = {min, max}` | scalar, fixed source |
| `keff` | k_eff (collision estimator) and the generation estimator, with every generation | none | both, eigenvalue mode |

```toml
[tallies.flux]                       # a flux spectrum in two cells
type        = "flux"
cells       = ["fuel", "water"]
energy_bins = ["1e-5 eV", "0.625 eV", "100 keV", "20 MeV"]

[tallies.flux_map]                   # a 3-D flux map
type = "flux"
mesh = { min = [-10, -10, -10], max = [10, 10, 10], dims = [20, 20, 20] }
```

- **Flux** is reported as track length per source neutron (cm), meaning
  integrated over each bin's volume, as OpenMC reports it. Divide by a
  cell's volume to get its average flux.
- **Flux estimators:**
  - `track-length` (the default) adds up flight-path lengths. It works in
    void cells and usually has the lower variance.
  - `collision` adds weight / Σ_t at each collision. It is zero in voids.
  - A `mesh` flux uses the collision estimator; track-length on a mesh isn't
    available yet.
- **Energies** are the neutron's energy before the event: before a
  collision, or as it leaks.
- **`energy_bins`** is a list of edges, `{ log = [min, max, n] }` or
  `{ lin = [min, max, n] }`, where `n` is the number of bins. Events outside
  the bins, the listed cells or the mesh are not counted in `mean`. The
  tally's `total` always counts every event.
- **Defaults:** a fixed-source model with no tallies gets a `leakage` tally,
  and an eigenvalue model always gets a `keff` tally.
- **Detector boxes** are only tested where a neutron stops, at collisions and
  surface crossings. Make the detector its own void cell, as in
  `examples/features/slab_detector.toml`, so neutrons can't fly through it
  unseen.

Each tally's result has a `total` and `total_sem`. Leakage also has
`mean_energy` and `mean_energy_sem`, and flux has its `estimator`. A filtered
tally also has:
- `axes`, for example `["cell", "energy"]`;
- `mean` and `sem` arrays with those axes;
- `cells`, `energy_bins` and `mesh` as applicable.

## Results

`pyneut run model.toml` prints a summary and writes `model_results/`. Use
`-o` to choose the folder.

- **`results.json`:**
  - every tally (bins, values, standard errors), and k_eff with each generation;
  - run information: engine, wall time, histories per second, event counters;
  - **provenance:** the input file's SHA-256, the input with every default
    filled in, the PyNeut version and git commit, and the nuclear-data
    directory with its fingerprint (the same one the validation suite
    records) and the AWR read for each nuclide.
- **`<tally>.vtk` and `<tally>.npz`** for each mesh tally. The VTK opens in
  ParaView and is summed over any cell and energy axes. The NPZ holds the full
  `mean` and `sem` arrays and their `axes`. In `results.json`, mesh tallies
  list these files instead of their arrays.

## Command line

```
pyneut check model.toml [--data DIR] [--samples N]
pyneut run   model.toml [-o DIR] [--data DIR] [--engine vector|scalar]
                        [--particles N] [--seed S] [--workers N]
pyneut --version
```

`check` exits with status 1 if there is any error, so it can be used in
scripts and CI.

## From Python

```python
from pyneut.model import Model

model = Model.from_toml("examples/validation/fixed_source/pb208_elastic.toml")
model.run.particles = 50_000          # override anything before solving
model.run.engine = "scalar"
results = model.solve()
print(results.summary())
print(results["leakage"]["total"], results["leakage"]["total_sem"])
results.save("out/")
```

`Model.from_dict(...)` takes the same structure as a Python dict. The
low-level API that `Model.solve()` drives (`simulate_single_particle`,
`vector_transport.run_transport`, `criticality.run_keff`) is still available
and documented in the [API reference](api-reference.md).

## Errors and checks

Loading a model reports **every** problem at once, each with its location and
a suggestion where one helps:

```
model.toml: 4 problem(s):
  - run.engin: unknown key (did you mean 'engine'?)
  - cells.shield.region: missing operator between '-shell' and '-top'; write '&' for intersection or '|' for union
  - cells.shield.material: unknown material 'leed' (did you mean 'lead'?)
  - source.energy: cannot read energy '2 Mev'; write a number in eV or a number with a unit (eV, keV, MeV), e.g. "2 MeV"
```

`pyneut check` also finds:
- missing data files;
- a source that starts outside every cell;
- overlapping cells, by sampling random points in the geometry's bounding box.
